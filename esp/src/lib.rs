//! Layered DP counting partitions of {1..n*k} into n unlabeled equal-sum blocks of size k.
//!
//! Semantics match `../equal_sum_partitions.py` exactly (same states, same pruning, same
//! labeled counts), but counts are kept mod a prime p < 2^63. Per-layer fingerprints follow
//! `../FINGERPRINT_SPEC.md`, so they can be compared with the Python reference layer by layer.
//!
//! Key packing: a canonical state (sorted (slots, sum) blocks) is stored without its largest
//! block, which is rebuilt from the fixed layer totals Σs = rest, Σr = rest(rest+1)/2.
//! Each stored block is a dense code in 1..F (0 = absent) enumerating the feasible pairs
//! (s, r) with s(s+1)/2 <= r <= s(2·rest+1-s)/2; fields are bits(F-1) wide.
//! Keys are u64 when they fit in 63 bits, else u128.
//!
//! Memory: each layer is held either in RAM or as checksummed record files on disk. The next
//! layer is built in flat hash tables, one per shard (bits 40.. of the key hash). When the
//! tables would not fit the memory budget, the target shards are split into ranges built in
//! successive passes over the input, and each finished range is written to disk.

pub mod sys;

use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

pub const PRIMES: [u64; 3] = [(1 << 63) - 25, (1 << 63) - 165, (1 << 63) - 259];
pub const MAX_BLOCKS: usize = 16;

/// splitmix64 finalizer.
#[inline]
pub fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

/// (slots_left, sum_left) of one unfinished block.
pub type Block = (u32, u32);

/// u64 hash of a canonical state, as defined in FINGERPRINT_SPEC.md.
pub fn state_hash(blocks: &[Block]) -> u64 {
    let mut h = 0x9E3779B97F4A7C15u64 ^ blocks.len() as u64;
    for &(s, r) in blocks {
        assert!(s < 256 && r < 65536);
        h = mix64(h ^ (((s as u64) << 16) | r as u64));
    }
    h
}

#[inline]
pub fn add_mod(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b; // a, b < p < 2^63
    if s >= p { s - p } else { s }
}

#[inline]
pub fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((a as u128 * b as u128) % p as u128) as u64
}

pub fn pow_mod(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1 % p;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = mul_mod(r, b, p);
        }
        b = mul_mod(b, b, p);
        e >>= 1;
    }
    r
}

#[inline]
fn lo(s: u32) -> u32 {
    s * (s + 1) / 2
}

fn bits(x: u32) -> u32 {
    32 - x.leading_zeros()
}

pub trait Key: Copy + Eq + Send + Sync + 'static {
    /// Never a valid key: layouts use at most BYTES*8-1 bits.
    const EMPTY: Self;
    const BYTES: usize;
    fn from_u128(x: u128) -> Self;
    fn to_u128(self) -> u128;
    fn fold(self) -> u64;
    fn write_le(self, out: &mut Vec<u8>);
    fn read_le(b: &[u8]) -> Self;
}

impl Key for u64 {
    const EMPTY: Self = u64::MAX;
    const BYTES: usize = 8;
    #[inline]
    fn from_u128(x: u128) -> Self {
        x as u64
    }
    #[inline]
    fn to_u128(self) -> u128 {
        self as u128
    }
    #[inline]
    fn fold(self) -> u64 {
        self
    }
    #[inline]
    fn write_le(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }
    #[inline]
    fn read_le(b: &[u8]) -> Self {
        u64::from_le_bytes(b.try_into().unwrap())
    }
}

impl Key for u128 {
    const EMPTY: Self = u128::MAX;
    const BYTES: usize = 16;
    #[inline]
    fn from_u128(x: u128) -> Self {
        x
    }
    #[inline]
    fn to_u128(self) -> u128 {
        self
    }
    #[inline]
    fn fold(self) -> u64 {
        (self as u64) ^ mix64((self >> 64) as u64)
    }
    #[inline]
    fn write_le(self, out: &mut Vec<u8>) {
        out.extend_from_slice(&self.to_le_bytes());
    }
    #[inline]
    fn read_le(b: &[u8]) -> Self {
        u128::from_le_bytes(b.try_into().unwrap())
    }
}

const KEY_SALT: u64 = 0x2545F4914F6CDD1D;

#[inline]
fn key_hash<K: Key>(k: K) -> u64 {
    mix64(k.fold() ^ KEY_SALT)
}

#[inline]
fn shard_of(h: u64, shard_bits: u32) -> usize {
    ((h >> 40) & ((1u64 << shard_bits) - 1)) as usize
}

/// Linear-probing table of (key, count mod p). Slot index comes from hash bits 0..40, which
/// are independent of the shard bits 40..52.
pub struct FlatTable<K: Key> {
    slots: Vec<(K, u64)>,
    len: usize,
}

impl<K: Key> FlatTable<K> {
    pub fn with_entries(n: usize) -> Self {
        FlatTable { slots: vec![(K::EMPTY, 0); (n + n / 4).max(16)], len: 0 }
    }

    #[inline]
    fn slot(h: u64, cap: usize) -> usize {
        (((h << 24) as u128 * cap as u128) >> 64) as usize
    }

    #[inline]
    pub fn add(&mut self, k: K, h: u64, c: u64, p: u64) {
        if (self.len + 1) * 10 > self.slots.len() * 9 {
            self.grow();
        }
        let cap = self.slots.len();
        let mut i = Self::slot(h, cap);
        loop {
            let e = &mut self.slots[i];
            if e.0 == k {
                e.1 = add_mod(e.1, c, p);
                return;
            }
            if e.0 == K::EMPTY {
                *e = (k, c);
                self.len += 1;
                return;
            }
            i += 1;
            if i == cap {
                i = 0;
            }
        }
    }

    fn grow(&mut self) {
        let old = std::mem::take(&mut self.slots);
        self.slots = vec![(K::EMPTY, 0); old.len() * 3 / 2 + 16];
        self.len = 0;
        for (k, c) in old {
            if k != K::EMPTY {
                self.add(k, key_hash(k), c, u64::MAX);
            }
        }
    }

    #[inline]
    pub fn get(&self, k: K, h: u64) -> Option<u64> {
        let cap = self.slots.len();
        let mut i = Self::slot(h, cap);
        loop {
            let e = self.slots[i];
            if e.0 == k {
                return Some(e.1);
            }
            if e.0 == K::EMPTY {
                return None;
            }
            i += 1;
            if i == cap {
                i = 0;
            }
        }
    }

    /// Occupied entries, compacted in place.
    pub fn into_compact(mut self) -> Vec<(K, u64)> {
        let mut w = 0;
        for r in 0..self.slots.len() {
            if self.slots[r].0 != K::EMPTY {
                self.slots[w] = self.slots[r];
                w += 1;
            }
        }
        self.slots.truncate(w);
        self.slots.shrink_to_fit();
        self.slots
    }
}

/// Bit layout of packed keys for the layer where 1..rest remain.
#[derive(Clone, Debug)]
pub struct Layout {
    pub n: usize,
    pub rest: u32,
    fb: u32,
    pub bits: u32,
    base: Vec<u32>,   // base[s] = code of (s, lo(s)); base[s] = F for s > min(k, rest)
    dec: Vec<Block>,  // code -> (s, r); dec[0] unused
}

impl Layout {
    pub fn new(n: usize, k: u32, rest: u32) -> Self {
        let smax = k.min(rest);
        let mut base = vec![0u32; k as usize + 2];
        let mut dec = vec![(0, 0)];
        for s in 1..=smax {
            base[s as usize] = dec.len() as u32;
            dec.extend((0..=s * (rest - s)).map(|off| (s, lo(s) + off)));
        }
        let f = dec.len() as u32;
        for b in &mut base[smax as usize + 1..] {
            *b = f;
        }
        let fb = bits(f - 1);
        let total = (n as u32 - 1) * fb;
        assert!(total <= 127, "key needs {total} bits (n={n}, k={k}, rest={rest})");
        Layout { n, rest, fb, bits: total, base, dec }
    }

    pub fn wide(&self) -> bool {
        self.bits > 63
    }

    /// Pack a canonical state (sorted, non-empty unless rest == 0).
    #[inline]
    pub fn encode(&self, blocks: &[Block]) -> u128 {
        let mut key = 0u128;
        if let Some((_, stored)) = blocks.split_last() {
            for (i, &(s, r)) in stored.iter().enumerate() {
                let code = self.base[s as usize] + (r - lo(s));
                assert!(s >= 1 && code < self.base[s as usize + 1], "block ({s},{r}) out of range");
                key |= (code as u128) << (i as u32 * self.fb);
            }
        }
        key
    }

    /// Like `encode`, but None if a stored block is outside this layer's feasible range
    /// (such a state cannot occur in the layer).
    pub fn try_encode(&self, blocks: &[Block]) -> Option<u128> {
        let mut key = 0u128;
        if let Some((_, stored)) = blocks.split_last() {
            for (i, &(s, r)) in stored.iter().enumerate() {
                if s == 0 || s as usize + 1 >= self.base.len() || r < lo(s) {
                    return None;
                }
                let code = self.base[s as usize] + (r - lo(s));
                if code >= self.base[s as usize + 1] {
                    return None;
                }
                key |= (code as u128) << (i as u32 * self.fb);
            }
        }
        Some(key)
    }

    /// Unpack into `out`, rebuilding the dropped largest block; returns the block count.
    #[inline]
    pub fn decode(&self, key: u128, out: &mut [Block; MAX_BLOCKS]) -> usize {
        let mask = (1u128 << self.fb) - 1;
        let (mut len, mut ss, mut rs) = (0, 0, 0);
        for i in 0..self.n - 1 {
            let code = ((key >> (i as u32 * self.fb)) & mask) as usize;
            if code == 0 {
                break;
            }
            let (s, r) = self.dec[code];
            out[len] = (s, r);
            len += 1;
            ss += s;
            rs += r;
        }
        if ss < self.rest {
            out[len] = (self.rest - ss, self.rest * (self.rest + 1) / 2 - rs);
            len += 1;
        }
        len
    }
}

// ---------------------------------------------------------------------------------------
// Layer storage

/// Records per disk read chunk.
const CHUNK: usize = 1 << 22;

/// One record file: (key LE, count LE) records, with an order-independent checksum.
pub struct Part {
    path: PathBuf,
    records: u64,
    checksum: u64,
    bytes: u64,
}

pub enum Data<K: Key> {
    Mem(Vec<Vec<(K, u64)>>),
    Disk(Vec<Part>),
}

impl<K: Key> Data<K> {
    fn mem_bytes(&self, states: u64) -> f64 {
        match self {
            Data::Mem(_) => states as f64 * size_of::<(K, u64)>() as f64,
            Data::Disk(_) => 3.0 * (CHUNK * (K::BYTES + 8)) as f64,
        }
    }
}

pub enum AnyData {
    Narrow(Data<u64>),
    Wide(Data<u128>),
}

#[inline]
fn rec_check<K: Key>(k: K, c: u64) -> u64 {
    mix64(k.fold() ^ mix64(c ^ 0x5851F42D4C957F2D))
}

/// Tracks bytes of live part files and the peak.
#[derive(Default)]
pub struct DiskUse {
    pub now: u64,
    pub peak: u64,
}

fn write_part<K: Key>(path: &Path, shards: &[Vec<(K, u64)>], disk: &mut DiskUse) -> Part {
    let mut w = BufWriter::with_capacity(1 << 23, File::create(path).expect("create part file"));
    let (mut records, mut checksum) = (0u64, 0u64);
    let mut buf = Vec::with_capacity(1 << 20);
    for sh in shards {
        for &(k, c) in sh {
            k.write_le(&mut buf);
            buf.extend_from_slice(&c.to_le_bytes());
            checksum = checksum.wrapping_add(rec_check(k, c));
            if buf.len() >= 1 << 20 {
                w.write_all(&buf).expect("write part file");
                buf.clear();
            }
        }
        records += sh.len() as u64;
    }
    w.write_all(&buf).expect("write part file");
    w.flush().expect("flush part file");
    let bytes = records * (K::BYTES as u64 + 8);
    disk.now += bytes;
    disk.peak = disk.peak.max(disk.now);
    Part { path: path.to_path_buf(), records, checksum, bytes }
}

/// Stream all parts in chunks; a reader thread overlaps IO with `f`. Verifies record counts,
/// file lengths and checksums.
fn read_parts<K: Key>(parts: &[Part], mut f: impl FnMut(&[(K, u64)])) {
    let rec = K::BYTES + 8;
    std::thread::scope(|s| {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<(K, u64)>>(2);
        s.spawn(move || {
            let mut bytes = vec![0u8; CHUNK * rec];
            for part in parts {
                let mut r = BufReader::with_capacity(1 << 23, File::open(&part.path).expect("open part"));
                let (mut left, mut sum) = (part.records, 0u64);
                while left > 0 {
                    let nrec = left.min(CHUNK as u64) as usize;
                    r.read_exact(&mut bytes[..nrec * rec]).expect("read part");
                    let v: Vec<(K, u64)> = bytes[..nrec * rec]
                        .chunks_exact(rec)
                        .map(|b| (K::read_le(&b[..K::BYTES]), u64::from_le_bytes(b[K::BYTES..].try_into().unwrap())))
                        .collect();
                    for &(k, c) in &v {
                        sum = sum.wrapping_add(rec_check(k, c));
                    }
                    left -= nrec as u64;
                    if tx.send(v).is_err() {
                        return;
                    }
                }
                assert_eq!(r.read(&mut [0u8; 1]).unwrap(), 0, "trailing bytes in {:?}", part.path);
                assert_eq!(sum, part.checksum, "checksum mismatch in {:?}", part.path);
            }
        });
        for v in rx {
            f(&v);
        }
    });
}

fn delete_parts(parts: Vec<Part>, disk: &mut DiskUse) {
    for p in parts {
        std::fs::remove_file(&p.path).expect("remove part file");
        disk.now -= p.bytes;
    }
}

// ---------------------------------------------------------------------------------------
// Transitions

struct StepCtx<'a> {
    lin: &'a Layout,
    lout: &'a Layout,
    m: u32,
    lo_t: Vec<u32>,
    hi_t: Vec<i64>,
    p: u64,
    shard_bits: u32,
}

/// Per-thread transition buffers for the shards of one pass, flushed into the locked tables.
struct Buffers<'a, K: Key> {
    tables: &'a [Mutex<FlatTable<K>>],
    first: usize,
    bufs: Vec<Vec<(K, u64, u64)>>,
    shard_bits: u32,
    p: u64,
}

const FLUSH: usize = 256;

impl<'a, K: Key> Buffers<'a, K> {
    fn new(tables: &'a [Mutex<FlatTable<K>>], first: usize, shard_bits: u32, p: u64) -> Self {
        Buffers { tables, first, bufs: (0..tables.len()).map(|_| Vec::new()).collect(), shard_bits, p }
    }

    #[inline]
    fn push(&mut self, k: K, c: u64) {
        let h = key_hash(k);
        let sh = shard_of(h, self.shard_bits).wrapping_sub(self.first);
        if sh >= self.tables.len() {
            return; // another pass builds this shard
        }
        self.bufs[sh].push((k, h, c));
        if self.bufs[sh].len() >= FLUSH {
            self.flush(sh);
        }
    }

    fn flush(&mut self, sh: usize) {
        let p = self.p;
        let mut t = self.tables[sh].lock().unwrap();
        for (k, h, c) in self.bufs[sh].drain(..) {
            t.add(k, h, c, p);
        }
    }
}

impl<K: Key> Drop for Buffers<'_, K> {
    fn drop(&mut self) {
        for sh in 0..self.bufs.len() {
            if !self.bufs[sh].is_empty() {
                self.flush(sh);
            }
        }
    }
}

/// Place m into every block type of every state in `chunk`; returns candidate transitions.
fn expand<KI: Key, KO: Key>(chunk: &[(KI, u64)], cx: &StepCtx, buf: &mut Buffers<KO>) -> u64 {
    let (m, p, rest) = (cx.m, cx.p, cx.m - 1);
    let mut b = [(0, 0); MAX_BLOCKS];
    let mut nb = [(0, 0); MAX_BLOCKS];
    let mut t = 0u64;
    for &(key, cnt) in chunk {
        let len = cx.lin.decode(key.to_u128(), &mut b);
        let mut i = 0;
        while i < len {
            let blk = b[i];
            let mut j = i + 1;
            while j < len && b[j] == blk {
                j += 1;
            }
            let mult = j - i;
            let (s, r) = blk;
            if r >= m {
                t += 1;
                let nblk = (s - 1, r - m);
                if nblk.0 != 0 || nblk.1 == 0 {
                    // b without index i, plus nblk unless it is (0, 0), kept sorted.
                    let mut nl = 0;
                    let mut placed = nblk.0 == 0;
                    for (idx, &x) in b[..len].iter().enumerate() {
                        if idx == i {
                            continue;
                        }
                        if !placed && nblk < x {
                            nb[nl] = nblk;
                            nl += 1;
                            placed = true;
                        }
                        nb[nl] = x;
                        nl += 1;
                    }
                    if !placed {
                        nb[nl] = nblk;
                        nl += 1;
                    }
                    let mut ok = true;
                    let mut slots = 0;
                    for &(bs, br) in &nb[..nl] {
                        if br < cx.lo_t[bs as usize] || br as i64 > cx.hi_t[bs as usize] {
                            ok = false;
                            break;
                        }
                        slots += bs;
                    }
                    if ok {
                        assert_eq!(slots, rest);
                        let mut w = 0;
                        for _ in 0..mult {
                            w = add_mod(w, cnt, p);
                        }
                        buf.push(KO::from_u128(cx.lout.encode(&nb[..nl])), w);
                    }
                }
            }
            i = j;
        }
    }
    t
}

/// Build the target shards in `range` from the whole input; returns compacted shards.
fn run_pass<KI: Key, KO: Key>(input: &Data<KI>, cx: &StepCtx, range: Range<usize>,
                              per_shard: usize) -> (Vec<Vec<(KO, u64)>>, u64) {
    let tables: Vec<Mutex<FlatTable<KO>>> =
        range.clone().map(|_| Mutex::new(FlatTable::with_entries(per_shard))).collect();
    let transitions = AtomicU64::new(0);
    let init = || Buffers::new(&tables, range.start, cx.shard_bits, cx.p);
    let work = |buf: &mut Buffers<KO>, chunk: &[(KI, u64)]| {
        transitions.fetch_add(expand(chunk, cx, buf), Ordering::Relaxed);
    };
    match input {
        Data::Mem(shards) => shards.par_iter().for_each_init(init, |b, sh| work(b, sh)),
        Data::Disk(parts) => read_parts(parts, |v| v.par_chunks(1 << 15).for_each_init(init, work)),
    }
    let out = tables.into_par_iter().map(|t| t.into_inner().unwrap().into_compact()).collect();
    (out, transitions.into_inner())
}

fn fingerprint_impl<K: Key>(shards: &[Vec<(K, u64)>], lay: &Layout, p: u64) -> (u64, u64) {
    shards
        .par_iter()
        .map(|sh| {
            let mut b = [(0, 0); MAX_BLOCKS];
            let (mut cs, mut hs) = (0, 0);
            for &(key, c) in sh {
                let len = lay.decode(key.to_u128(), &mut b);
                cs = add_mod(cs, c, p);
                hs = add_mod(hs, mul_mod(state_hash(&b[..len]) % p, c, p), p);
            }
            (cs, hs)
        })
        .reduce(|| (0, 0), |a, b| (add_mod(a.0, b.0, p), add_mod(a.1, b.1, p)))
}

// ---------------------------------------------------------------------------------------
// Driver

pub fn is_feasible(n: usize, k: u32) -> bool {
    n >= 1 && k >= 1 && (k as u64 * (n as u64 * k as u64 + 1)) % 2 == 0
}

#[derive(Clone, Debug)]
pub struct RunOpts {
    pub n: usize,
    pub k: u32,
    pub p: u64,
    pub fingerprint: bool,
    pub max_states: Option<u64>,
    pub shard_bits: u32,
    pub verbose: bool,
    /// Target bytes for in-RAM layers plus hash tables.
    pub mem_budget: u64,
    /// Parent directory for spilled layers (a per-run subdirectory is created and removed).
    pub work_dir: PathBuf,
}

impl RunOpts {
    pub fn new(n: usize, k: u32, p: u64) -> Self {
        RunOpts { n, k, p, fingerprint: false, max_states: None, shard_bits: 12, verbose: false,
                  mem_budget: 4 << 30, work_dir: std::env::temp_dir() }
    }
}

#[derive(Clone, Debug)]
pub struct LayerStat {
    pub m: u32,
    pub states: u64,
    pub transitions: u64,
    pub seconds: f64,
    pub key_bits: u32,
    pub passes: u32,
    pub input_on_disk: bool,
    pub output_on_disk: bool,
    pub fp: Option<(u64, u64)>,
}

#[derive(Clone, Debug, Default)]
pub struct RunResult {
    pub labeled_mod_p: u64,
    pub peak_states: u64,
    pub total_seconds: f64,
    pub aborted: bool,
    pub peak_disk_bytes: u64,
    pub layers: Vec<LayerStat>,
}

struct Ctx<'a> {
    o: &'a RunOpts,
    dir: PathBuf,
    disk: DiskUse,
}

impl Ctx<'_> {
    fn part_path(&mut self, m: u32, pass: u32) -> PathBuf {
        std::fs::create_dir_all(&self.dir).expect("create work dir");
        self.dir.join(format!("layer{m:03}_part{pass:04}.bin"))
    }
}

/// Re-reading the input more often than this costs more than going over budget.
const MAX_PASSES: usize = 64;

/// Builds the layer after placing m. Returns (data, stats, fingerprint seconds, kept input).
/// With `keep_input` the input is moved to disk (if needed) and its parts are returned
/// instead of deleted.
fn layer_step<KI: Key, KO: Key>(input: Data<KI>, in_states: u64, est_out: f64, cx: &StepCtx,
                                ctx: &mut Ctx, keep_input: bool)
                                -> (Data<KO>, LayerStat, f64, Option<Data<KI>>) {
    let budget = ctx.o.mem_budget as f64;
    let nsh = 1usize << cx.shard_bits;
    let table_bytes = size_of::<(KO, u64)>() as f64 * 1.3;
    let mut input = input;
    let mut avail = budget - input.mem_bytes(in_states);
    let single = est_out * table_bytes <= avail;
    if keep_input || (!single && avail < budget / 2.0) {
        if let Data::Mem(shards) = &input {
            // Not enough room left for useful passes: move the input to disk first.
            let path = ctx.part_path(cx.m, 9999);
            let part = write_part(&path, shards, &mut ctx.disk);
            input = Data::Disk(vec![part]);
            avail = budget - input.mem_bytes(in_states);
        }
    }

    let mut st = LayerStat { m: cx.m, states: 0, transitions: 0, seconds: 0.0, key_bits: cx.lout.bits,
                             passes: 0, input_on_disk: matches!(input, Data::Disk(_)),
                             output_on_disk: !single, fp: None };
    let (mut fp, mut fp_secs) = ((0u64, 0u64), 0.0);
    let mut per_shard = est_out / nsh as f64;
    let mut parts = Vec::new();
    let mut mem_out = Vec::new();
    let mut start = 0;
    while start < nsh {
        let span = if single {
            nsh
        } else {
            let fit = (avail.max(0.0) / (table_bytes * per_shard)) as usize;
            if fit < nsh / MAX_PASSES && st.passes == 0 {
                eprintln!("  warning: m={} needs more than {MAX_PASSES} passes for the memory budget; \
                           running over budget", cx.m);
            }
            fit.max(nsh / MAX_PASSES).clamp(1, nsh - start)
        };
        let (shards, t) = run_pass::<KI, KO>(&input, cx, start..start + span, (per_shard * 1.05) as usize + 16);
        if st.passes == 0 {
            st.transitions = t;
        }
        st.states += shards.iter().map(|s| s.len() as u64).sum::<u64>();
        if ctx.o.fingerprint {
            let tf = Instant::now();
            let f = fingerprint_impl(&shards, cx.lout, cx.p);
            fp = (add_mod(fp.0, f.0, cx.p), add_mod(fp.1, f.1, cx.p));
            fp_secs += tf.elapsed().as_secs_f64();
        }
        if single {
            mem_out = shards;
        } else {
            let path = ctx.part_path(cx.m - 1, st.passes);
            parts.push(write_part(&path, &shards, &mut ctx.disk));
        }
        start += span;
        st.passes += 1;
        per_shard = st.states as f64 / start as f64 * 1.02;
    }
    let mut kept = None;
    if let Data::Disk(old) = input {
        if keep_input {
            kept = Some(Data::Disk(old));
        } else {
            delete_parts(old, &mut ctx.disk);
        }
    }
    st.fp = ctx.o.fingerprint.then_some(fp);
    (if single { Data::Mem(mem_out) } else { Data::Disk(parts) }, st, fp_secs, kept)
}

fn step_any(input: AnyData, in_states: u64, est_out: f64, cx: &StepCtx, ctx: &mut Ctx,
            keep: bool) -> (AnyData, LayerStat, f64, Option<AnyData>) {
    use AnyData::*;
    macro_rules! go {
        ($d:expr, $wi:ident, $wo:ident) => {{
            let (o, s, f, kept) = layer_step($d, in_states, est_out, cx, ctx, keep);
            ($wo(o), s, f, kept.map($wi))
        }};
    }
    match (input, cx.lout.wide()) {
        (Narrow(d), false) => go!(d, Narrow, Narrow),
        (Narrow(d), true) => go!(d, Narrow, Wide),
        (Wide(d), false) => go!(d, Wide, Narrow),
        (Wide(d), true) => go!(d, Wide, Wide),
    }
}

/// Count of the first state (the only one in the final layer), or 0 for an empty layer.
fn first_count<K: Key>(d: Data<K>, disk: &mut DiskUse) -> u64 {
    match d {
        Data::Mem(v) => v.into_iter().flatten().map(|e| e.1).next().unwrap_or(0),
        Data::Disk(parts) => {
            let mut c = None;
            read_parts::<K>(&parts, |v| c = c.or(v.first().map(|e| e.1)));
            delete_parts(parts, disk);
            c.unwrap_or(0)
        }
    }
}

/// Final layer of a forward run, plus the layer kept for the meet-in-the-middle join.
struct Forward {
    data: AnyData,
    lay: Layout,
    kept: Option<(AnyData, Layout, u64)>,
}

/// Forward DP placing m = N down to `last_m`. With `keep_rest`, the layer where 1..keep_rest
/// remain is kept on disk.
fn forward(ctx: &mut Ctx, res: &mut RunResult, last_m: u32, keep_rest: Option<u32>) -> Forward {
    let o = ctx.o;
    let (n, k) = (o.n, o.k);
    let big_n = n as u32 * k;
    let s_total = k * (big_n + 1) / 2;
    let mut lay = Layout::new(n, k, big_n);
    let key0 = lay.encode(&vec![(k, s_total); n]);
    let mut data = if lay.wide() {
        AnyData::Wide(Data::Mem(vec![vec![(key0, 1 % o.p)]]))
    } else {
        AnyData::Narrow(Data::Mem(vec![vec![(key0 as u64, 1 % o.p)]]))
    };
    let (mut prev, mut cur) = (1u64, 1u64);
    let mut kept = None;
    res.peak_states = 1;

    for m in (last_m..=big_n).rev() {
        let tl = Instant::now();
        let lout = Layout::new(n, k, m - 1);
        let rest = m - 1;
        let cx = StepCtx {
            lin: &lay,
            lout: &lout,
            m,
            lo_t: (0..=k).map(lo).collect(),
            hi_t: (0..=k).map(|s| if s <= rest { (s * (2 * rest + 1 - s) / 2) as i64 } else { -1 }).collect(),
            p: o.p,
            shard_bits: o.shard_bits,
        };
        // Growth is falling layer to layer, so the last ratio is a safe-side estimate.
        let est_out = cur as f64 * (cur as f64 / prev as f64).max(0.5) * 1.05 + 1024.0;
        let keep = keep_rest == Some(m);
        let (next, mut st, fp_secs, kept_in) = step_any(data, cur, est_out, &cx, ctx, keep);
        if let Some(d) = kept_in {
            kept = Some((d, lay.clone(), cur));
        }
        data = next;
        st.seconds = tl.elapsed().as_secs_f64() - fp_secs;
        (prev, cur) = (cur, st.states);
        lay = lout;
        res.peak_states = res.peak_states.max(st.states);
        if o.verbose {
            eprintln!("  m={m:3}  states={:>14}  transitions={:>16}  bits={:3}  passes={:3}{}{}  {:9.3}s  \
                       rss={:.0} MB  disk={:.1} GB",
                      st.states, st.transitions, st.key_bits, st.passes,
                      if st.input_on_disk { " in:disk" } else { "" }, if st.output_on_disk { " out:disk" } else { "" },
                      st.seconds, sys::peak_rss_bytes().unwrap_or(0) as f64 / 1048576.0, ctx.disk.now as f64 / 1e9);
        }
        let states = st.states;
        res.layers.push(st);
        if states == 0 {
            break;
        }
        if o.max_states.is_some_and(|cap| states > cap) {
            res.aborted = true;
            break;
        }
    }
    Forward { data, lay, kept }
}

fn new_ctx(o: &RunOpts) -> Ctx<'_> {
    // Unique per run, so concurrent runs in one process never share part files.
    static RUN_ID: AtomicU64 = AtomicU64::new(0);
    assert!(o.n <= MAX_BLOCKS && o.p < 1 << 63 && o.shard_bits <= 20);
    Ctx {
        o,
        dir: o.work_dir.join(format!("esp_{}x{}_{}_{}_{}", o.n, o.k, o.p, std::process::id(),
                                     RUN_ID.fetch_add(1, Ordering::Relaxed))),
        disk: DiskUse::default(),
    }
}

fn finish(ctx: Ctx, res: &mut RunResult, t0: Instant) {
    // Remaining part files (aborted runs) are removed with the directory.
    let _ = std::fs::remove_dir_all(&ctx.dir);
    res.peak_disk_bytes = ctx.disk.peak;
    res.total_seconds = t0.elapsed().as_secs_f64();
}

/// Full forward DP.
pub fn run(o: &RunOpts) -> RunResult {
    let mut res = RunResult::default();
    if !is_feasible(o.n, o.k) {
        return res;
    }
    let t0 = Instant::now();
    let mut ctx = new_ctx(o);
    let f = forward(&mut ctx, &mut res, 1, None);
    let last_ok = !res.aborted && res.layers.last().is_some_and(|l| l.m == 1 && l.states == 1);
    let c = match f.data {
        AnyData::Narrow(d) => first_count(d, &mut ctx.disk),
        AnyData::Wide(d) => first_count(d, &mut ctx.disk),
    };
    if last_ok {
        res.labeled_mod_p = c;
    }
    finish(ctx, &mut res, t0);
    res
}

#[derive(Clone, Debug, Default)]
pub struct MitmInfo {
    pub mid_states: u64,
    pub upper_states: u64,
    pub dead_mid_states: u64,
    pub join_passes: u32,
    pub join_seconds: f64,
}

/// Visit every chunk of a layer in parallel.
fn par_chunks_of<K: Key>(d: &Data<K>, f: impl Fn(&[(K, u64)]) + Sync) {
    match d {
        Data::Mem(shards) => shards.par_iter().for_each(|s| f(s)),
        Data::Disk(parts) => read_parts(parts, |v| v.par_chunks(1 << 15).for_each(&f)),
    }
}

/// Σ c_mid(P) · c_upper(P̄) · den(P̄) mod p, where den = (n-len)! ∏ mult! = n! / arrangements.
/// The upper layer is loaded one shard range at a time.
fn join_impl<KM: Key, KU: Key>(mid: &Data<KM>, lmid: &Layout, mid_states: u64, upper: &Data<KU>,
                               lup: &Layout, up_states: u64, o: &RunOpts, info: &mut MitmInfo) -> u64 {
    let (n, k, p) = (o.n, o.k, o.p);
    let big_n = n as u32 * k;
    let s_total = (k * (big_n + 1) / 2) as i64;
    let nsh = 1usize << o.shard_bits;
    let table_bytes = size_of::<(KU, u64)>() as f64 * 1.3;
    let avail = (o.mem_budget as f64 - mid.mem_bytes(mid_states) - upper.mem_bytes(up_states)).max(1e6);
    let passes = ((up_states as f64 * table_bytes / avail).ceil() as usize).clamp(1, nsh);
    let span = nsh.div_ceil(passes);
    let fact: Vec<u64> = (0..=n as u64).scan(1u64, |a, i| { *a *= i.max(1); Some(*a) }).collect();
    let total = Mutex::new(0u64);
    let found = AtomicU64::new(0);

    let mut start = 0;
    while start < nsh {
        let range = start..(start + span).min(nsh);
        let per_shard = (up_states as f64 / nsh as f64 * 1.1) as usize + 16;
        let tables: Vec<Mutex<FlatTable<KU>>> =
            range.clone().map(|_| Mutex::new(FlatTable::with_entries(per_shard))).collect();
        let init = || Buffers::new(&tables, range.start, o.shard_bits, p);
        let load = |buf: &mut Buffers<KU>, chunk: &[(KU, u64)]| {
            for &(key, c) in chunk {
                buf.push(key, c);
            }
        };
        match upper {
            Data::Mem(shards) => shards.par_iter().for_each_init(init, |b, sh| load(b, sh)),
            Data::Disk(parts) => read_parts(parts, |v| v.par_chunks(1 << 15).for_each_init(init, load)),
        }
        let tables: Vec<FlatTable<KU>> = tables.into_iter().map(|t| t.into_inner().unwrap()).collect();

        par_chunks_of(mid, |chunk| {
            let mut b = [(0, 0); MAX_BLOCKS];
            let mut comp = [(0, 0); MAX_BLOCKS];
            let (mut acc, mut hits) = (0u64, 0u64);
            'state: for &(key, c) in chunk {
                let len = lmid.decode(key.to_u128(), &mut b);
                let mut cl = 0;
                for &(s, r) in &b[..len] {
                    let cr = s_total - s as i64 * (big_n as i64 + 1) + r as i64;
                    let cs = k - s;
                    if cs == 0 && cr == 0 {
                        continue;
                    }
                    if cs == 0 || cr < 0 {
                        continue 'state;
                    }
                    comp[cl] = (cs, cr as u32);
                    cl += 1;
                }
                for _ in len..n {
                    comp[cl] = (k, s_total as u32);
                    cl += 1;
                }
                comp[..cl].sort_unstable();
                let Some(ck) = lup.try_encode(&comp[..cl]) else { continue };
                let ck = KU::from_u128(ck);
                let h = key_hash(ck);
                let sh = shard_of(h, o.shard_bits);
                if !range.contains(&sh) {
                    continue;
                }
                let Some(cu) = tables[sh - range.start].get(ck, h) else { continue };
                let mut den = fact[n - cl];
                let mut i = 0;
                while i < cl {
                    let mut j = i + 1;
                    while j < cl && comp[j] == comp[i] {
                        j += 1;
                    }
                    den *= fact[j - i];
                    i = j;
                }
                acc = add_mod(acc, mul_mod(mul_mod(c, cu, p), den % p, p), p);
                hits += 1;
            }
            let mut t = total.lock().unwrap();
            *t = add_mod(*t, acc, p);
            found.fetch_add(hits, Ordering::Relaxed);
        });
        start = range.end;
        info.join_passes += 1;
    }
    info.dead_mid_states = mid_states - found.into_inner();
    total.into_inner().unwrap()
}

fn join_any(mid: &AnyData, lmid: &Layout, mid_states: u64, upper: &AnyData, lup: &Layout, up_states: u64,
            o: &RunOpts, info: &mut MitmInfo) -> u64 {
    use AnyData::*;
    match (mid, upper) {
        (Narrow(a), Narrow(b)) => join_impl(a, lmid, mid_states, b, lup, up_states, o, info),
        (Narrow(a), Wide(b)) => join_impl(a, lmid, mid_states, b, lup, up_states, o, info),
        (Wide(a), Narrow(b)) => join_impl(a, lmid, mid_states, b, lup, up_states, o, info),
        (Wide(a), Wide(b)) => join_impl(a, lmid, mid_states, b, lup, up_states, o, info),
    }
}

fn drop_data(d: AnyData, disk: &mut DiskUse) {
    match d {
        AnyData::Narrow(Data::Disk(p)) | AnyData::Wide(Data::Disk(p)) => delete_parts(p, disk),
        _ => {}
    }
}

/// Meet in the middle with the complement x -> N+1-x: forward DP until 1..h remain
/// (h = N/2), then pair each mid state with its complement in the layer where 1..N-h remain.
pub fn run_mitm(o: &RunOpts) -> (RunResult, MitmInfo) {
    let mut res = RunResult::default();
    let mut info = MitmInfo::default();
    if !is_feasible(o.n, o.k) {
        return (res, info);
    }
    let big_n = o.n as u32 * o.k;
    let h = big_n / 2;
    if h == 0 {
        return (run(o), info);
    }
    let t0 = Instant::now();
    let mut ctx = new_ctx(o);
    let upper_rest = big_n - h;
    let f = forward(&mut ctx, &mut res, h + 1, (upper_rest != h).then_some(upper_rest));
    if !res.aborted {
        let mid_states = res.layers.last().map_or(1, |l| l.states);
        info.mid_states = mid_states;
        let tj = Instant::now();
        let t = match &f.kept {
            Some((up, lup, up_states)) => {
                info.upper_states = *up_states;
                join_any(&f.data, &f.lay, mid_states, up, lup, *up_states, o, &mut info)
            }
            None => {
                info.upper_states = mid_states;
                join_any(&f.data, &f.lay, mid_states, &f.data, &f.lay, mid_states, o, &mut info)
            }
        };
        info.join_seconds = tj.elapsed().as_secs_f64();
        res.labeled_mod_p = mul_mod(t, pow_mod((factorial(o.n) % o.p as u128) as u64, o.p - 2, o.p), o.p);
    }
    drop_data(f.data, &mut ctx.disk);
    if let Some((up, _, _)) = f.kept {
        drop_data(up, &mut ctx.disk);
    }
    finish(ctx, &mut res, t0);
    (res, info)
}

/// Chinese remainder for two residues; the result is below p0·p1 < 2^126.
pub fn crt2(a0: u64, p0: u64, a1: u64, p1: u64) -> u128 {
    let inv = pow_mod(p0 % p1, p1 - 2, p1);
    let diff = add_mod(a1, p1 - a0 % p1, p1);
    let t = mul_mod(diff, inv, p1);
    a0 as u128 + p0 as u128 * t as u128
}

pub fn factorial(n: usize) -> u128 {
    (1..=n as u128).product()
}
