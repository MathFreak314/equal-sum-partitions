use esp::*;

// OEIS A203986; A(n, k) = n blocks of size k.
const KNOWN: &[(usize, u32, u128)] = &[
    (2, 2, 1), (3, 3, 2), (4, 4, 392), (5, 5, 3245664), (2, 4, 4), (3, 4, 32), (5, 4, 6883),
    (3, 5, 305), (2, 6, 29), (3, 6, 4331), (4, 6, 2097719), (5, 6, 2549091482),
];

fn opts(n: usize, k: u32, p: u64) -> RunOpts {
    RunOpts { shard_bits: 4, ..RunOpts::new(n, k, p) }
}

#[test]
fn disk_and_passes_match_in_memory() {
    // A tiny budget forces multi-pass layers, input spills and on-disk outputs.
    let dir = std::env::temp_dir().join("esp_test_spill");
    let mut spilled = false;
    for (n, k) in [(5, 5), (5, 6), (4, 6)] {
        let mem = run(&RunOpts { fingerprint: true, ..opts(n, k, PRIMES[1]) });
        let disk = run(&RunOpts { fingerprint: true, shard_bits: 4, mem_budget: 256 << 10, work_dir: dir.clone(),
                                  ..opts(n, k, PRIMES[1]) });
        spilled |= disk.layers.iter().any(|l| l.passes > 1 && l.input_on_disk && l.output_on_disk);
        assert_eq!(mem.labeled_mod_p, disk.labeled_mod_p);
        for (a, b) in mem.layers.iter().zip(&disk.layers) {
            assert_eq!((a.m, a.states, a.transitions, a.fp), (b.m, b.states, b.transitions, b.fp));
        }
    }
    assert!(spilled, "budget too large to exercise multi-pass disk layers");
    assert!(std::fs::read_dir(&dir).map_or(true, |mut d| d.next().is_none()), "work dir not cleaned");
}

#[test]
fn known_values_all_primes() {
    for &(n, k, want) in KNOWN {
        let labeled = want * factorial(n);
        for p in PRIMES {
            assert_eq!(run(&opts(n, k, p)).labeled_mod_p as u128, labeled % p as u128, "A({n},{k}) mod {p}");
        }
    }
}

#[test]
fn mitm_known_values() {
    let dir = std::env::temp_dir().join("esp_test_mitm");
    for &(n, k, want) in KNOWN {
        let labeled = want * factorial(n);
        for p in [PRIMES[0], PRIMES[2]] {
            let (r, _) = run_mitm(&opts(n, k, p));
            assert_eq!(r.labeled_mod_p as u128, labeled % p as u128, "MITM A({n},{k}) mod {p}");
        }
        // Tiny budget: kept upper layer on disk, join in several passes.
        let (r, info) = run_mitm(&RunOpts { mem_budget: 256 << 10, work_dir: dir.clone(), ..opts(n, k, PRIMES[1]) });
        assert_eq!(r.labeled_mod_p as u128, labeled % PRIMES[1] as u128, "MITM/disk A({n},{k})");
        if (n, k) == (5, 6) {
            assert!(info.join_passes > 1);
        }
    }
    assert!(std::fs::read_dir(&dir).map_or(true, |mut d| d.next().is_none()), "work dir not cleaned");
}

#[test]
fn mitm_info_matches_python() {
    // From count_mitm() in equal_sum_partitions.py: (mid, dead, upper) states.
    for (n, k, mid, dead, upper) in [(5, 5, 10155, 1574, 10025), (5, 6, 97631, 14902, 97631)] {
        let (_, info) = run_mitm(&opts(n, k, PRIMES[0]));
        assert_eq!((info.mid_states, info.dead_mid_states, info.upper_states), (mid, dead, upper), "{n}x{k}");
    }
}

#[test]
fn parity_zeros_and_trivial() {
    for n in 2..7 {
        for k in 1..8u32 {
            if (k * (n as u32 - 1)) % 2 == 1 {
                assert_eq!(run(&opts(n, k, PRIMES[0])).labeled_mod_p, 0, "A({n},{k})");
            }
        }
    }
    for k in 1..6 {
        assert_eq!(run(&opts(1, k, PRIMES[0])).labeled_mod_p, 1);
    }
    assert_eq!(run(&opts(3, 1, PRIMES[0])).labeled_mod_p, 0);
}

#[test]
fn crt_reconstructs() {
    let x: u128 = 4704531536541600 * 1_000_000_007;
    let r = |p: u64| (x % p as u128) as u64;
    assert_eq!(crt2(r(PRIMES[0]), PRIMES[0], r(PRIMES[1]), PRIMES[1]), x);
}

#[test]
fn spec_vectors() {
    assert_eq!(mix64(0), 0);
    assert_eq!(mix64(1), 0x5692161D100B05E5);
    assert_eq!(mix64(0xFFFFFFFFFFFFFFFF), 0xB4D055FCF2CBBD7B);
    assert_eq!(mix64(0x9E3779B97F4A7C15), 0xE220A8397B1DCDAF);
    assert_eq!(state_hash(&[]), 0x9E3779B97F4A7C15);
    assert_eq!(state_hash(&[(7, 175); 7]), 0xDF5241AB2DEC08F1);
    assert_eq!(state_hash(&[(1, 3), (2, 10), (2, 10), (7, 175)]), 0xBFA99AF96964800C);
}

#[test]
fn layout_roundtrip() {
    let lay = Layout::new(7, 7, 25);
    let st = [(1, 3), (2, 10), (2, 10), (3, 30), (3, 40), (6, 132)];
    // Σs = 17 != 25 here, so append the block that makes the totals consistent.
    let mut full = st.to_vec();
    let (ss, rs): (u32, u32) = (full.iter().map(|b| b.0).sum(), full.iter().map(|b| b.1).sum());
    full.push((25 - ss, 25 * 26 / 2 - rs));
    let key = lay.encode(&full);
    let mut out = [(0, 0); MAX_BLOCKS];
    let len = lay.decode(key, &mut out);
    assert_eq!(&out[..len], &full[..]);
    assert!(!lay.wide(), "7x7 at m=26 should use u64 keys ({} bits)", lay.bits);
}

#[test]
fn fingerprint_matches_python_4x4() {
    // From test_equal_sum_partitions.py::test_fingerprint_layers.
    let r = run(&RunOpts { fingerprint: true, ..opts(4, 4, PRIMES[0]) });
    let got: Vec<_> = r.layers[..3].iter().map(|l| (l.m, l.fp.unwrap())).collect();
    assert_eq!(got, vec![(16, (4, 91071815559765149)), (15, (16, 1639661008511990903)),
                         (14, (60, 7279894382750824717))]);
}
