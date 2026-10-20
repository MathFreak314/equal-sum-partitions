# Per-layer fingerprint and run-record spec (schema 1)

Every implementation (Python reference, Rust port, …) must produce identical
per-layer fingerprints. If two independent implementations agree on every
fingerprint of every layer, they agree on the DP itself, not only on the final answer.

## Canonical state

After number `m` is placed (so `1..m-1` remain), a state is the list of all
**unfinished** blocks as pairs `(s, r)`:
- `s` is the number of slots left, with 1 ≤ s ≤ k.
- `r` is the sum still needed.

Rules:
- Blocks that reached `(0, 0)` are removed.
- Untouched blocks `(k, S)` are **included**.
- The list is sorted ascending, lexicographically by `(s, r)`, with duplicates kept.

The Rust port drops the largest block internally to pack keys. It must rebuild
the full canonical state before hashing.

## Block packing and state hash

All arithmetic is on unsigned 64-bit integers, wrapping mod 2⁶⁴.

```
mix64(z):                       # splitmix64 finalizer
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
    z = (z ^ (z >> 27)) * 0x94D049BB133111EB
    return z ^ (z >> 31)

pack(s, r) = (s << 16) | r      # requires s < 256, r < 65536

state_hash(blocks):             # blocks in canonical order
    h = 0x9E3779B97F4A7C15 ^ len(blocks)
    for (s, r) in blocks:
        h = mix64(h ^ pack(s, r))
    return h
```

Test vectors:

| input | value |
|---|---|
| `mix64(0)` | `0x0` |
| `mix64(1)` | `0x5692161d100b05e5` |
| `mix64(0xFFFFFFFFFFFFFFFF)` | `0xb4d055fcf2cbbd7b` |
| `mix64(0x9E3779B97F4A7C15)` | `0xe220a8397b1dcdaf` (first standard splitmix64 output from seed 0) |
| `state_hash([])` | `0x9e3779b97f4a7c15` |
| `state_hash([(7,175)] * 7)` | `0xdf5241ab2dec08f1` |
| `state_hash([(1,3),(2,10),(2,10),(7,175)])` | `0xbfa99af96964800c` |

## Layer fingerprint

For a layer `{state: count}`, where counts are labeled-assignment counts (the same
semantics as the reference DP), and for each prime p:

```
count_sum_mod_p = Σ count                    mod p
hash_sum_mod_p  = Σ state_hash(state)·count  mod p
```

`state_hash` is the full u64 value; reducing it mod p first gives the same result.
Both sums are independent of iteration order.

Primes, in this order:

```
p0 = 2^63 - 25  = 9223372036854775783
p1 = 2^63 - 165 = 9223372036854775643
p2 = 2^63 - 259 = 9223372036854775549
```

An implementation that keeps counts mod a single prime reports that prime's
column only.

Layers are reported for m = N, N−1, …, 1. The last layer is the state `[]`. Its
`count_sum` is the labeled total, which equals `result · n!`.

## Run record (JSON)

The record is a single JSON object. Any integer that can exceed 2⁵³ is stored as
a **decimal string**.

| key | contents |
|---|---|
| `schema` | `1` |
| `implementation` | e.g. `"python-reference"` or `"rust"` |
| `code` | `file`, `sha256` of the source, `git_commit`, `git_dirty` |
| `runtime` | interpreter or compiler version |
| `machine` | `os`, `cpu`, `logical_cpus`, `ram_bytes` |
| `params` | `n_blocks`, `block_size`, `N`, `S`, `max_states`, `algorithm` |
| `primes` | the primes used, as strings |
| `started_utc`, `finished_utc` | ISO-8601 timestamps |
| `aborted` | bool |
| `result`, `labeled` | string, or null if unknown (e.g. mod-p-only runs) |
| `labeled_mod_p` | list of strings, one per prime |
| `peak_states`, `total_seconds`, `peak_rss_bytes`, `disk_bytes` | measurements |
| `layers` | list of `{m, states, transitions, seconds, count_sum_mod_p, hash_sum_mod_p}` |

Reference records are in `runs/record_{n}x{k}.json`.
