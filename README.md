# equal-sum-partitions

Counts the ways to partition {1, 2, …, nk} into n blocks of k numbers that all have the same sum,
k(nk+1)/2. Blocks are unlabeled.

- The full table is [OEIS A203986](https://oeis.org/A203986).
- The square case n = k is [OEIS A321230](https://oeis.org/A321230): split 1..n² into n blocks of n with equal sums.

## Result

**A321230(7) = 4302062977337507872775**

That is the number of ways to split the numbers 1..49 into 7 groups of 7 that each sum to 175.

| n | A321230(n) |
|---|---|
| 0 | 1 |
| 1 | 1 |
| 2 | 1 |
| 3 | 2 |
| 4 | 392 |
| 5 | 3245664 |
| 6 | 6534071578530 |
| **7** | **4302062977337507872775** |

Rows n ≤ 6 are the terms already in the OEIS, and this code reproduces all of them. Row 7 is new.

## Method

The count is a layered dynamic program. The numbers m = N, N−1, …, 1 (where N = nk) are placed one at a time.

- **States.** A state is the sorted multiset of (slots left, sum still needed) over the unfinished blocks.
- **Transitions.** Placing m into a block type of multiplicity μ adds μ times the state's count. So the DP counts *labeled* assignments, and the answer is that total divided by n!.
- **Pruning.** After m is placed, only 1..m−1 remain. A block with s slots and remaining sum r must then satisfy s(s+1)/2 ≤ r ≤ s(2m−1−s)/2.

For n = 7 the largest layer has 1,522,654,175 states (m = 25). Two implementations exist:

- **`equal_sum_partitions.py`: the Python reference.**
  - Exact integer counts.
  - Includes a meet-in-the-middle variant, `count_mitm`.
  - Practical up to 6x6.
- **`esp/`: the Rust implementation** that produced the n = 7 result.
  - **Counts mod primes.** Counts are kept mod three primes just below 2⁶³: 2⁶³−25, 2⁶³−165 and 2⁶³−259.
    - Any two of the primes multiply to more than 49!/(7!)⁷, which bounds the labeled count, so two primes reconstruct the exact value by the Chinese remainder theorem.
    - The third prime is a check.
  - **Packed keys.**
    - A state is packed into one 64-bit key.
    - The largest block is dropped and rebuilt from the fixed per-layer totals Σs = m−1 and Σr = m(m−1)/2.
    - Every other block becomes a dense code over its feasible (s, r) pairs.
  - **Out-of-core layers.**
    - A layer that doesn't fit the memory budget is built in several passes, each over a range of hash shards.
    - The input and output layers stream through checksummed files on disk.
    - The 7x7 run needs about 12 GB of RAM and about 50 GB of temporary disk.
  - **Meet in the middle (`--mitm`).**
    - Runs the DP only until 1..⌊N/2⌋ remain.
    - Then pairs each midpoint state with its complement under x ↦ N+1−x in the layer where 1..N−⌊N/2⌋ remain.
    - The pairing is an out-of-core hash join.

## Validation

1. **Known values.**
   - Both implementations reproduce every value of A203986 they were tested on, including A321230(6) = 6534071578530.
   - The Python reference is also checked against brute force on 8 small cases.
2. **Layer-by-layer agreement between two implementations.**
   - [`FINGERPRINT_SPEC.md`](FINGERPRINT_SPEC.md) defines, for every layer, two order-independent fingerprints mod each prime:
     - Σ count
     - Σ hash(state)·count
   - Rust and Python match on every layer of 3x3, 4x4, 5x4, 4x6, 5x5, 5x6 and 6x6, under all three primes.
   - Rust matches both with layers held in RAM and with layers forced through disk in up to 64 passes.
   - The records are in `runs/python_*.json` and `runs/rust_*.json`.
3. **Three primes.** The 7x7 residues for the three primes reconstruct one value: every pair gives the same CRT result.
4. **Divisibility.** The reconstructed labeled count, 21682397405781039678786000, is divisible by 7!. Nothing in the computation forces this.
5. **Second algorithm.** The meet-in-the-middle run gives the same residue as the forward DP.
   - It shares only the first half of the layers with the forward DP.
   - It replaces the second half with the complement pairing.
6. **Independent parity check** ([`parity_check.py`](parity_check.py)).
   - The complement x ↦ N+1−x is an involution on the set of partitions, so A ≡ F (mod 2), where F counts the partitions it fixes.
   - F reduces to a ranked subset convolution over the ⌊N/2⌋ complementary pairs, a calculation entirely separate from the DP.
   - F matches brute force on 11 small cases.
   - F has the right parity for all 13 known values.
   - For 7x7, F = 36518238043, which predicts A321230(7) is odd. It is.
   - That value of F is exact: `python parity_check.py --exact` recomputes it in big integers.
7. **Consistent layer statistics.**
   - Per-layer state and transition counts are identical across the three prime runs and the meet-in-the-middle run.
   - For m ≥ 35 they also equal the counts from an earlier Python-reference probe, run on a different machine (`runs/python_7x7_partial_layers.csv`).

The 7x7 run records and logs are in [`runs/7x7_dev/`](runs/7x7_dev/), together with a note on which build produced each one. This command re-runs checks 3–7 from them:

```
python check_7x7.py runs/7x7_dev/forward_p0.json runs/7x7_dev/forward_p12.json runs/7x7_dev/mitm_p0.json \
    runs/7x7_dev/forward_p0.log runs/7x7_dev/forward_p12.log runs/7x7_dev/mitm_p0.log \
    runs/python_7x7_partial_layers.csv
```

## Reproducing

Requirements:
- Python 3.12+ with `pytest` and `numpy`
- Rust 1.85+ (edition 2024)

```
python -m pytest -q                      # Python reference: 55 tests
cd esp && cargo test --release           # Rust: 9 tests, incl. disk spill and MITM
cargo build --release && cd ..

# A(n,k) with three primes, CRT, and a run record:
esp/target/release/esp 6 6 --primes 0,1,2 --record out.json
python compare_records.py runs/python_6x6.json out.json

# The n = 7 result. Pick a memory budget below free RAM and a work directory with ~50 GB free:
esp/target/release/esp 7 7 --primes 0,1,2 --mem-gb 11 --work-dir /path/to/scratch --record 7x7.json
esp/target/release/esp 7 7 --mitm --primes 0 --mem-gb 11 --work-dir /path/to/scratch --record 7x7_mitm.json
python check_7x7.py 7x7.json 7x7_mitm.json
```

The code was developed and run on Windows 11. The Rust code has portable fallbacks for other systems, but it has not been tested there.

On an Intel Core Ultra 9 290HX (24 threads) with an 11 GB budget:
- each prime of the forward DP takes about 17–20 minutes;
- the meet-in-the-middle run takes about 11 minutes.

Each run record states its code commit, compiler, machine and timings.

## Layout

| path | contents |
|---|---|
| `equal_sum_partitions.py` | Python reference DP, meet in the middle, run records (`--record`) |
| `test_equal_sum_partitions.py` | Python tests |
| `esp/` | Rust implementation (`src/lib.rs` DP and storage, `src/main.rs` CLI, `tests/known.rs`) |
| `FINGERPRINT_SPEC.md` | canonical state, hash, per-layer fingerprints, run-record schema |
| `compare_records.py` | compares two run records layer by layer |
| `parity_check.py` | the independent parity check |
| `check_7x7.py` | 7x7 cross-checks and final value |
| `runs/` | reference run records for 3x3..6x6 (`python_*`, `rust_*`), the 7x7 Python probe, and the 7x7 runs (`7x7_dev/`) |

## How AI was used

This project was built with substantial help from an AI coding assistant, Anthropic's **Claude**, used through **Claude Code**:
- The Python reference implementation, its tests and the first scaling study came from an earlier Claude Code session. The model version for that session was not recorded.
- Everything else came from a Claude Code session on 2026-10-20 using **Claude Opus 5.5**.

**What the AI did:**
- **Algorithms.** It proposed and wrote all the code in this repository: the layered DP, its pruning bounds, the meet-in-the-middle pairing, counting modulo primes with CRT, the key packing, and the out-of-core multi-pass storage.
- **Validation.** It designed the validation: the per-layer fingerprint spec, the comparison and check scripts, and the parity argument behind `parity_check.py` along with its subset-convolution evaluation.
- **Runs.** It ran the computations and analysed their logs.
- **Documentation.** It wrote this README, `FINGERPRINT_SPEC.md`, `runs/7x7_dev/README.md`, and the code comments and docstrings.

**What Daniel Branscombe did:**
- Chose the problem and the goal.
- Set the requirements: correctness before speed, strong independent validation, and the available hardware.
- Directed and approved each stage.
- Provided the machines and chose when to run the computations.

At the time of publication **he had not yet independently verified the mathematics or reviewed the code line by line**, and he is now studying both.

**What this means for the result.** A321230(7) = 4302062977337507872775 should be treated as a computational result produced by AI-written software. Its credibility rests on the checks under [Validation](#validation), not on human review of the code:
- It reproduces every known value.
- Two implementations agree on every intermediate layer.
- Three primes agree.
- The labeled count is divisible by 7!.
- A second algorithm gives the same residue.
- An independent parity argument agrees.

An **independent recomputation**, especially one written from scratch, would be the strongest further confirmation, and is very welcome.

No AI tool is an author of this work. Because of the [OEIS policy on AI use](https://oeis.org/wiki/Use_of_AI_for_OEIS_Submissions), this value has **not** been submitted to the OEIS. It will only be submitted once the human author can explain and defend the methods himself.

## License

MIT; see [LICENSE](LICENSE). Author: Daniel Branscombe.
