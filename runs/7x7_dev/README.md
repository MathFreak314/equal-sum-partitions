# 7x7 runs from the development build (2026-10-20)

These are the runs that first established A321230(7) = 4302062977337507872775. All three ran on
one laptop: Intel Core Ultra 9 290HX, 24 threads, 32 GB RAM (about 16 GB free), with an 11.5 GB
memory budget.

| files | algorithm | primes | wall time |
|---|---|---|---|
| `forward_p0.*` | forward DP | p0 | 19.2 min |
| `forward_p12.*` | forward DP | p1, p2 | 16.9 + 18.4 min |
| `mitm_p0.*` | meet in the middle | p0 | 11.0 min |

Each run has these files:
- `*.json`: the run record ([FINGERPRINT_SPEC.md](../../FINGERPRINT_SPEC.md))
- `*.log`: per-layer progress on stderr
- `*.out`: the summary on stdout
- `*_layers.csv`: per-layer stats

**Provenance.** The commit hashes in these records come from the development history, before the
code was moved into this repository, so they don't appear in this repository's history.

- `forward_p12` and `mitm_p0` were built from a clean development commit.
- `forward_p0` was built from an earlier development commit plus one uncommitted change, hence
  `git_dirty: true`. That change added the call that keeps Windows from sleeping during a run.

Between those builds and this repository's first commit, the code changed in these ways:
- The forward driver was refactored to share code with meet in the middle.
- Each run got its own work-directory name.
- The default work directory became the system temp directory.
- Wording and crate metadata were updated.

None of these changes touch the DP transitions, pruning, key packing or modular arithmetic. Every
reference case was re-validated against the Python reference at this repository's commits
(`runs/python_*.json`, `runs/rust_*.json`).

To re-check:

```
python check_7x7.py runs/7x7_dev/forward_p0.json runs/7x7_dev/forward_p12.json runs/7x7_dev/mitm_p0.json \
    runs/7x7_dev/forward_p0.log runs/7x7_dev/forward_p12.log runs/7x7_dev/mitm_p0.log \
    runs/python_7x7_partial_layers.csv
```
