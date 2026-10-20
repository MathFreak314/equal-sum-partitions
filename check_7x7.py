"""Cross-check 7x7 run records and print A(7,7).

Checks:
  * the forward residues for the three primes reconstruct one value by CRT (every pair agrees);
  * that value is below 49!/(7!)^7 and divisible by 7!;
  * A(7,7) has the parity predicted by parity_check.py;
  * every meet-in-the-middle residue equals the forward residue for the same prime;
  * per-layer states and transitions agree on the layers they share, across all records, every
    prime run in the given logs (esp stderr) and any per-layer CSVs (e.g. the Python probe).

    python check_7x7.py runs/7x7_forward.json runs/7x7_mitm_p0.json runs/7x7_forward.log \
        runs/7x7_mitm_p0.log runs/python_7x7_partial_layers.csv
"""

import csv
import json
import math
import re
import sys

from equal_sum_partitions import PRIMES
from parity_check import fixed_count


def crt(a: int, pa: int, b: int, pb: int) -> int:
    return a + pa * (((b - a) * pow(pa, -1, pb)) % pb)


LOG_LINE = re.compile(r"m=\s*(\d+)\s+states=\s*(\d+)\s+transitions=\s*(\d+)")


def log_runs(path: str) -> list[dict]:
    """Per-layer {m: (states, transitions)} for each run in an esp log; a new run starts at m=49."""
    runs = []
    for line in open(path, encoding="utf-8", errors="replace"):
        g = LOG_LINE.search(line)
        if g:
            m, s, t = map(int, g.groups())
            if m == 49 or not runs:
                runs.append({})
            runs[-1][m] = (s, t)
    return runs


def main(paths) -> int:
    rec_paths = [p for p in paths if p.endswith(".json")]
    recs = [json.load(open(p)) for p in rec_paths]
    tables = [(f"{p} (record)", {l["m"]: (l["states"], l["transitions"]) for l in r["layers"]})
              for p, r in zip(rec_paths, recs)]
    for p in paths:
        if p.endswith(".log"):
            tables += [(f"{p} run {i + 1}", t) for i, t in enumerate(log_runs(p))]
        elif p.endswith(".csv"):
            rows = csv.DictReader(open(p))
            tables.append((p, {int(r["m"]): (int(r["states"]), int(r["transitions"])) for r in rows}))
    forward, mitm = {}, {}
    ok = True
    for path, r in zip(rec_paths, recs):
        dest = mitm if r["params"]["algorithm"] == "mitm" else forward
        for p, v in zip(r["primes"], r["labeled_mod_p"]):
            p, v = int(p), int(v)
            if p in dest and dest[p] != v:
                ok = False
                print(f"{path}: residue mod {p} = {v} disagrees with another record ({dest[p]})")
            dest[p] = v
    if set(forward) != set(PRIMES):
        print(f"forward residues missing for primes {sorted(set(PRIMES) - set(forward))}")
        return 1

    p0, p1, p2 = PRIMES
    labeled = crt(forward[p1], p1, forward[p2], p2)
    for pa, pb in [(p0, p1), (p0, p2)]:
        same = crt(forward[pa], pa, forward[pb], pb) == labeled
        ok &= same
        print(f"CRT(p={pa}, p={pb}) gives the same value: {same}")
    bound = math.factorial(49) // math.factorial(7) ** 7
    ok &= labeled < bound
    print(f"labeled < 49!/(7!)^7: {labeled < bound}")
    div = labeled % math.factorial(7) == 0
    ok &= div
    print(f"divisible by 7!: {div}")
    result = labeled // math.factorial(7)
    f, _ = fixed_count(7, 7)
    par = result % 2 == f % 2
    ok &= par
    print(f"parity matches the complement-fixed count F = {f}: {par}")
    for p, v in sorted(mitm.items()):
        same = forward[p] == v
        ok &= same
        print(f"meet-in-the-middle residue mod {p} equals forward: {same}")

    base_name, base = tables[0]
    for name, lay in tables[1:]:
        common = set(base) & set(lay)
        same = bool(common) and all(base[m] == lay[m] for m in common)
        ok &= same
        print(f"{name}: {len(common)} shared layers identical to {base_name}: {same}")

    print(f"labeled = {labeled}")
    print(f"A(7,7) = {result}")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
