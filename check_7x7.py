"""Cross-check the 7x7 run records: CRT result vs. the third prime, 7! divisibility, parity,
and identical per-layer statistics across runs.

    py check_7x7.py runs/7x7_forward_p0.json runs/7x7_forward_p12.json
"""

import json
import math
import sys

from equal_sum_partitions import PRIMES
from parity_check import fixed_count


def main(paths) -> int:
    recs = [json.load(open(p)) for p in paths]
    residues = {}
    for r in recs:
        for p, v in zip(r["primes"], r["labeled_mod_p"]):
            residues[int(p)] = int(v)
    assert set(residues) == set(PRIMES), residues
    p0, p1, p2 = PRIMES
    # CRT on (p1, p2), check against p0 (and also every other pair).
    def crt(a, pa, b, pb):
        t = ((b - a) * pow(pa, -1, pb)) % pb
        return a + pa * t
    labeled = crt(residues[p1], p1, residues[p2], p2)
    ok = True
    for pa, pb in [(p0, p1), (p0, p2)]:
        same = crt(residues[pa], pa, residues[pb], pb) == labeled
        ok &= same
        print(f"CRT({pa}, {pb}) agrees: {same}")
    third = labeled % p0 == residues[p0]
    ok &= third
    print(f"labeled mod p0 matches the p0 run: {third}")
    bound = math.factorial(49) // math.factorial(7) ** 7
    ok &= labeled < bound
    print(f"labeled < 49!/(7!)^7: {labeled < bound}")
    div = labeled % math.factorial(7) == 0
    ok &= div
    result = labeled // math.factorial(7)
    f, _ = fixed_count(7, 7)
    par = result % 2 == f % 2
    ok &= par
    print(f"divisible by 7!: {div}")
    print(f"parity matches complement-fixed count ({f}): {par}")

    base = recs[0]["layers"]
    for r in recs[1:]:
        same = len(r["layers"]) == len(base) and all(
            (x["m"], x["states"], x["transitions"]) == (y["m"], y["states"], y["transitions"])
            for x, y in zip(base, r["layers"]))
        ok &= same
        print(f"per-layer states/transitions identical to first record: {same}")
    print(f"labeled = {labeled}")
    print(f"A(7,7) = {result}")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
