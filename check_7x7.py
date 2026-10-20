"""Cross-check 7x7 run records and print A(7,7).

Checks:
  * the forward residues for the three primes reconstruct one value by CRT (every pair agrees);
  * that value is below 49!/(7!)^7 and divisible by 7!;
  * A(7,7) has the parity predicted by parity_check.py;
  * every meet-in-the-middle residue equals the forward residue for the same prime;
  * per-layer states and transitions agree across all records on the layers they share.

    py check_7x7.py runs/7x7_forward.json runs/7x7_mitm_p0.json
"""

import json
import math
import sys

from equal_sum_partitions import PRIMES
from parity_check import fixed_count


def crt(a: int, pa: int, b: int, pb: int) -> int:
    return a + pa * (((b - a) * pow(pa, -1, pb)) % pb)


def main(paths) -> int:
    recs = [json.load(open(p)) for p in paths]
    forward, mitm = {}, {}
    for r in recs:
        dest = mitm if r["params"]["algorithm"] == "mitm" else forward
        for p, v in zip(r["primes"], r["labeled_mod_p"]):
            dest[int(p)] = int(v)
    ok = True
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

    layers = [{l["m"]: (l["states"], l["transitions"]) for l in r["layers"]} for r in recs]
    for i, lay in enumerate(layers[1:], 1):
        common = set(layers[0]) & set(lay)
        same = all(layers[0][m] == lay[m] for m in common)
        ok &= same
        print(f"{paths[i]}: {len(common)} shared layers identical to {paths[0]}: {same}")

    print(f"labeled = {labeled}")
    print(f"A(7,7) = {result}")
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
