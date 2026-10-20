"""Compare two run records (FINGERPRINT_SPEC.md) layer by layer.

Checks result, labeled residues, and per layer: states, transitions, and both fingerprints for
every prime the two records share. Exit code 0 only if everything shared matches.

    py compare_records.py runs/record_6x6.json runs/rust_6x6.json
"""

import json
import sys


def compare(a: dict, b: dict) -> list[str]:
    errs = []
    for key in ("n_blocks", "block_size"):
        if a["params"][key] != b["params"][key]:
            errs.append(f"params.{key}: {a['params'][key]} != {b['params'][key]}")
    if errs:
        return errs
    pa, pb = a["primes"], b["primes"]
    shared = [(pa.index(p), pb.index(p)) for p in pa if p in pb]
    if not shared:
        errs.append("no shared primes")
    for x in ("result", "labeled"):
        if a.get(x) is not None and b.get(x) is not None and a[x] != b[x]:
            errs.append(f"{x}: {a[x]} != {b[x]}")
    for ia, ib in shared:
        if a["labeled_mod_p"][ia] != b["labeled_mod_p"][ib]:
            errs.append(f"labeled mod {pa[ia]}: {a['labeled_mod_p'][ia]} != {b['labeled_mod_p'][ib]}")
    la, lb = a["layers"], b["layers"]
    if len(la) != len(lb):
        errs.append(f"layer count: {len(la)} != {len(lb)}")
    for x, y in zip(la, lb):
        m = x["m"]
        if m != y["m"]:
            errs.append(f"layer order differs at m={m} vs {y['m']}")
            break
        for key in ("states", "transitions"):
            if x[key] != y[key]:
                errs.append(f"m={m} {key}: {x[key]} != {y[key]}")
        for ia, ib in shared:
            for key in ("count_sum_mod_p", "hash_sum_mod_p"):
                if key in x and key in y and x[key][ia] != y[key][ib]:
                    errs.append(f"m={m} {key} (p={pa[ia]}): {x[key][ia]} != {y[key][ib]}")
    return errs


def main(argv) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    a, b = (json.load(open(p)) for p in argv)
    errs = compare(a, b)
    n_fp = sum(1 for x in a["layers"] if "hash_sum_mod_p" in x)
    shared = [p for p in a["primes"] if p in b["primes"]]
    if errs:
        print(f"MISMATCH ({len(errs)}):")
        for e in errs[:50]:
            print("  " + e)
        return 1
    print(f"OK: {a['params']['n_blocks']}x{a['params']['block_size']}, {len(a['layers'])} layers, "
          f"{n_fp} fingerprinted, {len(shared)} shared primes, result {a.get('result')} / {b.get('result')}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
