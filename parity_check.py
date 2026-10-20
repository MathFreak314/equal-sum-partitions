"""Independent parity check for A203986 / A321230 via the complement x -> N+1-x.

The complement map is an involution on the set of equal-sum partitions, so
A(n,k) ≡ F(n,k) (mod 2), where F counts the partitions it fixes.

In a fixed partition each block is one of:
  * type A: self-complementary. That is k/2 complementary pairs {i, N+1-i} (k even), or
    (k-1)/2 pairs plus the centre (N+1)/2 (k odd; then there is exactly one such block).
    Its sum is automatically S.
  * type B: swapped with its partner block B' = N+1-B. Together the pair of blocks covers k
    complementary pairs, with B taking one element of each; the weight of a k-set of pairs is
    (number of choices with sum S) / 2.

F = Σ over (a type-A blocks, t type-B block pairs) of (ordered covers of the P = ⌊N/2⌋ pairs)/(a! t!).
Ordered covers come from a ranked subset convolution:
    Σ_X (-1)^{P-|X|} ĝA(X)^a ĝB(X)^t,
where ĝ is the zeta transform over subsets X of the pairs. Everything is computed mod 2^64
with numpy uint64, which leaves F modulo 2^(64 - v2(a! t!)).

    python parity_check.py            # validates on small cases and known values, then reports 7x7
    python parity_check.py --exact    # also computes F(7,7) exactly (a few minutes)
"""

import itertools
import math
import sys

import numpy as np

from test_equal_sum_partitions import KNOWN

M64 = 1 << 64


def _zeta(a: np.ndarray, P: int) -> None:
    for i in range(P):
        v = a.reshape(-1, 2, 1 << i)
        v[:, 1, :] += v[:, 0, :]


def fixed_count(n: int, k: int, exact: bool = False) -> tuple[int, int]:
    """(F mod 2^bits, bits) for partitions of 1..nk fixed by x -> nk+1-x.

    With exact=True the final sums use Python integers and F is returned exactly (bits = 0).
    That mode is slower and needs the zeta values to fit in int64, which they do for n, k <= 7.
    """
    N = n * k
    if k * (N + 1) % 2:
        return 0, (0 if exact else 64)
    S = k * (N + 1) // 2
    P = N // 2
    size = 1 << P

    # Type-B weights on k-subsets of pairs 1..P.
    gB = np.zeros(size, dtype=np.uint64)
    if k <= P:
        combos = np.array(list(itertools.combinations(range(1, P + 1), k)), dtype=np.int64)
        signs = np.array(list(itertools.product((0, 1), repeat=k)), dtype=np.int64)
        sums = combos.sum(1)[:, None] + ((N + 1) - 2 * combos) @ signs.T
        f = (sums == S).sum(1)
        assert (f % 2 == 0).all()
        masks = (np.int64(1) << (combos - 1)).sum(1)
        gB[masks] = (f // 2).astype(np.uint64)
        del sums
    _zeta(gB, P)

    pc = np.zeros(size, dtype=np.int64)
    for i in range(P):
        pc[1 << i: 1 << (i + 1)] = pc[: 1 << i] + 1
    alpha = k // 2
    gA = np.array([math.comb(j, alpha) % M64 for j in range(P + 1)], dtype=np.uint64)[pc]
    neg = (P - pc) % 2 == 1

    if k % 2:
        cases = [(1, (n - 1) // 2)]          # the centre forces exactly one type-A block
    else:
        cases = [(n - 2 * t, t) for t in range(n // 2 + 1)]

    if exact:
        # Σ_X (-1)^{P-|X|} C(|X|,α)^a ĝB(X)^t, grouped by |X|, in exact integers.
        g = gB.astype(np.int64)
        total = 0
        for a, t in cases:
            if a * alpha + t * k != P:
                continue
            ordered = 0
            for j in range(P + 1):
                ca = math.comb(j, alpha) ** a
                if ca == 0:
                    continue
                powsum = sum(int(x) ** t for x in g[pc == j].tolist()) if t else math.comb(P, j)
                ordered += (-1) ** (P - j) * ca * powsum
            d = math.factorial(a) * math.factorial(t)
            assert ordered % d == 0
            total += ordered // d
        return total, 0

    total, bits = 0, 64
    with np.errstate(over="ignore"):
        for a, t in cases:
            if a * alpha + t * k != P:
                continue
            term = np.ones(size, dtype=np.uint64)
            for _ in range(a):
                term *= gA
            for _ in range(t):
                term *= gB
            ordered = (int(term[~neg].sum(dtype=np.uint64)) - int(term[neg].sum(dtype=np.uint64))) % M64
            d = math.factorial(a) * math.factorial(t)
            v = (d & -d).bit_length() - 1
            assert ordered % (1 << v) == 0
            b = 64 - v
            val = ((ordered >> v) * pow(d >> v, -1, 1 << b)) % (1 << b)
            bits = min(bits, b)
            total += val
    return total % (1 << bits), bits


def brute_fixed(n: int, k: int) -> int:
    """Count fixed partitions by enumerating all partitions (small cases only)."""
    N = n * k
    S = k * (N + 1) // 2
    if k * (N + 1) % 2:
        return 0
    found = 0

    def rec(remaining, blocks):
        nonlocal found
        if not remaining:
            parts = {frozenset(b) for b in blocks}
            if {frozenset(N + 1 - x for x in b) for b in parts} == parts:
                found += 1
            return
        first, rest = remaining[0], remaining[1:]
        for combo in itertools.combinations(rest, k - 1):
            if first + sum(combo) == S:
                left = [x for x in rest if x not in combo]
                rec(left, blocks + [(first,) + combo])

    rec(list(range(1, N + 1)), [])
    return found


def main() -> int:
    ok = True
    print("exact fixed counts vs brute force:")
    for n, k in [(2, 2), (2, 4), (3, 3), (3, 4), (4, 3), (2, 6), (4, 4), (3, 5), (2, 8), (3, 6), (5, 3)]:
        f, bits = fixed_count(n, k)
        b = brute_fixed(n, k)
        good = f == b % (1 << bits)
        ok &= good
        print(f"  {n}x{k}: F = {f}  brute = {b}  {'OK' if good else 'MISMATCH'}")

    print("parity of known values:")
    known = dict(KNOWN)
    known[(6, 6)] = 6534071578530
    for (n, k), a in sorted(known.items()):
        f, _ = fixed_count(n, k)
        good = f % 2 == a % 2
        ok &= good
        print(f"  A({n},{k}) = {a}  F = {f}  parity {'OK' if good else 'MISMATCH'}")

    print("exact mode vs mod-2^64 mode:")
    for n, k in [(3, 5), (4, 4), (5, 5), (4, 6), (6, 6)]:
        fe, _ = fixed_count(n, k, exact=True)
        f, bits = fixed_count(n, k)
        good = fe % (1 << bits) == f
        ok &= good
        print(f"  {n}x{k}: exact F = {fe}  {'OK' if good else 'MISMATCH'}")

    f, bits = fixed_count(7, 7)
    print(f"7x7: F mod 2^{bits} = {f}  =>  A(7,7) is {'odd' if f % 2 else 'even'}")
    if "--exact" in sys.argv:
        fe, _ = fixed_count(7, 7, exact=True)
        good = fe % (1 << bits) == f
        ok &= good
        print(f"7x7: exact F = {fe}  {'(agrees)' if good else 'MISMATCH'}")
    print("all checks passed" if ok else "CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
