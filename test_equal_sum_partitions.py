import itertools

import pytest

from equal_sum_partitions import count, count_with_stats

# OEIS A203986; A(n, k) = n blocks of size k.
KNOWN = {
    (2, 2): 1,
    (3, 3): 2,
    (4, 4): 392,
    (5, 5): 3245664,
    (2, 4): 4,
    (3, 4): 32,
    (5, 4): 6883,
    (3, 5): 305,
    (2, 6): 29,
    (3, 6): 4331,
    (4, 6): 2097719,
    (5, 6): 2549091482,
}


@pytest.mark.parametrize("nk,expected", sorted(KNOWN.items()))
def test_known_values(nk, expected):
    assert count(*nk) == expected


@pytest.mark.parametrize("n,k", [(n, k) for n in range(2, 7) for k in range(1, 8)
                                 if (k * (n - 1)) % 2 == 1])
def test_parity_zero(n, k):
    assert count(n, k) == 0


@pytest.mark.parametrize("k", range(1, 6))
def test_single_block(k):
    assert count(1, k) == 1


def test_block_size_one():
    # n blocks of one element each can't share a sum unless n == 1.
    assert count(3, 1) == 0


def brute(n, k):
    """Count unlabeled equal-sum partitions by canonical enumeration."""
    N = n * k
    S = k * (N + 1) // 2
    if k * (N + 1) % 2:
        return 0

    def rec(remaining):
        if not remaining:
            return 1
        first, rest = remaining[0], remaining[1:]
        total = 0
        for combo in itertools.combinations(rest, k - 1):
            if first + sum(combo) == S:
                left = [x for x in rest if x not in combo]
                total += rec(left)
        return total

    return rec(list(range(1, N + 1)))


@pytest.mark.parametrize("n,k", [(2, 2), (2, 4), (3, 3), (3, 4), (4, 3), (2, 6), (4, 4), (3, 5)])
def test_matches_brute_force(n, k):
    assert count(n, k) == brute(n, k)


def test_stats_recorded():
    st = count_with_stats(4, 4)
    assert st.result == 392
    assert st.peak_states == max(l.states for l in st.layers)
    assert st.layers[-1].m == 1 and st.layers[-1].states == 1


@pytest.mark.parametrize("nk,expected", sorted(KNOWN.items()))
def test_mitm_matches_known(nk, expected):
    from equal_sum_partitions import count_mitm
    assert count_mitm(*nk)[0] == expected


# Test vectors from FINGERPRINT_SPEC.md; a port must reproduce these exactly.
def test_mix64_vectors():
    from equal_sum_partitions import mix64
    assert mix64(0) == 0
    assert mix64(1) == 0x5692161D100B05E5
    assert mix64(0xFFFFFFFFFFFFFFFF) == 0xB4D055FCF2CBBD7B
    assert mix64(0x9E3779B97F4A7C15) == 0xE220A8397B1DCDAF   # standard splitmix64, seed 0


def test_state_hash_vectors():
    from equal_sum_partitions import state_hash
    assert state_hash(()) == 0x9E3779B97F4A7C15
    assert state_hash(((7, 175),) * 7) == 0xDF5241AB2DEC08F1
    assert state_hash(((1, 3), (2, 10), (2, 10), (7, 175))) == 0xBFA99AF96964800C


def test_fingerprint_layers():
    import json
    from equal_sum_partitions import PRIMES, run_record
    st = count_with_stats(4, 4, fingerprint=True)
    assert st.labeled == 392 * 24
    assert st.layers[-1].count_sum == st.labeled
    assert [(l.m, l.count_sum % PRIMES[0], l.hash_sum % PRIMES[0]) for l in st.layers[:3]] == [
        (16, 4, 91071815559765149), (15, 16, 1639661008511990903), (14, 60, 7279894382750824717)]
    rec = json.loads(json.dumps(run_record(st, "start", "end")))
    assert rec["result"] == "392" and rec["labeled_mod_p"][0] == str(9408)
    assert len(rec["layers"]) == 16 and len(rec["layers"][0]["hash_sum_mod_p"]) == 3


def test_fingerprint_does_not_change_counts():
    a = count_with_stats(5, 4)
    b = count_with_stats(5, 4, fingerprint=True)
    assert a.result == b.result == 6883
    assert [l.states for l in a.layers] == [l.states for l in b.layers]
