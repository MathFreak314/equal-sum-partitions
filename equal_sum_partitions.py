"""Count partitions of {1..n*k} into n unlabeled blocks of k elements with equal sums.

Each block must sum to S = k(nk+1)/2. This is the OEIS A203986 table; its n = k
diagonal is A321230 (n=7: the numbers 1..49 in 7 blocks of 7, each summing to 175).

Algorithm: layered (forward) memoized DP.
  * Numbers are placed from N = n*k down to 1.
  * A state is the sorted tuple of (slots_remaining, sum_remaining) for each
    unfinished block; blocks reaching (0, 0) are dropped. Sorting merges blocks
    whose remaining needs are identical.
  * Only one layer {state: count} is held at a time.
  * Placing m into a distinct block type is weighted by that type's
    multiplicity, so the DP counts labeled assignments; the final answer is
    divided by n!.
  * After placing m (so 1..m-1 remain), a block with s slots and r remaining sum
    must satisfy s(s+1)/2 <= r <= s(2m-1-s)/2 (sum of the s largest of 1..m-1).
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import tracemalloc
from dataclasses import dataclass, field

# Validation fingerprints; byte-exact definition in FINGERPRINT_SPEC.md.
PRIMES = (2**63 - 25, 2**63 - 165, 2**63 - 259)
_M64 = (1 << 64) - 1


def mix64(z: int) -> int:
    """splitmix64 finalizer on a u64."""
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _M64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _M64
    return z ^ (z >> 31)


def state_hash(state: tuple) -> int:
    """u64 hash of a canonical (sorted, (0,0)-free) state."""
    h = 0x9E3779B97F4A7C15 ^ len(state)
    for s, r in state:
        assert 0 <= s < 256 and 0 <= r < 65536, (s, r)
        h = mix64(h ^ ((s << 16) | r))
    return h


def layer_fingerprint(layer: dict[tuple, int]) -> tuple[int, int]:
    """Exact (sum of counts, sum of state_hash*count); reduce mod p as needed."""
    total = 0
    weighted = 0
    for state, cnt in layer.items():
        total += cnt
        weighted += state_hash(state) * cnt
    return total, weighted


@dataclass
class LayerStat:
    m: int              # number just placed
    states: int         # live states after placing m
    transitions: int    # candidate transitions generated (before pruning)
    seconds: float      # wall time for this layer (excluding fingerprinting)
    count_sum: int | None = None     # exact sum of counts, if fingerprinted
    hash_sum: int | None = None      # exact sum of state_hash*count, if fingerprinted


@dataclass
class RunStats:
    n_blocks: int
    block_size: int
    result: int = 0
    labeled: int = 0                         # result * n_blocks!
    peak_states: int = 0
    total_seconds: float = 0.0
    peak_rss_bytes: int | None = None        # process peak working set / maxrss
    peak_traced_bytes: int | None = None     # tracemalloc peak, if enabled
    layers: list[LayerStat] = field(default_factory=list)
    aborted: bool = False                    # stopped early by max_states


def peak_rss_bytes() -> int | None:
    """Peak resident memory of this process, or None if unavailable."""
    if sys.platform == "win32":
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        k32 = ctypes.WinDLL("kernel32")
        psapi = ctypes.WinDLL("psapi")
        k32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        if psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb):
            return pmc.PeakWorkingSetSize
        return None
    try:
        import resource
    except ImportError:
        return None
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss if sys.platform == "darwin" else rss * 1024


def _is_feasible(n_blocks: int, block_size: int) -> bool:
    if n_blocks < 1 or block_size < 1:
        return False
    return (block_size * (n_blocks * block_size + 1)) % 2 == 0


def _step(layer: dict[tuple, int], m: int, block_size: int) -> tuple[dict[tuple, int], int]:
    """Place number m into every block type of every state; return (next layer, transitions)."""
    rest = m - 1
    # Feasibility bounds for a block with s slots once only 1..rest remain.
    lo = [s * (s + 1) // 2 for s in range(block_size + 1)]
    hi = [s * (2 * rest + 1 - s) // 2 if s <= rest else -1 for s in range(block_size + 1)]

    nxt: dict[tuple, int] = {}
    get = nxt.get
    transitions = 0
    for state, cnt in layer.items():
        L = len(state)
        i = 0
        while i < L:
            blk = state[i]
            j = i + 1
            while j < L and state[j] == blk:
                j += 1
            mult = j - i
            s, r = blk
            if r >= m:
                transitions += 1
                s1, r1 = s - 1, r - m
                if s1 == 0:
                    if r1 == 0:
                        new = state[:i] + state[i + 1:]
                    else:
                        new = None
                else:
                    new = tuple(sorted(state[:i] + ((s1, r1),) + state[i + 1:]))
                if new is not None:
                    ok = True
                    slots = 0
                    for bs, br in new:
                        if br < lo[bs] or br > hi[bs]:
                            ok = False
                            break
                        slots += bs
                    if ok:
                        assert slots == rest, (m, new)
                        nxt[new] = get(new, 0) + cnt * mult
            i = j
    return nxt, transitions


def count_with_stats(n_blocks: int, block_size: int, verbose: bool = False,
                     trace_memory: bool = False, max_states: int | None = None,
                     fingerprint: bool = False, out=sys.stdout) -> RunStats:
    stats = RunStats(n_blocks, block_size)
    if not _is_feasible(n_blocks, block_size):
        return stats

    N = n_blocks * block_size
    S = block_size * (N + 1) // 2

    if trace_memory:
        tracemalloc.start()
    t_start = time.perf_counter()

    layer: dict[tuple, int] = {((block_size, S),) * n_blocks: 1}
    stats.peak_states = 1

    for m in range(N, 0, -1):
        t_layer = time.perf_counter()
        nxt, transitions = _step(layer, m, block_size)
        layer = nxt
        stats.peak_states = max(stats.peak_states, len(layer))
        ls = LayerStat(m, len(layer), transitions, time.perf_counter() - t_layer)
        if fingerprint:
            ls.count_sum, ls.hash_sum = layer_fingerprint(layer)
        stats.layers.append(ls)
        if verbose:
            print(f"  m={m:3d}  states={ls.states:>12,}  transitions={ls.transitions:>14,}"
                  f"  {ls.seconds:9.3f}s", file=out, flush=True)
        if not layer:
            break
        if max_states is not None and len(layer) > max_states:
            stats.aborted = True
            break

    if not stats.aborted:
        labeled = layer.get((), 0)
        fact = math.factorial(n_blocks)
        assert labeled % fact == 0, (labeled, fact)
        stats.labeled = labeled
        stats.result = labeled // fact
    stats.total_seconds = time.perf_counter() - t_start
    if trace_memory:
        stats.peak_traced_bytes = tracemalloc.get_traced_memory()[1]
        tracemalloc.stop()
    stats.peak_rss_bytes = peak_rss_bytes()
    return stats


def count_mitm(n_blocks: int, block_size: int, verbose: bool = False,
               out=sys.stdout) -> tuple[int, dict]:
    """Meet-in-the-middle count using the complement map x -> N+1-x.

    Runs the forward DP only until 1..h remain (h = N // 2). A way to fill the
    bottom numbers 1..h with per-block needs (s_i, r_i) maps under the
    complement to a fill of the top h numbers with (s_i, s_i(N+1) - r_i), whose
    count is already in the forward layer where 1..N-h remain. Returns
    (count, info) where info reports how many midpoint states have no completion.
    """
    if not _is_feasible(n_blocks, block_size):
        return 0, {}
    n, k = n_blocks, block_size
    N = n * k
    S = k * (N + 1) // 2
    h = N // 2
    fact = math.factorial(n)

    layer: dict[tuple, int] = {((k, S),) * n: 1}
    upper = layer if N - h == N else None     # layer with 1..N-h remaining
    for m in range(N, h, -1):
        layer, _ = _step(layer, m, k)
        if verbose:
            print(f"  m={m:3d}  states={len(layer):>12,}", file=out, flush=True)
        if m - 1 == N - h:
            upper = layer
    mid = layer                               # 1..h remaining

    def arrangements(state: tuple) -> int:
        a = fact // math.factorial(n - len(state))
        i = 0
        while i < len(state):
            j = i
            while j < len(state) and state[j] == state[i]:
                j += 1
            a //= math.factorial(j - i)
            i = j
        return a

    labeled = 0
    dead = 0
    for state, cnt in mid.items():
        comp = [(k - s, S - s * (N + 1) + r) for s, r in state]
        comp += [(k, S)] * (n - len(state))
        comp = tuple(sorted(b for b in comp if b != (0, 0)))
        c = upper.get(comp, 0)
        if c == 0:
            dead += 1
            continue
        arr = arrangements(comp)
        assert c % arr == 0, (comp, c, arr)
        labeled += cnt * (c // arr)
    assert labeled % fact == 0, (labeled, fact)
    return labeled // fact, {"mid_states": len(mid), "dead_mid_states": dead,
                             "upper_states": len(upper)}


def count(n_blocks: int, block_size: int) -> int:
    """Number of partitions of {1..n_blocks*block_size} into n_blocks unlabeled
    blocks of block_size elements, all with the same sum."""
    return count_with_stats(n_blocks, block_size).result


def _mb(b: int | None) -> str:
    return "n/a" if b is None else f"{b / 2**20:,.1f} MB"


def _cpu_name() -> str:
    if sys.platform == "win32":
        try:
            import winreg
            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                                r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as k:
                return winreg.QueryValueEx(k, "ProcessorNameString")[0].strip()
        except OSError:
            pass
    return platform.processor() or "unknown"


def _total_ram_bytes() -> int | None:
    if sys.platform == "win32":
        import ctypes

        class MEMSTAT(ctypes.Structure):
            _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                        ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                        ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                        ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                        ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

        ms = MEMSTAT()
        ms.dwLength = ctypes.sizeof(MEMSTAT)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(ms)):
            return ms.ullTotalPhys
        return None
    try:
        return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None


def _code_version() -> dict:
    path = os.path.abspath(__file__)
    with open(path, "rb") as f:
        info = {"file": os.path.basename(path), "sha256": hashlib.sha256(f.read()).hexdigest(),
                "git_commit": None, "git_dirty": None}
    try:
        cwd = os.path.dirname(path)
        info["git_commit"] = subprocess.run(["git", "rev-parse", "HEAD"], cwd=cwd, capture_output=True,
                                            text=True, check=True).stdout.strip()
        info["git_dirty"] = bool(subprocess.run(["git", "status", "--porcelain", "--", info["file"]],
                                                cwd=cwd, capture_output=True, text=True,
                                                check=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        pass
    return info


def _utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def run_record(st: RunStats, started_utc: str, finished_utc: str,
               max_states: int | None = None) -> dict:
    """JSON-ready run record (schema in FINGERPRINT_SPEC.md). Big values are decimal strings."""
    n, k = st.n_blocks, st.block_size
    N = n * k
    layers = []
    for ls in st.layers:
        d = {"m": ls.m, "states": ls.states, "transitions": ls.transitions,
             "seconds": round(ls.seconds, 6)}
        if ls.count_sum is not None:
            d["count_sum_mod_p"] = [str(ls.count_sum % p) for p in PRIMES]
            d["hash_sum_mod_p"] = [str(ls.hash_sum % p) for p in PRIMES]
        layers.append(d)
    return {
        "schema": 1,
        "implementation": "python-reference",
        "code": _code_version(),
        "runtime": {"python": sys.version, "python_implementation": platform.python_implementation()},
        "machine": {"os": platform.platform(), "cpu": _cpu_name(), "logical_cpus": os.cpu_count(),
                    "ram_bytes": _total_ram_bytes()},
        "params": {"n_blocks": n, "block_size": k, "N": N, "S": k * (N + 1) // 2,
                   "max_states": max_states, "algorithm": "forward"},
        "primes": [str(p) for p in PRIMES],
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "aborted": st.aborted,
        "result": None if st.aborted else str(st.result),
        "labeled": None if st.aborted else str(st.labeled),
        "labeled_mod_p": None if st.aborted else [str(st.labeled % p) for p in PRIMES],
        "peak_states": st.peak_states,
        "total_seconds": round(st.total_seconds, 3),
        "peak_rss_bytes": st.peak_rss_bytes,
        "disk_bytes": None,
        "layers": layers,
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("n_blocks", type=int)
    p.add_argument("block_size", type=int, nargs="?", help="defaults to n_blocks")
    p.add_argument("-q", "--quiet", action="store_true", help="no per-layer output")
    p.add_argument("--tracemalloc", action="store_true",
                   help="also measure Python heap peak (slows the run noticeably)")
    p.add_argument("--max-states", type=int, metavar="K",
                   help="abort (no result) once a layer exceeds K states; for probing growth")
    p.add_argument("--csv", metavar="PATH", help="write per-layer stats to a CSV file")
    p.add_argument("--record", metavar="PATH",
                   help="write a JSON run record with per-layer fingerprints (slows the run)")
    a = p.parse_args(argv)
    k = a.block_size if a.block_size is not None else a.n_blocks

    print(f"A({a.n_blocks},{k}): numbers 1..{a.n_blocks * k}, block sum "
          f"{k * (a.n_blocks * k + 1) / 2:g}", flush=True)
    started = _utc_now()
    st = count_with_stats(a.n_blocks, k, verbose=not a.quiet, trace_memory=a.tracemalloc,
                          max_states=a.max_states, fingerprint=a.record is not None)
    finished = _utc_now()
    print(f"result        = {'ABORTED (max-states)' if st.aborted else st.result}")
    print(f"peak states   = {st.peak_states:,}")
    print(f"total time    = {st.total_seconds:.2f}s")
    print(f"peak RSS      = {_mb(st.peak_rss_bytes)}")
    if a.tracemalloc:
        print(f"peak traced   = {_mb(st.peak_traced_bytes)}")
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("m,states,transitions,seconds\n")
            for ls in st.layers:
                f.write(f"{ls.m},{ls.states},{ls.transitions},{ls.seconds:.6f}\n")
    if a.record:
        with open(a.record, "w") as f:
            json.dump(run_record(st, started, finished, a.max_states), f, indent=1)
            f.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
