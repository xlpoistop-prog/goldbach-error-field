"""gcfar2 -- same physics, better scheduling.

Differences from gcfar.py:
  * ONE worker pool for the whole run (no re-fork per k);
  * fine-grained tasks (8 chunks per worker) dispatched with imap_unordered, so the pool
    stays fed right up to the end of each k (work stealing), instead of 100 equal chunks
    whose slowest member idles the other 99;
  * progress and utilisation printed per k.
"""
import argparse
import os
import sys
import os
import time
import numpy as np
from numba import njit

CPHI = 0.125
TWOC2 = 1.3203236316937392


def sieve_lam(X, chunk=1 << 27, verbose=True):
    t0 = time.time()
    base = int(X ** 0.5) + 1
    sv = np.ones(base + 1, bool)
    sv[:2] = False
    for i in range(2, int(base ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    bp = np.nonzero(sv)[0].astype(np.int64)
    lam = np.zeros(X + 1, dtype=np.float32)
    parts = []
    for c in range(0, X + 1, chunk):
        hi = min(c + chunk, X + 1)
        seg = np.ones(hi - c, bool)
        if c == 0:
            seg[0] = False
            if hi > 1:
                seg[1] = False
        for p in bp:
            p = int(p)
            start = max(p * p, ((c + p - 1) // p) * p)
            if start >= hi:
                continue
            seg[start - c::p] = False
        idx = np.nonzero(seg)[0] + c
        parts.append(idx)
        lam[idx] = np.log(idx.astype(np.float32))
        if verbose and (c // chunk) % 25 == 0:
            print("  sieve %5.1f%%  %.0fs" % (100.0 * hi / X, time.time() - t0), flush=True)
    primes = np.concatenate(parts).astype(np.int64)
    for p in bp:
        pp = int(p) * int(p)
        lp = np.float32(np.log(np.float32(p)))
        while pp <= X:
            lam[pp] = lp
            pp *= int(p)
    if verbose:
        print("  sieve done: %d primes, %.1fs" % (primes.size, time.time() - t0), flush=True)
    return lam, primes


@njit(cache=True, fastmath=True, nogil=True)
def G_at(Ns, lam, primes, k, out):
    for i in range(Ns.size):
        N = Ns[i] * k
        lo = N // 3 + 1
        hi = N // 2
        j0 = np.searchsorted(primes, lo, side="right")
        j1 = np.searchsorted(primes, hi, side="left")
        s = 0.0
        for j in range(j0, j1):
            n = primes[j]
            m = N - n
            lm = lam[m]
            if lm != 0.0:
                u = n / N
                c = np.cos(9.42477796076938 * (u - 0.5))
                s += float(lam[n]) * float(lm) * (c * c * c * c)
        out[i] = 2.0 * s / N


_LAM = _PRIMES = _NS = None


def _task(t):
    k, a, b = t
    o = np.empty(b - a)
    G_at(_NS[a:b], _LAM, _PRIMES, k, o)
    return k, a, o


def sing(N):
    s = TWOC2
    n = N
    while n % 2 == 0:
        n //= 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            s *= (d - 1.0) / (d - 2.0)
            while n % d == 0:
                n //= d
        d += 2
    if n > 2:
        s *= (n - 1.0) / (n - 2.0)
    return s


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def block_se(a, b, width=40, reps=300, seed=1):
    B = a.size // width
    if B < 3:
        return float("nan")
    x = a[:B * width].reshape(B, width)
    y = b[:B * width].reshape(B, width)
    x = x - x.mean(1, keepdims=True)
    y = y - y.mean(1, keepdims=True)
    v = (x * y).sum(1) / np.sqrt((x * x).sum(1) * (y * y).sum(1))
    rng = np.random.default_rng(seed)
    m = np.array([v[rng.integers(0, B, B)].mean() for _ in range(reps)])
    return float(m.std(ddof=1))


def run(n0, m, ks, workers, tag, per=8):
    global _LAM, _PRIMES, _NS
    import multiprocessing as mp
    t0 = time.time()
    Ns = (n0 + 2 * np.arange(m)).astype(np.int64)
    allk = [1] + [k for k in ks if k != 1]
    X = int(max(allk) * Ns[-1] * 2 // 3) + 10
    print("window N0=%d M=%d ks=%s  sieve to %.3e (%.0f GB)" % (n0, m, ks, X, 4.0 * X / 1e9), flush=True)
    _LAM, _PRIMES = sieve_lam(X)
    _NS = Ns

    nc = max(workers * per, 4)
    edges = np.linspace(0, m, nc + 1).astype(int)
    chunks = [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    tasks = [(k, a, b) for k in allk for a, b in chunks]
    print("  %d tasks (%d k x %d chunks), one pool of %d" % (len(tasks), len(allk), len(chunks), workers),
          flush=True)

    out = {k: np.empty(m) for k in allk}
    ctx = mp.get_context("fork")
    done = 0
    t1 = time.time()
    with ctx.Pool(workers) as pool:
        for k, a, o in pool.imap_unordered(_task, tasks, chunksize=1):
            out[k][a:a + o.size] = o
            done += 1
            if done % max(1, len(tasks) // 20) == 0:
                el = time.time() - t1
                print("    %3d/%d tasks  %.0fs  (%.1f%% of the run)" % (done, len(tasks), el,
                      (time.time() - t0) / max(1e-9, el) * 100), flush=True)
    print("  compute done in %.0fs" % (time.time() - t1), flush=True)

    def eps(k):
        S = np.array([sing(int(v)) for v in k * Ns])
        return out[k] / (CPHI * S) - 1.0

    e1 = eps(1)
    print("  k=1   sd(eps)=%.4f" % e1.std(), flush=True)
    rows = []
    for k in allk:
        if k == 1:
            continue
        ek = eps(k)
        r, se = corr(e1, ek), block_se(e1, ek)
        rows.append((k, r, se, e1.std(), ek.std()))
        print("  k=%-3d c=%+.4f +- %.4f   sd %.4f/%.4f" % (k, r, se, e1.std(), ek.std()), flush=True)
    with open("ck_%s.csv" % tag, "w") as f:
        f.write("N0,M,k,c,se,sd1,sd2\n")
        for k, r, se, s1, s2 in rows:
            f.write("%d,%d,%d,%.4f,%.4f,%.4f,%.4f\n" % (n0, m, k, r, se, s1, s2))
    print("written ck_%s.csv   total %.0fs" % (tag, time.time() - t0), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n0", type=float, default=1e8)
    ap.add_argument("--m", type=int, default=10000)
    ap.add_argument("--k", type=str, default="2,3,4,5,6,7,9,10,12,15")
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--per", type=int, default=8)
    ap.add_argument("--tag", type=str, default="run")
    a = ap.parse_args()
    run(int(a.n0), a.m, [int(x) for x in a.k.split(",")], a.workers or min(120, os.cpu_count() or 8),
        a.tag, a.per)
