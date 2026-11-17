"""Large-N ladder for c(k): N0 = 1e7, 3e7, 1e8.

Same estimator as everywhere else: eps(N) = G(N)/(c_phi S(N)) - 1, c(k) = corr(eps(N), eps(kN)),
window N = N0 + 2i, i < M; block-bootstrap standard error (block 40).

One worker per k, each sieving its own range up to k*N_max.
"""
import os
import sys
import os
import csv
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
CPHI = 0.125


def sieve32(nmax):
    sv = np.ones(nmax + 1, bool)
    sv[:2] = False
    for i in range(2, int(nmax ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    lam = np.zeros(nmax + 1, dtype=np.float32)
    pr = np.nonzero(sv)[0]
    lam[pr] = np.log(pr.astype(np.float32))
    for p in pr[pr * pr <= nmax]:
        pk = int(p) * int(p)
        lp = np.float32(np.log(np.float32(p)))
        while pk <= nmax:
            lam[pk] = lp
            pk *= int(p)
    pp = np.nonzero(lam > 0)[0].astype(np.int64)
    return lam, pp, pr


def sing_one(N, pr):
    s = 2.0 * float(np.prod([1 - 1.0 / (p - 1.0) ** 2 for p in pr[pr > 2][:3000]]))
    n = N
    while n % 2 == 0:
        n //= 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            if d > 2:
                s *= (d - 1.0) / (d - 2.0)
            while n % d == 0:
                n //= d
        d += 2
    if n > 2:
        s *= (n - 1.0) / (n - 2.0)
    return s


def G(Ns, lam, pp):
    out = np.empty(Ns.size)
    for i, N in enumerate(Ns):
        half = N // 2
        j0 = np.searchsorted(pp, N / 3.0, side="right")
        j1 = np.searchsorted(pp, half, side="left")
        n = pp[j0:j1]
        m = N - n
        u = n.astype(np.float64) / N
        c = np.cos(3.0 * np.pi * (u - 0.5))
        w = c * c * c * c
        out[i] = 2.0 * np.sum(lam[n].astype(np.float64) * lam[m].astype(np.float64) * w) / N
    return out


def job(args):
    k, N0, M, nmax = args
    lam, pp, pr = sieve32(nmax)
    Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
    g1 = G(Ns, lam, pp)
    g2 = G(k * Ns, lam, pp)
    s1 = np.array([sing_one(N, pr) for N in Ns])
    s2 = np.array([sing_one(N, pr) for N in k * Ns])
    e1 = g1 / (CPHI * s1) - 1
    e2 = g2 / (CPHI * s2) - 1
    se = tool.block_se(e1, e2, L=40, reps=200, seed=1)
    return k, N0, M, float(tool.corr(e1, e2)), float(se), float(e1.std()), float(e2.std())


if __name__ == "__main__":
    plans = [(10_000_000, 3000, [2, 3, 4, 5]),
             (30_000_000, 2000, [2, 3, 4]),
             (100_000_000, 1200, [2, 3])]
    rows = []
    for N0, M, ks in plans:
        nmax = int(max(ks) * (N0 + 2 * M) * 1.02) + 10
        jobs = [(k, N0, M, nmax) for k in ks]
        print("N0=%d M=%d k=%s  sieve to %.2e  (%d workers)" % (N0, M, ks, nmax, len(jobs)), flush=True)
        with Pool(len(jobs)) as pool:
            res = pool.map(job, jobs)
        for k, n0, m, c, se, sd1, sd2 in sorted(res):
            print("   k=%d  c=%+.4f +- %.4f   sd %.4f/%.4f" % (k, c, se, sd1, sd2), flush=True)
            rows.append([n0, m, k, "%.4f" % c, "%.4f" % se, "%.4f" % sd1])
    with open(os.path.join(OUT, "ck_big.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N0", "M", "k", "c_k", "se", "sd_eps"])
        w.writerows(rows)
    print("written lab/out/ck_big.csv")
