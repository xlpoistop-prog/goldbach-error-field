"""Weight-robustness table, with the correct lower cut for each weight support."""
import os, sys, csv
import numpy as np
from multiprocessing import Pool
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool, ngb

STATE = {}
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
WEIGHTS = [(0, "cos4", 1 / 3, 2 / 3), (1, "flat", 1 / 3, 2 / 3),
           (2, "narrow", 5 / 12, 7 / 12), (3, "wide", 1 / 4, 3 / 4),
           (4, "beta", 1e-9, 1.0)]
KS = [2, 3, 4, 6, 10]


def init_worker(nmax):
    lam, pp, primes = tool.sieve(nmax)
    STATE["lam"] = lam; STATE["pp"] = pp
    STATE["spf"] = tool.spf_sieve(nmax)
    STATE["twoC2"] = tool.c2_constant()


def job(args):
    k, Ns, wid, lo, hi = args
    lam = STATE["lam"]; pp = STATE["pp"]
    r1 = ngb.rho(Ns, lam, pp, STATE["spf"], STATE["twoC2"], wid, lo, hi)
    rk = ngb.rho(k * Ns, lam, pp, STATE["spf"], STATE["twoC2"], wid, lo, hi)
    se = tool.block_se(r1, rk, L=40, reps=300, seed=1)
    return k, tool.corr(r1, rk), se


if __name__ == "__main__":
    N0, M = 12000, 20000
    Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
    nmax = int(10 * Ns[-1]) + 10
    rows = []
    for wid, name, lo, hi in WEIGHTS:
        with Pool(10, initializer=init_worker, initargs=(nmax,)) as pool:
            res = pool.map(job, [(k, Ns, wid, lo, hi) for k in KS])
        res.sort()
        rows.append((name, res))
        print(name, " ".join("%d:%+.4f(%.4f)" % (k, c, se) for k, c, se in res), flush=True)
    with open(os.path.join(OUT, "ck_weight.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["weight"] + ["c%d" % k for k in KS] + ["se%d" % k for k in KS])
        for name, res in rows:
            w.writerow([name] + ["%.4f" % x[1] for x in res] + ["%.4f" % x[2] for x in res])
    print("ck_weight.csv rewritten")
