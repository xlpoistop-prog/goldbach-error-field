"""Parallel driver on top of the numba implementation."""
import argparse, sys, time
import numpy as np
from multiprocessing import Pool
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool
import ngb

STATE = {}


def init_worker(nmax):
    lam, pp, primes = tool.sieve(nmax)
    STATE["lam"] = lam; STATE["pp"] = pp
    STATE["spf"] = tool.spf_sieve(nmax)
    STATE["twoC2"] = tool.c2_constant()


def job(args):
    k, Ns, wid = args
    lam = STATE["lam"]; pp = STATE["pp"]
    r1 = ngb.rho(Ns, lam, pp, STATE["spf"], STATE["twoC2"], wid)
    rk = ngb.rho(k * Ns, lam, pp, STATE["spf"], STATE["twoC2"], wid)
    c = tool.corr(r1, rk)
    se = tool.block_se(r1, rk, L=40, reps=300, seed=1)
    h = len(r1) // 2
    return k, c, se, tool.corr(r1[:h], rk[:h]), tool.corr(r1[h:], rk[h:])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n0", type=int, default=12000)
    ap.add_argument("--m", type=int, default=40000)
    ap.add_argument("--step", type=int, default=2)
    ap.add_argument("--k", type=str, default="2,3,4,5,6")
    ap.add_argument("--wid", type=int, default=0)
    ap.add_argument("--nmax", type=int, default=0)
    ap.add_argument("--procs", type=int, default=11)
    a = ap.parse_args()
    ks = [int(x) for x in a.k.split(",") if x.strip()]
    Ns = (a.n0 + a.step * np.arange(a.m)).astype(np.int64)
    nmax = a.nmax if a.nmax else int(max(ks) * Ns[-1]) + 10
    t0 = time.time()
    with Pool(a.procs, initializer=init_worker, initargs=(nmax,)) as pool:
        res = pool.map(job, [(k, Ns, a.wid) for k in ks])
    print("M=%d N0=%d step=%d wid=%d nmax=%d  (%.1fs)" % (a.m, a.n0, a.step, a.wid, nmax, time.time() - t0))
    print(" k   c(k)      se     half1    half2")
    for k, c, se, h1, h2 in res:
        print("%3d %+8.4f %7.4f %+8.4f %+8.4f" % (k, c, se, h1, h2))
