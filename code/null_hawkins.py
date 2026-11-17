"""The missing null: a random-sieve model with the primes' variance structure.

Why the earlier nulls failed:
  * independent sampling (fourvar2/null_par) has the right local densities but is *mass*
    dominated, i.e. smooth in N, while the real field is white;
  * the first random-sieve attempt (randsieve.py) had the wrong density (flat instead of
    1/log n), hence the wrong number of pairs.

This model:
  1. keep the odd integers;  for every odd prime p <= P delete one random residue class
     (this reproduces the sieve's negative short-range correlations, hence sub-Poisson mass);
  2. thin with probability t(n) = c_t / log n, c_t chosen so that the density is exactly 1/log n;
  3. E[G] is analytic and includes the singular-series factor of the shift:
        P(n, m both alive) = Q2 * t(n) t(m) * prod_{p | m-n, p<=P} (p-1)/(p-2),
        Q2 = prod_{p<=P, p>2} (1 - 2/p).
Four-point correlations factorise by construction (all choices are independent).
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
X = 3_000_000
P = 1612
N0, M = 12000, 15000
NS = (N0 + 2 * np.arange(M)).astype(np.int64)
CPHI = 0.125
STATE = {}


def init_worker():
    # the model starts from ALL odd integers, not from the primes
    sv_odds = np.zeros(X + 1, bool)
    sv_odds[3::2] = True
    sv = np.ones(X + 1, bool)
    sv[:2] = False
    for i in range(2, int(X ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    primes = np.nonzero(sv)[0]
    sp = primes[(primes > 2) & (primes <= P)]
    dens = 0.5 * float(np.prod(1.0 - 1.0 / sp))
    ct = 1.0 / dens
    Q2 = float(np.prod(1.0 - 2.0 / sp))
    # singular-series factor prod_{p | h, p <= P} (p-1)/(p-2) for h <= X/3
    H = int(2 * NS[-1] // 3) + 10
    Sfac = np.ones(H + 1)
    for p in sp:
        if p > H:
            break
        Sfac[p::p] *= (p - 1.0) / (p - 2.0)
    logn = np.log(np.maximum(np.arange(X + 1), 2))
    t = np.minimum(ct / logn, 1.0)
    t[:3] = 0.0
    t[2::2] = 0.0                      # odds only
    lam, pp, primes_all = tool.sieve(X)
    STATE.update(dict(sv=sv, sv_odds=sv_odds, sp=sp, dens=dens, ct=ct, Q2=Q2, Sfac=Sfac,
                      logn=logn, t=t, lam=lam, pp=pp))
    STATE["EG1"] = EG(NS)
    STATE["EG2"] = EG(2 * NS)


def EG(Ns):
    """E[G](N) = (2/N) Q2 sum_{n odd in window} log n log m t(n)t(m) Sfac(|N-2n|) w(n/N)."""
    Q2 = STATE["Q2"]; Sfac = STATE["Sfac"]; t = STATE["t"]
    out = np.empty(Ns.size)
    for i, N in enumerate(Ns):
        n = np.arange(N // 3 + 1, N // 2 + 1)
        n = n[n % 2 == 1]
        m = N - n
        h = np.abs(N - 2 * n)
        w = tool.weight(n / N)
        out[i] = 2.0 * Q2 * np.sum(np.log(n) * np.log(m) * t[n] * t[m] * Sfac[h] * w) / N
    return out


def G_model(Ns, lam_a, pp):
    out = np.empty(Ns.size)
    for i, N in enumerate(Ns):
        half = N // 2
        j0 = np.searchsorted(pp, N / 3.0, side="right")
        j1 = np.searchsorted(pp, half, side="left")
        n = pp[j0:j1]
        if n.size == 0:
            out[i] = 0.0
            continue
        m = N - n
        w = tool.weight(n / N)
        out[i] = 2.0 * np.sum(lam_a[n] * lam_a[m] * w) / N
    return out


def blockcorr(a, b, width=40):
    B = len(a) // width
    v = np.empty(B)
    for i in range(B):
        x = a[i * width:(i + 1) * width]; y = b[i * width:(i + 1) * width]
        x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
        v[i] = np.mean(x * y)
    return float(v.mean()), float(v.std() / np.sqrt(B))


def job(seed):
    rng = np.random.default_rng(seed)
    sv = STATE["sv_odds"].copy()
    for p in STATE["sp"]:
        r = int(rng.integers(0, p))
        sv[r::p] = False
    sv[:3] = False
    sv[2::2] = False
    keep = rng.random(X + 1) < STATE["t"]
    alive = sv & keep
    lam_a = np.where(alive, STATE["logn"], 0.0)
    idx = np.flatnonzero(alive)          # the model's own support, NOT the real prime powers
    G1 = G_model(NS, lam_a, idx)
    G2 = G_model(2 * NS, lam_a, idx)
    e1 = G1 / STATE["EG1"] - 1.0
    e2 = G2 / STATE["EG2"] - 1.0
    loc, locse = blockcorr(e1, e2)
    return (seed, float(tool.corr(e1, e2)), loc, locse, float(e1.std()),
            int(alive.sum()), float(G1.mean()))


if __name__ == "__main__":
    init_worker()
    print("density of the sieved set: %.5f (1/log 3e6 = %.5f),  c_t = %.3f"
          % (STATE["dens"], 1.0 / np.log(X), STATE["ct"]), flush=True)
    lam, pp, _ = tool.sieve(X)
    print("E[G](12000) = %.5f   (the real G(12000) = 0.44936)" % STATE["EG1"][0], flush=True)
    print("E[G](12000) from the real primes: %.5f"
          % (2.0 * np.sum(lam[pp[(pp > 4000) & (pp < 6000)]] *
                          lam[12000 - pp[(pp > 4000) & (pp < 6000)]] *
                          tool.weight(pp[(pp > 4000) & (pp < 6000)] / 12000.0)) / 12000), flush=True)
    with Pool(10, initializer=init_worker) as pool:
        res = pool.map(job, range(40))
    a = np.array([[r[1] for r in res], [r[2] for r in res]])
    sd = np.array([r[4] for r in res])
    print("\nrandom-sieve null (40 realisations, 10 cores):")
    print("  sd(eps) = %.4f +- %.4f   (real: 0.0862)" % (sd.mean(), sd.std() / np.sqrt(len(sd))))
    print("  full-window   c(2) = %+.4f  (spread %.4f, se of mean %.4f)"
          % (a[0].mean(), a[0].std(ddof=1), a[0].std(ddof=1) / np.sqrt(len(res))))
    print("  within-block  c(2) = %+.4f  (spread %.4f, se of mean %.4f)"
          % (a[1].mean(), a[1].std(ddof=1), a[1].std(ddof=1) / np.sqrt(len(res))))
    print("  real primes:  full -0.1400,  within-block -0.1348 +- 0.0087")
    with open(os.path.join(OUT, "null_hawkins.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["seed", "c2_full", "c2_within", "se_within", "sd_eps", "survivors", "Gmean"])
        for r in res:
            w.writerow([r[0]] + ["%.4f" % x for x in r[1:5]] + [r[5], "%.4f" % r[6]])
    print("written lab/out/null_hawkins.csv")
