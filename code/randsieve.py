"""The decisive null: a random sieve.

The independent-sampling null of fourvar2.py has the right local densities but the wrong
*variance* structure: the primes are strongly sub-Poisson for smooth weighted counts (the
sieve's negative short-range correlations), and an independent model is not.

This script builds the natural null with the same variance structure: a random sieve.

    keep only odd n;  for every odd prime p <= sqrt(X) pick a random residue r_p and delete
    every n == r_p (mod p).

Survivors: ~ (1/2) prod_{p odd <= sqrt X}(1 - 1/p) ~ e^{-gamma}/log sqrt X * X, i.e. the same
density as the primes, the same local structure, the same short-range (sieve) correlations,
and - by construction - *factorised* long-range four-point correlations.

E[G] is computed analytically: for a pair (n, N-n) the probability that both survive is
Q2 * prod_{p | N-2n} (p-1)/(p-2),  Q2 = prod_{p odd <= sqrt X} (1 - 2/p).
"""
import os
import sys
import os
import csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
X = 2_600_000
N0, M = 12000, 15000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)

# ---------- primes up to sqrt(X) and the singular-series factor
sv = np.ones(int(np.sqrt(X)) + 1, bool); sv[:2] = False
for i in range(2, int(np.sqrt(X)) + 1):
    if sv[i]:
        sv[i * i::i] = False
psmall = np.nonzero(sv)[0]
psmall = psmall[psmall > 2]            # p=2 is handled by "odds only"
print("sieving with %d odd primes, largest %d; predicted density %.5f (1/log X = %.5f)"
      % (psmall.size, psmall[-1], 0.5 * np.prod(1 - 1.0 / psmall), 1.0 / np.log(X)))
Q2 = float(np.prod(1.0 - 2.0 / psmall))

H = 200_000
Sfac = np.ones(H + 1)
for p in psmall:
    if p * p > H:
        pass
    idx = np.arange(p, H + 1, p)
    Sfac[idx] *= (p - 1.0) / (p - 2.0)

NODD = np.arange(3, X, 2)


def expected(Nsx):
    out = np.empty(Nsx.shape)
    NODDl = NODD
    for i, N in enumerate(Nsx):
        n = NODDl[(NODDl > N / 3) & (NODDl < N // 2)]
        h = N - 2 * n
        w = tool.weight(n / N)
        out[i] = 2.0 * Q2 * np.sum(np.log(n) * np.log(N - n) * w * Sfac[np.abs(h)]) / N
    return out


def rand_sieve(seed):
    rng = np.random.default_rng(seed)
    alive = np.zeros(X + 1, bool)
    alive[3::2] = True                      # odd numbers only
    for p in psmall:
        r = int(rng.integers(0, p))
        alive[r::p] = False
    alive[0:3] = False
    return alive


def Gmodel(Nsx, logw):
    """G(N) = (2/N) sum_{n<N/2, n in M} log n log(N-n) w(n/N)."""
    out = np.empty(Nsx.shape)
    idx = np.nonzero(logw)[0]
    idx = idx[idx >= 3]
    for i, N in enumerate(Nsx):
        j = np.searchsorted(idx, N // 2)
        n = idx[:j]
        if n.size == 0:
            out[i] = 0.0
            continue
        m = N - n
        out[i] = 2.0 * np.sum(logw[n] * logw[m] * tool.weight(n / N)) / N
    return out


def local_corr(a, b, width=40):
    B = len(a) // width
    v = np.empty(B)
    for i in range(B):
        x = a[i * width:(i + 1) * width]; y = b[i * width:(i + 1) * width]
        x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
        v[i] = np.mean(x * y)
    return float(v.mean()), float(v.std() / np.sqrt(B))


EG1 = expected(Ns)
EG2 = expected(2 * Ns)
print("E[G] at N=12000: %.5f   (the real G at 12000 is 0.44936)" % EG1[0])

rows = []
for seed in range(14):
    alive = rand_sieve(seed)
    lg = np.where(alive, np.log(np.maximum(np.arange(X + 1), 2)), 0.0)
    G1 = Gmodel(Ns, lg)
    G2 = Gmodel(2 * Ns, lg)
    e1 = G1 / EG1 - 1
    e2 = G2 / EG2 - 1
    full = tool.corr(e1, e2)
    loc, locse = local_corr(e1, e2)
    # sub-Poisson check: relative fluctuation of the model's (weighted) density in the window
    rows.append((seed, full, loc, locse, e1.std(), e2.std(), int(alive.sum())))
    print("  seed %2d: full %+.4f  within %+.4f+-%.4f  sd %.3f/%.3f  survivors %d"
          % (seed, full, loc, locse, e1.std(), e2.std(), int(alive.sum())), flush=True)

arr = np.array([[r[1] for r in rows], [r[2] for r in rows]])
print("\n  E[G] analytic vs empirical (first 5 N): ")
alive0 = rand_sieve(99)
lg0 = np.where(alive0, np.log(np.maximum(np.arange(X + 1), 2)), 0.0)
G0 = Gmodel(Ns[:5], lg0)
print("    analytic " + " ".join("%.4f" % v for v in EG1[:5]))
print("    empirical " + " ".join("%.4f" % v for v in G0))
print("\n  random sieve: full-window mean %+.4f (se %.4f), within-block mean %+.4f (se %.4f)"
      % (arr[0].mean(), arr[0].std(ddof=1) / np.sqrt(arr.shape[1]),
         arr[1].mean(), arr[1].std(ddof=1) / np.sqrt(arr.shape[1])))
print("  real primes:  full-window -0.1400, within-block -0.1348 +- 0.0087")
with open(os.path.join(OUT, "null_randsieve.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["seed", "c2_full", "c2_within", "se_within", "sd_eps1", "sd_eps2", "survivors"])
    for r in rows:
        w.writerow([r[0]] + ["%.4f" % x for x in r[1:6]] + [r[6]])
print("written lab/out/null_randsieve.csv")
