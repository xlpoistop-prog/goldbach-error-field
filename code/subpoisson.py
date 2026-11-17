"""Is the prime mass in a short interval sub-Poisson?  Compare with a random set of the
same density, and with the naive Poisson prediction."""
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool

NMAX = 2_600_000
lam, pp, primes = tool.sieve(NMAX)
N0, M = 12000, 15000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
NODD = np.arange(3, 2_100_000, 2)

# weights and index sets per N (independent of the prime set)
W = []
for N in Ns:
    n = NODD[(NODD > N / 3) & (NODD < 2 * N / 3)]
    W.append((n, tool.weight(n / N)))

p_odd = 2.0 / np.log(np.arange(3, 2_100_000))   # density of primes among odd numbers


def ratio_sd(stat):
    num = np.array([np.sum(stat(n) * w) for n, w in W])
    den = np.array([np.sum((2.0 / np.log(n)) * w) for n, w in W])
    r = num / den - 1.0
    return r.std(), r.mean()


print("relative fluctuation of the weighted prime mass in [N/3,2N/3]:")
sd_real, mean_real = ratio_sd(lambda n: (lam[n] > 0).astype(float))
print("  real primes            : sd = %.4f   mean ratio = %.3f" % (sd_real, mean_real))
w0 = None
Ntest = 12000
n0 = NODD[(NODD > Ntest / 3) & (NODD < 2 * Ntest / 3)]
w0 = tool.weight(n0 / Ntest)
eff = (w0 * 2 / np.log(n0)).sum() ** 2 / ((w0 ** 2) * 2 / np.log(n0)).sum()
print("  naive Poisson at N=12000: 1/sqrt(%.0f effective terms) = %.4f" % (eff, 1 / np.sqrt(eff)))
for seed in range(4):
    rng = np.random.default_rng(seed)
    pfull = np.zeros(NMAX + 1)
    idx = np.arange(3, NMAX + 1, 2)
    pfull[idx] = np.minimum(2.0 / np.log(idx), 1.0)
    a = (rng.random(NMAX + 1) < pfull).astype(float)
    sd_m, mean_m = ratio_sd(lambda n: a[n])
    print("  random set, seed %d    : sd = %.4f   mean ratio = %.3f" % (seed, sd_m, mean_m))

print("\nthe same for the *sieve* structure (coprime to 2,3,5,7,11,13,17,19,23 then thinned):")
mask = np.ones(NMAX + 1, bool)
mask[:2] = False
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    mask[p::p] = False
frac = float(np.prod(1.0 - 1.0 / np.array([2, 3, 5, 7, 11, 13, 17, 19, 23])))
dens = np.zeros(NMAX + 1)
idx = np.arange(2, NMAX + 1)
dens[2:] = np.where(mask[2:], 1.0 / (np.log(idx) * frac), 0.0)
for seed in range(4):
    rng = np.random.default_rng(100 + seed)
    a = (rng.random(NMAX + 1) < dens).astype(float)
    sd_m, mean_m = ratio_sd(lambda n: a[n])
    print("  sieve model, seed %d   : sd = %.4f   mean ratio = %.3f" % (seed, sd_m, mean_m))

print("\nunweighted count, detrended by its expectation:")
for tag, stat in [("real primes", lambda n: (lam[n] > 0).astype(float))]:
    num = np.array([np.sum(stat(n)) for n, w in W])
    den = np.array([np.sum(2.0 / np.log(n)) for n, w in W])
    print("  %s: sd of count/expectation = %.4f   (Poisson %.4f)"
          % (tag, (num / den - 1).std(), 1 / np.sqrt(den.mean())))
