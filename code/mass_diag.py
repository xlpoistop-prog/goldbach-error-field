"""Diagnostic: what is the size of the weighted-prime-mass fluctuation?

Compares (a) the real primes, (b) a random set with the same density, (c) the Poisson
prediction, for the weighted mass in [N/3, 2N/3].
"""
import os
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool

NMAX = 2_600_000
lam, pp, primes = tool.sieve(NMAX)
N0, M = 12000, 15000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
NODD = np.arange(3, 2_100_000, 2)


def wmass(Ns, indicator, denom=None):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD[(NODD > N / 3) & (NODD < 2 * N / 3)]
        out[i] = np.sum(indicator(n) * tool.weight(n / N))
    return out


def expected(Ns):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD[(NODD > N / 3) & (NODD < 2 * N / 3)]
        out[i] = np.sum(tool.weight(n / N) * 2.0 / np.log(n))   # 2/log n: primes among odds
    return out


ind_real = lambda n: (lam[n] > 0).astype(float)
m_real = wmass(Ns, ind_real)
E = expected(Ns)
print("real primes, weighted mass in [N/3,2N/3]:")
print("  mean = %.2f   expected = %.2f   ratio = %.4f" % (m_real.mean(), E.mean(), m_real.mean() / E.mean()))
print("  relative sd = %.4f" % (m_real.std() / m_real.mean()))
print("  Poisson prediction (1/sqrt(effective terms)):")
w = None
Ntest = 12000
n = NODD[(NODD > Ntest / 3) & (NODD < 2 * Ntest / 3)]
ww = tool.weight(n / Ntest)
p = 2.0 / np.log(n)
eff = (ww * p).sum() ** 2 / ((ww ** 2) * p).sum()
print("    effective terms at N=12000: %.0f  ->  1/sqrt = %.4f" % (eff, 1 / np.sqrt(eff)))

rng = np.random.default_rng(0)
for seed in range(3):
    rng = np.random.default_rng(seed)
    pfull = np.zeros(NMAX + 1)
    idx = np.arange(3, NMAX + 1, 2)
    pfull[idx] = np.minimum(2.0 / np.log(idx), 1.0)
    a = (rng.random(NMAX + 1) < pfull).astype(float)
    m = wmass(Ns, lambda x: a[x])
    print("  model seed %d: mean = %.2f  relative sd = %.4f" % (seed, m.mean(), m.std() / m.mean()))

print("\nautocorrelation of the real weighted mass (index lags, step 2):")
mz = (m_real - m_real.mean()) / m_real.std()
for lag in [1, 2, 5, 10, 50, 100, 500, 1000]:
    print("  lag %5d: %+.4f" % (lag, np.mean(mz[:-lag] * mz[lag:])))

print("\nsd of the *unweighted* prime count in [N/3,2N/3]:")
cnt = np.array([np.sum((lam[NODD[(NODD > N / 3) & (NODD < 2 * N / 3)]] > 0)) for N in Ns[:2000]])
print("  mean %.1f  sd %.1f  relative %.4f   (Poisson 1/sqrt(mean) = %.4f)"
      % (cnt.mean(), cnt.std(), cnt.std() / cnt.mean(), 1 / np.sqrt(cnt.mean())))
