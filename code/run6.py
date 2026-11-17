import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tool import *

NMAX = 1_200_000
lam, pp, primes = sieve(NMAX)
spf = spf_sieve(NMAX)
twoC2 = c2_constant()
CPHI = 0.125
N0, M = 12000, 15000
Ns = N0 + 2 * np.arange(M)

NODD = np.arange(3, NMAX // 2 + 1, 2)
logden = np.zeros(NMAX + 1); logden[3:] = 1.0 / np.log(np.arange(3, NMAX + 1))


def G_gen(a, Ns):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD[NODD < N // 2]
        m = N - n
        out[i] = 2.0 * np.sum(a[n] * a[m] * weight(n / N)) / N
    return out


def E_gen(p, Ns):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD[NODD < N // 2]
        m = N - n
        out[i] = 2.0 * np.sum(p[n] * p[m] * weight(n / N)) / N
    return out


p_odd = logden * (np.arange(NMAX + 1) % 2)
E1 = E_gen(p_odd, Ns)
E2 = E_gen(p_odd, 2 * Ns)
vals = []
print("== iid Bernoulli(1/log n) on odd n, 8 independent draws ==")
for seed in range(8):
    rng = np.random.default_rng(seed)
    a = (rng.random(NMAX + 1) < p_odd).astype(float)
    s1 = G_gen(a, Ns) / E1 - 1
    s2 = G_gen(a, 2 * Ns) / E2 - 1
    c = corr(s1, s2)
    vals.append(c)
    print("  seed %d: c(2)=%+.4f   sd(rho)=%.3f  ac1=%+.3f" % (seed, c, s1.std(), corr(s1[:-1], s1[1:])))
print("  mean %+.4f  sd %.4f  (expected sd of mean %.4f)" % (np.mean(vals), np.std(vals), np.std(vals) / np.sqrt(len(vals))))

print("\n== real data: weighted-mass fluctuation vs rho ==")
NODD_all = np.arange(3, NMAX + 1, 2)
weighted = np.zeros(NMAX + 1)
# weighted prime mass in the support, weight w(n/N), normalised by its expectation
def mass(Ns):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD_all[(NODD_all > N / 3) & (NODD_all < 2 * N / 3)]
        lm = lam[n]
        lm2 = lam[N - n]
        out[i] = 2.0 * np.sum((lm > 0) * weight(n / N))
    return out


def massE(Ns):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD_all[(NODD_all > N / 3) & (NODD_all < 2 * N / 3)]
        out[i] = 2.0 * np.sum(weight(n / N) / np.log(n))
    return out


G1 = G_arr(Ns, lam, pp); S1 = sing_arr(Ns, spf, twoC2)
G2 = G_arr(2 * Ns, lam, pp); S2 = sing_arr(2 * Ns, spf, twoC2)
rho1 = G1 / (CPHI * S1); rho2 = G2 / (CPHI * S2)
m1 = mass(Ns) / massE(Ns) - 1
m2 = mass(2 * Ns) / massE(2 * Ns) - 1
print("  sd(rho)=%.4f  sd(mass)=%.4f" % (rho1.std(), m1.std()))
print("  corr(rho(N), mass(N))    = %+.4f" % corr(rho1, m1))
print("  corr(rho(2N), mass(2N))  = %+.4f" % corr(rho2, m2))
print("  corr(mass(N), mass(2N))  = %+.4f" % corr(m1, m2))
print("  corr(rho(N), mass(2N))   = %+.4f" % corr(rho1, m2))
