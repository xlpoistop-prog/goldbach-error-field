"""Is the multiplicative correlation a generic property of 'autocorrelation-like'
counting, or specific to the real primes?  Surrogate models."""
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tool import *

NMAX = 2_600_000
lam, pp, primes = sieve(NMAX)
spf = spf_sieve(NMAX)
twoC2 = c2_constant()
CPHI = 0.125
N0, M = 12000, 15000
Ns = N0 + 2 * np.arange(M)

NODD = np.arange(3, NMAX // 2 + 1, 2)


def G_gen(a, Ns):
    """G(N) = (2/N) sum_{n<N/2, n odd} a[n] a[N-n] w(n/N);  only odd n matter."""
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        n = NODD[NODD < N // 2]
        m = N - n
        out[i] = 2.0 * np.sum(a[n] * a[m] * weight(n / N)) / N
    return out


logden = np.zeros(NMAX + 1)
logden[3:] = 1.0 / np.log(np.arange(3, NMAX + 1))

print("== surrogate A: iid Bernoulli(1/log n) on odd n ==")
for seed in range(3):
    rng = np.random.default_rng(seed)
    a = (rng.random(NMAX + 1) < logden).astype(np.float64)
    a[::2] = 0.0
    G1 = G_gen(a, Ns); G2 = G_gen(a, 2 * Ns); G3 = G_gen(a, 3 * Ns)
    r1 = G1 / G1.mean(); r2 = G2 / G2.mean(); r3 = G3 / G3.mean()
    print("  seed %d: c(2)=%+.4f  c(3)=%+.4f  mean #pairs~%.0f"
          % (seed, corr(r1, r2), corr(r1, r3), a.sum()))

print("\n== surrogate B: local-structure random sieve (exclude 2,3,5,7,11,13) ==")
mask = np.zeros(NMAX + 1, dtype=bool)
mask[3:] = True
for p in [3, 5, 7, 11, 13]:
    mask[p::p] = False
mask[::2] = False
dens = np.where(mask, logden, 0.0)
for seed in range(3):
    rng = np.random.default_rng(100 + seed)
    a = (rng.random(NMAX + 1) < np.minimum(dens, 1.0)).astype(np.float64)
    G1 = G_gen(a, Ns); G2 = G_gen(a, 2 * Ns); G3 = G_gen(a, 3 * Ns)
    r1 = G1 / G1.mean(); r2 = G2 / G2.mean(); r3 = G3 / G3.mean()
    print("  seed %d: c(2)=%+.4f  c(3)=%+.4f  #terms=%d"
          % (seed, corr(r1, r2), corr(r1, r3), int(np.sum(a))))

print("\n== real primes, stability across sub-windows ==")
G1 = G_arr(Ns, lam, pp); S1 = sing_arr(Ns, spf, twoC2)
G2 = G_arr(2 * Ns, lam, pp); S2 = sing_arr(2 * Ns, spf, twoC2)
rho1 = G1 / (CPHI * S1); rho2 = G2 / (CPHI * S2)
print("  full window      c(2)=%+.4f" % corr(rho1, rho2))
for b in range(5):
    sl = slice(b * M // 5, (b + 1) * M // 5)
    print("  sub-window %d (N=%d..%d)  c(2)=%+.4f" % (b, Ns[sl][0], Ns[sl][-1], corr(rho1[sl], rho2[sl])))

print("\n== real primes, decile table (R6 replication) ==")
qs = np.quantile(rho1, np.arange(0, 1.01, 0.1))
dec = np.clip(np.searchsorted(qs, rho1, side="right") - 1, 0, 9)
print("   q      E[rho(N)]   E[rho(2N)]")
for d in range(10):
    m = dec == d
    print("  %2d    %+8.4f    %+8.4f   (n=%d)" % (d + 1, rho1[m].mean(), rho2[m].mean(), m.sum()))
print("  regression slope rho(2N) on rho(N) = %+.4f" % np.polyfit(rho1, rho2, 1)[0])
