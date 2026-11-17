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
Ns_int = Ns.astype(np.int64)

G1 = G_arr(Ns, lam, pp); S1 = sing_arr(Ns, spf, twoC2)
G2 = G_arr(2 * Ns, lam, pp); S2 = sing_arr(2 * Ns, spf, twoC2)
rho1 = G1 / (CPHI * S1); rho2 = G2 / (CPHI * S2)

print("== (1) how local is the effect?  pooled corr within blocks ==")
for B in [30, 75, 150, 300, 750]:
    w = M // B
    prod = []
    for i in range(B):
        a = rho1[i * w:(i + 1) * w]; b = rho2[i * w:(i + 1) * w]
        a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
        prod.append(np.mean(a * b))
    print("  blocks=%3d (width %5d in N): c(2) = %+.4f +- %.4f"
          % (B, 2 * w, np.mean(prod), np.std(prod) / np.sqrt(B)))

print("\n== (2) raw unweighted count of prime pairs R(N), no log weights ==")


def Rcount(Ns, lam, pp):
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        half = N // 2
        j = np.searchsorted(pp, half)
        n = pp[:j]
        m = N - n
        out[i] = 2.0 * np.sum((lam[n] > 0) * (lam[m] > 0)) / N
    return out


R1 = Rcount(Ns, lam, pp); R2 = Rcount(2 * Ns, lam, pp)
print("  c(2) for R/S = %+.4f" % corr(R1 / S1, R2 / S2))
print("  c(2) for log R   = %+.4f" % corr(np.log(R1), np.log(R2)))

print("\n== (3) is there a mod-q structure?  R^2 of rho on N mod q, and c(2) after removing it ==")
for q in [6, 12, 30, 60, 210, 420, 2310, 4620, 30030]:
    idx = (Ns_int // 2) % q
    lab, inv = np.unique(idx, return_inverse=True)
    cnt = np.bincount(inv)
    means = np.bincount(inv, weights=rho1) / cnt
    pred = means[inv]
    r2 = 1 - np.sum((rho1 - pred) ** 2) / np.sum((rho1 - rho1.mean()) ** 2)
    r1r = rho1 - pred
    idx2 = (2 * Ns_int // 2) % q
    lab2, inv2 = np.unique(idx2, return_inverse=True)
    cnt2 = np.bincount(inv2)
    means2 = np.bincount(inv2, weights=rho2) / cnt2
    r2r = rho2 - means2[inv2]
    print("  q=%6d (#classes %4d): R^2=%7.4f   c(2)|residual = %+.4f" % (q, lab.size, r2, corr(r1r, r2r)))

print("\n== (4) surrogates with the theoretical mean (kills smooth trend) ==")
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


for tag, p in [("iid 1/log n", np.where(np.arange(NMAX + 1) % 2 == 1, logden, 0.0))]:
    rng = np.random.default_rng(7)
    a = (rng.random(NMAX + 1) < p).astype(float)
    Eg1 = E_gen(p, Ns); Eg2 = E_gen(p, 2 * Ns)
    for t in range(2):
        s1 = G_gen(a, Ns) / Eg1 - 1
        s2 = G_gen(a, 2 * Ns) / Eg2 - 1
        print("  %s seed %d: c(2)=%+.4f  sd=%.3f" % (tag, t, corr(s1, s2), s1.std()))
        rng = np.random.default_rng(1000 + t)
        a = (rng.random(NMAX + 1) < p).astype(float)
