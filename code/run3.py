import sys, time
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tool import *

NMAX = 1_200_000
lam, pp, primes = sieve(NMAX)
spf = spf_sieve(NMAX)
twoC2 = c2_constant()
CPHI = 0.125
N0, M = 12000, 15000
Ns = N0 + 2 * np.arange(M)
ks = [2, 3, 4, 6]

Gd, Sd = {}, {}
Gd[1] = G_arr(Ns, lam, pp); Sd[1] = sing_arr(Ns, spf, twoC2)
for k in ks:
    Gd[k] = G_arr(k * Ns, lam, pp); Sd[k] = sing_arr(k * Ns, spf, twoC2)
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "d1.npz"),
         Ns=Ns, **{f"G{k}": Gd[k] for k in Gd}, **{f"S{k}": Sd[k] for k in Sd})

rho = {k: Gd[k] / (CPHI * Sd[k]) for k in [1] + ks}
print("mean rho1=%.6f sd=%.6f" % (rho[1].mean(), rho[1].std()))

print("\n== binned means of rho(1) (30 bins of 500 consecutive N, i.e. width 1000) ==")
B = 30
w = M // B
binmean = np.array([rho[1][i * w:(i + 1) * w].mean() for i in range(B)])
binc = np.array([Ns[i * w:(i + 1) * w].mean() for i in range(B)])
print("  bin centers:", " ".join("%d" % c for c in binc[:10]), "...")
print("  bin means  :", " ".join("%+.4f" % v for v in binmean))
print("  sd of bin means = %.4f   (expected from pure noise: %.4f)"
      % (binmean.std(), rho[1].std() / np.sqrt(w)))

print("\n== correlation of the BINNED field (kills all structure below width %d) ==" % (2 * w))
for k in ks:
    bk = np.array([rho[k][i * w:(i + 1) * w].mean() for i in range(B)])
    print("  k=%d  binned corr = %+.4f" % (k, corr(binmean, bk)))

print("\n== correlation within narrow blocks (kills all structure above block width) ==")
for B2 in [5, 10, 30]:
    w2 = M // B2
    for k in ks:
        prod = []
        for i in range(B2):
            a = rho[1][i * w2:(i + 1) * w2]
            b = rho[k][i * w2:(i + 1) * w2]
            a = (a - a.mean()) / a.std()
            b = (b - b.mean()) / b.std()
            prod.append(np.mean(a * b))
        print("  B=%2d k=%d  local corr = %+.4f" % (B2, k, np.mean(prod)))

print("\n== cross-correlation rho(N) vs rho(N+Delta) ==")
for D in [2, 20, 200, 1000, 5000, 10000]:
    j = D // 2
    print("  Delta=%5d  corr = %+.4f" % (D, corr(rho[1][:-j], rho[1][j:])))

print("\n== correlation of rho(1) with simple arithmetic features ==")
Ns_int = Ns.astype(np.int64)
lamN = np.zeros(M)
mu2 = np.zeros(M)
om = np.zeros(M)
for i, N in enumerate(Ns_int):
    n = int(N); f = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        f.append(n)
    om[i] = len(f)
    s = 0.0
    for p in f:
        s += np.log(p) / p
    lamN[i] = s
    mu2[i] = (-1) ** len(f)
print("  corr(rho, omega(N))      = %+.4f" % corr(rho[1], om))
print("  corr(rho, (-1)^omega(N)) = %+.4f" % corr(rho[1], mu2))
print("  corr(rho, sum log p/p)   = %+.4f" % corr(rho[1], lamN))
print("  corr(rho, log N)         = %+.4f" % corr(rho[1], np.log(Ns)))
print("  corr(rho, N)             = %+.4f" % corr(rho[1], Ns.astype(float)))
# local model R^2: indicators p|N for p<=29 plus omega
ps = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
X = [np.ones(M)]
for p in ps:
    X.append((Ns_int % p == 0).astype(float))
X.append(om)
X = np.array(X).T
beta, *_ = np.linalg.lstsq(X, rho[1], rcond=None)
pred = X @ beta
print("  local model (p|N for p<=29 + omega): R^2 = %.4f" % (1 - np.sum((rho[1] - pred) ** 2) / np.sum((rho[1] - rho[1].mean()) ** 2)))
