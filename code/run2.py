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

Gd = {}
Sd = {}
Gd[1] = G_arr(Ns, lam, pp)
Sd[1] = sing_arr(Ns, spf, twoC2)
for k in ks:
    Gd[k] = G_arr(k * Ns, lam, pp)
    Sd[k] = sing_arr(k * Ns, spf, twoC2)

e = {k: CPHI * Sd[k] - Gd[k] for k in [1] + ks}
rho = {k: Gd[k] / (CPHI * Sd[k]) for k in [1] + ks}          # ratio G/main
rel = {k: e[k] / Sd[k] for k in [1] + ks}                    # error in units of S

print("sd(e1)=%.6f  mean ratio rho(1)=%.6f  sd(rho)=%.6f" % (e[1].std(), rho[1].mean(), rho[1].std()))

print("\n== which normalization carries the effect? (k: c) ==")
norms = {
    "Y (theirs: sqrt(S)*logN/sqrtN)": lambda k, N: e[k] / (np.sqrt(Sd[k]) * np.log(N) / np.sqrt(N.astype(float))),
    "const sigma=1               ": lambda k, N: e[k],
    "sqrt(S), no N-dep          ": lambda k, N: e[k] / np.sqrt(Sd[k]),
    "logN/sqrtN, no S           ": lambda k, N: e[k] * np.sqrt(N.astype(float)) / np.log(N),
    "S*logN/sqrtN               ": lambda k, N: e[k] / (Sd[k] * np.log(N) / np.sqrt(N.astype(float))),
    "error/S  (relative error)  ": lambda k, N: e[k] / Sd[k],
    "raw G                      ": lambda k, N: Gd[k],
    "G/S                        ": lambda k, N: rho[k],
    "log G                      ": lambda k, N: np.log(Gd[k]),
}
for name, f in norms.items():
    row = []
    for k in ks:
        row.append("%+.4f" % corr(f(1, Ns), f(k, k * Ns)))
    print("  %s  " % name + "  ".join("k=%d:%s" % (k, r) for k, r in zip(ks, row)))
