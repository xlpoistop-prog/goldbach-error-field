"""T1: variance stabilisation and inverse-variance weighting.

The field eps(N) is heteroskedastic: sd(eps(N)) ~ log N / sqrt(N).  A Pearson correlation of
two such vectors weights the window unevenly.  Here the correlation is recomputed after
(i) dividing by the known sd profile, (ii) inverse-variance weighting, (iii) on a window of
constant relative width.
"""
import os
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool
import ngb

NMAX = 2_600_000
lam, pp, primes = tool.sieve(NMAX)
spf = tool.spf_sieve(NMAX)
twoC2 = tool.c2_constant()
CPHI = 0.125


def fields(N0, M, k=2):
    Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
    e1 = ngb.Garr(Ns, lam, pp) / (CPHI * tool.sing_arr(Ns, spf, twoC2)) - 1
    e2 = ngb.Garr(k * Ns, lam, pp) / (CPHI * tool.sing_arr(k * Ns, spf, twoC2)) - 1
    return Ns, e1, e2


def wcorr(x, y, w):
    mx = np.average(x, weights=w); my = np.average(y, weights=w)
    cxy = np.average((x - mx) * (y - my), weights=w)
    return cxy / np.sqrt(np.average((x - mx) ** 2, weights=w) * np.average((y - my) ** 2, weights=w))


Ns, e1, e2 = fields(12000, 15000)
g1 = np.log(Ns) / np.sqrt(Ns)
g2 = np.log(2 * Ns) / np.sqrt(2 * Ns)
print("window N0=12000, M=15000, k=2")
print("  raw Pearson                      : %+.4f" % tool.corr(e1, e2))
print("  after dividing by logN/sqrtN     : %+.4f" % tool.corr(e1 / g1, e2 / g2))
print("  inverse-variance weighted Pearson: %+.4f" % wcorr(e1, e2, 1.0 / (g1 * g2) ** 2))
print("  Spearman (rank)                  : %+.4f"
      % np.corrcoef(np.argsort(np.argsort(e1)), np.argsort(np.argsort(e2)))[0, 1])
m = (np.abs(e1 - e1.mean()) / e1.std() < 3) & (np.abs(e2 - e2.mean()) / e2.std() < 3)
print("  trimmed at |z|<3                 : %+.4f  (keeps %.0f%%)" % (tool.corr(e1[m], e2[m]), 100 * m.mean()))

print("\nsame statistic on a window of constant relative width [N0, 1.3 N0]:")
for N0, M in [(12000, 1800), (50000, 7500), (200000, 30000)]:
    Ns_, a, b = fields(N0, M)
    print("  N0=%7d  M=%5d: raw %+.4f   stabilised %+.4f"
          % (N0, M, tool.corr(a, b), tool.corr(a / (np.log(Ns_) / np.sqrt(Ns_)), b / (np.log(2 * Ns_) / np.sqrt(2 * Ns_)))))
