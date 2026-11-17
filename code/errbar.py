"""How large is the error bar on c(2)?  Blocks vs realization-level coherence.

The within-block estimator pools many short-block correlations.  If those are independent,
se ~ sd/sqrt(B); if they are coherent over stretches, the true error is larger.  This script
measures the block-correlation autocorrelation for (a) the real primes, (b) the local-density
random model, and reports a Newey-West style inflated error.
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
N0, M = 12000, 15000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
CPHI = 0.125


def block_corrs(a, b, width):
    B = len(a) // width
    out = np.empty(B)
    for i in range(B):
        x = a[i * width:(i + 1) * width]
        y = b[i * width:(i + 1) * width]
        x = (x - x.mean()) / x.std()
        y = (y - y.mean()) / y.std()
        out[i] = np.mean(x * y)
    return out


def nw_se(v):
    """Newey-West standard error of the mean of v, with a simple automatic lag."""
    n = v.size
    m = v.mean()
    c = v - m
    v0 = np.mean(c * c) / n
    lag = int(4 * (n / 100.0) ** (2.0 / 9.0))
    s = v0
    for L in range(1, lag + 1):
        cov = np.mean(c[:-L] * c[L:])
        s += 2 * (1 - L / (lag + 1)) * cov / n
    return float(np.sqrt(max(s, 1e-12)))


def report(tag, a, b):
    print(tag)
    print("  full-window corr = %+.4f" % tool.corr(a, b))
    for width in [20, 40, 100]:
        v = block_corrs(a, b, width)
        print("  width %3d, %3d blocks: mean %+.4f  naive se %.4f  NW se %.4f  (sd %.3f, ac1 %+.3f)"
              % (width, v.size, v.mean(), v.std() / np.sqrt(v.size), nw_se(v), v.std(),
                 np.corrcoef(v[:-1], v[1:])[0, 1]))


# real
G1 = ngb.Garr(Ns, lam, pp)
G2 = ngb.Garr(2 * Ns, lam, pp)
S1 = tool.sing_arr(Ns, spf, twoC2)
S2 = tool.sing_arr(2 * Ns, spf, twoC2)
e1 = G1 / (CPHI * S1) - 1
e2 = G2 / (CPHI * S2) - 1
report("REAL primes", e1, e2)

# model
P = 97
sv = np.ones(P + 1, bool); sv[:2] = False
for i in range(2, int(P ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = False
sp = np.nonzero(sv)[0]
cop = np.ones(NMAX + 1, bool); cop[:2] = False
for p in sp:
    cop[p::p] = False
frac = float(np.prod(1.0 - 1.0 / sp))
idx = np.arange(2, NMAX + 1)
prob = np.zeros(NMAX + 1)
prob[2:] = np.where(cop[2:], np.minimum(1.0 / (np.log(idx) * frac), 1.0), 0.0)


def EG(Nsx):
    out = np.empty(Nsx.shape)
    for i, N in enumerate(Nsx):
        j = np.searchsorted(pp, N // 2)
        n = pp[:j]
        m = N - n
        out[i] = 2.0 * np.sum(np.log(n) * np.log(m) * prob[n] * prob[m] * tool.weight(n / N)) / N
    return out


EG1 = EG(Ns); EG2 = EG(2 * Ns)
for seed in [0, 3]:
    rng = np.random.default_rng(seed)
    a = (rng.random(NMAX + 1) < prob)
    lam_a = np.where(a, np.log(np.maximum(np.arange(NMAX + 1), 2)), 0.0)
    f1 = ngb.Garr(Ns, lam_a, pp) / EG1 - 1
    f2 = ngb.Garr(2 * Ns, lam_a, pp) / EG2 - 1
    report("MODEL seed %d" % seed, f1, f2)
