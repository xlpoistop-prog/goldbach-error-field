"""T2, final version.

(1) four-point factorisation: exact (ratio 1 to machine precision) - see fourvar.py.
(2) measurable covariance of the two pair counts (with the correct normalisation).
(3) null model with factorised four-point correlations and the exact local densities,
    evaluated with the *same* within-block estimator that the real data is evaluated with.
"""
import os
import sys
import os
import csv
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool
import ngb

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)

NMAX = 2_600_000
lam, pp, primes = tool.sieve(NMAX)
spf = tool.spf_sieve(NMAX)
twoC2 = tool.c2_constant()
N0, M = 12000, 15000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
CPHI = 0.125


def local_corr(a, b, blocks):
    w = len(a) // blocks
    vals = []
    for i in range(blocks):
        x = a[i * w:(i + 1) * w]
        y = b[i * w:(i + 1) * w]
        x = (x - x.mean()) / x.std()
        y = (y - y.mean()) / y.std()
        vals.append(np.mean(x * y))
    vals = np.array(vals)
    return float(vals.mean()), float(vals.std() / np.sqrt(blocks))


# ---------- real data
G1 = ngb.Garr(Ns, lam, pp)
G2 = ngb.Garr(2 * Ns, lam, pp)
S1 = tool.sing_arr(Ns, spf, twoC2)
S2 = tool.sing_arr(2 * Ns, spf, twoC2)
e1 = G1 / (CPHI * S1) - 1.0
e2 = G2 / (CPHI * S2) - 1.0
print("real data (M=%d):" % M)
print("  sd(eps(N)) = %.4f   sd(eps(2N)) = %.4f" % (e1.std(), e2.std()))
print("  full-window corr = %+.4f" % tool.corr(e1, e2))
for B in [30, 150, 375, 750]:
    m, se = local_corr(e1, e2, B)
    print("  within-block corr, %3d blocks of width %4d: %+.4f +- %.4f" % (B, M // B, m, se))
cov = float(np.mean((e1 - e1.mean()) * (e2 - e2.mean())))
print("  Cov = %+.3e   (14%% of sd*sd = %.3e)" % (cov, 0.14 * e1.std() * e2.std()))

# ---------- null model
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
        w = tool.weight(n / N)
        out[i] = 2.0 * np.sum(np.log(n) * np.log(m) * prob[n] * prob[m] * w) / N
    return out


EG1 = EG(Ns)
EG2 = EG(2 * Ns)
print("\nnull model: random set, local densities exact mod every p <= %d, density 1/log n" % P)
print("  mean number of model primes below 2.6e6: %.0f (the primes: %d)" % (prob.sum(), primes.size))
rows = []
for seed in range(16):
    rng = np.random.default_rng(seed)
    a = (rng.random(NMAX + 1) < prob)
    lam_a = np.where(a, np.log(np.maximum(idx0 := np.arange(NMAX + 1), 2)), 0.0)
    ga = ngb.Garr(Ns, lam_a, pp)
    gb = ngb.Garr(2 * Ns, lam_a, pp)
    f1 = ga / EG1 - 1.0
    f2 = gb / EG2 - 1.0
    full = tool.corr(f1, f2)
    loc, locse = local_corr(f1, f2, 375)
    rows.append((seed, full, loc, locse, f1.std(), f2.std()))
    print("  seed %2d: full %+.4f   within-block(375) %+.4f +- %.4f   sd %.3f/%.3f"
          % (seed, full, loc, locse, f1.std(), f2.std()), flush=True)

arr = np.array([[r[1] for r in rows], [r[2] for r in rows]])
print("\n  null full-window:  mean %+.4f  sd %.4f  se %.4f" % (arr[0].mean(), arr[0].std(ddof=1), arr[0].std(ddof=1) / 4))
print("  null within-block: mean %+.4f  sd %.4f  se %.4f" % (arr[1].mean(), arr[1].std(ddof=1), arr[1].std(ddof=1) / 4))
print("  real within-block: %+.4f +- 0.009" % local_corr(e1, e2, 750)[0])
print("  => deviation in units of the null spread: %.1f sigma" % (abs(local_corr(e1, e2, 750)[0] - arr[1].mean()) / (arr[1].std(ddof=1) / 2.0)))
with open(os.path.join(OUT, "null_fourvar.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["seed", "c2_full", "c2_within", "se_within", "sd_eps1", "sd_eps2"])
    for r in rows:
        w.writerow([r[0]] + ["%.4f" % x for x in r[1:]])
print("written lab/out/null_fourvar.csv")
