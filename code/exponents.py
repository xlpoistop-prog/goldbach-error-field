"""Is the error field a function of the exponents v_p(N)?  (The decisive test between
"a property of N's arithmetic" and "a property of the prime configuration".)

If eps(N) carried a multiplicative-noise component m(N) = prod_p (1 + xi_p(v_p(N))) with
independent xi_p, then (i) its correlation kernel would be multiplicative in k -- exactly the
observed pattern -- and (ii) it would be a function of the exponent vector.  Test (ii) by
regressing eps on the cells of (v_2, v_3, v_5, v_7, ...) and computing the (adjusted) R^2.
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
N0, M = 12000, 30000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
e = ngb.Garr(Ns, lam, pp) / (0.125 * tool.sing_arr(Ns, spf, twoC2)) - 1


def vp_arr(Ns, p):
    v = np.zeros(Ns.size, dtype=np.int64)
    n = Ns.copy()
    while True:
        m = (n % p == 0)
        if not m.any():
            break
        v[m] += 1
        n[m] //= p
    return v


def cells_r2(ps, cap):
    idx = np.zeros(Ns.size, dtype=np.int64)
    for p in ps:
        v = np.minimum(vp_arr(Ns, p), cap)
        idx = idx * (cap + 1) + v
    lab, inv = np.unique(idx, return_inverse=True)
    cnt = np.bincount(inv)
    means = np.bincount(inv, weights=e) / cnt
    pred = means[inv]
    sse = np.sum((e - pred) ** 2)
    sst = np.sum((e - e.mean()) ** 2)
    r2 = 1 - sse / sst
    k = lab.size
    r2a = 1 - (1 - r2) * (Ns.size - 1) / (Ns.size - k)
    return r2, r2a, k, cnt.min()


print("regressing eps(N) on the cells of the exponent vector (M=%d, sd(eps)=%.4f)" % (M, e.std()))
for ps, cap in [([2], 3), ([2, 3], 3), ([2, 3, 5], 2), ([2, 3, 5, 7], 2), ([2, 3, 5, 7, 11], 2)]:
    r2, r2a, k, cmin = cells_r2(ps, cap)
    print("  primes %-18s cells %4d (min %3d pts):  R^2 = %.4f   adjusted = %+.4f"
          % (str(ps), k, cmin, r2, r2a))

# what R^2 would a multiplicative-noise component of amplitude a produce?
# if eps = a*m + noise with m unit variance and m a function of the cells, R^2 ~ a^2/var(eps)
print("\n  a multiplicative component of amplitude a would give R^2 ~ (a/%.3f)^2 * (fraction of m"
      % e.std())
print("  resolved by the cells). For the observed c(p)=-0.135 one needs a ~ 0.135*sd(eps) = %.4f"
      % (0.135 * e.std()))
print("  -> R^2 of order %.3f if the cells resolve m; observe the table above." % (0.135 ** 2))
# sanity: a pure multiplicative surrogate should be caught
rng = np.random.default_rng(0)
sp = [2, 3, 5, 7]
xi = {p: rng.normal(0, 0.14, 4) for p in sp}
m = np.ones(Ns.size)
for p in sp:
    m *= 1.0 + xi[p][np.minimum(vp_arr(Ns, p), 3)]
m = (m - m.mean()) / m.std()
sur = 0.05 * m + 0.05 * rng.standard_normal(Ns.size)
e_backup = e
e = (sur - sur.mean()) / sur.std() * e_backup.std()
r2, r2a, k, cmin = cells_r2([2, 3, 5, 7], 2)
print("\n  control: a synthetic field with a 5%% multiplicative component on p<=7 gives")
print("           R^2 = %.4f (adjusted %+.4f) on the same cells" % (r2, r2a))
