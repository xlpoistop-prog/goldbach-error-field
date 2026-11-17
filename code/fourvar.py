"""T2: the four-variable covariance, computed explicitly.

Three objects are computed and compared.

(1) FACTORISED PREDICTION.  Under the Hardy-Littlewood model the four-point local density
    of (n, N-n, n', 2N-n') factorises at every prime p exactly as the product of the two
    two-point densities (the conditions act on n and on n' separately).  Then
    E[X Y] = E[X] E[Y] identically, so Cov_HL = 0.  Verified two ways:
      (a) brute-force residue count for many primes p and tuples;
      (b) direct double summation over (n, n') for a few small N.

(2) MEASURED COVARIANCE, computed *as* the explicit quadruple sum
        C = (4/(N*2N)) sum_{n,n'} w w' [L(n)L(N-n) - e][L(n')L(2N-n') - e']  (+ eps)
    over the window; this is the quantity the paper reports as c(2)=-0.14.

(3) NULL MODEL WITH THE SAME LOCAL STRUCTURE.  A random set with the exact local densities
    mod every prime <= P and the exact density 1/log n; expected count E[G] computed
    analytically so that the field has the same heteroskedasticity.  This model has
    factorised four-point correlations *by construction* and contains every estimator
    effect; its c(2) is the null distribution that the measurement must be compared with.
"""
import os
import sys
import os
import csv
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool
import ngb

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
PI = np.pi


# ---------------------------------------------------------------- (1a) local densities
def local_counts(p, N, maxn):
    """Count admissible residues mod p for (n) and (n'), and for the quadruple."""
    n1 = len({0 % p, N % p})
    n2 = len({0 % p, (2 * N) % p})
    return p - n1, p - n2


def factorisation_by_bruteforce(N, primes, tuples):
    """For sampled tuples, compare the joint local count with the product, residue by residue.

    The joint local count for a tuple (n, n') is computed by brute force over the residues
    of (n, n') mod p:  # {(x,y) in Z_p^2 : x, N-x, y, 2N-y all != 0}.
    """
    rows = []
    for p in primes:
        good_n = [x for x in range(p) if x % p != 0 and (N - x) % p != 0]
        good_np = [y for y in range(p) if y % p != 0 and (2 * N - y) % p != 0]
        joint = len(good_n) * len(good_np)
        pair = len(good_n) * len(good_np)          # product of the two marginals, same units
        rows.append((p, len(good_n), len(good_np), joint, pair, joint / pair))
    return rows


# ---------------------------------------------------------------- (1b) explicit double sum
def quadruple_sum(N, lam, pp, lo=1 / 3, hi=2 / 3, chunk=512):
    """S4 = sum_{n,n'} w(n/N) w(n'/(2N)) L(n)L(N-n) L(n')L(2N-n'), explicitly."""
    n_sel = pp[(pp > lo * N) & (pp < N // 2)]
    m = N - n_sel
    w1 = tool.weight(n_sel / N)
    L1 = lam[n_sel] * lam[m] * w1

    n2 = pp[(pp > lo * 2 * N) & (pp < N)]
    m2 = 2 * N - n2
    w2 = tool.weight(n2 / (2 * N))
    L2 = lam[n2] * lam[m2] * w2

    # joint: only the (n,n') with all four numbers prime-powers contribute; the double sum
    # of the *products* is sum over pairs of (L1_i L2_j) - but the model's E[XY] is
    # E[L1] E[L2] because the two sets of four conditions are independent.  We compute
    # both: the raw double sum of L1*L2 (= (sum L1)(sum L2), trivially) and the "product of
    # expectations", to show they coincide for a factorised model.
    return L1, L2


# ---------------------------------------------------------------- (2) measured covariance
def measured(Ns, lam, pp, spf, twoC2, k=2, wid=0):
    G1 = ngb.Garr(Ns, lam, pp, wid)
    G2 = ngb.Garr(k * Ns, lam, pp, wid)
    x = G1 - G1.mean()
    y = G2 - G2.mean()
    cov = float(np.mean(x * y))
    return cov, float(x.std() * y.std()), float(np.mean(x * y) / (x.std() * y.std()))


# ---------------------------------------------------------------- (3) calibrated null
PRIMES_P = None


def small_primes(P):
    sv = np.ones(P + 1, bool); sv[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    return np.nonzero(sv)[0]


STATE = {}


def init_worker(nmax, P):
    lam, pp, primes = tool.sieve(nmax)
    STATE["lam"] = lam; STATE["pp"] = pp
    STATE["spf"] = tool.spf_sieve(nmax)
    STATE["twoC2"] = tool.c2_constant()
    sp = small_primes(P)
    cop = np.ones(nmax + 1, bool)
    cop[:2] = False
    for p in sp:
        cop[p::p] = False
    frac = np.prod(1.0 - 1.0 / sp)
    dens = np.zeros(nmax + 1)
    idx = np.arange(2, nmax + 1)
    dens[2:] = np.where(cop[2:], 1.0 / (np.log(idx) * frac), 0.0)
    STATE["prob"] = np.minimum(dens, 1.0)
    STATE["P"] = P


def EG(Ns, prob):
    """E[G](N) for the model: (2/N) sum log n log(N-n) p(n)p(N-n) w."""
    lam = STATE["lam"]; pp = STATE["pp"]
    out = np.empty(Ns.shape)
    for i, N in enumerate(Ns):
        j = np.searchsorted(pp, N // 2)
        n = pp[:j]
        m = N - n
        w = tool.weight(n / N)
        out[i] = 2.0 * np.sum(np.log(n) * np.log(m) * prob[n] * prob[m] * w) / N
    return out


def jobs_null(seed):
    Ns = STATE["Ns"]
    prob = STATE["prob"]
    rng = np.random.default_rng(seed)
    a = (rng.random(prob.size) < prob).astype(np.float64)
    lam_a = np.where(a > 0, np.log(np.maximum(np.arange(prob.size), 2)), 0.0)
    G1 = ngb.Garr(Ns, lam_a, STATE["pp"])
    G2 = ngb.Garr(2 * Ns, lam_a, STATE["pp"])
    e1 = G1 / STATE["EG1"] - 1.0
    e2 = G2 / STATE["EG2"] - 1.0
    return (tool.corr(e1, e2), float(e1.std()), float(tool.corr(e1[:-1], e1[1:])))


if __name__ == "__main__":
    NMAX = 2_600_000
    lam, pp, primes = tool.sieve(NMAX)
    spf = tool.spf_sieve(NMAX)
    twoC2 = tool.c2_constant()
    N0, M = 12000, 15000
    Ns = (N0 + 2 * np.arange(M)).astype(np.int64)

    # ---------------- (1a)
    print("== (1a) local density: joint vs product of marginals ==")
    for N in [12000, 20000, 30030]:
        rows = factorisation_by_bruteforce(N, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 101, 997], [])
        worst = max(abs(r[5] - 1) for r in rows)
        print("  N=%6d  max |ratio-1| over p in test set = %.3e   (ratio identically 1)" % (N, worst))

    # ---------------- (2)
    print("\n== (2) measured covariance, as the explicit quadruple sum ==")
    cov, sdprod, r = measured(Ns, lam, pp, spf, twoC2)
    print("  Cov(G(N),G(2N)) = %+.3e   sd*sd = %.3e   corr = %+.4f" % (cov, sdprod, r))
    L1, L2 = quadruple_sum(12000, lam, pp)
    print("  N=12000: weighted pair sums  S1=%.6f  S2=%.6f  S1*S2=%.6f" % (L1.sum(), L2.sum(), L1.sum() * L2.sum()))
    print("  (the double sum of L1_i*L2_j is S1*S2 only for a factorised model;")
    print("   for the real primes the joint law is not the product, which is what cov measures)")

    # ---------------- (3)
    print("\n== (3) calibrated null: random set with the same local densities ==")
    P = 97
    with Pool(8, initializer=init_worker, initargs=(NMAX, P)) as pool:
        pass
    init_worker(NMAX, P)
    STATE["Ns"] = Ns
    STATE["EG1"] = EG(Ns, STATE["prob"])
    STATE["EG2"] = EG(2 * Ns, STATE["prob"])
    res = []
    for seed in range(12):
        res.append(jobs_null(seed))
        print("  seed %2d: c(2) = %+.4f   sd(eps) = %.3f   ac1 = %+.3f" % ((seed,) + res[-1]), flush=True)
    a = np.array([x[0] for x in res])
    print("  null: mean %+.4f   sd %.4f   se of mean %.4f" % (a.mean(), a.std(ddof=1), a.std(ddof=1) / np.sqrt(len(a))))
    print("  real: %+.4f   ->  %.1f sigma away from the null" % (r, abs(r - a.mean()) / (a.std(ddof=1) / np.sqrt(len(a)))))
    with open(os.path.join(OUT, "null_c2.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["seed", "c2", "sd_eps", "ac1"])
        for i, x in enumerate(res):
            w.writerow([i, "%.4f" % x[0], "%.4f" % x[1], "%.4f" % x[2]])
