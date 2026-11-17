"""Is the amplitude modulation of c(k) a local-density effect?

For each big window (N0, M) we split the index range into S sub-windows and measure, inside each:
    c2 = within-sub correlation of eps(N), eps(2N)
    c3 = the same for k=3
    mu = the mean weighted prime mass deviation m(N) in that sub-window

Three tests, in increasing specificity:
  (i)  scatter of c2 across sub-windows against its own bootstrap error  (is there extra scatter?)
  (ii) corr(c2, c3) across sub-windows   (is the modulation COMMON to all k?)
  (iii) corr(c2, mu)  and corr(c3, mu)   (is it the LOCAL PRIME DENSITY?)
"""
import sys
import time
import numpy as np
from numba import njit
import gcfar2 as g2

CPHI = g2.CPHI
WITH_MASS = False
PHI2_INT = 0.125            # int_0^1 cos^4 = 0.125
UG = None


@njit(cache=True, fastmath=True, nogil=True)
def mass_at(Ns, lam, primes, out):
    for i in range(Ns.size):
        N = Ns[i]
        lo = N // 3 + 1
        hi = 2 * N // 3
        j0 = np.searchsorted(primes, lo, side="right")
        j1 = np.searchsorted(primes, hi, side="left")
        s = 0.0
        for j in range(j0, j1):
            u = primes[j] / N
            c = np.cos(9.42477796076938 * (u - 0.5))
            s += c * c * c * c
        out[i] = s


def expected_mass(Ns, nu=2001):
    """sum_{n in (N/3,2N/3)} phi^2(n/N)/log n, by quadrature."""
    u = np.linspace(1.0 / 3, 2.0 / 3, nu)
    c = np.cos(3 * np.pi * (u - 0.5))
    w = c ** 4
    out = np.empty(Ns.size)
    for i, N in enumerate(Ns):
        lg = np.log(N * u)
        out[i] = N * np.trapezoid(w / lg, u)
    return out


def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return float(np.dot(a, b) / d) if d > 0 else float("nan")


def sub_corr(x, y, S):
    """within-sub-window correlation for each of S contiguous blocks, + bootstrap se"""
    B = x.size // S
    vals, ses = [], []
    for s in range(S):
        a = x[s * B:(s + 1) * B]
        b = y[s * B:(s + 1) * B]
        vals.append(corr(a, b))
        w = max(10, B // 40)
        n = a.size // w
        if n >= 3:
            xa = a[:n * w].reshape(n, w); ya = b[:n * w].reshape(n, w)
            xa = xa - xa.mean(1, keepdims=True); ya = ya - ya.mean(1, keepdims=True)
            v = (xa * ya).sum(1) / np.sqrt((xa * xa).sum(1) * (ya * ya).sum(1))
            rng = np.random.default_rng(s)
            ses.append(float(np.array([v[rng.integers(0, n, n)].mean() for _ in range(200)]).std()))
        else:
            ses.append(float("nan"))
    return np.array(vals), np.array(ses)


def one(n0, m, S=10, workers=120, tag=""):
    import multiprocessing as mp
    t0 = time.time()
    Ns = (n0 + 2 * np.arange(m)).astype(np.int64)
    X = int(2 * Ns[-1]) + 10
    print("\n=== N0=%.0e  M=%d  S=%d  sieve to %.2e ===" % (n0, m, S, X), flush=True)
    lam, primes = g2.sieve_lam(X, verbose=False)
    print("  sieve: %d primes  %.0fs" % (primes.size, time.time() - t0), flush=True)
    g2._LAM, g2._PRIMES, g2._NS = lam, primes, Ns

    nc = max(workers * 8, 8)
    edges = np.linspace(0, m, nc + 1).astype(int)
    chunks = [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:]) if b > a]

    eps = {}
    ctx = mp.get_context("fork")
    for k in (1, 2, 3):
        out = np.empty(m)
        tasks = [(k, a, b) for a, b in chunks]
        with ctx.Pool(workers) as pool:
            for kk, a, o in pool.imap_unordered(g2._task, tasks, chunksize=1):
                out[a:a + o.size] = o
        Sk = np.array([g2.sing(int(v)) for v in k * Ns])
        eps[k] = out / (CPHI * Sk) - 1.0
        print("  k=%d  sd(eps)=%.5f  %.0fs" % (k, eps[k].std(), time.time() - t0), flush=True)

    # local mass (test iii only -- it is the slow single-threaded step)
    if WITH_MASS:
        cnt = np.empty(m)
        mass_at(Ns, lam, primes, cnt)
        mu = cnt - expected_mass(Ns)
        mu = mu / mu.std()
        print("  mass done %.0fs" % (time.time() - t0), flush=True)
    else:
        mu = np.zeros(m)

    c2, se2 = sub_corr(eps[1], eps[2], S)
    c3, se3 = sub_corr(eps[1], eps[3], S)
    # the mass of each sub-window
    B = m // S
    msub = np.array([mu[s * B:(s + 1) * B].mean() for s in range(S)])

    print("  ---- sub-window table ----")
    print("   s     c2      +-      c3      +-      mu")
    for s in range(S):
        print("  %2d  %+.3f  %.3f   %+.3f  %.3f   %+.2f" % (s, c2[s], se2[s], c3[s], se3[s], msub[s]))

    chi2_c2 = float(np.sum(((c2 - np.average(c2, weights=1 / se2**2)) / se2) ** 2))
    print("  (i)   scatter of c2: chi2 = %.1f / %d dof (se-weighted mean %+.3f)"
          % (chi2_c2, S - 1, np.average(c2, weights=1 / se2**2)))
    print("  (ii)  corr(c2, c3) = %+.3f  over %d sub-windows" % (corr(c2, c3), S))
    print("  (iii) corr(c2, mu) = %+.3f    corr(c3, mu) = %+.3f" % (corr(c2, msub), corr(c3, msub)))
    np.savez("sub_%s.npz" % tag, c2=c2, se2=se2, c3=c3, se3=se3, mu=msub)
    return dict(c2=c2, c3=c3, mu=msub, se2=se2, se3=se3)


if __name__ == "__main__":
    import json
    out = {}
    for n0, m, tag in [(3e8, 20000, "n3e8"), (1e9, 10000, "n1e9")]:
        out[tag] = one(n0, m, S=10, tag=tag)
