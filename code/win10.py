"""Ten windows at the SAME scale: is the amplitude of c(2) stable, or modulated?

Each window: N0_j = N0 + j*step, M consecutive even N.  Per window we measure c(2), c(3) and the
mean weighted prime mass mu.  Then:
   chi2 of the ten c(2) against a common value   -> is there excess scatter at all?
   corr(c2, mu) over the ten windows             -> is the scatter the local prime density?
"""
import sys
import time
import numpy as np
from numba import njit
import gcfar2 as g2

CPHI = g2.CPHI


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


_LAM = _PRIMES = None


def _mass_task(t):
    a, b, Ns, out = t
    mass_at(Ns[a:b], _LAM, _PRIMES, out[a:b])
    return a, out[a:b].copy()


def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return float(np.dot(a, b) / d) if d > 0 else float("nan")


def block_se(a, b, width=40, reps=300, seed=1):
    B = a.size // width
    if B < 3:
        return float("nan")
    x = a[:B * width].reshape(B, width); y = b[:B * width].reshape(B, width)
    x = x - x.mean(1, keepdims=True); y = y - y.mean(1, keepdims=True)
    v = (x * y).sum(1) / np.sqrt((x * x).sum(1) * (y * y).sum(1))
    rng = np.random.default_rng(seed)
    return float(np.array([v[rng.integers(0, B, B)].mean() for _ in range(reps)]).std())


def expected_mass(Ns, nu=1501):
    u = np.linspace(1.0 / 3, 2.0 / 3, nu)
    c = np.cos(3 * np.pi * (u - 0.5)); w = c ** 4
    out = np.empty(Ns.size)
    for i, N in enumerate(Ns):
        out[i] = N * np.trapezoid(w / np.log(N * u), u)
    return out


if __name__ == "__main__":
    import multiprocessing as mp
    global_n0 = float(sys.argv[1]) if len(sys.argv) > 1 else 1e9
    M = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    NW = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    step = int(sys.argv[4]) if len(sys.argv) > 4 else 2000000
    workers = 120

    t0 = time.time()
    starts = [int(global_n0) + j * step for j in range(NW)]
    X = int(2 * (max(starts) + 2 * M)) + 10
    print("scale %.0e   %d windows, M=%d, step=%d   sieve to %.2e" % (global_n0, NW, M, step, X), flush=True)
    lam, primes = g2.sieve_lam(X, verbose=True)
    _LAM, _PRIMES = lam, primes
    g2._LAM, g2._PRIMES = lam, primes
    ctx = mp.get_context("fork")
    rows = []
    for j, s0 in enumerate(starts):
        Ns = (s0 + 2 * np.arange(M)).astype(np.int64)
        g2._NS = Ns
        nc = max(workers * 4, 8)
        e = np.linspace(0, M, nc + 1).astype(int)
        chunks = [(int(a), int(b)) for a, b in zip(e[:-1], e[1:]) if b > a]
        eps = {}
        for k in (1, 2, 3):
            out = np.empty(M)
            with ctx.Pool(workers) as pool:
                for kk, a, o in pool.imap_unordered(g2._task, [(k, a, b) for a, b in chunks], chunksize=1):
                    out[a:a + o.size] = o
            Sk = np.array([g2.sing(int(v)) for v in k * Ns])
            eps[k] = out / (CPHI * Sk) - 1.0
        mc = np.empty(M)
        with ctx.Pool(workers) as pool:
            for a, o in pool.imap_unordered(_mass_task, [(a, b, Ns, mc) for a, b in chunks], chunksize=1):
                mc[a:a + o.size] = o
        mu = float(np.mean(mc - expected_mass(Ns)))
        c2, c3 = corr(eps[1], eps[2]), corr(eps[1], eps[3])
        se2, se3 = block_se(eps[1], eps[2]), block_se(eps[1], eps[3])
        rows.append((j, s0, c2, se2, c3, se3, mu))
        print("  win %2d  N0=%.6e  c2=%+.4f(%.4f)  c3=%+.4f(%.4f)  mu=%+.3e   [%.0fs]"
              % (j, s0, c2, se2, c3, se3, mu, time.time() - t0), flush=True)

    c2 = np.array([r[2] for r in rows]); se2 = np.array([r[3] for r in rows])
    c3 = np.array([r[4] for r in rows]); mu = np.array([r[6] for r in rows])
    w = 1 / se2**2
    mean = np.average(c2, weights=w)
    chi2 = float(np.sum(((c2 - mean) / se2) ** 2))
    print("\n  c2: weighted mean %+.4f +- %.4f   scatter %.4f   chi2 = %.1f / %d dof"
          % (mean, 1 / np.sqrt(w.sum()), c2.std(ddof=1), chi2, NW - 1))
    print("  c2 vs mu: corr = %+.3f      c2 vs c3: corr = %+.3f" % (corr(c2, mu), corr(c2, c3)))
    np.savez("win10_%g.npz" % global_n0, c2=c2, se2=se2, c3=c3, mu=mu, starts=np.array(starts))
    print("saved win10_%g.npz   total %.0fs" % (global_n0, time.time() - t0), flush=True)
