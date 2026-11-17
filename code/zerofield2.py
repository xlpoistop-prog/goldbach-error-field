"""Zero-driven field test -- parallel version.

For the window N = 12000 + 2i (i < M) we build
    L(N) = sum_j N^{i g_j} / g_j                        (one zero: density-like)
    P(N) = |L(N)|^2 = sum_jk N^{i(g_j+g_k)}/(g_j g_k)   (two zeros: pair-count-like)
with real positive coefficients (the most favourable choice: no cancellation, so whatever
correlation appears is an upper bound over phase choices).

L is computed once per N-grid and shared by both fields; the N-grids are chunked and processed
by a fork pool, so the work is spread over all cores.
"""
import sys
import time
import numpy as np
import multiprocessing as mp

N0, M = 12000, 15000
_G = None
_GRID = None


def z_of_t(t):
    n = np.floor(np.sqrt(np.maximum(t, 2.0) / (2 * np.pi)))
    th = t / 2 * np.log(np.maximum(t, 2.0) / (2 * np.pi)) - t / 2 - np.pi / 8
    Z = np.zeros_like(t)
    for k in range(1, int(n.max()) + 1):
        m = n >= k
        if m.any():
            Z[m] += 2.0 * np.cos(th[m] - t[m] * np.log(k)) / np.sqrt(k)
    return Z


def zeros_up_to(tmax, dt=0.02, t0=10.0, iters=55):
    t = np.arange(t0, tmax + dt, dt)
    Z = z_of_t(t)
    idx = np.nonzero(np.sign(Z[:-1]) * np.sign(Z[1:]) < 0)[0]
    a, b = t[idx].copy(), t[idx + 1].copy()
    Za = Z[idx].copy()
    for _ in range(iters):
        c = 0.5 * (a + b)
        Zc = z_of_t(c)
        same = np.sign(Zc) == np.sign(Za)
        a = np.where(same, c, a)
        b = np.where(same, b, c)
        Za = np.where(same, Zc, Za)
    return 0.5 * (a + b)


def _L_chunk(bounds):
    a, b = bounds
    lg = np.log(_GRID[a:b])
    out = np.empty(b - a, dtype=np.complex128)
    step = 400
    for s in range(0, b - a, step):
        e = min(s + step, b - a)
        out[s:e] = np.exp(1j * lg[s:e, None] * _G[None, :]) @ (1.0 / _G)
    return a, out


def L_of(Ns, workers=120):
    """L(N) for a grid, parallel over chunks."""
    global _GRID
    _GRID = Ns
    n = Ns.size
    edges = np.linspace(0, n, max(workers * 8, 8) + 1).astype(int)
    bounds = [(int(a), int(b)) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    out = np.empty(n, dtype=np.complex128)
    with mp.get_context("fork").Pool(workers) as pool:
        for a, o in pool.imap_unordered(_L_chunk, bounds, chunksize=1):
            out[a:a + o.size] = o
    return out


def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = np.sqrt(np.dot(a, a) * np.dot(b, b))
    return float(np.dot(a, b) / d) if d > 0 else float("nan")


def rvm(T):
    """Riemann-von Mangoldt main term for the number of zeros with 0 < gamma < T."""
    import math
    return (T / (2 * math.pi)) * (math.log(T / (2 * math.pi)) - 1.0) + 7.0 / 8.0

if __name__ == "__main__":
    TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 40000.0
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    t0 = time.time()
    g = zeros_up_to(TMAX)
    _G = g
    print("%d zeros up to t=%.0f  (%.1fs)" % (g.size, TMAX, time.time() - t0), flush=True)

    base = (N0 + 2 * np.arange(M)).astype(np.float64)
    data = {2: -0.1401, 3: -0.1227, 4: -0.0632, 6: 0.0235, 10: 0.0357}

    t1 = time.time()
    L1 = L_of(base, workers)
    print("L(base) done %.1fs" % (time.time() - t1), flush=True)
    for kind in ("linear", "pair"):
        A = L1.real if kind == "linear" else (L1.real ** 2 + L1.imag ** 2)
        print("\n--- %s field: lag1 = %+.3f (data -0.027), lag10 = %+.3f (data -0.031), sd = %.4g ---"
              % (kind, corr(A[:-1], A[1:]), corr(A[:-10], A[10:]), A.std()), flush=True)
        print("   k     c_sim      c_data")
        for k in (2, 3, 4, 6, 10):
            t2 = time.time()
            Lk = L_of(k * base, workers)
            B = Lk.real if kind == "linear" else (Lk.real ** 2 + Lk.imag ** 2)
            print("   %-3d  %+.4f    %+.4f      (%.1fs)" % (k, corr(A, B), data[k], time.time() - t2),
                  flush=True)
    print("\ntotal %.1fs" % (time.time() - t0), flush=True)
