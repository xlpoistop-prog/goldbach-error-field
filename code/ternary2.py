"""Ternary Goldbach error field, numerically stable version.

    T(N) = (1/N^2) sum_{n1+n2+n3=N, N/6 <= ni <= N/2} Lambda(n1)Lambda(n2)Lambda(n3)

a *balanced box*: each of the three primes lies between N/6 and N/2, so no tiny element can
dominate (the analogue of the binary support [1/3,2/3]).  Flat weight, so T is computed
directly; the box is N-dependent, which rules out FFT but keeps everything exact and stable.

Normalisation: S3(N) = base * prod_{p|N,p>2} (p-2)(p-1)/(p^2-3p+3), plus an empirical smooth
N-dependence fitted in log N (a smooth factor cannot create fine-scale correlation).

Check: mean(T/S3 | 3|N) / mean(T | 3 does not divide) = 2/3.
"""
import os
import sys
import os
import numpy as np
from numba import njit
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool

X = 1_500_000
N0, M = 30_001, 6_000   # ternary needs ODD N
STATE = {}


@njit(cache=True, fastmath=True)
def G3(Ns, pp, lam):
    """T_box(N) = (1/N^2) sum_{n1+n2+n3=N} Lambda Lambda Lambda over the balanced box."""
    out = np.empty(Ns.size)
    npp = pp.size
    for i in range(Ns.size):
        N = Ns[i]
        lo = N / 6.0
        hi = N / 2.0
        j0 = np.searchsorted(pp, lo, side="right")
        j1 = np.searchsorted(pp, hi, side="left")
        s = 0.0
        for a in range(j0, j1):
            n1 = pp[a]
            l1 = lam[n1]
            # n2 in [max(lo, N-n1-hi), min(hi, N-n1-lo)]
            a0 = N - n1 - hi
            a1 = N - n1 - lo
            if a0 < lo:
                a0 = lo
            if a1 > hi:
                a1 = hi
            k0 = np.searchsorted(pp, a0, side="right")
            k1 = np.searchsorted(pp, a1, side="left")
            for b in range(k0, k1):
                n2 = pp[b]
                n3 = N - n1 - n2
                l3 = lam[n3]
                if l3 != 0.0:
                    s += l1 * lam[n2] * l3
        out[i] = s / (N * N)
    return out


def init_worker(nmax):
    lam, pp, primes = tool.sieve(nmax)
    STATE["lam"] = lam
    STATE["pp"] = pp
    STATE["primes"] = primes


def job(chunk):
    Ns = chunk
    return G3(Ns, STATE["pp"], STATE["lam"])


def s3rel(N):
    s = 1.0
    nn = N
    while nn % 2 == 0:
        nn //= 2
    d = 3
    while d * d <= nn:
        if nn % d == 0:
            s *= (d - 2.0) * (d - 1.0) / (d * d - 3.0 * d + 3.0)
            while nn % d == 0:
                nn //= d
        d += 2
    if nn > 2:
        s *= (nn - 2.0) * (nn - 1.0) / (nn * nn - 3.0 * nn + 3.0)
    return s


if __name__ == "__main__":
    Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
    print("ternary box count, N in [%d,%d], M=%d" % (Ns[0], Ns[-1], M), flush=True)
    with Pool(10, initializer=init_worker, initargs=(X,)) as pool:
        chunks = np.array_split(Ns, 10)
        parts = pool.map(job, chunks)
    T = np.concatenate(parts)
    print("  T(%.0f) = %.6f   mean %.6f" % (Ns[0], T[0], T.mean()), flush=True)

    # validation against a direct double loop
    lam, pp, primes = tool.sieve(X)
    for N in [30000, 30006]:
        s = 0.0
        lo, hi = N / 6.0, N / 2.0
        for n1 in pp[(pp > lo) & (pp < hi)]:
            for n2 in pp[(pp > max(lo, N - n1 - hi)) & (pp < min(hi, N - n1 - lo))]:
                n3 = N - n1 - n2
                if lam[n3]:
                    s += lam[n1] * lam[n2] * lam[n3]
        print("  check N=%d: numba %.6f  direct %.6f" % (N, T[np.searchsorted(Ns, N)], s / N ** 2), flush=True)

    S3 = np.array([s3rel(int(N)) for N in Ns])
    U = T / S3
    m3 = (Ns % 3 == 0)
    print("\nsingular-series check: mean(U | 3|N)/mean(U | 3∤N) = %.5f   (theory 0.66667)"
          % (U[m3].mean() / U[~m3].mean()), flush=True)

    lg = np.log(Ns.astype(float))
    cf = np.polyfit(lg, np.log(U), 3)
    eps = U / np.exp(np.polyval(cf, lg)) - 1.0
    print("field: sd %.4f  ac1 %+.3f" % (eps.std(), np.corrcoef(eps[:-1], eps[1:])[0, 1]), flush=True)

    print("\n k   c_ternary(k)   se      (binary at N0=12000 for comparison)")
    binref = {2: -0.140, 3: -0.123, 4: -0.063, 5: -0.151, 6: 0.024, 7: -0.158, 9: -0.039, 10: 0.036, 12: 0.004, 15: 0.019}
    for k in [2, 3, 4, 5, 6, 7, 9, 10, 12, 15]:
        Nk = k * Ns
        ok = Nk <= X
        Uk = T[ok] / np.array([s3rel(int(v)) for v in Nk[ok]])
        smk = np.exp(np.polyval(cf, np.log(Nk[ok].astype(float))))
        e2 = Uk / smk - 1.0
        c = tool.corr(eps[ok], e2)
        se = tool.block_se(eps[ok], e2, L=40, reps=200, seed=1)
        print(" %2d   %+.4f        %.4f    (%+.3f)" % (k, c, se, binref[k]), flush=True)
