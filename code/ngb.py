"""numba-accelerated Goldbach error field.

G(N) = (2/N) sum_{N/3<n<N/2} Lam(n) Lam(N-n) w(n/N)
with w(u) = cos^4(3 pi (u-1/2)) (default) or flat on (1/3,2/3), etc.
"""
import sys
import os
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool

PI = 3.14159265358979323846


@njit(cache=True, fastmath=True)
def _G(Ns, j0, j1, pp, lam, wid):
    out = np.empty(Ns.size)
    for i in range(Ns.size):
        N = Ns[i]
        s = 0.0
        inv = 1.0 / N
        for j in range(j0[i], j1[i]):
            n = pp[j]
            m = N - n
            lm = lam[m]
            if lm != 0.0:
                u = n * inv
                if wid == 0:            # cos4
                    if u > 0.3333333333333333 and u < 0.6666666666666666:
                        c = np.cos(3.0 * PI * (u - 0.5))
                        s += lam[n] * lm * (c * c) * (c * c)
                elif wid == 1:          # flat
                    s += lam[n] * lm
                elif wid == 2:          # narrow cos4, support (5/12,7/12)
                    if u > 0.4166666666666667 and u < 0.5833333333333334:
                        c = np.cos(9.0 * PI * (u - 0.5))
                        s += lam[n] * lm * (c * c) * (c * c)
                elif wid == 3:          # wide cos4, support (1/4,3/4)
                    if u > 0.25 and u < 0.75:
                        c = np.cos(1.5 * PI * (u - 0.5))
                        s += lam[n] * lm * (c * c) * (c * c)
                elif wid == 4:          # beta^2
                    if u > 0.0 and u < 1.0:
                        s += lam[n] * lm * 16.0 * (u * (1.0 - u)) ** 2
        out[i] = 2.0 * s * inv
    return out


def Garr(Ns, lam, pp, wid=0, lo=1 / 3, hi=2 / 3):
    Ns = np.asarray(Ns, dtype=np.int64)
    j0 = np.searchsorted(pp, lo * Ns.astype(np.float64), side="right")
    j1 = np.searchsorted(pp, Ns.astype(np.float64) / 2.0, side="left")
    return _G(Ns, j0, j1, pp, lam, wid)


def rho(Ns, lam, pp, spf, twoC2, wid=0, lo=1 / 3, hi=2 / 3):
    G = Garr(Ns, lam, pp, wid, lo, hi)
    S = tool.sing_arr(Ns, spf, twoC2)
    return G / S
