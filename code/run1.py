import sys, time
import numpy as np
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tool import *

t0 = time.time()
NMAX = 2_600_000
lam, pp, primes = sieve(NMAX)
spf = spf_sieve(NMAX)
twoC2 = c2_constant()
print("sieve done %.1fs  #pp=%d  2C2=%.6f" % (time.time() - t0, pp.size, twoC2))

CPHI = 0.125  # int_{1/3}^{2/3} cos^4(3 pi (u-1/2)) du
N0, M = 12000, 15000
Ns = N0 + 2 * np.arange(M)

t0 = time.time()
GN = G_arr(Ns, lam, pp)
print("G done %.1fs" % (time.time() - t0))
SN = sing_arr(Ns, spf, twoC2)
e = CPHI * SN - GN
print("N=%d  G=%.6f  cphi*S=%.6f  ratio=%.6f" % (Ns[0], GN[0], CPHI * SN[0], GN[0] / (CPHI * SN[0])))

sig = lambda N: 0.7 * np.sqrt(sing_arr(np.asarray(N, dtype=np.int64), spf, twoC2) * 0.091146) * np.log(N) / np.sqrt(np.asarray(N, dtype=np.float64))
Y1 = e / sig(Ns)
print("sd(e)=%.6f   sd(Y)=%.4f  mean(Y)=%.4f" % (e.std(), Y1.std(), Y1.mean()))

print("\n== reproduce c(k), M=%d ==" % M)
ks = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 16, 21, 25, 27]
t0 = time.time()
res = {}
for k in ks:
    Nk = k * Ns
    Gk = G_arr(Nk, lam, pp)
    Sk = sing_arr(Nk, spf, twoC2)
    ek = CPHI * Sk - Gk
    Yk = ek / sig(Nk)
    res[k] = corr(Y1, Yk)
    print("  k=%2d  c(k)=%+.4f" % (k, res[k]), flush=True)
print("time %.1fs" % (time.time() - t0))

print("\n== c(2) split-half ==")
h = M // 2
G2 = G_arr(2 * Ns, lam, pp); S2 = sing_arr(2 * Ns, spf, twoC2); e2 = CPHI * S2 - G2; Y2 = e2 / sig(2 * Ns)
print("  full  c(2)=%+.4f" % corr(Y1, Y2))
print("  first half  c(2)=%+.4f" % corr(Y1[:h], Y2[:h]))
print("  second half c(2)=%+.4f" % corr(Y1[h:], Y2[h:]))
print("  block-bootstrap se (L=50) = %.4f" % block_se(Y1, Y2, L=50))
print("  block-bootstrap se (L=200)= %.4f" % block_se(Y1, Y2, L=200))

print("\n== autocorrelation of Y(N) in N (lags in index, step 2) ==")
for lag in [1, 2, 5, 10, 20, 50, 100]:
    print("  lag %4d: %+.4f" % (lag, corr(Y1[:-lag], Y1[lag:])))

np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_M15000.npz"), Ns=Ns, GN=GN, SN=SN, e=e, Y1=Y1)
