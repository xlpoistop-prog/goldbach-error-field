"""Clean, independent re-implementation of the 'Goldbach error field' measurement.

G(N) = (2/N) * sum_{n<N/2} Lam(n) Lam(N-n) w(n/N),   w(u)=cos^4(3*pi*(u-1/2)) on [1/3,2/3]
(so w = phi^2 with phi = cos^2(3 pi (u-1/2)), c_phi = int phi^2 = 1/8 = 0.125)

S(N) = 2 C2 prod_{p|N, p>2} (p-1)/(p-2)
e(N) = c_phi S(N) - G(N)

Nothing here uses the original scripts; this is a from-scratch implementation
used to test the reported claims.
"""
import numpy as np

C2 = None


def sieve(NMAX):
    """Lam(p^k) = log p, else 0, as float64 array; plus sorted prime-power list."""
    sv = np.ones(NMAX + 1, dtype=bool)
    sv[:2] = False
    for i in range(2, int(NMAX ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    primes = np.nonzero(sv)[0]
    lam = np.zeros(NMAX + 1, dtype=np.float64)
    for p in primes:
        lp = np.log(float(p))
        pk = int(p)
        while pk <= NMAX:
            lam[pk] = lp
            pk *= p
    pp = np.nonzero(lam > 0)[0]
    return lam, pp, primes


def spf_sieve(NMAX):
    spf = np.zeros(NMAX + 1, dtype=np.int64)
    for i in range(2, NMAX + 1):
        if spf[i] == 0:
            spf[i::i] = np.where(spf[i::i] == 0, i, spf[i::i])
    return spf


def c2_constant(npairs=200000):
    """2*C2, C2 = prod_{p>2}(1-1/(p-1)^2), computed by a long sieve."""
    N = 4 * npairs
    sv = np.ones(N + 1, dtype=bool)
    sv[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    ps = np.nonzero(sv)[0]
    ps = ps[ps > 2]
    return 2.0 * np.prod(1.0 - 1.0 / (ps - 1.0) ** 2)


def sing_arr(Ns, spf, twoC2):
    """S(N) for an int array Ns, using smallest-prime-factor table."""
    Ns = np.asarray(Ns, dtype=np.int64)
    out = np.full(Ns.shape, twoC2, dtype=np.float64)
    n = Ns.copy()
    for _ in range(40):
        act = n > 1
        if not act.any():
            break
        p = np.where(act, spf[n], 1)
        isodd = p > 2
        out = np.where(isodd, out * (p - 1.0) / np.maximum(p - 2.0, 1e-30), out)
        n = np.where(act, n // np.maximum(p, 1), 1)
        # strip remaining powers of the same prime
        while True:
            again = (n % np.maximum(p, 1) == 0) & act & (n > 1) & (p > 1)
            if not again.any():
                break
            n = np.where(again, n // np.maximum(p, 1), n)
    return out


def weight(u):
    x = 3.0 * np.pi * (u - 0.5)
    c = np.cos(x)
    return np.where((u > 1.0 / 3.0) & (u < 2.0 / 3.0), c * c * c * c, 0.0)


def G_arr(Ns, lam, pp):
    """G(N) for each N in Ns (float64), all n+m=N, no coprimality restriction."""
    Ns = np.asarray(Ns, dtype=np.int64)
    out = np.empty(Ns.shape, dtype=np.float64)
    L = lam.size - 1
    for i, N in enumerate(Ns):
        half = N // 2
        j = np.searchsorted(pp, half)
        n = pp[:j]
        if n.size and n[-1] >= L + 1:
            n = n[n < L]
        m = N - n
        out[i] = 2.0 * np.sum(lam[n] * lam[m] * weight(n / N)) / N
    return out


def zs(v):
    v = np.asarray(v, dtype=np.float64)
    return (v - v.mean()) / v.std()


def corr(a, b):
    a = zs(a); b = zs(b)
    return float(np.mean(a * b))


def block_se(a, b, L=50, reps=400, seed=0):
    """Moving-block bootstrap standard error of corr(a,b)."""
    rng = np.random.default_rng(seed)
    a = zs(a); b = zs(b)
    n = a.size
    m = n // L
    vals = np.empty(reps)
    for r in range(reps):
        st = rng.integers(0, n - L + 1, size=m)
        ia = (st[:, None] + np.arange(L)[None, :]).ravel()
        aa = a[ia]; bb = b[ia]
        aa = (aa - aa.mean()) / aa.std()
        bb = (bb - bb.mean()) / bb.std()
        vals[r] = np.mean(aa * bb)
    return float(vals.std())
