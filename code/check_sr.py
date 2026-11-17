import sys
sys.stdout.reconfigure(encoding="utf-8")


def odd_primes(n):
    s = list(range(3, n + 1, 2))
    out = []
    while s:
        p = s[0]
        out.append(p)
        s = [x for x in s if x % p]
    return out


P = odd_primes(300000)
C2 = 1.0
for p in P:
    C2 *= (1 - 1.0 / (p - 1) ** 2)
twoC2 = 2 * C2
print("C2 = %.9f   2C2 = %.9f   1/(2C2) = %.9f" % (C2, twoC2, 1 / twoC2))

for N in (2, 6, 30, 2310, 30030, 12000, 131072):
    fac = [p for p in P if N % p == 0]
    S = twoC2
    for p in fac:
        S *= (p - 1.0) / (p - 2.0)

    R1 = 1.0                      # prod over p | N  (includes p = 2)
    if N % 2 == 0:
        R1 *= 0.5
    for p in fac:
        R1 *= 1.0 - 1.0 / p

    R2 = 1.0                      # prod over p > 2, p does not divide N
    for p in P:
        if N % p != 0:
            R2 *= 1.0 + 1.0 / (p * (p - 2.0))

    print("N=%7d  fac=%-18s  S=%.6f  R1=%.6f  R2=%.6f  S*R1*R2=%.6f"
          % (N, str(fac), S, R1, R2, S * R1 * R2))

# --- gate ---
assert abs(S * R1 * R2 - 1.0) < 1e-6, "the identity S(N)R_N(0)=1 failed: %.10f" % (S * R1 * R2)
print("OK: S*R1*R2 = 1 to 1e-6")
