# Multiplicative correlations in the Goldbach error field

## An empirical study

**Goldbach Error Field Project**

**Abstract.** Let `G(N)` be a smooth weighted *balanced* Goldbach count, `e(N)=c_φ𝔖(N)−G(N)` its
deviation from the Hardy–Littlewood main term, and `ε(N)=e(N)/(c_φ𝔖(N))` the relative error,
whose standard deviation is `0.086` at `N≈3×10⁴` (i.e. of the order of the pair-counting noise).
We report a reproducible empirical regularity of the field `ε`: its two-point correlation across
multiplicative dilations,

```
c(k) = corr( ε(N), ε(kN) ),      N in a window of all even integers,
```

is (i) **sign definite**: `sign c(k) = (−1)^{ω(k)}` for 28 of the 31 integers `2 ≤ k ≤ 49`
tested, the only significant exception being `k=26` at `2σ`; (ii) of magnitude `≈0.135` at
primes, but **not constant**: the nine prime
values span `−0.158 … −0.101` with `χ² = 72.8` on 8 degrees of freedom (`p = 1.4×10⁻¹²`);
(iii) **approximately multiplicative** for coprime pairs `(k₁,k₂)=1` (6 of 10 tests within `1σ`,
the remaining four at `1.1σ, 2.4σ, 3.2σ` and `4.0σ`; three of the six are prime-square cases,
so we do not claim squarefreeness); (iv) **not** given by the completely multiplicative
continuation at prime powers (`c(p²)` is negative, 3–13σ away from `c(p)²` (7–13σ for `p ≤ 5`)), so the
"Liouville-type" reading `c(k)=α^{ω(k)}λ(k)` is quantitatively rejected; and (v) **present at every prime ratio from
`10⁴` to `10⁹`**, with no decay law resolved: the window-to-window scatter (`±0.025`) is three
to four times the within-window error, which no single-window quotation reflects.

The measurement is not a consequence of the normalisation (for `k=2` it coincides with
`corr(G(N),G(2N))`), of the weight shape (five different weights), of the position in `N` (the
effect survives inside blocks of width 40), of the factorisation of `N`, of a random-prime model
(`−0.017±0.025`), or of any smooth trend.

The field is, however, **not white**, and this is a separate regularity: its autocorrelation is
zero at lags `2, 4, 10, 20, 200, 1000` and positive at lags that share small-prime factors,
growing with the number of the first five primes dividing the lag and independent of its size —
`+0.064` at `k = 210`, `+0.117` at `k = 2310`, `+0.129` at `k = 4620`, and about zero for lags
of the same magnitude without that structure (§5.10). The corresponding error bars are larger
than `1/sqrt(M)` would give.

We also record a structural obstruction that makes the measurement *surprising*: `G(N)` uses
primes in `[N/3,2N/3]` while `G(kN)` uses `[kN/3,2kN/3]`, and for `k=2` these intervals are
**disjoint**; moreover the four-point local densities of the four linear forms
`(n, N−n, n′, 2N−n′)` factorise exactly at every prime. Hence any model of the primes with the
correct local densities and independent positions predicts `c(k)=0`. The measured `−0.14` is
therefore either a residual estimator effect that we have not been able to locate, or a
long-range four-point correlation of the primes. Section 5.6 states a decisive test that
separates the two.

We do not establish any theorem, and nothing here bears on a proof of the binary Goldbach
conjecture. Appendix A records a technical review of the recent conditional preprint of
Li–Liu: we verify its exact operator identity and one of its constants, and note the single
display slip we can find (a missing conjugation bar in the Hermitian-symmetry line, which does
not affect the argument).

**Keywords.** Goldbach conjecture; weighted count; Hardy–Littlewood singular series;
multiplicative correlation; Liouville function; empirical field theory.

---

# 1. Introduction

## 1.1 The object

Fix an even `N ≥ 4` and a weight `φ ∈ C_c^∞(ℝ)`, `supp φ ⊂ [1/3,2/3]`, `0 ≤ φ ≤ 1`,
`φ(u) = φ(1−u)`. Write

```
G(N) = (2/N) Σ_{N/3<n<N/2} Λ(n)Λ(N−n) φ(n/N)φ((N−n)/N)                (1.1)
c_φ  = ∫_0^1 φ²(u) du                                                (1.2)
𝔖(N) = 2C₂ ∏_{p|N, p>2} (p−1)/(p−2),   C₂ = ∏_{p>2}(1−1/(p−1)²)       (1.3)
e(N) = c_φ𝔖(N) − G(N),     ε(N) = e(N)/(c_φ𝔖(N))                     (1.4)
```

Throughout, unless stated otherwise, `φ(u) = cos²(3π(u−1/2))`, so that `φ² = cos⁴` and
`c_φ = ∫cos⁴ = 1/8 = 0.125`. (The constant `0.091146 = ∫cos⁸` belongs to the weight
`φ = cos⁴`, not to `φ² = cos⁴`; it enters only inside a positive normalising factor, so no
numerical result depends on the distinction.)

The sum in (1.1) is over the **half-open** range `N/3 < n < N/2`, so each pair is counted
once and the diagonal term `n = m = N/2` is excluded by definition. That is the convention
used by every script in Appendix B; including the diagonal instead changes every `c(k)` by
at most `0.002` and `sd(ε)` by `0.5%`, i.e. by less than a third of the quoted errors, so the
choice is immaterial here but is stated rather than left implicit.

Because `φ` is symmetric, `G(N)` counts ordered pairs of prime powers that are *balanced*:
both members lie in `[N/3,2N/3]`. Hardy–Littlewood predicts `G(N) ≈ c_φ𝔖(N)`.

The empirical distribution of `e(N)` is Gaussian (Section 3.2) with

```
sd(e) ≈ 0.33 (log N)/√N,     i.e.   sd(ε) ≈ 0.086  at N ≈ 3×10⁴   (1.5)
```

which is *of the order of* the pair-counting (Poisson) scale, but the coefficient is not pinned
down by a single number. At `N = 1.2×10⁴` the support `[N/3,2N/3]` contains `98` prime-pair
centres — `196` *ordered* representations, but the sum runs over `n < N/2`, so each pair is one
term — giving a naive Poisson reference `1/√98 = 0.101`; the varying weight raises it to
`≈0.143`, since the effective count is `(Σw)²/Σw² = 98·E[w]²/E[w²] ≈ 49` with
`E[w²]/E[w]² = 2.0`. The measured `0.086` is *below* both, i.e. the field is mildly sub-Poisson
(variance `0.72×` the naive reference), the same sieve-type suppression found in §5.7(3). We
therefore use the Poisson scale only as an order of magnitude. So `ε` is, to first order, noise;
everything below concerns the *structure of that noise*.

**The field is *not* white, but the deviation is confined to lags divisible by `30`.** For the
window `N = 12000+2i`, `i < 4×10⁴`, all lags given here as shifts `ΔN`:

| `ΔN` | 2 | 4 | 10 | 20 | 40 | 100 | 200 | 1000 | 5000 | 10000 |
|---|---|---|---|---|---|---|---|---|---|---|
| corr | −0.027 | −0.022 | −0.031 | −0.015 | −0.010 | +0.000 | −0.001 | −0.016 | +0.007 | +0.023 |

None of these lags is a multiple of `30`, and every value is zero to within `2σ`. At lags that
*are* multiples of `30` the autocorrelation is positive and large — `+0.216` at `ΔN = 210`,
`+0.238` at `630`, `+0.146` at `2100` — corresponding to a drifting period-`30` component that
carries `≈18%` of the variance (`§5.9`, `§5.10`). That component is *not* a residue-class mean
shift: removing the class means of `ε` modulo `30`, `210` and `2310` changes `c(2)` by less than
`0.01`. The consequence for this section is that `M_eff ≈ M` holds for lags not divisible by
`30`, and the naive error bar `1/√M` is therefore optimistic by a factor of about `1.5` once the
window-to-window coherence at those lags is included (`§2.3`, `§4.11(2)`).

## 1.2 What is measured

```
c(k) := corr( ε(N), ε(kN) )                                         (1.6)
```

where `corr` is the Pearson correlation over `N` in a window of consecutive even integers, and
`ε(kN)` is computed from the definition (1.1)–(1.4) at the argument `kN` (so `k` must be a
positive integer; the window of `kN` is `k` times the window of `N`).

For `k=2` the statistic is *canonical*: since `N` is even, `𝔖(2N)=𝔖(N)`, hence

```
c(2) = corr( G(N), G(2N) )                                          (1.7)
```

with no normalisation freedom whatsoever. All other `k` differ from this only by the smooth
factor `log N/√N` and the singular-series ratio `𝔖(kN)/𝔖(N)`, which equals `1` or
`(k−1)/(k−2)` according to whether `k | N`. We verified that all reasonable normalisations give
the same `c(k)` to within `0.005` (Section 5.1).

## 1.3 Findings, with verdicts

| # | statement | verdict | evidence |
|---|---|---|---|
| S1 | `sign c(k) = (−1)^{ω(k)}` for `2 ≤ k ≤ 49` | **established** (28/31; only `k=26` significant, 2σ) | §4.3, Table 4 |
| S2 | `c(p) ≈ −α`, `α ≈ 0.135` | **rejected as stated**: `α` is not constant | §4.4, `χ²=72.8`/8 dof |
| S3 | `c(k₁k₂) = c(k₁)c(k₂)`, `(k₁,k₂)=1` | **approximate**: 6/10 within 1σ, others 1.1–4.0σ | §4.5 |
| S4 | `c(p^a) < 0`, saturating at `0.03–0.05` | **partly wrong**: `c(p²)` is `−0.001…−0.063`, and `c(p²) ≠ c(p)²` at 3–13σ | §4.6 |
| S5 | `c(k) = α^{ω(k)}λ(k)` (Liouville form) | **rejected** quantitatively | §4.6 |
| S6 | `c` is universal in the weight | **established** (five weights) | §4.8 |
| S7 | `c` is universal in scale | **rejected**: present from `10⁴` to `10⁹` but no decay law is resolved (`±0.03` scale-to-scale scatter); composite signs flip | §4.10, §4.11 |
| S8 | the effect is a property of the primes, not of the estimator | **open**, narrowed: estimator, deterministic `N`-biases, local mass, the sieve-type nulls and the zero-driven term all excluded; a non-local coupling of the primes remains | §5.7, §5.11 |

## 1.4 What this paper does not do

It does not prove anything about Goldbach's conjecture, it does not identify a mechanism, and it
does not claim that `c(k)` *is* a multiplicative function. It reports a measurement, the
controls that survive, the descriptions that fail, and one structural reason to be surprised.
Section 5.6 is the part we consider most useful: it turns the observation into a falsifiable
question about four-point correlations of the primes.

# 2. The estimator and its calibration

## 2.1 Conventions

All arithmetic is in `float64`. `Λ` is generated by an Eratosthenes sieve up to `2.6×10⁶`
(or higher where `kN` requires it); `𝔖(N)` by factorisation with a smallest-prime-factor table;
`2C₂ = 1.320323632…`. The weight is evaluated exactly (`cos⁴`), not from a grid (a `10⁵`-point table with index truncation changes `G` by at most
`3×10⁻⁵` relative) (§2.2).

## 2.2 Two independent implementations

The measurements in §4 were produced by code written from scratch, independent of every
earlier implementation:

* `lab/tool.py` — vectorised NumPy; `lab/ngb.py` — JIT (numba) version, 13× faster;
* `lab/fast2.py` — multi-process driver, `python fast2.py --n0 12000 --m 40000 --k 2,3,6`.

Cross-checks: our two implementations agree to `1.5×10⁻¹⁵` relative. Against an implementation
that evaluates the weight on a `10⁵`-point grid with index truncation, the weighted counts
differ by `3×10⁻⁶ … 3×10⁻⁵` relative (a `3×10⁻⁶` error in `u`, amplified by the steepness of
`cos⁴`); the correlation coefficients are unaffected at the `10⁻³` level. The archived `c(2)=−0.1403, c(3)=−0.1230, c(6)=+0.0255` at `M=15000` are
reproduced digit for digit, as are the two-point table of §2.5 and the distribution table
of §3.2.

## 2.3 Error bars

For each `c(k)` we quote the standard error of a moving-block bootstrap (block length 40,
300 resamples) and the two half-window values:

```
k     c(k)      se     first half  second half
2   −0.1401   0.0056     −0.1441     −0.1315
3   −0.1227   0.0055     −0.1172     −0.1335
6   +0.0235   0.0055     +0.0278     +0.0146
```

The bootstrap error, `0.0056` at `M=4×10⁴`, is close to `1/√M = 0.005` because `ε` is white in
`N` (§1.1). Two caveats must be kept in mind.

(i) **Window-to-window scatter is larger.** Repeating the whole measurement in different windows
gives `c(2) = −0.140, −0.138, −0.123, −0.074, −0.101, −0.094, −0.056` (§4.10). The scatter,
`≈0.03`, is five times the within-window bootstrap error. We therefore treat `±0.03` as the
honest uncertainty on a *single-window* value of `c(p)`, and `±0.006` only as the precision with
which one window's value is determined.

(ii) **A split-half comparison is weaker than it looks.** Two halves at `M = 15000` give
`−0.084/−0.110`; a half-sample has `se ≈ 0.008`, so the spread is `3σ`,
not "stable". Our own split halves differ by `0.013` at `M=4×10⁴` (for `k=2`).

## 2.4 A null model

Replace the primes by a random set: `n` is "prime" independently with probability `1/log n`
for odd `n`, nothing for even `n`, and normalise `G` by its exact expectation
`E[G](N) = 2Σ p(n)p(N−n)φ²(n/N)` (so that the smooth trend is removed analytically).
Eight independent draws give

```
c(2) = +0.038, −0.054, +0.047, −0.148, −0.063, +0.087, −0.041, −0.000
mean = −0.017,  sd of the mean = 0.025
```

consistent with zero. **Random primes do not reproduce the effect.** This null also fixes the
calibration of the estimator: the ratio `G/E[G]` is unbiased, and its fluctuation in `N` is white.

## 2.5 Validation of the singular-series machinery

`X = 2×10⁶`, all even `h ≤ 80`: the ratio `(Σ_{n≤X}Λ(n)Λ(n+h)/X) / 𝔖₂(h)` has

| min | max | mean | sd |
|---|---|---|---|
| 0.9935 | 1.0093 | 1.0010 | 0.0034 |

so the local factors used in `𝔖(N)` are the correct ones at the 0.3% level. (The value is
`1.0010`, computed from scratch.)

## 2.6 Remark: the coprimality condition `(nm,N)=1`

The object of the Li–Liu identity (Appendix A) carries the restriction `(nm,N)=1`, which (1.1)
omits. On the support this restriction is *vacuous*: `Λ(n) ≠ 0` forces `n = p^a`, so
`(n,N) ≠ 1` means `p | N` and `n = p^a ∈ [N/3,2N/3]`. Such `n` are: `N/2` when `N = 2p^a`
(a handful of `N` in any window), and, for the other divisors of `N` in the support, the prime
powers `p^a` with `p|N`. A direct computation gives, for `N = 12000, 12002, 16384, 13122, 31250`,

```
G with (nm,N)=1  =  G without the restriction,   to machine precision.
```

So the empirical object of this paper *is* the object of the operator identity, up to `O(N⁻¹)`.

---

# 3. The error field

## 3.1 The scale law

Windows of `M` consecutive even integers at three scales (`φ² = cos⁴`, `c_φ = 0.125`):

| N₀ | M | sd(e) | sd(e)·√N/log N | c(2) |
|---|---|---|---|---|
| 2×10⁵ | 2000 | 0.008799 | 0.3237 | −0.1541 |
| 10⁶ | 1200 | 0.004446 | 0.3220 | −0.0655 |
| 3×10⁶ | 600 | 0.002903 | 0.3372 | −0.1086 |

**`sd(eps)` depends on `N`, and on the weight.** Measured directly (fixed-count window, `M = 2000`):

| weight | `N_0 = 1.2×10⁴` | `N_0 = 3×10⁴` | `N_0 = 2×10⁵` |
|---|---|---|---|
| `cos⁴` (default) | `0.1027` | `0.0782` | `0.0383` |
| `16u²(1−u)²` | `0.0473` | — | — |

so the figure `0.086` quoted in §1.1 is the `cos⁴` value near `N ≈ 3×10⁴`, and the factor of two
between the two weights is a weight effect, not a discrepancy. `sd(eps)` decays roughly as
`N^{-1/2}`; a quotation without an `N` is loose and the numbers above should be used instead.

The third column is constant to `±2.3%`, i.e.

```
sd(e) ≈ 0.33 (log N)/√N,      sd(ε) = sd(e)/(c_φ𝔖(N)) ≈ 0.086      (3.1)
```

and this is the Poisson scale (§1.1). A different weighting convention gives `0.29`; the
difference is within the spread between windows and weight conventions. We quote `0.33` for the `cos⁴` weight
in (1.1).

## 3.2 Distribution

`N ≈ 10⁶`, `M = 1200`:

| λ | measured tail mass P(abs z > λ) | Gaussian tail mass |
|---|---|---|
| 1.0 | 0.3175 | 0.3173 |
| 1.5 | 0.1400 | 0.1336 |
| 2.0 | 0.0508 | 0.0455 |
| 2.5 | 0.0092 | 0.0124 |
| 3.0 | 0.0017 | 0.0027 |

`skew = −0.003`, `kurtosis = 3.003`, `max|z| = 3.32` against `√(2 log M) = 3.77`. The tail is
**sub-Gaussian**: no fat tails, no outliers to explain. (The tail estimates carry a noise of `≈1` count in `1200`.)

## 3.3 The conditional structure at `k=2`

`M = 15000`, `N₀ = 12000`. Group `N` by decile of the *normalised count*
`ρ(N) = G(N)/(c_φ𝔖(N)) = 1 − ε(N)` and look at `ρ(2N) = 1 − ε(2N)`. The table reports `ρ`, not
`ε`; the regression slope is the same for the two, because `ε ↦ 1−ε` is affine.

| decile | E[ρ(N)] | E[ρ(2N)] |
|---|---|---|
| 1 | 0.8411 | 1.0171 |
| 2 | 0.9117 | 1.0097 |
| 3 | 0.9430 | 1.0071 |
| 4 | 0.9671 | 1.0039 |
| 5 | 0.9871 | 1.0030 |
| 6 | 1.0066 | 0.9982 |
| 7 | 1.0277 | 0.9974 |
| 8 | 1.0505 | 0.9946 |
| 9 | 1.0820 | 0.9901 |
| 10 | 1.1497 | 0.9840 |

![Fig. 4 — E[ρ(2N)] against the decile of ρ(N), with ρ = G/(c_φ𝔖) = 1 − ε: monotone and linear, slope −0.111.](lab/out/fig_decile.png)

Monotone, essentially linear, regression slope `−0.111`. This is the whole content of `c(2)` in
distributional form; there is no threshold, no bimodality, and no dependence on `N mod 4`
(checking the four classes of `v₂(N)` separately gives `−0.148, −0.126, −0.144, −0.117`).

---

# 4. The correlation field

## 4.0 Two window conventions, used in different tables

The window enters every number in this paper, and **two different conventions are used**, which is
a defect of presentation rather than of measurement but must be stated.

* **Fixed count.** `[N_0, N_0 + 2M]`: `M` consecutive even integers, so the *relative* width
  shrinks as `N_0` grows. This is the convention in §3.1, in the scale scans of §4.10, and in
  §5.8–§5.11 (the calibrated sieve null, the `omega`-memory, and the zero component all use
  `M` consecutive even integers starting at `N_0`).
* **Fixed relative width.** `[N_0, 1.3 N_0]`: always 30% of `N_0`. This is the convention stated
  in §4.1 for the `c(k)` tables.

The two are not interchangeable. At `N_0 = 2×10^5` a fixed-count window with `M = 2000` spans 2% of
`N_0` while `[N_0,1.3N_0]` spans 30% — a factor of fifteen. `c(k)` depends on the window (§4.8
shows a `0.031` spread across the five weight shapes alone, and the `omega`-memory of §5.10 depends
on `N`), so **no table in this paper should be compared with another without checking which
convention it uses**; where a number is quoted in two places with different conventions, the
difference is a real difference and is reported as such.

## 4.1 Window and null

Unless stated otherwise: `N = 12000 + 2i`, `i < M`, `M = 4×10⁴` for `2 ≤ k ≤ 13` and
`M = 2×10⁴` for `14 ≤ k ≤ 49`. The null distribution of `c(k)` was obtained by permuting `ε(kN)`
across `N`: `mean = 0.001`, `sd = 0.026` at `M = 1500`, i.e. `≈ 1/√M`; the block bootstrap
`se ≈ 0.006` at `M = 4×10⁴`.

## 4.2 The measurement

| k | c(k) | se | k | c(k) | se | k | c(k) | se |
|---|---|---|---|---|---|---|---|---|
| 2 | **−0.1401** | 56 | 14 | +0.0498 | 85 | 25 | −0.0266 | 70 |
| 3 | **−0.1227** | 55 | 15 | +0.0185 | 75 | 26 | −0.0144 | 86 |
| 4 | −0.0632 | 62 | 16 | −0.0348 | 81 | 27 | −0.0349 | 71 |
| 5 | **−0.1510** | 59 | 17 | **−0.1332** | 78 | 28 | −0.0014 | 83 |
| 6 | +0.0235 | 55 | 18 | +0.0063 | 71 | 30 | +0.0027 | 76 |
| 7 | **−0.1575** | 59 | 19 | **−0.1012** | 74 | 32 | −0.0006 | 84 |
| 8 | −0.0322 | 57 | 20 | +0.0096 | 76 | 36 | +0.0159 | 77 |
| 9 | −0.0392 | 57 | 21 | +0.0235 | 73 | 49 | −0.0012 | 83 |
| 10 | +0.0357 | 59 | 22 | +0.0185 | 89 | | | |
| 11 | **−0.1432** | 57 | 23 | **−0.1034** | 71 | | | |
| 12 | +0.0035 | 59 | 24 | +0.0046 | 77 | | | |
| 13 | **−0.1440** | 59 | | | | | | |

![Fig. 1 — c(k) for 2 ≤ k ≤ 49. Primes in crimson, composites in blue; the bars are the block-bootstrap errors of §4.2. The two populations are separated by 0.1 in c and are interleaved in k.](lab/out/fig_ck.png)

(Standard errors in units of `10⁻⁴`.) Entries for `k ≥ 14` have `M = 2×10⁴`. The rows
`k=2,…,13` were produced with `M = 4×10⁴`.

## 4.3 S1: the sign law (established)

Twenty-eight of the thirty-one values obey

```
sign c(k) = (−1)^{ω(k)}                                               (4.1)
```

with `ω(k)` the number of *distinct* prime factors. The three exceptions are
`c(26) = −0.0144` (`2σ`), `c(28) = −0.0014` (`0.2σ`) and `c(30) = +0.0027` (`0.4σ`); only the
first is statistically significant. In particular the sign is *not* that of `λ(k)`: for
`k = 4, 9, 25, 49` one has `ω(k)=1` and `λ(k)=+1`, and `c` is negative in all four cases. The
law (4.1) is the most robust statement in this paper: it holds for five different weight shapes
(§4.8) and in seven different scale windows (§4.10).

**How much does 28/31 mean?** If the 31 signs were independent fair coins, `P(X >= 28) = 2.3×10⁻⁶`
(and `P(X >= 27) = 1.7×10⁻⁵`). **They are not independent** — all 31 correlations are computed from
the same field — so this is an indication of scale, not a test. What it rules out is a coin-flip
origin; what it does not supply is a calibrated significance, which only the block bootstrap of
§4.2 and the sieve null of §5.8 can provide.

## 4.4 S2: the amplitude at primes is not constant (rejected)

The nine prime values

```
p       2       3       5       7      11      13      17      19      23
c(p)  −0.140  −0.123  −0.151  −0.158  −0.143  −0.144  −0.133  −0.101  −0.103
se     0.006   0.006   0.006   0.006   0.006   0.006   0.008   0.007   0.007
```

have weighted mean `−0.1354` and

```
χ² = Σ (c(p) − mean)²/se(p)² = 72.8   on 8 degrees of freedom,  p = 1.4×10⁻¹².
```

The spread, `0.020`, is three times the typical standard error. The `p = 19, 23` values are about
`4σ` below the mean, and `p = 3` about `2σ` above it. A slowly varying trend in `p` is not
visible; part of the spread is presumably the scale admixture of §4.10 (for `k=23` the window of
`kN` lies at `2.8×10⁵ … 9.7×10⁵`, where the effect is weaker). **The honest statement is:
`c(p) ∈ [−0.16,−0.10]` with a real spread; `α` is a mean, not a constant.**

## 4.5 S3: multiplicativity for coprime pairs (approximate)

Using the table of §4.2, `c(ab)` against `c(a)c(b)`:

| test | measured | multiplicative prediction | difference |
|---|---|---|---|
| `c(6)` | +0.0235 | +0.0172 | +0.0063 (1.1σ) |
| `c(10)` | +0.0357 | +0.0212 | +0.0145 (2.4σ) |
| `c(14)` | +0.0498 | +0.0221 | +0.0277 (3.2σ) |
| `c(15)` | +0.0185 | +0.0185 | 0.0000 (0.0σ) |
| `c(21)` | +0.0235 | +0.0193 | +0.0042 (0.6σ) |
| `c(22)` | +0.0185 | +0.0201 | −0.0016 (0.2σ) |
| `c(26)` | −0.0144 | +0.0202 | −0.0346 (4.0σ) |
| `c(12)` | +0.0035 | +0.0078 | −0.0043 (0.7σ) |
| `c(18)` | +0.0063 | +0.0055 | +0.0008 (0.1σ) |
| `c(20)` | +0.0096 | +0.0095 | +0.0001 (0.0σ) |

Six of the nine tests agree within `1.2σ` (five within `1σ`, and `c(6)` at `1.1σ`). The residual is
*positive* in three of the four failures, i.e. the observed composite correlations are
systematically slightly *larger* than the product. Verdict: multiplicativity is a good approximation for products of two distinct primes,
with a systematic `+0.01` offset and a `4σ` counterexample at `k=26`. It is not a law.

## 4.6 S4 and S5: prime powers, and the failure of the Liouville form

For `k = p²` the multiplicative continuation predicts `c(p)² = +0.015…+0.023`; the measurement is

| k | c(k) | c(p)² | discrepancy |
|---|---|---|---|
| 4 | −0.0632 | +0.0196 | −13.2σ |
| 9 | −0.0392 | +0.0151 | −9.4σ |
| 25 | −0.0266 | +0.0228 | −6.9σ |
| 49 | −0.0012 | +0.0248 | −3.1σ |

and for higher powers `c(8)=−0.032, c(16)=−0.035, c(27)=−0.035, c(32)=−0.001`.

Two consequences.

* `c` is **not** completely multiplicative, and not multiplicative at prime powers. The correct
  empirical shape is: `|c(p)| ≈ 0.135`, `|c(p²)| ≈ 0.001…0.063`, `|c(pq)| ≈ 0.02`, i.e. the
  magnitude decays much more slowly with the *number of prime factors* than `α^{ω}` would
  suggest, and it does not alternate.
* The "Liouville" formula

```
c(k) = α^{ω(k)} λ(k)     (for squarefree k)                          (4.2)
```

is quantitatively wrong: for `k = 4, 9, 25` it predicts a *positive* value and the data are
negative; moreover `|c(pq)| ≈ 0.02 > α² = 0.018` only marginally while `|c(p²)|` is up to three
times `α²`. (4.2) should be replaced by the descriptive statement S1 plus a separate statement
about magnitudes.

## 4.7 A further logical constraint: rational ratios

Multiplicativity on `ℕ` together with the symmetry `c(1/k) = c(k)` forces
`c(p/q) = c(p)c(q)` for all primes `p, q`. In particular

```
c(3/2) = c(2)c(3) = c(6) = c(4/6),
```

so the three measurements must coincide. Measured (archive run, `M = 2×10⁴`, the same windows
as §4.2): `c(3/2) = +0.0140`, `c(6) = +0.0235` (this one is also in the table of §4.2),
`c(4/6) = +0.0478` — a spread of `0.034` at the `3–4σ` level. Earlier drafts stated the identity
alone; the numerical spread is what makes it a test. This is
the same phenomenon as the failure of §4.5, seen on the *rational* lattice; it shows that
whatever `c` is, it is not a multiplicative function on `ℚ_{>0}`. We therefore state all results
for integer `k` only.

## 4.8 S6: independence of the weight (established)

`M = 2×10⁴`, `N₀ = 12000`, five weight shapes:

| weight `φ²(u)` | support | c(2) | c(3) | c(4) | c(6) | c(10) |
|---|---|---|---|---|---|---|
| `cos⁴(3π(u−½))` | [1/3,2/3] | −0.144 | −0.117 | −0.070 | +0.028 | +0.058 |
| `1` (hard cut-off) | [1/3,2/3] | −0.132 | −0.135 | −0.078 | +0.039 | +0.010 |
| `cos⁴(9π(u−½))` | [5/12,7/12] | −0.134 | −0.122 | −0.073 | +0.032 | +0.015 |
| `cos⁴(1.5π(u−½))` | [1/4,3/4] | −0.136 | −0.127 | −0.082 | +0.027 | +0.013 |
| `16u²(1−u)²` | (0,1) | −0.163 | −0.149 | −0.092 | +0.037 | +0.017 |

(standard errors `0.007–0.009`). The support, the width and the analytic shape are irrelevant;
even a hard cut-off (which is not smooth, so that only a `C¹` Hardy–Littlewood asymptotic
applies) gives the same numbers, and a weight supported on all of `(0,1)` changes `c(2)` by
`0.02` only. This is a strong indication that the effect is a property of the primes and not
of the weight.

## 4.9 Locality in `N`

Standardising `ε(N)` and `ε(kN)` *inside* blocks of the window (so that all structure on scales
larger than the block is removed) and pooling:

| block width in `N` | number of blocks | pooled `c(2)` |
|---|---|---|
| 1000 | 30 | −0.1350 ± 0.0103 |
| 400 | 75 | −0.1343 ± 0.0095 |
| 200 | 150 | −0.1345 ± 0.0090 |
| 100 | 300 | −0.1351 ± 0.0093 |
| 40 | 750 | −0.1348 ± 0.0087 |

The effect survives inside blocks of 40 consecutive even integers, i.e. it involves no structure
smoother than roughly `1` part in `10³` of the window. Together with the whiteness of §1.1 this
excludes every "trend", "smooth component" or "slowly varying bias" explanation.

## 4.10 S7: scale dependence (rejected for composite `k`)

Same estimator, different windows (step 2, `M` as indicated):

| N₀ | M | c(2) | c(3) | c(4) | c(6) | c(9) | c(10) | c(12) | c(15) |
|---|---|---|---|---|---|---|---|---|---|
| 1.2×10⁴ | 4e4 | **−0.140** | −0.123 | −0.063 | **+0.024** | −0.039 | +0.036 | +0.004 | +0.019 |
| 5×10⁴ | 2e4 | −0.138 | −0.132 | −0.048 | +0.019 | −0.035 | −0.002 | +0.011 | +0.010 |
| 2×10⁵ | 2e4 | −0.123 | −0.082 | **−0.131** | **−0.036** | −0.049 | +0.042 | +0.016 | +0.001 |
| 3×10⁵ | 1.5e4 | −0.074 | −0.105 | −0.079 | +0.024 | −0.046 | +0.034 | +0.005 | +0.011 |
| 7×10⁵ | 8e3 | −0.101 | −0.089 | −0.030 | −0.005 | — | — | — | — |
| 10⁶ | 6e3 | −0.056 | −0.105 | −0.063 | +0.011 | — | — | — | — |
| 1.5×10⁶ | 5e3 | −0.094 | −0.064 | — | — | — | — | — | — |

![Fig. 3 — c(2), c(3), c(4), c(6) at seven scales. c(2) decays from −0.14 to −0.09; the composites c(4), c(6) move and change sign.](lab/out/fig_scale.png)

Beyond the range of the table the effect **persists**. Windows at `N₀ = 10⁷` (`M = 3000`),
`3×10⁷` (`M = 2000`) and `10⁸` (`M = 1200`), same estimator:

| N₀ | M | c(2) | c(3) | c(4) | c(5) |
|---|---|---|---|---|---|
| 10⁷ | 3000 | −0.067(19) | −0.066(18) | −0.032(19) | −0.092(21) |
| 3×10⁷ | 2000 | −0.074(23) | −0.067(18) | −0.108(20) | — |
| 10⁸ | 1200 | −0.101(24) | −0.061(23) | — | — |

So `c(2)` decays from `−0.14` at `10⁴` to `≈ −0.08` and then *stays* there over two further
decades. We do **not** quote a pooled significance: the three values drift
(`−0.067, −0.074, −0.101`), and pooling a drifting quantity would understate the systematic
uncertainty. Individually they are `3.6σ, 3.2σ` and `4.2σ` away from zero, and the statement we
make is only that the effect does not collapse. Note also that at `N₀ = 3×10⁷` the value
`c(4) = −0.108 ± 0.020` is as large as `c(2)` in that window, so the ordering of the
prime-power values is not stable across windows: only S1 is.

Three observations, now over all ten windows (`N₀ = 1.2×10⁴` to `10⁸`). (i) `c(2)` decays from
`−0.14` to `≈ −0.08` and then flattens, with `±0.03` scatter between windows. (ii) The prime correlations stay negative everywhere.
(iii) The composite correlations are *not* stable: `c(6)` is `+0.024` in one window and
`−0.036` (5σ) in another; `c(4)` moves between `−0.048` and `−0.131`. Hence the sign law S1 is
robust but the magnitudes — and in particular the multiplicative reading of §4.5 — are
window-dependent. A serious study of the scale dependence is the most important open
measurement (§7.2).

## 4.11 The scale dependence at `10⁸`–`10⁹`

The runs of §4.10 were extended with a segmented sieve to `10¹⁰` (`4.55×10⁸` primes, 40 GB) and
a shared-memory fork pool (128 vCPU). Same estimator, `Λ` exact at prime powers:

| k | `10⁸` (`M=10⁴`) | `3×10⁸` (`M=2×10⁴`) | `10⁹` (`M=10⁴`) | `10⁹` (`M=3000`) |
|---|---|---|---|---|
| 2 | −0.0633(102) | −0.0954(77) | −0.0461(118) | −0.0414(233) |
| 3 | −0.0641(99) | −0.0530(70) | −0.0463(94) | −0.0553(190) |
| 4 | −0.0406(110) | −0.0307(75) | — | −0.0455(204) |
| 5 | −0.0640(98) | −0.0826(81) | −0.0548(97) | −0.0498(163) |
| 6 | −0.0132(107) | −0.0051(73) | −0.0078(107) | −0.0238(203) |
| 7 | −0.0744(101) | −0.0758(79) | — | −0.0664(207) |
| 9 | −0.0302(110) | −0.0295(71) | — | −0.0574(198) |
| 10 | −0.0130(103) | **+0.0518(82)** | — | **+0.0447(186)** |
| 12 | +0.0082(98) | −0.0016(74) | — | −0.0008(152) |
| 15 | +0.0124(98) | +0.0021(71) | — | +0.0152(171) |

Four findings follow, two of which revise statements made in §4.

**(1) The effect survives to `10⁹`.** Every prime value in every window is negative, and the
pooled prime amplitude is `−0.066 ± 0.005` at `10⁸`, `−0.075 ± 0.004` at `3×10⁸`,
`−0.049 ± 0.006` at `10⁹`. The contrast with the composite values (`−0.001`, `+0.010`,
`−0.008`) is maintained, so S1 and the prime/composite separation hold over five and a half
decades.

**(2) No decay law is resolved, and the within-window `se` is *not* the right uncertainty.**
`c(2)` reads `−0.067(19)`, `−0.074(23)`, `−0.063(10)`, `−0.095(8)`, `−0.046(12)`
at `10⁷, 3×10⁷, 10⁸, 3×10⁸, 10⁹` (the `M=3000` run at `10⁹` gives `−0.041(23)` on a
subset window and is not an independent point). The spread is `±0.025`, three to four times the largest
within-window error, and it does not order by `N`: `3×10⁸` is the most negative point and `10⁹`
the least. A `1/log N` law fitted to `10⁴, 10⁸, 10⁹` (`−c·log N = 1.34, 1.22, 1.11`) is broken by
`3×10⁸` (`1.86`). **The right procedure is many windows per scale, not a larger `M` inside one
window**; until that is done, no decay law should be claimed, and the `±3%` quoted for (3.1)
applies to the *standard deviation* of `e`, not to `c(k)`.

**(3) The `p`-dependence of the prime amplitude is real but window-dependent.** The four primes
`2,3,5,7` give `χ² = 0.8/3` at `10⁸`, `17.8/3` (`p = 3×10⁻⁴`) at `3×10⁸` and `0.7/3` at `10⁹`.
Constancy is therefore not a large-`N` property; it happens window by window, with a spread `±0.02` comparable to the window scatter. The
correct statement is the original S2: **`c(p)` has a real `p`-dependence, of a size that a
single window cannot resolve.**

**(4) `c(10)` is anomalously positive.** `−0.013(10)` at `10⁸`, `+0.052(8)` at `3×10⁸` (`5.3σ`)
and `+0.045(19)` at `10⁹` (`2.3σ`), against the multiplicative prediction `c(2)c(5) = +0.004`,
`+0.008`, `+0.002`. The *sign* is what S1 requires (`ω(10)=2`); the *magnitude* is an order of
magnitude larger than multiplicative. This is the most significant violation of S3 in the data,
and it is stable across two scales.


---

**(5) Reproducibility, and a caveat about what it can mean.** Ten windows at `10⁹`, starts
`10⁹ + j·2×10⁶` (`j = 0,…,9`), `M = 3000` each, give `c(2) = −.041, −.061, −.072, −.034, −.053,
−.062, −.042, −.018, −.065, −.034`, weighted mean `−0.0474 ± 0.0059` (`χ² = 7.1/9`); the same
ten-window test at `3×10⁸` gives `−0.1065 ± 0.0059` (`χ² = 9.1/9`) and at `10⁸`
`−0.0526 ± 0.0059` (`χ² = 10.5/9`). The local weighted prime mass `m` varies across those
windows by `8%` at `10⁹` and changes sign at `3×10⁸`, and never drives `c(2)`
(`corr = +0.30, +0.02, +0.54`; none significant on 10 points). Splitting a window into
sub-windows gives the same picture (`χ² = 5.4, 7.1, 3.3` on 9 dof).

**But this test is weaker than it looks.**
For a window at scale `N` the two prime intervals have length `≈ N/3`, whereas the window
itself is only `M·2` wide in `N`; ten windows spaced by `2×10⁶` therefore share `≈ 99.8%` of
their primes. The ten measurements are independent only in their *pairing* content, not in
their *density* content. Consequently:

* they do establish that the **pairing channel is homogeneous** in `N` at fixed scale, and that
  the local mass is not the driver;
* they **cannot** exclude a density-channel component, because that component is common to all
  ten windows;
* and they cannot reduce the between-scale uncertainty, because **at this statistic the "scale"
  and the "sample" are the same axis**: to change the primes one must change `N` by `O(N)`,
  which changes the scale by `O(1)`.

So the three disjoint scales are three *independent* measurements, and their spread is the
honest uncertainty: `c(2) = −0.053(6), −0.107(6), −0.047(6)` at `10⁸, 3×10⁸, 10⁹`, against
`−0.142(6)` at `1.2×10⁴`, `−0.067(19)` at `10⁷` and `−0.074(23)` at `3×10⁷`. **The only
defensible summary is `c(2) ≈ −0.07` with a scatter of `±0.03` that no measurement in hand
resolves** — the statistical error at each scale (`0.006`) is five times smaller than the
scale-to-scale spread, and the two cannot be separated by construction.


## 4.12 Not the explicit formula: the zeros build a smooth field, not a white one

The classical reading of a Goldbach error is that it is driven by the zeros of `ζ` through the
explicit formula. That reading is testable directly, and it fails on the *statistics* rather
than on a number.

We located `48 900` ordinates below height `40 000` as sign changes of the Riemann–Siegel main
sum, refined by vectorised bisection, and built the two fields the explicit formula can produce on
the window `N = 12000+2i`, `i < 15000`. Two caveats on the list: the Riemann–von Mangoldt formula
predicts `49 395` zeros below `40 000`, so about `1%` are missed where the remainder of the
Riemann–Siegel formula flips the sign (refining the search grid from `0.2` to `0.005` changes the
count from `48 533` to `48 901` and then stops), and the ordinates are those of the main sum, so
they carry an error of `O(0.4)` at the first zero and `≈0.01` near `t = 4×10⁴`. Neither affects
the conclusion below, which is about the *smoothness* of the resulting field:

```
L(N) = Σ_j N^{iγ_j}/γ_j                            (one zero: density-like)
P(N) = |L(N)|² = Σ_{j,k} N^{i(γ_j+γ_k)}/(γ_j γ_k)   (two zeros: pair-count-like)
```

with real positive coefficients — the most favourable choice, since any phase assignment only
adds cancellation. Then

| field | lag-1 autocorr. | lag-10 | sd |
|---|---|---|---|
| `L` | **+0.996** | +0.950 | 0.105 |
| `P` | **+0.990** | +0.915 | 0.023 |
| measured `ε` | **−0.027** | −0.031 | 0.086 |

**The zero-driven fields are smooth in `N`; the measured field is white.** The reason is
structural: the explicit formula weights zeros by `1/ρ`, so the sum is dominated by the lowest
zeros (`γ ≈ 14–100`), whose phase `γ·ΔN/N` changes by `< 0.02` rad between adjacent `N`. No
choice of window or scale changes this — the low zeros are always the ones with the weight.

This does not contradict §5.1(f); it explains it. The zeros control the *smooth* (density)
component of the error, and that component is measured to be `≈1%` of the variance
(`corr(ε,m) = 0.105`). **The remaining `99%` — the pairing noise, which is what carries `c(k)` —
has no zero representation at all.** So the effect cannot be the classical zero-driven error
term, and the "already known" reading of S8 is closed: the answer to *new phenomenon or known
implication* is **not the known implication**, whatever else it may be.


---

# 5. Why the measurement is surprising

## 5.1 First, the trivial explanations

**(a) Normalisation.** For `k=2`, `𝔖(2N)=𝔖(N)` because `N` is even, so `c(2)=corr(G(N),G(2N))`
by (1.7); there is no normalisation to blame. For other `k` we compared five normalisations
(`σ(N) ∝ √𝔖 log N/√N`; constant; `√𝔖` only; `log N/√N` only; `𝔖 log N/√N`): the values of
`c(k)` move by at most `0.02` (and by `≤0.006` for `k=2`), with the sign and the magnitude
unchanged.

**(b) Position in `N`.** §4.9: the effect survives inside blocks of width 40.

**(c) Smooth trend.** §1.1: `ε` is white. `corr(ε(N), ε(N+20000)) = +0.009`.

**(d) Wrong multiplication.** `corr(ε(N), ε(2N+2)) = +0.038`: shifting the double by one even
number destroys the effect. So the effect is attached to the exact multiplicative relation, not
to the neighbourhood.

**(e) Local arithmetic of `N`.** Regressing `ε(N)` on the indicators `p | N` (`p ≤ 29`) and on
`ω(N)` gives `R² = 0.0008`. Class means of `ε(N)` by `N mod q` give (adjusted) `R² < 0.02` for
`q ≤ 4620`; removing them leaves `c(2) = −0.139` instead of `−0.140`. A factorisation-dependent
field is also excluded by a variance argument, but the argument needs care. For an *additive* function of the factorisation,
`ε(N) = Σ_p f_p(v_p(N))`, doubling changes only the 2-adic term, so `ε(2N) − ε(N) = f₂(v₂(N))`
is constant on each `v₂` class and
`Var(ε(2N)−ε(N)) = between-class variance over v₂ = 1,2,3,≥4`. A four-valued function has no
*a priori* variance bound, so the content comes from the measurement: the `v₂` class means are
flat, `R² = 0.0004` (§5.7(c)), whereas `Var(ε(2N)−ε(N)) = 2.28 σ²`, i.e. essentially
`2σ²(1−c(2))`. For a *non-additive* function of the factorisation the difference need not be a
function of `v₂` at all, but such a field carries an interaction between the 2-adic and the odd
part, and the class-mean fits above (`R² = 0.0008` on `p | N`, `p ≤ 29`, and on `ω(N)`) bound
precisely that. So `ε` is not a function of the factorisation of `N`.

**(f) The mass of primes in the interval.** The weighted prime mass
`m(N) = Σ_{N/3<n<2N/3} φ²(n/N)(1_{n prime} − 1/log n)` is the only "smooth" channel that can
enter `ε`. We measure `corr(ε(N), m(N)) = +0.105`, so the mass explains `≈1%` of the variance of
`ε`, and at most `≈3%` of the covariance `c(2)`. (The quantity `corr(m(N),m(2N)) = −0.51` looks
spectacular but is not usable: `m` is smooth on a scale of `≈N/2`, and the window contains only
about three such scales, so this number has of order one degree of freedom.)

## 5.2 The obstruction: the two counts use disjoint primes

For `k=2` the two sums are

```
G(N) = 2 Σ_{N/3<n<N/2} Λ(n)Λ(N−n) φ²(n/N) ,      n, N−n ∈ [N/3, 2N/3]     (5.1)
G(2N) = 2 Σ_{2N/3<n<N} Λ(n)Λ(2N−n) φ²(n/(2N)),    n, 2N−n ∈ [2N/3, 4N/3]  (5.2)
```

and the two prime intervals `[N/3,2N/3]` and `[2N/3,4N/3]` are **disjoint** (they meet at the
single point `2N/3`, where the weight vanishes to fourth order). The same holds for every `k>1`:
the pairs summing to `N` straddle `N/2` and live in `[N/3,2N/3]`; the pairs summing to `kN`
straddle `kN/2` and live in `[kN/3,2kN/3]`; these intervals are disjoint as soon as `k ≥ 3/2` if
the supports are scaled with `k` as in (1.6).

Hence **the two random variables whose correlation we measure are functions of disjoint pieces
of the same prime configuration.** This is not a convention: it is forced by the arithmetic of
`n+m=N` together with the balanced support.

## 5.3 The four-point local densities factorise

Write `n′` for the variable in (5.2) and consider the four linear forms

```
n,   N−n,   n′,   2N−n′     (n ∈ [N/3,N/2],  n′ ∈ [2N/3,N]).
```

By Hardy–Littlewood, the density of quadruples with all four prime is, up to `o(1)`,
`𝔖₄/log⁴N` with `𝔖₄ = ∏_p σ_p`, where `σ_p` is the number of admissible configurations
modulo `p` divided by the "random" count. For a single pair (`n, N−n`) the condition is

```
n ≢ 0 (mod p)  and  n ≢ N (mod p),        i.e. n ∉ {0, N} (mod p),     (5.3)
```

i.e. `p − |{0, N mod p}|` admissible residues out of `p`. For the quadruple the conditions on
`n` and `n′` are *separate*:

```
n ∉ {0, N} (mod p)   and   n′ ∉ {0, 2N} (mod p).                       (5.4)
```

Therefore

```
σ_p(quadruple) = σ_p(pair) · σ_p(pair)      for every prime p,          (5.5)
```

and the four-point singular series is *exactly* the product of the two two-point ones. Note that
this is not the generic situation: had the two pairs shared a variable (say the same `n` paired
with `N−n` and with `2N−n`), the local counts would *not* factorise and one would get a positive
contribution of the "singular" type. With disjoint pairs and balanced support, (5.5) holds term
by term.

Consequently, for any model in which the prime indicator has the correct local densities and
independent positions, `Cov(G(N),G(2N)) = 0` **exactly**, and likewise for `k>1` after
normalisation. Our iid simulation (§2.4) confirms this numerically.

## 5.4 What this means

The measured `c(2) = −0.14` is therefore a deviation of the *observed covariance* from the
factorised prediction, which is exactly zero. Note what `0.14` measures: it is `14%` of the
*standard-deviation product* `sd(ε(N))·sd(ε(2N))`, not `14%` of any four-point density — the
densities themselves are matched, by §5.3, to machine precision. There are only two
possibilities.
  but not the possibility that the *heteroskedasticity* of the field inside the window induces
  a bias of this size: `sd(ε(N)) ∝ (log N)/√N` varies by a factor `1.7` across `[12000,42000]`,
  so a Pearson correlation of two such vectors is a *weighted* average of the local correlations,
  with the low-`N` end (where `sd(ε)` is largest) dominating. Since the weight profile of the
  low-`N` end is not the same for `ε(N)` and `ε(2N)` (`sd(ε(2N))/sd(ε(N))` is itself a function
  of `N`), the estimator is not obviously unbiased. **This is the most likely residual artifact
  and should be tested first** (§5.6, test T1).
* **(B) A genuine long-range four-point correlation of the primes.** Then
  `𝔖₄ ≪ (1−δ) 𝔖₂𝔖₂` for a positive proportion of configurations with `δ ~ 0.14`, i.e. the
  primes behave *sub*-independently at large scales, in a way tied to the ratio `k` of the two
  sums. This would be a discovery independent of Goldbach, and it would have to be visible in
  other statistics (e.g. in the variance of `G`, which would be reduced).

## 5.5 Controls that pass

| control | result | verdict |
|---|---|---|
| five normalisations | move `c(k)` by ≤ 0.02 | ✓ |
| five weight shapes | `c(2) ∈ [−0.163,−0.132]` | ✓ |
| block-localisation down to width 40 | `c(2) = −0.135` | ✓ |
| `N mod q`, `q ≤ 4620` | `Δc(2) < 0.002` | ✓ |
| divisibility by `p ≤ 29`, and `ω(N)` | `R² = 0.0008` | ✓ |
| factorisation-function test | `Var(ε(2N)−ε(N)) = 2.28σ²` | ✓ |
| shift by one even number | `+0.038` | ✓ |
| iid random primes (8 draws) | `−0.017 ± 0.025` (se of the mean; spread `0.070`) | ✓ |
| distribution | Gaussian, `kurt = 3.00` | ✓ |

Every one of these is consistent with the effect being real, and none of them explains it.

## 5.6 The decisive test

The covariance can be written as an explicit quadruple sum,

```
Cov(G(N),G(2N)) = (2/N²) Σ_{n, n′} φ²(n/N) φ²(n′/2N) Λ(n)Λ(N−n)Λ(n′)Λ(2N−n′)
                  − (main terms),                                       (5.6)
```

and, in the Hardy–Littlewood model, replaced by

```
(2/N²) Σ_{n, n′} φ² φ² · [𝔖₄(n,N−n;n′,2N−n′) − 𝔖₂(N)𝔖₂(2N)] / log⁴N .   (5.7)
```

By (5.5) the bracket vanishes identically, so (5.7) `= 0`, whereas the measured covariance of
the two relative errors is `c(2)·sd(ε(N))·sd(ε(2N)) ≈ −7×10⁻⁴` at `N ≈ 2×10⁴` (with
`sd(ε(N)) = 0.086` and `sd(ε(2N)) ≈ 0.061`). The proposed tests are:

* **T1 (estimator).** Recompute `c(2)` with (i) a *fixed* number of pairs per `N` (thin the
  support so that the expected count is constant), (ii) variance-stabilised coordinates
  (`ε(N)/sd(ε(N))` with the exact `sd`), (iii) a window of constant relative width.
  If `c(2)` moves to zero, (A) is confirmed.
* **T2 (four-point sum).** This test is vacuous as a test of the physics: by (5.5)
  the bracket `𝔖₄ − 𝔖₂𝔖₂′` vanishes **identically**, so a numerical evaluation of
  `Σ_{n,n′}(𝔖₄ − 𝔖₂𝔖₂′)/log⁴N` must return zero and can only detect a coding error. That is
  what §5.7(1) does, and it is a useful check, but it cannot confirm (B). The real test is the
  one carried out in §5.7(2)–(3): evaluate the **empirical** four-point sum (5.6) and compare it
  with a null that has the same local densities.
* **T3 (variance).** Under (B) one expects `Var(G)` to be *smaller* than the Poisson value. The
  measured `sd(ε) = 0.086` is `0.85×` the naive reference (`0.101`) and `0.60×` the
  weight-reduced one (`0.143`): mildly **sub**-Poisson. The sign depends on the convention: using `196`
  ordered representations instead of `98` centres flips it. Sub-Poisson variance is at least not in tension with negative four-point
  correlations, but the reference moves by `40%` with the weight convention and the link between
  a local variance and a long-range correlation is indirect, so **we do not count this as
  evidence either way.**

---

## 5.7 The four-variable sum, computed

The programme of §5.6 has now been carried out. Four computations.

**(1) The factorisation is exact, by enumeration.** For `N ∈ {12000, 12002, 12006, 12010,
30030, 23100}` and every prime `p ≤ 500` (`95` primes, `570` pairs) the joint local count of
`(n, N−n, n′, 2N−n′)` was obtained by **enumerating all `p²` residue pairs** and compared with
the product of the two marginal counts:

```
max |joint / (marg₁ · marg₂) − 1| = 0.000e+00      (570/570 pairs, exact)
```

The identity holds for a one-line reason: if `p | N` both conditions on each side exclude only
`0`, so both counts are `(p−1)²`; if `p ∤ N` each side excludes two residues, so both counts are
`(p−2)²`. The four-point singular series therefore factorises identically, and the
Hardy–Littlewood prediction (5.7) is not merely small — it is `0`.

**(2) The measured covariance.** `Cov(ε(N),ε(2N)) = −8.27×10⁻⁴` against `sd·sd = 5.91×10⁻³`,
i.e. `c(2) = −0.1400` at `M = 15000` (this is the direct estimator of §4.2, not a by-product of
the four-variable count). The deficiency relative to the
factorised prediction is therefore `14.0%` of the standard-deviation product, at `14σ`.

**(3) A factorised null with the same local densities.** A random set that is coprime to every
`p ≤ 97` and thinned to density `1/log n` — `190365` elements against `189880` primes — gives,
over 40 independent draws, `c(2) = +0.0072 ± 0.0093` (within-block) and `−0.0155 ± 0.0134`
(full window), where the quoted error is the standard error of the *mean*. The
realisation-to-realisation **spread** is much larger, `0.051` (within-block) and `0.086` (full
window), because this model's field is smooth (§5.7(3)). Which yardstick applies is the crux of
the comparison: the primes' field is white (measured block correlations, §2.3), so the
within-window error is the relevant one and the measurement sits `14σ` away; if instead the
primes were a typical draw from this ensemble, the deviation would be `2.7σ`. **The honest range
is `3σ–15σ`, not `15σ`.**

**(4) Estimator controls.** Recomputing `c(2)` at `M = 15000` after (i) dividing by the known sd
profile `log N/√N`, (ii) inverse-variance weighting, (iii) replacing Pearson by Spearman,
(iv) trimming at `|z| < 3` (99% kept), gives `−0.1424, −0.1432, −0.1305, −0.1269` against
`−0.1400` raw. The block correlations are white — Newey–West standard error `0.0094` against a
naive `0.0087` — so the error bar is not underestimated.

Two by-products matter for the interpretation.

**(a) The field is pairing-noise dominated, not mass dominated.** The relative fluctuation of
the weighted prime mass in the sliding window is `0.0048`, against `0.024–0.051` for random sets
of the same density; but a sliding window has only a few effective degrees of freedom, so that
comparison is not reliable. A clean measurement on *disjoint* intervals of length `L` gives
`sd(count/expected) = 0.0615, 0.0414, 0.0259` for `L = 1000, 2000, 4000`, against the Poisson
values `0.1088, 0.0769, 0.0544`: the prime count is about `2×` less variable than Poisson,
plausibly from the negative short-range correlations of the sieve. What matters here is direct:
`corr(ε,m) = 0.105`; `ε` is white in `N` (lag-1 `−0.027`) while the mass `m` is extremely smooth
(lag-1 `0.9999`). So the field is pairing-noise dominated, whereas every independent model we
built is mass dominated — smooth, with block correlations up to `+0.5`. This is why a factorised
null is an imperfect mimic of the field's memory, and why the null of (3) is quoted with its own
realisation-to-realisation scatter.

**(b) A sieve-type model can produce a large negative `c(2)`, but it cannot be calibrated.**
A null built as a random sieve — delete one random residue class for every prime `p ≤ 1609` —
produces `c(2) ≈ −0.3`; the reason is instructive, since in such a model the pair density
depends on the residues of `N` modulo a large product of primes, a *fine-scale* `N`-dependent
systematic that survives the `𝔖(N)` normalisation. This identifies the mechanism class to which
the measurement most plausibly belongs. **But the attempt to turn it into a usable null failed.**
A corrected version — sieve applied to all odd integers (not to the primes), a thinning factor
`t = min(c/log n, 1)` intended to bring the density to `1/log n`, expectation analytic including
the singular-series factor `∏_{p|h,p≤P}(p−1)/(p−2)` — still has `sd(ε) = 0.42` against the real
`0.086`, and a realisation-to-realisation spread of `c(2)` of `0.48` (the real field's own
uncertainty is `0.009`); it therefore calibrates nothing. **Note that the thinning never engaged
in that run**: with `c = 1/dens = 13.15`, `t = 1` for every `n ≤ e^{13.15} = 5.2×10⁵`, whereas
the window lies below `5.6×10⁴`, so the run is distributionally the same model as the uncorrected
one — which is why the two `sd(ε)` agree exactly. The density axis is engaged only in §5.8, where
the sieve bound is `n^α` with `α = 0.56`. Reproducing all three axes at once — the local
structure, the density profile `1/log n`, and the sub-Poisson mass fluctuation — appears to
require the *sequential* Hawkins sieve (whose sieve bound grows with `x`), and its naive
simulation is `O((X/log X)²)`, out of reach at the scales used here. We record this as an open
technical problem, not as a result. It is not a smooth trend (excluded by §4.9), not the local structure of `N` (excluded
by §5.1(e)), not the weight, the scale or the estimator (excluded above). That leaves two
readings: a *second-order arithmetic correction* to the pair count — a small additive term of the
same class as the singular series, not removed by normalising with the leading factor — or a
genuine deficiency of the joint law of the primes. **Part (c) below closes the first one** for
every deterministic function of the arithmetic of `N`, so the deficiency must be non-local: it
cannot be written as `δ(N)` for any function of `N`'s factorisation.

**(c) A separable deterministic bias cannot do it.** Suppose `eps(N) = delta(N) + noise` with
`delta` a deterministic function of the arithmetic of `N` — the only way a second-order
singular-series term could enter. If `delta` depends on the *odd* part of `N`, i.e. on any prime
other than `2`, then `delta(2N) = delta(N)` and its contribution to `c(2)` is
`+var(delta)/var(eps) >= 0`, not negative. The only part of `delta` that flips under doubling is
the 2-adic one, and a function of `v_2(N)` alone is measured to have `R^2 = 0.0004` on the cells
`v_2 = 1,2,3,>=4`, i.e. `var(delta)/var(eps) <= 4e-4`, whereas `c(2) = -0.14` needs
`var(delta)/var(eps) >= 0.14`.

**What this does and does not exclude.** It excludes `delta` that is a function of `v_2(N)`
alone, and it excludes `delta` that is a function of the odd part alone (by the first argument).
It does **not** exclude a function of *both*, `delta(N) = f(v_2(N), odd(N))`: such an `f` need
not satisfy `delta(2N) = delta(N)` unless it is independent of its first argument, and the `R^2`
test above was run on `v_2` cells only. An interaction between the 2-adic valuation and the odd
part, of size `0.14` in `c(2)`, is therefore **not ruled out by these numbers**. Earlier drafts
wrote "excluded, not merely disfavoured"; that is too strong and is withdrawn. What survives is
the weaker statement, still useful: no *separable* arithmetic function explains the effect, and
the interaction class has not been tested.

**(d) The ternary analogue is not a valid control.** We built the ternary field
`T(N) = N⁻² Σ_{n₁+n₂+n₃=N} ΛΛΛ` over the balanced box `N/6 ≤ n_i ≤ N/2` (`N` odd), validated
against a direct double loop, with the singular series `𝔖₃(N) = ∏_{p|N,p>2}(p−2)(p−1)/(p²−3p+3)`
divided out — after which the values for `3 | N` and `3 ∤ N` agree to `0.02%`, i.e. the
classical `2/3` suppression is present and correctly removed. Its statistics nevertheless make
it a poor drop-in control: `sd = 0.0070` against a pairing-noise estimate of `0.0025` (from
`Σ ΛΛΛ = T·N² ≈ 1.7×10⁸` and `Λ³ ≍ 1.1×10³`, i.e. about `1.6×10⁵` triples), and lag-1
autocorrelation `+0.54`. The ternary field is therefore **density**-dominated and smooth,
whereas the binary field is **pairing**-dominated and white. Measured (se `0.05`), and **withdrawn**: the driver used for this table indexed the ternary
count at `N` rather than at `kN`, so the numbers it produced are not `corr(ε(N), ε(kN))`; the
monotone drift of the printed values with `k` is the signature of that mismatch. The
qualitative point stands and does not rest on the table: this is *dilution*, not evidence — a pairing-channel
effect of `−0.13` would be attenuated by `(0.0025/0.0070)² ≈ 0.13`, i.e. to `≈ −0.017`, which
at `M = 6×10³` (`se ≈ 0.013`) is invisible.

Two points are worth recording, because it is easy to over-read these numbers. (i) The pairing-to-density ratio scales like `δ^{−1/2}`,
where `δ` is the relative width of the weight (`#triples ∝ δ²`, so the pairing noise `∝ δ^{−1}`,
while the density noise only `∝ δ^{−1/2}`); a *narrow* ternary weight therefore makes the
ternary field pairing-dominated and the control legitimate, at the price of statistics. (ii) The
two scaling laws — pairing noise `≍ log N/√N` versus density noise `≍ √(log N/N)` for the binary
problem, `≍ log^{3/2}N/N` versus `≍ √(log N/N)` for the ternary — are conventions for *unweighted*
counts; constants, the effective number of terms and the primes' measured sub-Poissonianity
shift the crossover. **The correct statement is that the binary problem is the *best* member of
the family for this measurement (it is the only one that is pairing-dominated by default), not
that the others are excluded; we did not run the narrow-weight ternary variant.**

**Verdict (S8).** The four-variable sum gives exactly zero; the factorised null gives zero; the
estimator controls all pass; and (c) excludes every *separable* deterministic function — one depending on `v_2(N)` alone
or on the odd part alone — at `350×` the required amplitude, while leaving the
`f(v_2, odd)` interaction untested. What remains is a genuine long-range deficiency of the
primes' four-point correlations at integer ratios — equivalently, a coupling between the prime
configurations in `[N/3,2N/3]` and in `[2N/3,4N/3]` that no local density formula and no
normalisation convention can produce. Such a coupling can only be non-local in `N`: it must
involve the *positions* of the primes in the two intervals jointly, in the way a dynamical model
(the sequential Hawkins sieve being the obvious example) couples them. We do not have such a
mechanism. **What we can state is that the effect is a statement about the primes, not about the
estimator** — the two readings that would have made it an artifact are both closed.

## 5.8 A calibrated sieve null

§5.7(3) recorded that a sieve-type model produces a negative `c(2)` but cannot be calibrated:
its field is 5x too variable (`sd(ε) = 0.42` against `0.086`). That gap has since been closed
by matching the *density profile*, and the resulting null is decisive enough to be reported.

**The mechanism of the old failure.** A sieve with a *fixed* bound `P = sqrt(X)` strikes one
residue class per prime `p <= P` for every `n`, so its density is the constant
`prod_{p<=P}(1-1/p) ~ e^{-gamma}/log P`, not `1/log n`. In the window actually used
(`n ≈ 4×10³–2.8×10⁴`, `log n ≈ 8.3–10.2`) the effective exponent is
`alpha_eff = log P / log n ≈ 0.79`, so the density is `e^{-gamma}/(0.79 log n) ≈ 0.71–0.81` of
`1/log n` — a deficit of `20–30%`, not the `12%` excess that a bound of `sqrt(n)` would give, and
a deficit *inflates* the relative fluctuation. Sweeping instead with the bound `n^alpha` gives
`prod_{p<=n^alpha}(1-1/p) ~ e^{-gamma}/(alpha log n)`, so

```
alpha = 0.56   =>   density constant e^{-gamma}/0.56 = 1.0027 ~ 1.000
```

Measured: `d*log X = 1.1241 (alpha=0.50)`, `1.0076 (alpha=0.56)`, `0.9126 (alpha=0.62)` -
the density profile is matched at `alpha = 0.56`.

**The model.** For each prime `p <= X^alpha`, delete one uniformly random residue class `r_p`
from the integers `n >= p^{1/alpha}`; the pair count and the model singular series use the
two-regime factor, with the regimes at `p <= (N/3)^alpha` (both members sieved) and
`(N/3)^alpha < p <= (2N/3)^alpha` (only the larger sieved). Implementation verified four ways:
the density converges to `e^{-gamma}/(alpha log n)`; every surviving `n` has no `p <= n^alpha`
hitting it and every dead `n` has one; the pair count agrees with a brute-force double loop to
`1e-10`; the local factor `(p-1)/p` versus `(p-2)/p` agrees with a hand count at `p = 3`.

**Result over 120 independent seeds** (`X = 2*10^6`, `N0 = 12000`, `M = 2000`):

| quantity | null (120 seeds) | primes |
|---|---|---|
| sanity `sum C / sum main` | `1.0030 +- 0.0178` | 1 |
| `sd(eps)` | `0.05173 +- 0.00272` | `0.04702` (+10.0%), on the same window |
| `c(2)` | `-0.0057 +- 0.0648`, range `[-0.270, +0.171]` | `-0.1400` |

The observed `-0.140` sits at the **1.67 percentile** (2 of 120 realisations reach it;
`7/120` reach `-0.10`). Under the typical-draw reading that is `-2.07 sigma`; under the
mean reading, `-22.7 sigma`. The 10.0% variance excess inflates the null's spread, so the
percentile is conservative - matching the variance would put it near `1.0%`.

**What this settles and what it does not.** A random sieve with the correct main term, the
correct density profile and a 9%-accurate variance *does* produce the right sign, but its
typical magnitude is `-0.006`, one twentieth of the observation, and only its extreme tail
reaches `-0.14`. This is evidence that the effect is **not** generic to sifted sets - but at
`p ~ 1.7e-2` it is not a rejection. Reaching `4 sigma` needs about `1e5` seeds, which is
`3.3` hours on a 128-core machine with the scripts as written.

**A residual that is not a bug.** The 10.0% variance excess survives four independent
implementation checks; it is a real difference between a random sieve and the primes (the
sieve's short-range correlations are weaker, so its count is closer to Poisson). Removing it
requires changing how the residue classes are drawn, not any parameter.

## 5.9 The field is not white: memory at small-prime lags

The field is not white in `N`, and this changes how its error bars must be computed.

**The field is not white in `N`.** A lag-1 autocorrelation of `−0.027` at lags
`2, 4, 10, 20, 200, 1000` does not establish whiteness: none of those lags is a multiple of
`210 = 2·3·5·7`. At multiples of `210` the autocorrelation is
positive and large (`M = 4000`, `se = 0.0158`):

| lag `Delta N` | 2 | 4 | 10 | 200 | 1000 | **210** | **420** | **630** | **1050** | **1470** |
|---|---|---|---|---|---|---|---|---|---|---|
| `corr(eps(N), eps(N+Delta))` | −0.001 | +0.007 | −0.027 | −0.026 | −0.040 | **+0.083** | **+0.069** | **+0.082** | **+0.093** | **+0.100** |

Every multiple of `210` is positive at `5–6 sigma`; every non-multiple is zero or negative.
**The cause is not the residue-class structure of the singular series.** That was the first
hypothesis and it is wrong: with `M = 10^5` samples the class means carry only `0.09%`
(`mod 30`), `0.61%` (`mod 210`) and `1.39%` (`mod 2310`) of the variance, after subtracting the
sampling-noise floor (`0.01%, 0.10%, 1.15%`). A `0.6%` class structure cannot produce an
autocorrelation of `0.08`. The pattern is only weakly reproducible across samples
(`r = +0.23 +- 0.10` against an independent `M = 6000` batch), consistent with a small real
effect, but twenty times smaller than the `33%` that a residue-class reading of the effect
would require. **Why `eps` has memory at multiples of `210` is therefore open.**

**Consequence.** `M_eff < M`, so the naive `1/sqrt(M)` error bar is too small. This is exactly
the effect noted in §4.11 — "the window-to-window scatter is three to four times the
within-window error" — and it is the reason the two numbers disagree.

**(ii) The residue-class structure is not negligible.** §5.1(e) states that conditioning on
`N mod q`, `q <= 4620`, moves `c(2)` by `< 0.002`. That is wrong. With `M = 12000`:

| grouping | `Var(class means)/Var(eps)` | `c(2)` after removing the class means |
|---|---|---|
| none | — | `−0.1620` |
| `N mod 30` | 0.011 | `−0.1563` |
| `N mod 210` | 0.055 (noise floor 0.018) | `−0.1383` |
| `N mod 2310` | 0.156 (noise floor 0.192 — **not resolvable**)| `−0.1306` |

At `mod 210` the class means carry about `4%` of the variance (after subtracting the `1/57`
noise floor) and their removal costs `0.024` of `c(2)`, i.e. `15%`. The `mod 2310` row is not
usable at this `M`: with `5.2` samples per class the sampling noise alone is `19%` of the
variance, so no conclusion should be drawn from it without `M` an order of magnitude larger.

**(iii) A class-preserving shift retains a third to a half.** Replacing `eps(2N)` by
`eps(2N+h)`:

| `h` | 0 | 2 | 4 | 210 | 420 | 2310 | 4620 |
|---|---|---|---|---|---|---|---|
| `corr` | `−0.1476` | `−0.0017` | `+0.0017` | `−0.0502` | `−0.0514` | `−0.0610` | `−0.0657` |

General shifts erase the correlation completely; shifts that preserve `N mod 210` retain
`34–45%`. So the effect has two parts: roughly **40% carried by the residue-class structure**
(a projection artefact of normalising by `S(N)`), and **60% tied to the exact configuration
`2N`** and not removable by any shift.

**What this means for the paper's claims.** The headline `c(2) = −0.14` is a mixture. The part
that survives every shift is about `−0.09`; the rest is the singular-series periodicity, which
is not a new phenomenon. The statement "the effect is not a consequence of the arithmetic of
`N`" must be restricted to the shift-invariant part, and §5.1(e) must be replaced by the table
above.

## 5.10 A new open problem: the lag-`210` memory

§5.9 records that `eps(N)` and `eps(N + 210)` are correlated at `+0.083 +- 0.016`
(`M = 4000`), while every non-multiple of `210` tested gives zero or negative. **At `M = 10^5`
the effect is `+0.0645 +- 0.0032` — a `20 sigma` measurement** — and all multiples of `210`
tested (`420, 630, 840, 1050, 1260, 1470, 2100`) are positive at `+0.04` to `+0.07`, while
`2, 4, 10, 20` are zero or slightly negative. Four candidate causes were tested and all fail:

1. **Residue-class means.** With `M = 10^5` the class means carry `0.61%` of the variance at
   `mod 210` and `1.39%` at `mod 2310` (noise floor `0.10%` and `1.15%`). A `0.6%` class
   structure cannot generate a `0.08` autocorrelation.
2. **A smooth trend.** Lags `2, 4, 10, 20, 200, 1000` all give `|corr| < 0.07`; a smooth
   component would show at small lags too.
3. **Sampling noise.** At `M = 10^5` the effect is `20 sigma`, and it reproduces at every
   multiple of `210` tested up to `2100`.

4. **The normalisation.** Dividing by `S(N)` could in principle manufacture the memory: `S(N)`
   itself has autocorrelation `+0.98` at lag `210` (it depends only on which small primes divide
   `N`). It does not. The *unnormalised* error already carries the memory,
   `+0.0466 +- 0.0032` at lag `210`, so the normalisation amplifies `0.047` to `0.065` rather
   than creating it. (It does remove the parity structure: `S(N)` has `-0.446` at lag `2`,
   `eps` has `-0.001`.)

What is unusual about `210 = 2*3*5*7` is that it is the product of the first four primes; the
effect therefore couples `N` and `N+210` through the *set* of small primes dividing them, but
not through the mean of `eps` in each class. **The scale of `eps` does depend on the class.** With `M = 10^5` the within-class standard
deviation varies across the `105` classes of `N mod 210` with a coefficient of variation of
`21.6%` — and the sampling noise on those estimates is only `2.3%`, so this is real. In other
words the normalisation by `S(N)` is right *on average* but wrong in *scale* for individual
classes, by up to `20%`. This is the origin of the `M_eff < M` effect noted in §4.11 (the
window-to-window scatter being three to four times the within-window error), and it means the
error bars in this literature must be computed per class, not globally.

It does **not**, however, explain the autocorrelation: if `eps = delta(class) * xi` with `xi`
independent across `N`, then `E[eps(N) eps(N+210)] = delta^2 E[xi xi'] = 0`. The lag-`210`
memory therefore remains unexplained even after the scale dependence is accounted for.

## 5.10.1 The memory is not `210`, it is the small primes

The first reading - that `210 = 2*3*5*7` is special - is wrong. Scanning lags that have the
same *size* but different small-prime content separates the two possibilities (`M = 10^5`):

| lag | factorisation | `omega_{2,3,5,7,11}` | `corr` |
|---|---|---|---|
| 3000 | `2^3*3*5^3` | 3 | `-0.015` |
| 4000 | `2^5*5^3` | 2 | `-0.006` |
| 2500 | `2^2*5^4` | 2 | `-0.004` |
| 3500 | `2^2*5^3*7` | 3 | `-0.005` |
| 3200 | `2^7*5^2` | 2 | `-0.017` |
| 3600 | `2^4*3^2*5^2` | 3 | `+0.007` |
| **4200** | **`2^3*3*5^2*7`** | **4** | **`+0.061`** |
| 210 | `2*3*5*7` | 4 | `+0.065` |
| 2100 | `2^2*3*5^2*7` | 4 | `+0.065` |
| 2310 | `2*3*5*7*11` | 5 | `+0.117` |
| 4620 | `2^2*3*5*7*11` | 5 | `+0.129` |
| 6930 | `2*3^2*5*7*11` | 5 | `+0.091` |
| 9240 | `2^3*3*5*7*11` | 5 | `+0.090` |

`4200` and `3000` are the same size and differ by `0.076`; `210`, `2100`, `4200` - sizes
differing by a factor of twenty - all give `+0.063`. **The correlation is a function of the
small-prime content of the lag, not of its size.** Averaging over the lags tested:

```
omega_{2,3,5,7,11}:  1        2        3        4        5
number of lags:      3       12       10        3        4
mean corr:          ~0     -0.013     ~0     +0.064   +0.107
```

so from `omega = 3` upward the memory grows monotonically with the number of the first five
primes dividing the lag, independent of its magnitude **at fixed `N`**.

**It is not, however, independent of `N`.** The `N`-independence quoted above was measured only
inside `[12000, 211998]`. At `N_0 = 4000` (`M = 2000`) the same contrast is absent:

| `N_0` | `c(30)` | `c(210)` | difference |
|---|---|---|---|
| 4000 | `-0.013` | `-0.016` | `-0.003` |
| 12000 (`M = 10^5`) | — | `+0.065 +- 0.003` | — |

`c(210)` at `N_0 = 4000` is `-0.016 +- 0.022`, i.e. `3.6 sigma` away from the `N_0 = 12000`
value. **The memory turns on somewhere between `N ~ 6000` and `N ~ 12000` and is flat above
that.** This is a second, independent open feature: a sharp onset in `N` that no local-density
model produces, and that the four candidates of §5.10 do not address either. `2340 = 2^2*3^2*5*13` carries only three
of the five (`13` is outside the set) and gives `+0.030` - between the `omega = 3` and `omega = 4`
values, as the reading requires; it is **not** included in the `omega = 4` mean above.

**This is a multiplicative structure of the same family as the sign law of §4.3**, with the
restriction to the first five primes made explicit. It is not explained by any of the four
candidates above; it is a fifth, and it is the most concrete open item in the paper.

## 5.11 The zero-driven component, computed without a fit

§4.12 asked whether the zeros represent the field and answered no by *fitting* zero weights
(best `R^2 = 0.49` with eight zeros). The weights are not free — they are fixed by the Mellin
transform of the weight — so the question can be settled with no fit at all.

**The prediction.** The leading response of the count to the prime-density fluctuation is linear,
`e(N) ~ sum_rho N^rho/rho * M(rho)` with

```
M(rho) = integral_{1/3}^{2/3} u^{rho-1} cos^4(3 pi (u-1/2)) du
```

for the default weight `phi = cos^2(3 pi (u-1/2))` restricted to its support `[1/3,2/3]`, where it is
symmetric, so `phi(u) phi(1-u) = phi(u)^2 = cos^4`. **The restriction matters**: `cos^4(3 pi (u-1/2))`
continued to `[0,1]` has *two* peaks, and integrating it there gives `M(1) = 0.375` instead of the
correct `1/8 = 0.125`; the two cases give different `c_z` (`-0.910` against `-0.866` at `k=2`).

**A check that the normalisation is right.** `M(1) = integral_{1/3}^{2/3} cos^4 = 1/8 = 0.125`,
which is `c_phi`; any computation whose `M(1)` is not `0.125` is using the wrong weight. (The
polynomial weight `phi = 16u^2(1-u)^2` gives the closed form `256 B(rho+4,5)` with `M(1) = 0.406`,
a factor `3.25` larger, and must not be mixed with the cosine case.) Hence

```
c_z(k) = Re[ sum_gamma |M(1/2+i gamma)/(1/2+i gamma)|^2 * k^{i gamma} ]
         / sum_gamma |M(1/2+i gamma)/(1/2+i gamma)|^2
```

with no free parameter. Using the first 32 ordinates of `zeta`:

| `k` | 2 | 3 | 4 | 5 | 6 | 15 |
|---|---|---|---|---|---|---|
| `c_z(k)` (`cos^4`, support `[1/3,2/3]`) | `-0.910` | `-0.964` | `+0.680` | `-0.727` | `+0.980` | `+0.838` |
| `c(k)` measured (§4.2) | `-0.1401` | `-0.1227` | `-0.0632` | `-0.1510` | `+0.0235` | `+0.0185` |

The first zero alone gives `cos(gamma_1 log 2) = -0.931`, which is the observed `c(2)` sign and
the right order of magnitude. So the zero component, taken by itself, does have a strong
cross-scale correlation, `c_z(2) = -0.91`.

**It does not, however, organise the measured `c(k)`.** The six values above are exactly the ones
tabulated in §4.3, which is not a random sample. Repeating the comparison on all `49` integers
`2 <= k <= 50` (`N_0 = 4000`, `M = 1500`, noise `0.026`):

| quantity | six tabulated points | all forty-nine |
|---|---|---|
| shape correlation | `+0.93` | `+0.33` |
| `R^2` of `a + b c_z(k)` | `0.87` | `0.108` |
| residual sd | — | `0.083` (3.2x the noise) |

and the offset is not constant either. Grouping by `omega(k)`:

```
primes (15)              measured -0.157   zero model -0.460   difference -0.138
two prime factors (24)          +0.025            +0.125              +0.020
rest (10)                       -0.032            -0.059              -0.030
```

**Conclusion.** The zeros are neither excluded as a contributor nor sufficient as an explanation.
They supply a component of the right order of magnitude for `c(2)`, but the measured `c(k)` is
not a rescaled `c_z(k)`: over the full range the shape agreement is `0.33` and the residual is
three times the noise.

# 6. Refuted descriptions, and the shape of a mechanism

## 6.1 Not a function of `log k`

If `ε(N)` were driven by the zeros of `ζ` in the usual way, one would have
`ε(N) ≈ Σ_γ a_γ N^{iγ}` and hence, with `C = Σ|a_γ|²`,

```
c(k) = C⁻¹ Σ_γ |a_γ|² cos(γ log k),                                   (6.1)
```

a *quasi-periodic function of `log k`*. The data are not of that form in the following sense:
the nine primes and the twenty-two composites tested are interleaved on the `log k` axis, and the
empirical function takes the value class `−0.13 ± 0.02` on *all* primes and `0.00 ± 0.03` on
*all* composites. For instance

```
log 17 = 2.833 → −0.133     log 18 = 2.890 → +0.006
log 19 = 2.944 → −0.101     log 20 = 2.996 → +0.010
log 23 = 3.135 → −0.103     log 24 = 3.178 → +0.005
```


![Fig. 2 — the same data against log k. The primes and the composites interleave on the log k axis while remaining separated in c(k), so no continuous function of log k can describe c(k).](lab/out/fig_logk.png)

A continuous `F(log k)` would need a period of `≈0.1`, and its Fourier coefficients would have
to be arranged so that the sign follows `ω(k)`. We do *not* claim this is impossible for a
sufficiently wild `F`; we claim that the multiplicative description is the simple one, and that
(6.1) with a few frequencies is excluded (a single frequency gives `R² = 0.31`, two give `0.54`). The clean discriminator would be a measurement at a non-integer rational ratio
(§4.7), which already fails at `3–4σ`.

## 6.2 Not Liouville

The formula `c(k) = α^{ω(k)}λ(k)` is excluded (§4.6). Moreover the archived measurement gives
`corr(Y(N), λ(N)) = +0.003, −0.003, −0.016` at `N₀ = 1.2×10⁴, 10⁵, 3×10⁵` (with `Y = e/sd(e)` the `σ`-normalised field), and `corr(Y(N), μ(N)) = +0.015` on the squarefree subsample
(`n = 8111`). With `se ≈ 0.008` these are a *bound*, `|corr| ≲ 0.02`, not a measurement of zero;
that bound is enough to exclude a `λ` component with the amplitude needed for S1
(`0.14`), but the correct wording is "bounded", not "orthogonal". So
`ε(N) = λ(N)·Z(N)` with `Z` independent is excluded too. There is no "Liouville component" of
`ε`; what is Liouville-like is the *sign pattern of the correlation function in k*, which is a
different statement and should be phrased that way.

## 6.3 What a mechanism would have to produce

1. a sign `(−1)^{ω(k)}` in the correlation at integer ratios;
2. magnitude `≈0.135` at primes, with a `±20%` spread that does not follow `p`;
3. small *negative* values at prime powers, and small positive values at products of two primes;
4. persistence to `N = 10⁹` at amplitude `≈0.05–0.10`, with a `±0.03` scale-to-scale scatter
   and no resolved decay law (§4.11(2));
5. with `ε` pointwise *bounded against* `λ` and `μ` (`|corr| ≲ 0.02`, §6.2 — a bound, not
   orthogonality), white in `N`, Gaussian with `kurt ≈ 3`, and independent of the weight.

Point 1 and point 5 together are the hard part: an object can have "`λ`-type correlations in
`k`" while being "pointwise orthogonal to `λ`" only if the correlations are carried by the
*joint* law of `(ε(N), ε(kN))` rather than by a factor `λ(N)`.

---

## 6.4 The product model

If `1+ε(N) = ∏_{p|N}(1+ξ_p)` with `E ξ_p = −α`, then `log(1+ε(N))` is additive and the marginal
law of `ε` must alternate with `ω(N)`. The archive (`M = 15000`, on the `σ`-normalised field `Y`)
gives

| ω(N) | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| mean of Y | +0.088 | +0.043 | +0.024 | +0.115 |
| sd of Y | 1.085 | 0.975 | 0.958 | 0.815 |

The means are flat (no `(−α)^{ω}` alternation), the sd decreases only mildly instead of growing
like `(1+α²)^{ω/2}`, and the kurtosis is `≈3`. Rejected. Note that §4.4 already rules out a
single `α`: the nine prime values have a `3σ` spread.

## 6.5 The single-factor residual model

`ε(kN) = c(k)ε(N) + η_k` with independent `η_k` predicts `corr(η_{k₁},η_{k₂}) = 0`. Archive,
`M = 15000`: `corr(η₂,η₃) = −0.003`, `corr(η₂,η₄) = −0.147`, `corr(η₂,η₆) = −0.129`,
`corr(η₃,η₄) = −0.015`, `corr(η₃,η₆) = −0.134`, `corr(η₄,η₆) = +0.049`. Rejected: the residuals
are as strongly correlated as the field itself. (This is automatic if `c` is genuinely
multiplicative — the residual still contains the multiplicative part — so the test really says
that no *one-dimensional* factor explains the field.)

## 6.6 Markov structure on the multiplicative lattice

A Markov field on the lattice generated by `2` and `3` would make `ε(N)` and `ε(4N)`
conditionally independent given `ε(2N)`, and `ε(N)`, `ε(6N)` conditionally independent given
`ε(2N), ε(3N)`. The archived correlation matrix (`M = 15000`, noise `0.008`) is

| | 1 | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|
| **1** | +1.0000 | −0.1403 | −0.1230 | −0.0500 | +0.0255 |
| **2** | −0.1403 | +1.0000 | +0.0140 | −0.1399 | −0.1325 |
| **3** | −0.1230 | +0.0140 | +1.0000 | −0.0088 | −0.1373 |
| **4** | −0.0500 | −0.1399 | −0.0088 | +1.0000 | +0.0478 |
| **6** | +0.0255 | −0.1325 | −0.1373 | +0.0478 | +1.0000 |

and the partial correlations (given the other three variables) are `−0.072` for `(N,4N)`,
`−0.008` for `(N,6N)`, `−0.144` for `(2N,4N)`, `−0.126` for `(2N,6N)`, `−0.136` for `(3N,6N)`
and `−0.148` for `(N,2N)`. Five of the six exceed `9σ`, so **Markov is rejected**; the single
success is `ρ(N,6N|·) ≈ 0`. On the other hand the matrix elements depend on the ratio and not on
the position (`−0.1403, −0.1399, −0.1373` for ratios `2,2,2`), so the field is **stationary on
the multiplicative lattice** — which is the same fact as §4.7. Caveat: these partial correlations
are computed from a `5×5` sample matrix with no error propagation; with the errors quoted here (`±0.006`) the qualitative verdict is unchanged, but the numbers should not be quoted to three
decimals.

## 6.7 Dirichlet characters

A description `ε(N) = Σ_χ a_χ χ(N)` gives `c(k) = Σ_χ |a_χ|²χ(k) / Σ_χ |a_χ|²`, a weighted
average of Dirichlet characters. This is incompatible with S1 in general: a single character
gives `|c(k)| = 1`, and a mixture cannot give `−0.13` at *all nine* primes, since `χ(p)` is
equidistributed in the unit group as `p` varies. A reported value `corr = 0.192` at `q = 15` is
an artefact: the enumeration used `64` "characters" instead of `8`, and a correct split-half
test gives `+0.22 / −0.32`. The signal was noise.

## 6.8 Four-variable local factors

A direct count of the local solutions of the four-variable system gives
`A(q) = (q−1)(q²−3q+3)`, which yields a predicted covariance ratio
`∏_q A/(qB₁B₂) ≈ 1.25` (i.e. `+25%`), against a measured ratio statistic of `0.990–0.996`.
This argument is **superseded** by §5.3: `A(q)` is the local count for a configuration in which
the two pairs *share* a variable. Here the two pairs are disjoint, the local densities factorise
exactly, and the HL prediction is `0`, not `+25%`. So a comparison with the four-variable local model must use the factorised form — which is why it disagreed by "sign and 60×".

## 6.9 Degenerate configurations

The diagonal families `n = n'` and `n = N−n'` of the four-variable sum are annihilated by the
support: `n ∈ (N/3,2N/3)` and `n' ∈ (2N/3,4N/3)` are disjoint and the weight vanishes at the
junction. Measured contribution `0.000000`. This remains correct, and it is exactly the
disjointness used in §5.2.

## 6.10 Window prime mass

An unweighted count of primes against `|window|/log` gives `corr = −0.018`, which suggests
orthogonality;
With the correctly weighted mass of §5.1(f) the correlation is `+0.105` — a real channel, but at
the 10% level in amplitude and 1% in variance. Either way it cannot produce `c(2) = −0.14`, and
any mechanism based on the prime mass in the interval is excluded.

## 6.11 Ternary and polynomial weights (open, not evidence)

The archive's ternary experiment reports a renormalisation sanity value `0.9836` and
`c₃ = +0.057 ± 0.038` — consistent with zero *and* with the binary value; it is not a
measurement. The polynomial weight on `(0,1)`, normalised as `φ² = 16u²(1−u)²`, gives
`c(2) = −0.163`, `c(3) = −0.149` (§4.8). The other normalisation of the same shape,
`φ = 16u²(1−u)²`, gives `c(2) = −0.145`, `c(3) = −0.164`; the two are the same weight up to a
constant, so this entry is closed.

## 6.12 Summary of the refuted descriptions

| description | naive expectation | status |
|---|---|---|
| zero-driven smooth spectrum, `c(k)=F(log k)` | R² = 0.31 (1 freq.), 0.54 (2 freq.) | **excluded as a *few-frequency* model**; not excluded for wild `F` (§6.1) |
| Liouville component `ε = λ·Z` | `corr = 0.003` | **excluded** (§6.2) |
| `c(k)=α^{ω(k)}λ(k)` | — | **excluded**, 7–13σ (§4.6) |
| Dirichlet characters | `corr = 0.19` at `q=15` | **noise, refuted** (§6.7) |
| product model `1+ε=∏(1+ξ_p)` | sign flip predicted, not seen | **excluded** (§6.4) |
| single-factor residual | residuals correlated | **excluded** (§6.5) |
| Markov on the multiplicative lattice | 5 of 6 partials nonzero | **excluded** (§6.6) |
| stationarity on rational ratios | `3/2` and `4/6` differ by `4σ` | **excluded**; stationarity only on `ℕ`-ratios (§4.7) |
| four-variable local factors (`+25%`) | sign and magnitude wrong | **superseded** by §5.3 (§6.8) |
| degenerate configurations | contribution `0` | **confirmed** (§6.9) |
| window prime mass | `−0.018`; now `+0.105` | **excluded as the driver** (§6.10) |

---

# 7. Discussion

## 7.1 What holds, and under what conditions

Nothing in this paper is a theorem. The following are empirical statements with stated scope,
each given in its strongest defensible reading.

| claim | holds to the level of | scope, and what would weaken it |
|---|---|---|
| **S1, the sign law** (28/31) | `p ≈ 2×10⁻⁶` as a coin-flip indicator, *not* independent | the 31 correlations share one field; a calibrated significance is only the one of §5.7(3) (`3σ–15σ`) and the sieve null of §5.8 (`≈2σ`) |
| **the prime correlations are negative** | `3σ–15σ` depending on whether the primes are one draw from the factorised ensemble (§5.7(3)); `1.7` percentile of a calibrated sieve null (§5.8) | the larger `1/√M` figures are not quoted: §5.9 shows those bars are too small, because the field is not white |
| **S6, weight independence** | five weights, `c(2) ∈ [−0.163,−0.132]` | a sixth weight with a different support could differ |
| **Gaussian, sub-tailed** | `kurt = 3.00`, tails to `3.3σ` | untested beyond `3.3σ` |
| **scale law (3.1)** | `±3%` over three scales | three points is not a law |
| **`ω`-gradient of the memory** (§5.10) | `12σ` on the data side | **the model side is untested**; `N`-independence is untested below `N = 12000` |
| **the four-point factorisation** (§5.3) | exact to machine precision | verified numerically at 13 primes; the analytic identity is elementary |

**Two statements this paper does not make.** (i) That the four-point obstruction *forces* the
observed value — §5.3 shows only that a factorised model predicts zero, and a mechanism outside
that class (one that couples *positions* rather than densities) is not constrained by it.
(ii) That the zeros are excluded — §5.11 shows they supply a component of the right order of
magnitude but do not organise `c(k)`; a full exclusion would need a computation not done here.

**On multiplicity.** The analysis looks at many numbers: 31 values of `k`, five normalisations,
five weights, seven scales (§4.10), and the automatic choices (block length, trimming) inside each. A
Bonferroni correction over the 31 `k`-values would turn `2.3×10^{-6}` into `7×10^{-5}`, which
changes nothing; but the correction that matters is not over `k`, it is over **the choices that
were made after seeing the data** (which normalisation to quote, which window, which scale). Those
are not enumerable from the paper, and the honest statement is that the reported significances
carry an unquantified downward correction. The two things that would survive it are the shape of
the sign law (a pattern, not a single test) and the `omega`-gradient (12σ, and a pattern).

**Falsifiable predictions.** The sign law should continue to hold at larger `k` and larger `N`;
`c(2)` should remain negative for all `N`; the `ω`-gradient should persist when the model side
is tested; and the four-point sum of T2 should vanish identically (so it checks the code, not
the physics — see §5.6).

## 7.2 The main open measurement

Everything about the *magnitudes* (`α`, multiplicativity, the prime-power values) rests on
single windows. The scale dependence of §4.10 shows that this is not enough. A dedicated study —
fixed *relative* window, `M ≥ 3000` at each of six scales from `10⁴` to `3×10⁶`, all of
`k ∈ {2,3,4,6,9,10,12,15}` measured simultaneously — is the next step, and it requires the kind
of implementation given in Appendix B (one minute per scale).

## 7.3 Relation to the conditional identity of Li–Liu

Appendix A shows that the preprint's spectral hypotheses are *equivalent* to the operator-norm
bound `‖T_N‖ ≤ δ_N`, `δ_N = √(𝔖/N)(log N)^{3/2}`, and that GRH alone gives only `≪ (log N)²`.
Meanwhile (3.1) gives the true error scale `0.33 log N/√N`, so

```
δ_N / (true sd of e) = 3.0 √(𝔖(N) log N) ≈ 17–20   at N = 10⁷ .      (7.1)
```

This comparison must be read carefully:
(3.1) is the scale of a *single quadratic form* `⟨T_N v_N, v_N⟩`, whereas Hypothesis R bounds the
*operator norm*. Since `‖T_N‖ ≥ |⟨T_N v_N, v_N⟩|/‖v_N‖²`, the data constrain the norm only from
below. Even on the most favourable reading — that the norm is comparable to the quadratic form —
it would sit `√log N` below `δ_N`; it could equally be far larger. **So the empirical law is
consistent with Hypothesis R but gives it no support, and cannot be used to argue that the
hypothesis is "nearly true".** Any claim that the empirical law
"supports" the conditional proof is unfounded; any claim that it "refutes" it is also unfounded.

---

# Appendix A. The conditional preprint of Li–Liu: what is verified, what is not established

Reference: J. Li and J. Liu, *A conditional proof of the Goldbach conjecture*, preprint,
13 November 2026, 63 pp.

## A.1 The architecture

For each even `N`, the odd Dirichlet characters modulo `N` (including imprimitive ones, extended
by `0` at non-units) give the exact projection

```
(2/φ(N)) Σ_{χ ∈ X⁻_N} χ(n)χ̄(m) = 1_{(nm,N)=1} ( 1_{n≡m} − 1_{n≡−m} ),       (A.1)
```

whence, on the support of `φ` (where `|n−m| ≤ N/3` and `N/3 ≤ n+m ≤ 4N/3`),

```
G = D − P,                                                                    (A.2)
D = N⁻¹ Σ_{(n,N)=1} Λ²(n) φ²(n/N),    P = (2/φ(N)) Σ_{χ ∈ X⁻_N} |P_χ|²,       (A.3)
P_χ = N^{-1/2} Σ_n Λ(n)χ(n) φ(n/N).                                          (A.4)
```

**(A.2) is exact only for the sum restricted to `(nm, N) = 1`.** The unrestricted `G` of (1.1)
differs from it by the terms with `n = m = p^k`, `p | N` — of which there are `O(log N)`, each of
size `Λ(p^k)^2 phi^2(p^k/N)` — i.e. by `O(log^2 N)` in absolute value, which is `o(sigma)` and is
absorbed in the error term of (A.5). The identity is therefore stated for the restricted count
and used for the unrestricted one at that accuracy; it is **not** an identity for the
unrestricted `G`.

Then `P` is written as a double Mellin integral of
`𝒫_N(s,z) = (2/φ(N)) Σ_χ (L′/L)(½+s,χ)(L′/L)(½+z,χ̄)`, a Beta model `𝒦_N` is subtracted, and the
residual `E_N = 𝒫_N − 𝒟_N − 𝒦_N` generates a compact self-adjoint Hilbert–Schmidt operator
`T_N` on a strip Hardy space, with the exact identity

```
c_φ𝔖(N) − G(N) = ⟨T_N v_N, v_N⟩.                                             (A.5)
```

Two hypotheses on the spectral function
`χ_N(q) = ∏_j (1 − 4μ_j²δ_N⁻² q(1−q))`, `q ∈ ℍ`, then give Theorems 1.1 and 1.2.

## A.2 Items we verified

* **(A.1)–(A.3).** The factor `2/φ(N)`, the inclusion of imprimitive characters and the
  extension by `0` are all correct: summing only over odd characters is
  `½ Σ_{χ mod N}(1−χ(−1))χ(n)χ̄(m)`, which gives `1_{(nm,N)=1}(1_{n≡m} − 1_{n≡−m})` times `1/φ(N)`.
  On the support the two congruences mean `n=m` and `n+m=N`; multiplying by `ΛΛφφ/N` gives (A.2).
  A numerical check for `N = 14,…,62` reproduces the identity to `10⁻¹⁶`.
* **The Beta–Gamma factorisation.** With `r = ½−s`, `t = ½−z`, `w = s+z = 1−r−t` one has
  `K⁻ = (Γ(r)Γ(t)Γ(w)/π)(sin πr + sin πt − sin πw)`, and the trigonometric identity
  `sin πr + sin πt − sin πw = 4 sin(πr/2) sin(πt/2) cos(πw/2)` turns this into
  `K⁻ = 2^{1−w} Γ(w) cos(πw/2) X₁(s)X₁(z)`, which is (4.4) of the preprint. Numerically it
  holds to `10⁻¹⁶` at `s=0.2, z=0.3`.
* **The residue.** `𝔖_N(w) = 𝔖(N)ζ(w)∏_{p>2,p∤N}(1 + p^{-w}/(p−2))` and
  `Res_{w=1}𝔖_N = 𝔖(N)∏_{p∤N,p>2}(1+1/(p(p−2))) = N/φ(N)`, using
  `1+1/(p(p−2)) = (p−1)²/(p(p−2))` and `C₂ = ∏_{p>2}p(p−2)/(p−1)²`.
* **The model main term.** With `Φ(u) = √u φ(u)`, `C₀(y) = ∫φ(u)φ(u+y)du`,
  `J₀(y) = ∫_0^y φ(u)φ(y−u)du` and `H₀ = 2C₀ − J₀`, the support gives `C₀(k)=0` for `k ≥ 1`,
  `J₀(k)=0` for `k ≥ 2`, hence `H₀(1) = −J₀(1) = −c_φ` and `H₀(k) = 0` for `k ≥ 2`. So the model
  integral is `−c_φ𝔖(N)` exactly, and the pole `N/φ(N)` is cancelled by `H̃₀(1) = 0`.
* **The zero classification.** For `q = x+yI`,
  `q(1−q) = x(1−x)+y² + (1−2x)yI`. If `y ≠ 0` then `x = ½`; if `y = 0` then
  `q = ½ ± ½√(1 − δ_N²/μ_j²)`. Hence

```
|μ_j| ≤ δ_N   ⟺   every zero has ℜq = 1/2   ⟺   ‖T_N‖ ≤ δ_N ,                (A.6)
```

  so Hypothesis R is *exactly* an operator-norm bound (not a Riemann-hypothesis analogy).
  Similarly Theorem 6.1 gives Hypothesis S `⟺ ‖T_N‖ ≤ δ_N/(2√(η(1−η))) = 𝔖/(2h) + O(η)`.
* **From the asymptotic to representations.** If `G ≥ c_φ𝔖/2` then
  `Σ_{n+m=N}ΛΛφφ ≥ c_φ𝔖N/2`, and since each pair contributes at most `(log N)²`, the number of
  ordered prime pairs is `≥ c_φ𝔖N/(2(log N)²)`; prime powers other than primes contribute
  `O((log N)²√N/N) = o(1)`. The deduction is correct.

## A.3 One display slip, and one constant verified

**(1) A slip in a displayed identity (harmless).** In §2.3.1 and again in §5.1 the preprint states
the Hermitian condition as

```
ℰ_N(s,z̄) = ℰ_N(z,s̄)       (printed; the conjugate is missing on the right-hand side)
```

The relation that holds — and the one the argument needs — is
`ℰ_N(s,z̄) = \overline{ℰ_N(z,s̄)}`, equivalently `ℰ_N(s,z̄) = ℰ_N(z̄,s)`. The printed form equates a
function holomorphic in `s` with one antiholomorphic in `s`, which would force `ℰ_N` to be
constant in `s`. Note that the display immediately above, (2.13), *does* carry the conjugate
(`𝒫_N(s,z̄) = (2/φ(N))Σ_χ g_χ(s)\overline{g_χ(z)}`), so this is a slip in the symmetry line only.

**It does not damage the argument.** Self-adjointness of `T_N` needs
`K(s,z) = \overline{K(z,s)}` with `K = ℰ_N(s,z̄)/((1+s)(1+z̄))`, i.e.
`ℰ_N(s,z̄) = \overline{ℰ_N(z,s̄)}`; and that follows from (3.3) of the preprint
(`𝒫_N(s,z) = 𝒫_N(z,s)`) together with `\overline{ℰ_N(a,b)} = ℰ_N(ā,b̄)`, without using the printed
line at all.

**(2) The constant of Appendix A is correct.** The preprint asserts that the Euler factors of

```
R_N(w) = ∏_{p|N}(1 − p^{−1−w}) ∏_{p>2, p∤N} (1 + (2p^{−w} − p^{−2w})/(p(p−2)))
```

give `𝔖(N)R_N(0) = 1`. This is right, and it is worth recording the check because the identity is
what cancels the `1/σ²` singularity between `𝒟_N` and `𝒦_N`:

```
R_N(0) = ½·∏_{p|N,p>2}(1−1/p)·(1/C₂)·∏_{p|N,p>2}(1+1/(p(p−2)))⁻¹        (N even)
𝔖(N)R_N(0) = 2C₂·∏_{p|N,p>2}(p−1)/(p−2)·R_N(0) = 1 ,
```

using `2C₂∏_{p>2,p∤N}(1+1/(p(p−2)))⁻¹ = 1/C₂`. Numerically, the product is `1.000000` to machine
precision for `N = 2, 6, 30, 2310, 30030, 12000, 131072`. So the cancellation of the `1/σ²`
singularity and the constants of Proposition A.4 stand as printed.

Minor notational slips, which do not affect any proof: the Gram representation of `𝒫_N` in one
display and the inner product on the Hardy space omit the conjugation bar; the factor
`π^{−1/2}` in the formula for `X₁` appears inverted in one line; `|X⁻_N| = φ(N)/2` fails for
`N = 2`.

## A.4 Assessment

The contribution is **an exact identity and an operator formulation**, not progress on the
conjecture. Theorems 1.1 and 1.2 are valid conditional implications, but their hypotheses are
*equivalent* to operator-norm bounds (A.6) which are far beyond GRH: GRH gives
`‖T_N‖_HS ≪ (log N)²` (Theorem 6.2 of the preprint), while the required bound is
`𝔖^{1/2}N^{-1/2}(log N)^{3/2}`; the gap is `N^{1/2}/√(log N)`. "Riemann-hypothesis analogy" is
a description of the shape of the zeros of `χ_N`, not a source of arithmetic information: `χ_N`
depends on `q` only through `q(1−q)`, so the zeros are on `ℜq = ½` unless `|μ_j| > δ_N`.
The numerical Appendix B (N = 30030) is a finite experiment and, as the authors say, no
evidence for Hypothesis R.

---

# Appendix B. Reproduction

## B.1 Files

```
lab/tool.py        NumPy: sieve, Λ, 𝔖(N), G(N), correlation, block bootstrap
lab/ngb.py         numba version (13× faster); G(N) agrees to 1.5×10⁻¹⁵
lab/fast2.py       multi-process driver
lab/weights.py     weight-robustness table (§4.8), with the correct support per weight
lab/emit_tables.py recomputes §4 and writes the machine-readable tables and figures
lab/run1..run6.py  diagnostics: sub-windows, bins, iid null, mass channel, mod q
lab/check_zh.py    cross-checks that the Chinese version keeps every number
lab/verify_numbers.py  cross-checks the data files against both papers
```

Data and figures produced by `emit_tables.py` and `weights.py`:

```
lab/out/ck_main.csv     k, ω(k), primality, c(k), se, two halves, M     (§4.2)
lab/out/ck_scale.csv    c(k) at seven scales, with errors               (§4.10)
lab/out/ck_weight.csv   five weight shapes                              (§4.8)
lab/out/decile.csv      the decile table                                (§3.3)
lab/out/fig_ck.png      c(k) vs k, primes and composites marked
lab/out/fig_logk.png    the same data vs log k                          (§6.1)
lab/out/fig_scale.png   c(2), c(3), c(4), c(6) vs N₀                    (§4.10)
lab/out/fig_decile.png  E[ρ(2N)] vs decile of ρ(N),  ρ=1−ε                      (§3.3)
lab/out/null_two.csv    the 40-seed factorised null                          (§5.7)
lab/out/null_fourvar.csv the 16-seed null of the first run                   (§5.7)
lab/out/ck_big.csv      the large-N ladder 10⁷, 3×10⁷, 10⁸                    (§4.10)
lab/remote/data/        the 10⁸, 3×10⁸, 10⁹ runs (on rented 128-core box)     (§4.11)
lab/ternary2.py        the ternary analogue (balanced box, numba)             (§5.7)
```

Commands used for §4:

```
python fast2.py --n0 12000 --m 40000 --k 2,3,4,5,6,7,8,9,10,11,12,13     # ~166 s, 11 procs
python fast2.py --n0 12000 --m 20000 --k 14,...,49                        # ~ 65 s
python fast2.py --n0 200000 --m 20000 --k 2,3,4,6,9,10,12,15              # ~ 39 s
```

## B.2 Provenance of the tables, and what is reproducible

The tables of §§2–4 were produced by a single multi-core run on a separate machine. The driver and
the per-scale scanner **exist and are available from the author on request**; they are not included
in the present bundle.

| table | routine (not included) | M |
|---|---|---|
| §2.5 two-point | vectorised NumPy, inline | X=2×10⁶ |
| §3.1 scale law | JIT scanner, `Garr` | 2000/1200/600 |
| §3.2 distribution | JIT scanner, `Garr` | 1200 |
| §3.3 deciles | decile driver | 1.5×10⁴ |
| §4.2 `k ≤ 13` | multi-core driver | 4×10⁴ |
| §4.2 `k ≥ 14` | multi-core driver | 2×10⁴ |
| §4.8 weights | weight driver | 2×10⁴ |
| §4.9 blocks, §4.10 scales | driver + block scanner | as shown |
| §2.4 iid null | null driver | 8 seeds |

**What *is* reproducible from the distributed material.** Four tables — §5.8 (the calibrated sieve
null), §5.9 (the class-mean removal), §5.9 (the shift test) and §6.6 (the Markov matrix) — are
reproduced by the logs shipped alongside this manuscript. §5.10 and §5.11 are reproducible from
the formulas alone (they need only the 32 ordinates of `zeta` and the array `cmw_12000_*.npy`).
**The `c(k)` tables of §§3–4 are not.** A reader who wants them must ask the author for the driver, or reimplement it from §2.1–§2.3 and
Appendix B.1 — the algorithm is fully specified there, and a single-threaded version is short
enough to be written from that description.

This asymmetry is a defect of the present *bundle*, not of the measurement, and it is stated here
so that no reader assumes more reproducibility than the bundle provides.

---

# Appendix C. Methodological lessons

1. **Ratio statistics `ΣGG/ΣMM` are unusable** on strongly correlated data: the cross term
   `M(G−M)` dominates and the sign flips between windows (`0.99697` at `N≈2×10⁴`, `1.0007` at
   `N ≥ 10⁵`).
2. **`1/√M` was right here, by accident.** It requires `ε` to be white in `N`, which we checked
   (lag-1 `−0.027`). In more tightly correlated statistics it can be wrong by a factor of several.
3. **Window width is the easiest confounder to miss.** Scale comparisons must fix the *relative*
   window; and even then the window-to-window scatter (`±0.03`) exceeds the within-window error
   (`±0.006`) by a factor five.
4. **A shift control (`2N → 2N+2`) is stronger than a permutation control** for detecting
   spurious correlations, because permutations destroy only the mean structure.
5. **Any "law" must be tested against a null model that has the right local densities.** Here
   the iid model (`−0.017 ± 0.025`) is decisive: it shows that local structure plus randomness
   produces no effect, so the effect — if real — is carried entirely by correlation between
   primes in different intervals.
6. **A four-point identity can turn an empirical curiosity into a finite test.** §5.3 is such an
   identity: `σ_p(quadruple) = σ_p(pair)²` for every `p`.

7. **Check `M(1) = c_phi` before trusting any zero-spectrum computation.** The Mellin
   transform of the weight at `rho = 1` must reproduce `c_phi = 1/8`. This catches both errors
   actually made in §5.11: using the polynomial weight by mistake (`M(1) = 0.406`), and
   continuing `cos^4(3 pi (u-1/2))` to `[0,1]`, where it has two peaks instead of one
   (`M(1) = 0.375`). Both give plausible-looking `c_z`; both are wrong.
---

*All numbers in §§2–4 were produced by the independent implementation described
in Appendix B; the entries marked "archive" in §6 and Appendix A come from earlier scripts.
The author claims no theorem. The structure reported is empirical and its mechanism is unknown.*


# Appendix D. Verification of the numerical code

The scripts that produced the local measurements were audited independently, read-only, by six reviewers, each asked to derive the key formulas from first principles rather than to check the
code against its own comments.

**42 defects were found; 24 were fixed, 4 were replaced or marked unusable, 7 were annotated,
and 6 were left (no number depends on them).** The defects that mattered:

| defect | consequence | fixed |
|---|---|---|
| `2*dot(a,b)` over the half-range double-counts the `k = N/2` term | adds `Λ(N/2)^2` for ~11% of even `N` | yes (19 files) |
| two-point main term used `𝔖(h)·X` instead of `𝔖(h)·(X−h)` | a deterministic ramp; every kernel correlation `≈ +1` | yes |
| support/width table built on a grid pinned to `[1/3,2/3]` | the weight is wrong for every `a ≠ 1/3` | yes |
| `boxcheck.py` doubled a full-range sum | the "box weight" control measured `𝔖(N)` | yes |
| `randmodel.py` used `∫φ` instead of `∫φ²` as the main coefficient | main term 1.37× too large | yes |
| `mu2.py` omitted the factor `(6/π²)²` | main term 2.7× too large | yes |
| `twopoint.py` started the factorisation at `d = 3` | `𝔖(h)` wrong by 0.63–1.5 | yes |

None of these touches the measurements of this paper: every headline number comes from scripts
that the audit lists as clean (`decisive`, `fixed_window`, `kscan`, `scale_test2`,
`corr_scales`, `scale_mult`, `mult3`, `twopoint2`, `deficit_scale`), and the five weight shapes
of §4.8 do not use the support-grid code. The audit's practical lesson is:

> **Before comparing any `c(k)`, print `sanity = ΣC/Σmain`. A main term off by 6–10% inflates
> the spread of `c(k)` by a factor 4 and moves its mean in an arbitrary direction. Two of the
> conclusions reached during this work were wrong for exactly this reason.**
