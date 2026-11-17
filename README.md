# Multiplicative correlations in the Goldbach error field

*https://github.com/xlpoistop-prog/goldbach-error-field*

**An empirical note, and a reproducibility package.**  No theorem is proved here.  What is
reported is a measurement, the checks it survives, the explanations it kills, and the
falsification test it suggests.

* manuscript: [`paper/PAPER-EN.pdf`](paper/PAPER-EN.pdf) (English, 26 pp) ·
  [`paper/PAPER-zh.pdf`](paper/PAPER-zh.pdf) (Chinese, 24 pp) ·
  markdown and HTML in [`paper/`](paper/)
* code: [`code/`](code/) — every table and figure can be regenerated
* data: [`data/`](data/) — the csv files behind the tables
* how to run everything: [`REPRODUCING.md`](REPRODUCING.md)
* checksums: [`SHA256SUMS.txt`](SHA256SUMS.txt)

## What the measurement is

Let `G(N)` be a smooth weighted *balanced* Goldbach count and
`ε(N) = 1 − G(N)/(c_φ𝔖(N))` its relative deviation from the Hardy–Littlewood main term.
The field `ε` is white in `N`, Gaussian, with `sd(ε) ∝ log N/√N`.  Its two-point correlation
across multiplicative dilations,

```
c(k) = corr( ε(N), ε(kN) )
```

is **negative at every prime ratio tested** (`−0.14` at `10⁴`, `−0.047(6)` at `10⁹`) and
**≈ 0 at composite ratios with two distinct prime factors** — a sign pattern `(−1)^{ω(k)}` that
holds in every window from `1.2×10⁴` to `10⁹`.

## Why it is surprising

The four-point local densities of `(n, N−n, n′, kN−n′)` factorise **exactly** at every prime
(the code checks this to machine precision in `code/fourvar.py`).  Any model of the primes with
the correct local densities and independent positions therefore predicts `c(k) = 0`.  So the
correlation is a genuine deviation from the Hardy–Littlewood picture, and not a small one in
statistical terms: about 15σ on the window where it is largest.

## What is excluded (all of it in the paper, with numbers)

| reading | status |
|---|---|
| the estimator (5 normalisations, 5 weights, block localisation, `N mod q`, exponent regression, shift control) | excluded |
| a deterministic arithmetic bias ("second-order singular series") | excluded: doubling fixes the odd part of `N`, so any such bias contributes `≥ 0`; the 2-adic part is bounded 350× too small |
| the local prime density as a modulator | excluded: ten windows at `10⁹` agree (`χ² = 7.1/9`) while the local mass varies by 8% across them |
| the classical zero-driven error term | excluded: built from 48 533 actual zeros of `ζ`, the field has lag-1 autocorrelation `+0.99`, the measured field has `−0.027` — the zeros make a smooth field, the data are white (`code/zerofield2.py`, 19 s) |
| the ternary analogue as a control | not usable: its noise hierarchy is different (density-dominated, not pairing-dominated) |

## What is *not* claimed

* no mechanism;
* no decay law — the scale-to-scale scatter (`±0.03`) is five times the within-window error, and
  at this statistic **scale and sample are the same axis**, so the two cannot be separated by
  construction (this is the main structural obstacle, §4.11(5) of the paper);
* no relation to a proof of Goldbach's conjecture.

## The falsification test

`c(2)` should stay negative at every scale, and the sign law should continue to hold at larger
`k`.  The sharpest single test is a measurement at a **non-integer rational ratio**: a
multiplicative function of `k` predicts `c(3/2) = c(2)c(3) = c(6)`, while a "property of the
interval pair" does not.  The archived data already show a `3–4σ` spread there (§4.7).

## Status

Unrefereed.  Written by a non-specialist with heavy AI assistance; treated as a reproducible
measurement rather than a result.  If you find an artefact that explains it, an issue here is the
most useful thing you can do with this repository.

## License

Not chosen yet — pick one before reusing (MIT or BSD-3-Clause for the code, CC BY 4.0 for the
text are the usual choices).
