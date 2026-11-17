# Code for *Multiplicative correlations in the Goldbach error field*

This directory reproduces every number, table and figure of the preprint.

```
code/     the scripts (flat; they import each other, so run them from this directory)
data/     the machine-readable data behind the tables, and the four figures
paper/    a snapshot of the manuscript (both languages) at release time
```

## Requirements

Python >= 3.10 with `numpy`, `numba` and `matplotlib`.  The numbers in this release were
produced with Python 3.13.5, numpy 2.5.4, numba 0.68.1, matplotlib 3.11.2 on Linux
(Debian 13, x86-64) and on Windows.  No GPU is used, and nothing needs a network connection.

`numba` is optional in the sense that `tool.py` works without it, but `ngb.py` (and therefore
`fast2.py`, `emit_tables.py`, `bigscale.py`) needs it; it is 13x faster.

## What to run

All commands are from `code/`.

| script | paper | what it does | cost |
|---|---|---|---|
| `tool.py` | §2.1 | the estimator: sieve, Λ, 𝔖(N), G(N), correlation, block bootstrap | library |
| `ngb.py` | §2.2 | numba re-implementation of G(N) (13x faster) | library |
| `errbar.py` | §2.3 | how large is the error bar on c(2) | seconds |
| `run6.py` | §2.4 | the iid null model (8 seeds) | ~1 min |
| `run1.py`, `run2.py`, `run5.py` | §5.1 | local-arithmetic, mass-channel and weight diagnostics | 1–3 min each |
| `run4.py` | §3.3 | the decile table (ρ vs 1−ε) | ~1 min |
| `run3.py` | §4.9, §4.10 | block localisation and the scale ladder | minutes |
| `fast2.py` | §4.2 | the main k-table (the paper's headline numbers) | ~3 min for k ≤ 13 at M = 4e4 on 11 cores |
| `weights.py` | §4.8 | the five weight shapes | ~5 min |
| `emit_tables.py` | §3, §4 | regenerates `data/*.csv` and the four figures | ~20 min |
| `bigscale.py` | §4.10 | the 10⁷ / 3×10⁷ / 10⁸ ladder | tens of minutes |
| `stab.py` | §5.7(4) | variance stabilisation and inverse-variance weighting | seconds |
| `subpoisson.py` | §5.7(a) | is the prime mass sub-Poisson | seconds |
| `mass_diag.py` | §5.1(f) | size of the weighted-mass fluctuation | seconds |
| `exponents.py` | §5.7(c) | regress ε on the exponent vector v_p(N) | seconds |
| `fourvar.py`, `fourvar2.py` | §5.7(1)(2) | the four-variable sum, computed explicitly | minutes |
| `null_par.py` | §5.7(3) | the 40-seed factorised null, all cores | ~2 min |
| `randsieve.py`, `null_hawkins.py` | §5.7(b) | the sieve-type nulls (the first is miscalibrated; see the paper) | ~10 min |
| `ternary2.py` | §5.7(d) | the ternary analogue on a balanced box | ~2 min |
| `gcfar2.py` | §4.11 | the 10⁸ / 3×10⁸ / 10⁹ runs (segmented sieve, fork pool) | see below |
| `win10.py` | §4.11(5) | ten windows at one scale; local-mass correlation | ~6 min at 10⁹ |
| `subwin2.py` | §4.11(5) | ten sub-windows inside one window | ~4 min |
| `zerofield2.py` | §4.12 | the zero-driven field from 48k actual ζ zeros | **19 s** |
| `check_sr.py` | §A.3 | verifies 𝔖(N)R_N(0) = 1 numerically | <1 s |
| `verify_numbers.py`, `check_zh.py`, `check_html.py` | — | self-checks: data vs papers, English vs Chinese, HTML tables | <1 s |
| `md2html.py` | — | builds the HTML from the markdown | <1 s |

## The large-N runs

`gcfar2.py`, `win10.py`, `subwin2.py` and `zerofield2.py` are written for a many-core Linux box
and were run on a rented 128-vCPU machine (AMD EPYC, 246 GB RAM).  They build one `float32`
array of `Λ` up to `10¹⁰` (40 GB) and share it with a fork pool, so they need

* for `gcfar2.py --n0 1e8`: ~4 GB and ~10 min on 100 cores;
* for `gcfar2.py --n0 1e9`: ~44 GB and ~30 min on 100 cores;
* for `win10.py 1e9 3000 10 2000000`: ~8 GB and ~6 min on 120 cores;
* for `zerofield2.py`: a laptop; it computes the zeros and the fields in 19 s.

On Windows the fork-based parallelism is unavailable (`multiprocessing` uses spawn), so those
four scripts are Linux-only as written.

## Notes

* `data/` holds exactly the files listed in Appendix B of the paper, with the same names.  The
  `*.csv` files carry four decimals; the paper prints two to four, so a literal comparison needs
  rounding — `verify_numbers.py` does that and reports zero mismatches.
* `win10.py` and `subwin2.py` share roughly 99.8% of their primes between windows at a fixed
  scale; §4.11(5) of the paper explains what that does and does not test.
* `randsieve.py` and `null_hawkins.py` are published *as failures*: neither reproduces the
  primes' variance structure well enough to calibrate a null.  The paper says so explicitly.
* The manuscript snapshot in `paper/` is a convenience for the self-checks; cite the published
  version, not this copy.

## License

No license file is included.  Add one before publishing (MIT or BSD-3-Clause are the usual
choices for research code); the author should decide.
