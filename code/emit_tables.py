"""Recompute the paper's main tables and figures, write machine-readable output.

Outputs (in lab/out/):
  ck_main.csv     c(k), se, halves for the main window (M=4e4 for k<=13, 2e4 for 14..49)
  ck_scale.csv    c(k) at several scales
  ck_weight.csv   weight-shape robustness
  fig_ck.png      c(k) vs k with error bars (primes marked)
  fig_logk.png    c(k) vs log k (same data)
  fig_scale.png   c(2),c(3),c(4),c(6) vs N0
  fig_decile.png  E[eps(2N)] vs decile of eps(N)
"""
import os, sys, csv
import numpy as np
from multiprocessing import Pool
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool, ngb

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
STATE = {}


def init_worker(nmax):
    lam, pp, primes = tool.sieve(nmax)
    STATE["lam"] = lam; STATE["pp"] = pp
    STATE["spf"] = tool.spf_sieve(nmax)
    STATE["twoC2"] = tool.c2_constant()


def job(args):
    k, Ns, wid = args
    lam = STATE["lam"]; pp = STATE["pp"]
    r1 = ngb.rho(Ns, lam, pp, STATE["spf"], STATE["twoC2"], wid)
    rk = ngb.rho(k * Ns, lam, pp, STATE["spf"], STATE["twoC2"], wid)
    c = tool.corr(r1, rk)
    se = tool.block_se(r1, rk, L=40, reps=300, seed=1)
    h = len(r1) // 2
    return k, c, se, tool.corr(r1[:h], rk[:h]), tool.corr(r1[h:], rk[h:])


def run(ks, n0, m, wid=0, procs=11):
    Ns = (n0 + 2 * np.arange(m)).astype(np.int64)
    nmax = int(max(ks) * Ns[-1]) + 10
    with Pool(procs, initializer=init_worker, initargs=(nmax,)) as pool:
        return pool.map(job, [(k, Ns, wid) for k in ks])


def isprime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def omega(n):
    c = 0; d = 2
    while d * d <= n:
        if n % d == 0:
            c += 1
            while n % d == 0:
                n //= d
        d += 1
    return c + (1 if n > 1 else 0)


if __name__ == "__main__":
    # ---- main table ----
    ks_small = list(range(2, 14))
    ks_big = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 32, 36, 49]
    rows = run(ks_small, 12000, 40000) + run(ks_big, 12000, 20000)
    rows.sort()
    with open(os.path.join(OUT, "ck_main.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["k", "omega", "is_prime", "c_k", "se", "half1", "half2", "M"])
        for k, c, se, h1, h2 in rows:
            w.writerow([k, omega(k), int(isprime(k)), "%.4f" % c, "%.4f" % se,
                        "%.4f" % h1, "%.4f" % h2, 40000 if k <= 13 else 20000])
    print("ck_main.csv written")

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    kk = np.array([r[0] for r in rows]); cc = np.array([r[1] for r in rows]); ss = np.array([r[2] for r in rows])
    pr = np.array([isprime(k) for k in kk], dtype=bool)
    ax.errorbar(kk[pr], cc[pr], yerr=ss[pr], fmt="o", ms=4, lw=1, color="crimson", label="k prime")
    ax.errorbar(kk[~pr], cc[~pr], yerr=ss[~pr], fmt="s", ms=3.5, lw=1, color="steelblue", label="k composite")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("k"); ax.set_ylabel("c(k) = corr(eps(N), eps(kN))")
    ax.legend(fontsize=8); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_ck.png"), dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    lg = np.log(kk)
    ax.errorbar(lg[pr], cc[pr], yerr=ss[pr], fmt="o", ms=4, lw=1, color="crimson", label="k prime")
    ax.errorbar(lg[~pr], cc[~pr], yerr=ss[~pr], fmt="s", ms=3.5, lw=1, color="steelblue", label="k composite")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("log k"); ax.set_ylabel("c(k)")
    ax.legend(fontsize=8); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_logk.png"), dpi=160); plt.close(fig)

    # ---- scale table ----
    scales = [(12000, 40000), (50000, 20000), (200000, 20000), (300000, 15000),
              (700000, 8000), (1000000, 6000), (1500000, 5000)]
    ks_s = [2, 3, 4, 6, 9, 10, 12, 15]
    srows = []
    for (n0, m) in scales:
        r = {x[0]: x for x in run(ks_s, n0, m)}
        srows.append((n0, m, r))
        print("scale", n0, "done", flush=True)
    with open(os.path.join(OUT, "ck_scale.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["N0", "M"] + ["c%d" % k for k in ks_s] + ["se%d" % k for k in ks_s])
        for n0, m, r in srows:
            w.writerow([n0, m] + ["%.4f" % r[k][1] for k in ks_s] + ["%.4f" % r[k][2] for k in ks_s])
    print("ck_scale.csv written")

    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    n0s = np.array([s[0] for s in srows], dtype=float)
    for k, col in zip([2, 3, 4, 6], ["crimson", "steelblue", "seagreen", "darkorange"]):
        v = np.array([s[2][k][1] for s in srows]); e = np.array([s[2][k][2] for s in srows])
        ax.errorbar(n0s, v, yerr=e, fmt="o-", ms=4, lw=1, label="k=%d" % k, color=col)
    ax.set_xscale("log"); ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("window start N0"); ax.set_ylabel("c(k)")
    ax.legend(fontsize=8); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_scale.png"), dpi=160); plt.close(fig)

    # ---- weight robustness ----
    wrows = []
    for wid, name in [(0, "cos4"), (1, "flat"), (2, "narrow"), (3, "wide")]:
        r = {x[0]: x for x in run([2, 3, 4, 6, 10], 12000, 20000, wid=wid)}
        wrows.append((name, r))
        print("weight", name, "done", flush=True)
    with open(os.path.join(OUT, "ck_weight.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["weight"] + ["c%d" % k for k in [2, 3, 4, 6, 10]] + ["se%d" % k for k in [2, 3, 4, 6, 10]])
        for name, r in wrows:
            w.writerow([name] + ["%.4f" % r[k][1] for k in [2, 3, 4, 6, 10]] + ["%.4f" % r[k][2] for k in [2, 3, 4, 6, 10]])
    print("ck_weight.csv written")

    # ---- decile figure ----
    init_worker(2_600_000)
    Ns = (12000 + 2 * np.arange(15000)).astype(np.int64)
    CPHI = 0.125      # ngb.rho returns G/S; the paper tabulates rho = G/(c_phi S) = 1 - eps
    r1 = ngb.rho(Ns, STATE["lam"], STATE["pp"], STATE["spf"], STATE["twoC2"], 0) / CPHI
    r2 = ngb.rho(2 * Ns, STATE["lam"], STATE["pp"], STATE["spf"], STATE["twoC2"], 0) / CPHI
    qs = np.quantile(r1, np.arange(0, 1.01, 0.1))
    dec = np.clip(np.searchsorted(qs, r1, side="right") - 1, 0, 9)
    m1 = np.array([r1[dec == d].mean() for d in range(10)])
    m2 = np.array([r2[dec == d].mean() for d in range(10)])
    np.savetxt(os.path.join(OUT, "decile.csv"), np.c_[np.arange(1, 11), m1, m2],
               delimiter=",", header="decile,E_rho_N,E_rho_2N", comments="")
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(m1, m2, "o-")
    ax.set_xlabel(r"$E[\rho(N)]$ by decile,  $\rho=G/(c_\varphi\mathfrak{S}(N))=1-\varepsilon$")
    ax.set_ylabel(r"$E[\rho(2N)]$")
    ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_decile.png"), dpi=160); plt.close(fig)
    print("decile.csv, fig_decile.png written")
