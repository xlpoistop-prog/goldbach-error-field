"""Two nulls, many seeds, all cores.

NULL A (independent): random set, coprime to every p <= 97, density 1/log n, independent
draws.  E[G] analytic.  Its field is mass-dominated (smooth) - not a perfect mimic.

NULL B (mass-pinned): same, but in every block of B integers the *number* of survivors is
pinned to its expectation.  This removes the mass channel by construction, so the field is
pairing-noise dominated like the real one, while the four-point correlations still factorise.
E[G] is analytic with q(n) = survivors_in_block / candidates_in_block.
"""
import os
import sys
import os
import csv
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tool
import ngb

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
X = 2_600_000
N0, M = 12000, 15000
Ns = (N0 + 2 * np.arange(M)).astype(np.int64)
CPHI = 0.125
P = 97
B = 2000
STATE = {}


def build_common(nmax):
    lam, pp, primes = tool.sieve(nmax)
    spf = tool.spf_sieve(nmax)
    twoC2 = tool.c2_constant()
    sv = np.ones(P + 1, bool); sv[:2] = False
    for i in range(2, int(P ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = False
    sp = np.nonzero(sv)[0]
    cop = np.ones(nmax + 1, bool); cop[:2] = False
    for p in sp:
        cop[p::p] = False
    frac = float(np.prod(1.0 - 1.0 / sp))
    cand = np.nonzero(cop)[0]
    cand = cand[cand >= 2]
    logn = np.log(np.maximum(np.arange(nmax + 1), 2))
    prob = np.zeros(nmax + 1)
    prob[cand] = np.minimum(1.0 / (logn[cand] * frac), 1.0)
    # block structure
    blk = cand // B
    nb = blk.max() + 1
    cnt_cand = np.bincount(blk, minlength=nb)
    exp_cnt = np.bincount(blk, weights=prob[cand], minlength=nb)
    ktar = np.floor(exp_cnt).astype(np.int64)
    q_blk = np.where(cnt_cand > 0, ktar / np.maximum(cnt_cand, 1), 0.0)
    STATE.update(dict(lam=lam, pp=pp, spf=spf, twoC2=twoC2, cand=cand, logn=logn,
                      prob=prob, blk=blk, ktar=ktar, nb=nb, q_blk=q_blk, cnt_cand=cnt_cand))


def init_worker(nmax):
    build_common(nmax)


def EG(Nsx, q):
    pp = STATE["pp"]; logn = STATE["logn"]
    out = np.empty(Nsx.shape)
    for i, N in enumerate(Nsx):
        j = np.searchsorted(pp, N // 2)
        n = pp[:j]
        m = N - n
        w = tool.weight(n / N)
        out[i] = 2.0 * np.sum(logn[n] * logn[m] * q[n] * q[m] * w) / N
    return out


def draw(seed, pinned):
    rng = np.random.default_rng(seed)
    cand = STATE["cand"]; blk = STATE["blk"]
    if not pinned:
        keep = rng.random(cand.size) < STATE["prob"][cand]
        return _mk(cand[keep])
    # pinned: within each block keep the ktar smallest uniforms
    u = rng.random(cand.size)
    order = np.lexsort((u, blk))
    cb = cand[order]; kb = blk[order]
    # rank within block
    starts = np.searchsorted(kb, np.arange(STATE["nb"]))
    rank = np.arange(cb.size) - starts[kb]
    keep = rank < STATE["ktar"][kb]
    return _mk(cb[keep])


def _mk(idx):
    a = np.zeros(X + 1, bool)
    a[idx] = True
    return a


def job(args):
    kind, seed = args
    alive = draw(seed, kind == "pinned")
    lg = np.where(alive, STATE["logn"], 0.0)
    q = STATE["q_blk"][STATE["blk"]] if kind == "pinned" else STATE["prob"]
    eg1 = EG(Ns, q); eg2 = EG(2 * Ns, q)
    G1 = ngb.Garr(Ns, lg, STATE["pp"])
    G2 = ngb.Garr(2 * Ns, lg, STATE["pp"])
    e1 = G1 / eg1 - 1; e2 = G2 / eg2 - 1
    W = 40
    nb = len(e1) // W
    v = np.empty(nb)
    for i in range(nb):
        x = e1[i * W:(i + 1) * W]; y = e2[i * W:(i + 1) * W]
        x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
        v[i] = np.mean(x * y)
    return kind, seed, float(tool.corr(e1, e2)), float(v.mean()), float(v.std() / np.sqrt(nb)), float(e1.std())


if __name__ == "__main__":
    init_worker(X)
    print("candidates: %d, blocks: %d, mean survivors per block: %.1f"
          % (STATE["cand"].size, STATE["nb"], STATE["ktar"].mean()))
    jobs = [("indep", s) for s in range(40)] + [("pinned", s) for s in range(24)]
    with Pool(10, initializer=init_worker, initargs=(X,)) as pool:
        res = pool.map(job, jobs)
    rows = []
    for kind in ["indep", "pinned"]:
        rr = [r for r in res if r[0] == kind]
        full = np.array([r[2] for r in rr]); loc = np.array([r[3] for r in rr])
        sds = np.array([r[5] for r in rr])
        print("%-6s  n=%2d  sd(eps)=%.3f   full-window %+.4f +- %.4f   within-block %+.4f +- %.4f"
              % (kind, len(rr), sds.mean(), full.mean(), full.std(ddof=1) / np.sqrt(len(rr)),
                 loc.mean(), loc.std(ddof=1) / np.sqrt(len(rr))))
        rows += [[r[0], r[1], "%.4f" % r[2], "%.4f" % r[3], "%.4f" % r[4], "%.4f" % r[5]] for r in rr]
    print("real primes: sd(eps)=0.086  full -0.1400   within-block -0.1348 +- 0.0087")
    with open(os.path.join(OUT, "null_two.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["kind", "seed", "c2_full", "c2_within", "se_within", "sd_eps"])
        w.writerows(rows)
    print("written lab/out/null_two.csv")
