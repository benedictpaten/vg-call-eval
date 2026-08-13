#!/usr/bin/env python3
"""Does any unused per-site signal rank calls better than GQ already does?

Motivated by the observation that the genotype likelihood uses only *which* alleles a read
fits: the mixture assumes each haplotype contributed exactly 1/|G| of the reads, with no
sampling model and no length adjustment. Several quantities the caller computes are
therefore never used to rank a call:

  - `best_ln`   how well reads fit their *best* allele in absolute terms -- the row divisor,
                which §5.3.1 flagged as a realignment trigger and never tested;
  - `share`     the fraction of reads whose best-fitting allele is one of the called ones,
                i.e. goodness of fit of the called genotype;
  - `balance`   how reads divide between the two alleles of a heterozygous call;
  - `depth`     reads at the site.

The question is not whether these correlate with correctness -- weak signals usually do --
but whether they add anything **beyond GQ**, which is already emitted and already ranks
calls well. So each is scored alone, and then against a logistic model fitted on held-out
data. An equal-weight rank-sum was tried first and is not reported: it dilutes a strong
signal with a weak one, and gave opposite answers on two subsets of the same data.

Logistic regression is hand-rolled because this environment has no numpy or sklearn. It is
a plain gradient ascent on standardised features -- adequate for asking whether an extra
feature moves a held-out AUC, which is all that is being asked.
"""

from __future__ import annotations

import argparse
import gzip
import math
import random
import re
import subprocess
from pathlib import Path

STEP = re.compile(r"([><])(\d+)")


def id_to_site(vid: str) -> str | None:
    m = STEP.findall(vid)
    if len(m) < 2:
        return None
    (d1, n1), (d2, n2) = m[0], m[-1]
    return f"{n1}{'+' if d1 == '>' else '-'}_{n2}{'+' if d2 == '>' else '-'}"


def auc(pos: list[float], neg: list[float]) -> float:
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    n1, n0 = len(pos), len(neg)
    if not n1 or not n0:
        return float("nan")
    i, rank, s = 0, 1, 0.0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        avg = (rank + rank + (j - i) - 1) / 2
        for k in range(i, j):
            if allv[k][1] == 1:
                s += avg
        rank += j - i
        i = j
    return (s - n1 * (n1 + 1) / 2) / (n1 * n0)


def standardise(rows: list[list[float]]) -> list[list[float]]:
    n_feat = len(rows[0])
    out = [r[:] for r in rows]
    for f in range(n_feat):
        vals = [r[f] for r in rows]
        m = sum(vals) / len(vals)
        sd = math.sqrt(sum((v - m) ** 2 for v in vals) / len(vals)) or 1.0
        for r in out:
            r[f] = (r[f] - m) / sd
    return out


def fit_logistic(X: list[list[float]], y: list[int], epochs: int = 300,
                 lr: float = 0.5) -> list[float]:
    n_feat = len(X[0])
    w = [0.0] * (n_feat + 1)
    n = len(X)
    for _ in range(epochs):
        g = [0.0] * (n_feat + 1)
        for xi, yi in zip(X, y):
            z = w[0] + sum(w[k + 1] * xi[k] for k in range(n_feat))
            p = 1.0 / (1.0 + math.exp(-max(-30.0, min(30.0, z))))
            d = yi - p
            g[0] += d
            for k in range(n_feat):
                g[k + 1] += d * xi[k]
        for k in range(n_feat + 1):
            w[k] += lr * g[k] / n
    return w


def score(w: list[float], xi: list[float]) -> float:
    return w[0] + sum(w[k + 1] * xi[k] for k in range(len(xi)))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work", required=True)
    p.add_argument("--support", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--min-svlen", type=int, default=0,
                   help="restrict to calls with an allele at least this many bp from REF")
    p.add_argument("--seed", type=int, default=3)
    args = p.parse_args()

    W = Path(args.work)
    sup = {}
    with open(args.support) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            sup[f[0]] = (int(f[2]), float(f[3]), [float(x) for x in f[4].split(",")])

    bd = {}
    with gzip.open(W / "results/aardvark-readlik/query.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            v = dict(zip(f[8].split(":"), f[9].split(":"))).get("BD")
            if v in ("TP", "FP"):
                bd[int(f[1])] = v

    q = subprocess.run(["bcftools", "query", "-f", "%POS\t%ID\t%REF\t%ALT[\t%GT\t%GQ\t%DP]\n",
                        str(W / "results/readlik.vcf.gz")], capture_output=True, text=True)
    feats, ys = [], []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 7:
            continue
        lab = bd.get(int(f[0]))
        if lab is None:
            continue
        alts = [a for a in f[3].split(",") if not a.startswith("<") and a != "*"]
        if args.min_svlen and not any(abs(len(a) - len(f[2])) >= args.min_svlen for a in alts):
            continue
        v = sup.get(id_to_site(f[1]))
        if not v:
            continue
        n_reads, mbl, counts = v
        idx = sorted({int(x) for x in f[4].replace("|", "/").split("/") if x.isdigit()})
        tot = sum(counts) or 1.0
        share = sum(counts[i] for i in idx if i < len(counts)) / tot
        balance = 0.5
        if len(idx) == 2 and idx[0] != idx[1] and max(idx) < len(counts):
            a, b = counts[idx[0]], counts[idx[1]]
            if a + b > 0:
                balance = min(a, b) / (a + b)
        gq = int(f[5]) if f[5].isdigit() else 0
        dp = int(f[6]) if f[6].isdigit() else 0
        feats.append([float(gq), share, mbl, balance, float(dp)])
        ys.append(1 if lab == "TP" else 0)

    names = ["GQ", "share", "best_ln", "balance", "depth"]
    npos, nneg = sum(ys), len(ys) - sum(ys)
    print(f"=== {args.label}"
          + (f", calls with a >={args.min_svlen} bp allele" if args.min_svlen else ", all calls")
          + f": {npos:,} TP / {nneg:,} FP ===")
    print(f"  {'signal':<10}{'AUC alone':>11}")
    for k, nm in enumerate(names):
        pos = [feats[i][k] for i in range(len(ys)) if ys[i] == 1]
        neg = [feats[i][k] for i in range(len(ys)) if ys[i] == 0]
        print(f"  {nm:<10}{auc(pos, neg):>11.4f}")

    rng = random.Random(args.seed)
    order = list(range(len(ys)))
    rng.shuffle(order)
    cut = len(order) // 2
    tr, te = order[:cut], order[cut:]
    Xs = standardise(feats)
    print(f"  {'model (held-out)':<26}{'AUC':>8}")
    for subset, nm in [([0], "GQ"), ([0, 2], "GQ+best_ln"), ([0, 1], "GQ+share"),
                       ([0, 3], "GQ+balance"), ([0, 4], "GQ+depth"),
                       ([0, 1, 2, 3, 4], "all five")]:
        Xtr = [[Xs[i][k] for k in subset] for i in tr]
        ytr = [ys[i] for i in tr]
        w = fit_logistic(Xtr, ytr)
        s_pos = [score(w, [Xs[i][k] for k in subset]) for i in te if ys[i] == 1]
        s_neg = [score(w, [Xs[i][k] for k in subset]) for i in te if ys[i] == 0]
        print(f"  {nm:<26}{auc(s_pos, s_neg):>8.4f}")


if __name__ == "__main__":
    main()
