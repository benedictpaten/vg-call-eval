#!/usr/bin/env python3
"""Does read *coverage* -- how many reads a call has, and how they split -- rank calls
better than GQ alone?

The genotype likelihood asks only which alleles a read fits. Its mixture assumes each
haplotype in the genotype contributed exactly 1/|G| of the reads; there is no sampling
model behind that split and no adjustment for how long the alleles are. So a het call whose
reads divide 40/2 scores the same as one that divides 21/21, and a call with 8 reads in a
50x neighbourhood scores the same as one with 50. Both look like defects the caller could
in principle notice.

Two related questions, kept separate on purpose:

  1. **Skew.** Under a binomial, how surprising is the observed split? The expectation is
     0.5, flat, with no length term -- and that is a measured result, not an assumption.
     The obvious geometric prior says otherwise: with 151 bp reads an allele carrying a
     120 bp insertion has almost no read that can span it and anchor on both sides, so
     `w_i = max(1, R - len_i - 2F)` and `p_i = w_i / sum_j w_j`. Scored that way, true het
     SV calls look impossible, and the resulting statistic ranks *worse* than chance
     (AUC 0.46 on the chr20 SV set). `allele_balance_by_length.py` shows why: among true
     ref/het calls the ALT share is 0.499-0.545 in every length class from 0 bp to over
     1 kb. The caller never required a read to span an allele -- a read covering one
     breakpoint still fits one allele better than the other -- so the geometry does not
     apply. The flat null is the right one.

  2. **Explained fraction.** What the same measurement *does* show moving with length is
     `sum(AD)/DP`: 1.00 at 0 bp, 0.83 by 256-1000 bp. At a site with long alleles more
     reads fit neither called allele best. So raw `share` is not comparable between a SNV
     and a 300 bp insertion, and is used here both raw and as a residual against the
     length-class mean among training-set true calls.

  3. **Depth.** Is this site's read count low for its neighbourhood? Compared against a
     rolling median over nearby calls rather than a global mean, because coverage varies
     smoothly along a chromosome and a global comparison mostly measures position.

Both are scored against GQ, not in isolation -- a weak signal usually correlates with
correctness, and the only interesting question is whether it adds anything to the ranking
already available. The final table is the one that matters: precision at matched recall,
which is what "a better PR curve" actually means.

Reads AD/BL straight out of the VCF (vg emits them as FORMAT fields), so no likelihood dump
is needed.
"""

from __future__ import annotations

import argparse
import gzip
import math
import random
import subprocess
from pathlib import Path

READ_LEN = 151.0   # measured from the GAF; these are Illumina reads
FLANK = 10.0       # bp of anchor each side before an alignment is informative


# ---------------------------------------------------------------- statistics


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


def log_binom_cdf(k: int, n: int, p: float) -> float:
    """log P(X <= k) for X ~ Binom(n, p). Direct summation; n here is read depth, so
    a few hundred at most."""
    if n <= 0:
        return 0.0
    p = min(max(p, 1e-9), 1 - 1e-9)
    k = min(max(k, 0), n)
    terms = []
    for i in range(k + 1):
        terms.append(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                     + i * math.log(p) + (n - i) * math.log1p(-p))
    m = max(terms)
    return m + math.log(sum(math.exp(t - m) for t in terms))


def skew_logp(counts: list[int], ps: list[float]) -> float:
    """How unlikely is the *least* supported called allele under its own expectation.

    One-sided on purpose: an allele with more reads than expected is not evidence against
    the call, but one with fewer is exactly the dropout case of interest.
    """
    n = sum(counts)
    if n <= 0 or len(counts) < 2:
        return 0.0
    return min(log_binom_cdf(c, n, p) for c, p in zip(counts, ps))


def rolling_median_depth(pos: list[int], dp: list[float], window: int = 201) -> list[float]:
    """Median depth over the `window` nearest calls by position. Calls arrive sorted."""
    n = len(dp)
    half = window // 2
    out = []
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        w = sorted(dp[lo:hi])
        out.append(w[len(w) // 2] if w else 0.0)
    return out


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


def fit_logistic(X, y, epochs: int = 400, lr: float = 0.5) -> list[float]:
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


def precision_at_recall(scores: list[float], ys: list[int], targets) -> dict:
    """Rank by score descending, walk down, report precision when recall of the labelled
    positives first reaches each target. Recall is over the calls in this set, so 1.0 means
    "keep everything" -- the point is what precision costs at each cut, not absolute recall."""
    order = sorted(range(len(ys)), key=lambda i: -scores[i])
    total_pos = sum(ys)
    out = {}
    tp = 0
    ti = 0
    tg = sorted(targets)  # ascending: each target is recorded the first time it is met
    for rank, i in enumerate(order, 1):
        tp += ys[i]
        while ti < len(tg) and total_pos and tp / total_pos >= tg[ti]:
            out[tg[ti]] = (tp / rank, rank - tp)
            ti += 1
        if ti >= len(tg):
            break
    return out


# ---------------------------------------------------------------- data


LEN_BUCKETS = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 1024]


def len_bucket(d: int) -> int:
    """Index of the octave-ish length class a call belongs to. Coarse on purpose: the
    quantity being estimated per class is a mean over calls, and the SV classes have only
    tens of members even on a whole chromosome."""
    b = 0
    for i, lo in enumerate(LEN_BUCKETS):
        if d >= lo:
            b = i
    return b


def aardvark_labels(work: Path) -> dict:
    bd = {}
    with gzip.open(work / "results/aardvark-readlik-z/query.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            v = dict(zip(f[8].split(":"), f[9].split(":"))).get("BD")
            if v in ("TP", "FP"):
                bd[int(f[1])] = v
    return bd


def truvari_labels(work: Path) -> dict:
    """TP/FP for the >=50 bp calls, from truvari's own query-side output.

    Necessary because aardvark is scored against the small-variant benchmark, which has no
    record at all above 50 bp; restricting *its* labels by size leaves a handful of calls
    whose labels mean nothing. Truvari is scored against the structural benchmark, so its
    labels are the only meaningful ones for SVs.
    """
    bd = {}
    d = work / "results/truvari-readlik-z"
    for fname, lab in (("tp-comp.vcf.gz", "TP"), ("fp.vcf.gz", "FP")):
        path = d / fname
        if not path.exists():
            raise SystemExit(f"missing {path} -- run truvari_sv.py first")
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                bd[int(line.split("\t", 2)[1])] = lab
    return bd


def load(work: Path, min_svlen: int, labels: str = "aardvark"):
    bd = aardvark_labels(work) if labels == "aardvark" else truvari_labels(work)

    q = subprocess.run(
        ["bcftools", "query",
         "-f", "%POS\t%REF\t%ALT[\t%GT\t%GQ\t%DP\t%AD\t%BL]\n",
         str(work / "results/readlik-z.vcf.gz")],
        capture_output=True, text=True)
    if q.returncode != 0:
        raise SystemExit(q.stderr.strip()[:400])

    rows = []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue
        pos = int(f[0])
        lab = bd.get(pos)
        if lab is None:
            continue
        ref, alt = f[1], f[2]
        alts = alt.split(",")
        if min_svlen and not any(a[0] not in "<*" and abs(len(a) - len(ref)) >= min_svlen
                                 for a in alts):
            continue
        gt = [int(x) for x in f[3].replace("|", "/").split("/") if x.isdigit()]
        if not gt:
            continue
        gq = float(f[4]) if f[4].replace(".", "").isdigit() else 0.0
        dp = float(f[5]) if f[5].replace(".", "").isdigit() else 0.0
        ad = [int(x) for x in f[6].split(",") if x.lstrip("-").isdigit()]
        try:
            bl = float(f[7])
        except ValueError:
            bl = 0.0
        if not ad or dp <= 0:
            continue
        lens = [len(ref)] + [len(a) if a[0] not in "<*" else len(ref) for a in alts]
        rows.append((pos, lab, gt, gq, dp, ad, bl, lens))
    return rows


def featurise(rows, train_idx=None):
    """Feature matrix, plus labels.

    `train_idx` selects the calls allowed to inform the length-class share baseline. It
    must be the training half: the baseline is estimated from true calls, so fitting it on
    everything would leak the label into a feature and inflate the held-out numbers.
    """
    rows.sort(key=lambda r: r[0])
    med = rolling_median_depth([r[0] for r in rows], [r[4] for r in rows])

    raw = []
    for (pos, lab, gt, gq, dp, ad, bl, lens), m in zip(rows, med):
        called = sorted({i for i in gt if i < len(ad)})
        share = sum(ad[i] for i in called) / dp
        counts = [ad[i] for i in called]
        # Flat null: measured to be right, see the module docstring.
        skew = skew_logp(counts, [1.0 / len(counts)] * len(counts)) if len(called) == 2 \
            else 0.0
        depth_ratio = dp / m if m > 0 else 1.0
        d = max((abs(lens[i] - lens[0]) for i in called), default=0)
        raw.append(([gq, share, bl, skew, depth_ratio], 1 if lab == "TP" else 0,
                    len_bucket(d)))

    # Length-class baseline for `share`, from training true calls only.
    allowed = set(train_idx) if train_idx is not None else set(range(len(raw)))
    acc = {}
    for i, (fv, y, b) in enumerate(raw):
        if y == 1 and i in allowed:
            s, n = acc.get(b, (0.0, 0))
            acc[b] = (s + fv[1], n + 1)
    base = {b: s / n for b, (s, n) in acc.items() if n >= 20}
    overall = (sum(s for s, _ in acc.values()) / sum(n for _, n in acc.values())
               if acc else 1.0)

    feats = [fv + [fv[1] - base.get(b, overall)] for fv, _, b in raw]
    ys = [y for _, y, _ in raw]
    return feats, ys


NAMES = ["GQ", "share", "best_ln", "skew", "depth_ratio", "share_resid"]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--min-svlen", type=int, default=0)
    p.add_argument("--labels", choices=["aardvark", "truvari"], default="aardvark",
                   help="aardvark scores against the small-variant benchmark; use truvari "
                        "for anything size-restricted above 50 bp, where aardvark has no "
                        "truth to score against")
    p.add_argument("--seed", type=int, default=3)
    args = p.parse_args()

    rows = load(Path(args.work), args.min_svlen, args.labels)
    if len(rows) < 100:
        raise SystemExit(f"only {len(rows)} labelled calls -- nothing to fit")

    # Split first: `share_resid` is calibrated on training true calls, so the split has to
    # exist before the features do. Sort here rather than inside featurise, so the indices
    # the split hands over refer to the same rows the features come back in.
    rows.sort(key=lambda r: r[0])
    rng = random.Random(args.seed)
    order = list(range(len(rows)))
    rng.shuffle(order)
    cut = len(order) // 2
    tr, te = order[:cut], order[cut:]
    feats, ys = featurise(rows, train_idx=tr)

    npos, nneg = sum(ys), len(ys) - sum(ys)
    scope = f"calls with a >={args.min_svlen} bp allele" if args.min_svlen else "all calls"
    print(f"=== {args.label}, {scope}, {args.labels} labels: "
          f"{npos:,} TP / {nneg:,} FP ===")
    print(f"  {'signal':<13}{'AUC alone':>11}")
    for k, nm in enumerate(NAMES):
        pos = [feats[i][k] for i in range(len(ys)) if ys[i] == 1]
        neg = [feats[i][k] for i in range(len(ys)) if ys[i] == 0]
        print(f"  {nm:<13}{auc(pos, neg):>11.4f}")

    Xs = standardise(feats)

    models = [([0], "GQ"),
              ([0, 3], "GQ+skew"),
              ([0, 4], "GQ+depth"),
              ([0, 1], "GQ+share"),
              ([0, 5], "GQ+share_resid"),
              ([0, 2], "GQ+best_ln"),
              ([0, 2, 5], "GQ+best_ln+share_resid"),
              ([0, 1, 2, 3, 4, 5], "all six")]

    # Precision alone hides the size of the effect when FPs are 3% of calls, so the FP
    # count that survives each cut is reported beside it -- that is the number a filter
    # would actually be judged on.
    # 0.70 and 0.80 are here because the pile-up guard (pileup_guard.py) lands near 0.72
    # recall on the SV sets, and a guard can only be judged against ranking at the recall
    # it actually achieves.
    targets = [0.70, 0.80, 0.90, 0.95, 0.99]
    print(f"\n  {'model (held-out)':<22}{'AUC':>8}" +
          "".join(f"{'P@R'+str(t):>9}{'FP':>7}" for t in targets))
    for subset, nm in models:
        Xtr = [[Xs[i][k] for k in subset] for i in tr]
        ytr = [ys[i] for i in tr]
        w = fit_logistic(Xtr, ytr)
        s_te = [score(w, [Xs[i][k] for k in subset]) for i in te]
        y_te = [ys[i] for i in te]
        a = auc([s for s, y in zip(s_te, y_te) if y == 1],
                [s for s, y in zip(s_te, y_te) if y == 0])
        pr = precision_at_recall(s_te, y_te, targets)
        cells = ""
        for t in targets:
            p_, fp_ = pr.get(t, (float("nan"), -1))
            cells += f"{p_:>9.4f}{fp_:>7d}"
        print(f"  {nm:<22}{a:>8.4f}" + cells)


if __name__ == "__main__":
    main()
