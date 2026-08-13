#!/usr/bin/env python3
"""Can the explained-read fraction go into GQ itself?

`share` = the fraction of reads whose best-fitting allele is one of the called ones, from
`AD` and `DP`. It is the only signal measured that improved on GQ in all eight cells --
two chromosomes, two graphs, two benchmarks -- so it is the one candidate for changing the
emitted quality rather than leaving it to a downstream filter.

The earlier evidence came from logistic fits with per-dataset weights, which is not
something that can ship. This asks the shippable question instead: is there a **single
fixed formula**, with no fitted parameter, that improves the ranking everywhere? Three
families, of which only the third has a probabilistic reading:

    scale     GQ' = GQ * share^a
    linear    GQ' = GQ - k * (1 - share) * 100
    phred     GQ' = min(GQ, -10 log10(1 - share))

The `phred` form treats the unexplained fraction as an error probability: if a tenth of the
reads at a site fit an allele the call does not contain, the call is not a Q60 call
whatever the likelihood ratio between the top two genotypes says. It needs no tuning, and
it can only ever lower GQ, which matters because GQ is a published field and raising it on
the strength of a heuristic would be worse than leaving it alone.

Scored as a ranking: precision and surviving false calls at matched recall, with recall on
the base side (see filter_lib.truth_counts). No train/test split, because these formulas
have nothing fitted to split over.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from filter_lib import DATASETS, collect, truth_counts  # noqa: E402


def auc(pos, neg):
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
        s += sum(avg for k in range(i, j) if allv[k][1] == 1)
        rank += j - i
        i = j
    return (s - n1 * (n1 + 1) / 2) / (n1 * n0)


MAX_Q = 60.0


def variants(emitted: bool = False):
    """Candidate qualities to rank by.

    With --emitted, compare what the caller actually wrote -- GQI (raw likelihood ratio)
    against GQ (discounted in the caller from the exact fractional support). That is the
    end-to-end check on the shipped implementation, and it is not quite the same
    experiment as the offline forms below, which reconstruct share from the *rounded*
    integer AD and so lose a little at low depth.
    """
    if emitted:
        return [("GQI (raw ratio)", lambda g, s, i: i),
                ("GQ (as emitted)", lambda g, s, i: g)]
    out = [("GQ (unchanged)", lambda g, s: g)]
    for a in (1.0, 2.0, 4.0):
        out.append((f"scale  a={a:g}", lambda g, s, a=a: g * (s ** a)))
    for k in (0.1, 0.25, 0.5):
        out.append((f"linear k={k:g}", lambda g, s, k=k: g - k * (1.0 - s) * 100.0))
    out.append(("phred", lambda g, s: min(g, MAX_Q if s >= 1.0
                                          else -10.0 * math.log10(1.0 - s))))

    # Plain `phred` is too harsh on small variants and the arithmetic says why: at DP 30 a
    # single stray read gives share 0.967, which caps a perfectly good SNV at Q14.8. Allow
    # a tolerance t of unexplained reads before the cap engages, so only the excess counts
    # as evidence against the call.
    def tol(g, s, t):
        excess = (1.0 - s - t) / (1.0 - t)
        if excess <= 0.0:
            return g
        return min(g, -10.0 * math.log10(min(1.0, excess)))

    for t in (0.05, 0.1, 0.2):
        out.append((f"phred t={t:g}", lambda g, s, t=t: tol(g, s, t)))
    return out


def at_recall(scores, ys, base_tp, total, targets):
    """Precision and surviving FP at each base-side recall target."""
    order = sorted(range(len(ys)), key=lambda i: -scores[i])
    tp_all = sum(ys)
    out = {}
    tp = 0
    ti = 0
    tg = sorted(targets)
    for rank, i in enumerate(order, 1):
        tp += ys[i]
        rec = base_tp * (tp / tp_all) / total if tp_all and total else 0.0
        while ti < len(tg) and rec >= tg[ti]:
            out[tg[ti]] = (tp / rank, rank - tp)
            ti += 1
        if ti >= len(tg):
            break
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--only", nargs="*")
    p.add_argument("--emitted", action="store_true",
                   help="compare the caller's own GQ against its GQI, rather than "
                        "reconstructing candidate formulas offline")
    args = p.parse_args()

    for label, w in DATASETS:
        if args.only and label not in args.only:
            continue
        W = Path(w)
        for kind, targets in (("aardvark", [0.90, 0.93]), ("truvari", [0.35, 0.42])):
            if not (W / f"results/{kind if kind == 'truvari' else 'aardvark'}"
                    f"-readlik").exists():
                continue
            rows = collect(W, kind)
            base_tp, total = truth_counts(W, kind)
            ys = [1 if r["label"] == "TP" else 0 for r in rows]
            bench = "small variants" if kind == "aardvark" else "SVs"
            print(f"\n=== {label}, {bench}: {sum(ys):,} TP / {len(ys) - sum(ys):,} FP ===")
            print(f"  {'form':<16}{'AUC':>8}" +
                  "".join(f"{'P@R' + str(t):>10}{'FP':>7}" for t in targets))
            base = None
            for name, fn in variants(args.emitted):
                sc = [fn(r["gq"], r["share"], r["gqi"]) if args.emitted
                      else fn(r["gq"], r["share"]) for r in rows]
                a = auc([s for s, y in zip(sc, ys) if y == 1],
                        [s for s, y in zip(sc, ys) if y == 0])
                pr = at_recall(sc, ys, base_tp, total, targets)
                cells = ""
                for t in targets:
                    v, fp = pr.get(t, (float("nan"), -1))
                    cells += f"{v:>10.4f}{fp:>7d}"
                mark = ""
                if base is None:
                    base = a
                else:
                    mark = f"  {a - base:+.4f}"
                print(f"  {name:<16}{a:>8.4f}" + cells + mark)


if __name__ == "__main__":
    main()
