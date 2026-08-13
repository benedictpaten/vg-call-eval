#!/usr/bin/env python3
"""Two-sided depth filter: drop calls with too few or too many reads.

Motivated by a sign reversal that a single linear depth term cannot express. In one
dataset, local depth ratio separates small-variant calls with AUC 0.65 (low depth means
false) and structural calls with AUC 0.37 (high depth means false). A min/max pair
expresses both at once, which a signed coefficient cannot.

Swept two ways, because it is not obvious which the threshold should be in:

  - **ratio**, DP over a rolling median of nearby calls. Robust to coverage varying along
    a chromosome, and to a different sequencing depth entirely.
  - **absolute** DP, which is what a user would reach for and what §5.3.3 tested.

Recall is computed against the truth total from the *unfiltered* run, so the denominator
does not move when the filter removes calls. On the truvari side, removing a query true
positive is assumed to remove one truth-side match; that is exact when matching is 1:1 and
slightly pessimistic when it is not.

Scored on both benchmarks at once. A depth cut that improves SVs by damaging small variants
is not a filter, it is a trade -- and small variants outnumber SVs several hundred to one.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from filter_lib import DATASETS, collect, prf, truth_counts


def evaluate(rows, counts, lo, hi, key):
    base_tp, total = counts
    kept = [r for r in rows if lo <= r[key] <= hi]
    tp = sum(1 for r in kept if r["label"] == "TP")
    fp = len(kept) - tp
    tp_all = sum(1 for r in rows if r["label"] == "TP")
    return prf(tp, fp, tp_all, base_tp, total), tp


def gq_precision_at(rows, tp_target):
    """Precision that plain GQ ranking reaches while retaining tp_target true calls.

    The only test a hard filter has to pass. Any filter that discards calls buys precision
    by giving up recall, and so does a GQ threshold; the question is never whether the
    filter improves precision but whether it improves it *more than simply lowering the GQ
    cut to the same recall would have*. The pile-up guard failed exactly here.
    """
    order = sorted(range(len(rows)), key=lambda i: -rows[i]["gq"])
    tp = 0
    for rank, i in enumerate(order, 1):
        if rows[i]["label"] == "TP":
            tp += 1
            if tp >= tp_target:
                return tp / rank
    return float("nan")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["ratio", "absolute"], default="ratio")
    p.add_argument("--only", nargs="*", help="restrict to these dataset labels")
    args = p.parse_args()

    key = "ratio" if args.mode == "ratio" else "dp"
    INF = float("inf")
    if args.mode == "ratio":
        # The loose end matters as much as the tight one: §5.3.3's worst pile-ups sit at
        # DP 7,873 against a chromosome median of 29, a ratio near 270. A grid that stops
        # at 3 never sees them.
        grid = [(0.0, INF), (0.5, INF), (0.2, INF), (0.1, INF),
                (0.0, 2.0), (0.5, 2.0), (0.5, 1.5), (0.5, 3.0), (0.0, 3.0),
                (0.0, 5.0), (0.0, 10.0), (0.0, 20.0), (0.0, 50.0),
                (0.2, 10.0), (0.2, 5.0)]
    else:
        grid = [(0, INF), (10, INF), (5, INF), (3, INF),
                (0, 100), (10, 100), (10, 60), (10, 200), (0, 200),
                (0, 150), (0, 300), (0, 600), (0, 1500),
                (5, 300), (5, 150)]

    for label, w in DATASETS:
        if args.only and label not in args.only:
            continue
        W = Path(w)
        print(f"\n=== {label} ({args.mode} depth cut) ===")
        cache = {}
        for kind in ("aardvark", "truvari"):
            # Skip a benchmark this dataset genuinely lacks, but say so -- a silent skip
            # here prints an empty table that looks like a result.
            if not (W / f"results/{'truvari' if kind == 'truvari' else 'aardvark'}"
                    f"-readlik").exists():
                print(f"  (no {kind} output in {W})")
                continue
            cache[kind] = (collect(W, kind), truth_counts(W, kind))
        head = f"  {'min':>6}{'max':>8}"
        for kind in cache:
            k = "smvar" if kind == "aardvark" else "SV"
            head += (f"{k + ' prec':>12}{k + ' rec':>11}{k + ' F1':>11}{'dF1':>9}"
                     f"{'GQ prec':>10}")
        print(head)
        base = {}
        for lo, hi in grid:
            fmt = "{:>6.2f}" if args.mode == "ratio" else "{:>6.0f}"
            line = "  " + fmt.format(lo) + ("{:>8}".format("-") if hi == float("inf")
                                            else fmt.format(hi).rjust(8))
            for kind, (rows, counts) in cache.items():
                (pr, rc, f1), tp_kept = evaluate(rows, counts, lo, hi, key)
                base.setdefault(kind, f1)
                gqp = gq_precision_at(rows, tp_kept)
                line += (f"{pr:>12.4f}{rc:>11.4f}{f1:>11.4f}{f1 - base[kind]:>+9.4f}"
                         f"{gqp:>10.4f}")
            print(line)


if __name__ == "__main__":
    main()
