#!/usr/bin/env python3
"""How many "missed" structural variants were actually called, spelled as smaller events?

A tandem repeat expansion is one event to the benchmark and often many events to the
graph: the region is built as a chain of small bubbles, so a 300 bp expansion of a 30 bp
period comes out as ten separate ~30 bp insertions. Every one of them is below the 50 bp
cut, so the structural comparison never sees them and the truth record is scored as a
miss. The same events are then false positives against the small-variant benchmark, whose
truth spells the locus as one large insertion. Lost twice, for a representation
disagreement rather than a calling error.

The test needs no re-calling. For each false negative, sum the net length change of every
call within a window and ask whether it matches the truth SVLEN. Reported with two
sensitivities and a chance control, because the loose settings will credit coincidences in
repeat-dense regions and the number is worthless without knowing how often that happens.

This is a **bound, not a metric**. Matching net length change is much weaker than matching
sequence; `truvari refine` does the latter properly by re-aligning the region with MAFFT,
and agrees closely (see docs/tier2-sv-errors.md). Do not quote this in place of refine.

It is also the mirror of the same-length-record finding in hap32_precision.py: there the
caller merges too much, one record for a handful of scattered SNVs; here it splits too
much, many records for one repeat expansion. Both are the same disagreement about
decomposition, and record-by-record matching punishes whichever direction it goes.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import gzip
import json
import random
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from hap32_precision import DATASETS, WORK, records  # noqa: E402


def load_calls(vcf: Path) -> dict:
    """chrom -> sorted [(pos, net length change of the called allele)], non-reference only."""
    by = collections.defaultdict(list)
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            alts = f[4].split(",")
            sample = dict(zip(f[8].split(":"), f[9].split(":")))
            idx = [int(t) for t in sample.get("GT", "").replace("|", "/").split("/")
                   if t.isdigit()]
            if not any(i > 0 for i in idx):
                continue
            best = 0
            for i in set(idx):
                if 0 < i <= len(alts):
                    d = len(alts[i - 1]) - len(f[3])
                    if abs(d) > abs(best):
                        best = d
            if best:
                by[f[0]].append((int(f[1]), best))
    for k in by:
        by[k].sort()
    return by


def near(calls: dict, chrom: str, pos: int, window: int) -> list:
    arr = calls.get(chrom, [])
    lo = bisect.bisect_left(arr, (pos - window,))
    hi = bisect.bisect_right(arr, (pos + window,))
    return arr[lo:hi]


def explained(deltas: list, svlen: int, tol: float) -> bool:
    """Do the same-direction calls sum to the truth length change, within tol?"""
    same = [d for _, d in deltas if (d > 0) == (svlen > 0)]
    if not same or svlen == 0:
        return False
    return abs(sum(same) - svlen) / abs(svlen) <= tol


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="readlik")
    ap.add_argument("--windows", nargs="+", type=int, default=[50, 200, 500])
    ap.add_argument("--tolerances", nargs="+", type=float, default=[0.1, 0.2, 0.3])
    args = ap.parse_args()

    out = {}
    for ds, sub in DATASETS:
        res = WORK / sub / "results"
        calls = load_calls(res / f"{args.arm}.vcf.gz")
        fn = records(res / f"truvari-{args.arm}" / "fn.vcf.gz")

        print(f"\n=== {ds}: {len(fn)} false negatives ===")
        print(f"  {'window':>7s} " + " ".join(f"{'tol ' + str(t):>9s}" for t in args.tolerances))
        grid = {}
        for w in args.windows:
            cells = []
            for tol in args.tolerances:
                n = sum(1 for r in fn
                        if explained(near(calls, r["chrom"], r["pos"], w), r["svlen"], tol))
                grid[f"w{w}_tol{tol}"] = n
                cells.append(f"{n / len(fn) * 100:8.1f}%")
            print(f"  {w:>7d} " + " ".join(cells))

        # Chance control: the same test at positions drawn from the call set itself, so the
        # local density of calls is preserved and only the truth position is wrong.
        rng = random.Random(0)
        w, tol = args.windows[len(args.windows) // 2], args.tolerances[-1]
        hits = 0
        for r in fn:
            arr = calls.get(r["chrom"], [])
            if not arr:
                continue
            p = rng.choice(arr)[0]
            if explained(near(calls, r["chrom"], p, w), r["svlen"], tol):
                hits += 1
        print(f"  chance rate at window {w}, tol {tol}: {hits / len(fn) * 100:.1f}% "
              f"against {grid[f'w{w}_tol{tol}'] / len(fn) * 100:.1f}% observed")

        # How many records does an explained event take?
        parts = []
        for r in fn:
            d = near(calls, r["chrom"], r["pos"], w)
            if explained(d, r["svlen"], tol):
                parts.append(sum(1 for _, x in d if (x > 0) == (r["svlen"] > 0)))
        if parts:
            print(f"  explained events use a median of {statistics.median(parts):.0f} records; "
                  f"{sum(1 for p in parts if p > 1) / len(parts) * 100:.0f}% use more than one")
        out[ds] = {"n_fn": len(fn), "grid": grid, "chance": hits,
                   "median_parts": statistics.median(parts) if parts else None}

    dest = WORK / "sv-atlas" / f"fn-decomposition-{args.arm}.json"
    dest.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
