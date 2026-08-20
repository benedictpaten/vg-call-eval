#!/usr/bin/env python3
"""Assert docs/pangenie-comparison.md still quotes the F1s its arm actually scored.

That page is written by hand, because it interleaves two tools' numbers with prose about what the
difference means -- so nothing regenerates it, and it has drifted twice. There is no sensible way to
generate it, so this checks the other direction: recompute the autosome F1s from the score directory
and require the page to contain each one. Called from scripts/test_harness.sh, which skips it when
the scored arm is not on disk.

Exit 0 if every figure is present, 1 otherwise, naming the ones that are not.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def small(rows, vtype: str) -> float:
    tp = fp = fn = 0
    for r in rows:
        if r["contig"] not in AUTOSOMES:
            continue
        for x in r.get("aardvark") or []:
            if (x.get("comparison") == "GT" and x.get("region_label") == "ALL"
                    and x.get("filter") == "ALL" and x.get("variant_type") == vtype):
                tp += int(x.get("truth_tp", 0) or 0)
                fp += int(x.get("query_fp", 0) or 0)
                fn += int(x.get("truth_fn", 0) or 0)
    return 2 * tp / (2 * tp + fp + fn) if tp else float("nan")


def sv(rows) -> float:
    tp = fp = fn = 0
    for r in rows:
        if r["contig"] not in AUTOSOMES:
            continue
        t = r.get("truvari") or {}
        tp += int(t.get("TP-base", 0) or 0)
        fp += int(t.get("FP", 0) or 0)
        fn += int(t.get("FN", 0) or 0)
    return 2 * tp / (2 * tp + fp + fn) if tp else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score", default="work/wgs-single/score/per-contig.json")
    ap.add_argument("--page", default="docs/pangenie-comparison.md")
    args = ap.parse_args()

    rows = json.loads(Path(args.score).read_text())
    doc = Path(args.page).read_text()
    want = (("ALL", small(rows, "ALL")), ("SNV", small(rows, "Snv")),
            ("Indel", small(rows, "JointIndel")), ("SV >=50 bp", sv(rows)))
    missing = [f"{name} {value:.4f}" for name, value in want if f"{value:.4f}" not in doc]
    for m in missing:
        print(f"  {args.page} does not quote {m}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
