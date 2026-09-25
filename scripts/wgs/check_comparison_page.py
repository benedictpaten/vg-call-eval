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

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import bench_metrics as bm  # noqa: E402

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def main() -> int:
    ap = argparse.ArgumentParser()
    # The arm the page's headline tables are drawn from. The page also quotes older arms as
    # history, so pointing this at one of those passes on a sentence rather than a table.
    ap.add_argument("--score", default="work/wgs-mm095/score/per-contig.json")
    ap.add_argument("--page", default="docs/pangenie-comparison.md")
    args = ap.parse_args()

    rows = json.loads(Path(args.score).read_text())
    doc = Path(args.page).read_text()
    want = (("ALL", bm.small(rows, AUTOSOMES, "ALL").f1),
            ("SNV", bm.small(rows, AUTOSOMES, "Snv").f1),
            ("Indel", bm.small(rows, AUTOSOMES, "JointIndel").f1),
            ("SV >=50 bp", bm.sv(rows, AUTOSOMES).f1))
    missing = [f"{name} {value:.4f}" for name, value in want if f"{value:.4f}" not in doc]
    for m in missing:
        print(f"  {args.page} does not quote {m}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
