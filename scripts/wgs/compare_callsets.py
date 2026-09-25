#!/usr/bin/env python3
"""Compare two call sets scored through the identical bench_wgs.py path.

Autosomes are the headline and chrX is reported apart from them, because the two are answering
different questions. On the autosomes both tools are given the same panel, the same reads and the
same truth, so a difference is a difference in how they weigh evidence. chrX is not like that: a
tool that calls a male chrX diploid throughout is wrong there by construction, and folding that
into a genome-wide F1 would report a ploidy-handling difference as though it were an accuracy one.

chrY is excluded from both, identically, for the reference mismatch documented in wgs-results.md.

Counts are summed and the rates recomputed. Averaging per-contig F1s would weight chr21 like chr1.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import bench_metrics as bm  # noqa: E402

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", required=True, help="NAME=path/to/per-contig.json")
    p.add_argument("--b", required=True, help="NAME=path/to/per-contig.json")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    sets = {}
    for spec in (args.a, args.b):
        name, path = spec.split("=", 1)
        sets[name] = json.loads(Path(path).read_text())
    names = list(sets)

    report = {}
    for scope, contigs in (("autosomes (chr1-22)", set(AUTOSOMES)), ("chrX", {"chrX"})):
        print(f"\n== {scope} ==\n")
        # Truth-side TP beside FN, query-side TP beside FP: each is the count its rate is over.
        print(f"{'':10} " + "".join(f"{n:>44}" for n in names))
        print(f"{'':10} " + "".join(f"{'truth TP':>10}{'FN':>8}{'query TP':>10}{'FP':>9}{'F1':>7}"
                                    for _ in names))
        for label, vtype in (("ALL", "ALL"), ("SNV", "Snv"), ("Indel", "JointIndel"), ("SV", None)):
            cells = ""
            for n in names:
                c = bm.small(sets[n], contigs, vtype) if vtype else bm.sv(sets[n], contigs)
                cells += f"{c.truth_tp:10,}{c.truth_fn:8,}{c.query_tp:10,}{c.query_fp:9,}{c.f1:7.4f}"
                report.setdefault(scope, {}).setdefault(label, {})[n] = {
                    "truth_tp": c.truth_tp, "fn": c.truth_fn, "query_tp": c.query_tp,
                    "fp": c.query_fp, "f1": c.f1, "recall": c.recall, "precision": c.precision}
            print(f"{'SV>=50bp' if label == 'SV' else label:10} " + cells)

    print("\nRecall and precision on the autosomes:\n")
    print(f"{'':10} " + "".join(f"{n:>26}" for n in names))
    print(f"{'':10} " + "".join(f"{'recall':>13}{'precision':>13}" for _ in names))
    for label in ("ALL", "SNV", "Indel", "SV"):
        cells = ""
        for n in names:
            d = report["autosomes (chr1-22)"][label][n]
            cells += f"{d['recall']:13.4f}{d['precision']:13.4f}"
        print(f"{label:10} " + cells)

    if args.out:
        Path(args.out).write_text(json.dumps(report, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
