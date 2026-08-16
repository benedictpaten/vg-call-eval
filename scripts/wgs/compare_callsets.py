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
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]


def f1(tp: int, fp: int, fn: int) -> float:
    return 2 * tp / (2 * tp + fp + fn) if tp else float("nan")


def pick(rows, comparison, vtype):
    for r in rows or []:
        if r.get("comparison", "").upper() == comparison and r.get("variant_type") == vtype:
            return r
    return None


def totals(results, contigs, vtype):
    tp = fp = fn = 0
    for r in results:
        if r["contig"] not in contigs:
            continue
        row = pick(r.get("aardvark"), "GT", vtype)
        if row:
            tp += int(row.get("truth_tp", 0) or 0)
            fp += int(row.get("query_fp", 0) or 0)
            fn += int(row.get("truth_fn", 0) or 0)
    return tp, fp, fn


def sv_totals(results, contigs):
    tp = fp = fn = 0
    for r in results:
        if r["contig"] not in contigs:
            continue
        s = r.get("truvari") or {}
        tp += int(s.get("TP-base", 0) or 0)
        fp += int(s.get("FP", 0) or 0)
        fn += int(s.get("FN", 0) or 0)
    return tp, fp, fn


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
        print(f"{'':10} " + "".join(f"{n:>34}" for n in names))
        print(f"{'':10} " + "".join(f"{'TP':>10}{'FP':>9}{'FN':>8}{'F1':>7}" for _ in names))
        for label, vtype in (("ALL", "ALL"), ("SNV", "Snv"), ("Indel", "JointIndel")):
            cells = ""
            for n in names:
                tp, fp, fn = totals(sets[n], contigs, vtype)
                cells += f"{tp:10,}{fp:9,}{fn:8,}{f1(tp,fp,fn):7.4f}"
                report.setdefault(scope, {}).setdefault(label, {})[n] = {
                    "tp": tp, "fp": fp, "fn": fn, "f1": f1(tp, fp, fn),
                    "recall": tp / (tp + fn) if tp + fn else None,
                    "precision": tp / (tp + fp) if tp + fp else None}
            print(f"{label:10} " + cells)
        cells = ""
        for n in names:
            tp, fp, fn = sv_totals(sets[n], contigs)
            cells += f"{tp:10,}{fp:9,}{fn:8,}{f1(tp,fp,fn):7.4f}"
            report.setdefault(scope, {}).setdefault("SV", {})[n] = {
                "tp": tp, "fp": fp, "fn": fn, "f1": f1(tp, fp, fn)}
        print(f"{'SV>=50bp':10} " + cells)

    print("\nRecall and precision on the autosomes:\n")
    print(f"{'':10} " + "".join(f"{n:>26}" for n in names))
    print(f"{'':10} " + "".join(f"{'recall':>13}{'precision':>13}" for _ in names))
    for label in ("ALL", "SNV", "Indel"):
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
