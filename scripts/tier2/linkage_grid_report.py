#!/usr/bin/env python3
"""Re-render a finished linkage grid from the cached per-point scores, broken out by class.

`linkage_grid.py` prints one summary line per point and writes a single JSON that the next run
overwrites. That is enough to pick a winner and not enough to see *what moved* -- a gain in
overall F1 can be insertions improving while deletions regress, and the summary line cannot show
it. Every point's full score is already on disk as `score-<tag>.json`, so this needs no re-runs.

What can and cannot be broken out. Small variants come from the benchmark's own stratification,
so genotype F1 is available for ALL / SNV / Insertion / Deletion directly. Structural variants
have one overall F1 and a `sv_by_class` table of (true positives, total) by type, size and
zygosity -- counts against the truth set, so **recall only**. There is no per-class false-positive
count, so there is no per-class SV precision or F1, and reporting the recall column as though it
were F1 would flatter any change that trades recall for precision, which is exactly what this
model does.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE.parent.parent / "work"
ATLAS = WORK / "sv-atlas"


def load(tag: str) -> dict | None:
    p = ATLAS / f"score-{tag}.json"
    return json.loads(p.read_text()) if p.exists() else None


def gt_f1(s: dict, vtype: str) -> float | None:
    for r in s.get("smallvar") or []:
        if (r.get("comparison") == "GT" and r.get("region_label") == "ALL"
                and r.get("filter") == "ALL" and r.get("variant_type") == vtype):
            try:
                return float(r["metric_f1"])
            except (KeyError, ValueError):
                return None
    return None


def sv_recall(s: dict, prefix: str) -> tuple[float | None, int]:
    """Pooled recall over every size and zygosity class of one SV type."""
    tp = n = 0
    for key, v in (s.get("sv_by_class") or {}).items():
        if key.startswith(prefix) and v:
            tp += v[0]
            n += v[1]
    return (tp / n if n else None), n


def fmt(v, width=7, prec=4):
    return f"{v:{width}.{prec}f}" if isinstance(v, float) else f"{'-':>{width}}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="chr20-34hap")
    ap.add_argument("--weights", nargs="+", default=["0", "2", "4", "6", "9", "13"])
    ap.add_argument("--freq-priors", nargs="+", default=["0", "0.25", "0.5", "1"])
    ap.add_argument("--scale", default="10000")
    args = ap.parse_args()

    ds = args.dataset
    points = []
    # The control is loaded once, outside the loops: it is inert on both axes, so emitting it per
    # `f` would repeat one row four times and make the table look like the model does nothing.
    if any(float(w) == 0 for w in args.weights):
        s = load(f"{ds}-lgrid-w0")
        if s is not None:
            points.append(("0", "-", s))
    for f in args.freq_priors:
        for w in args.weights:
            if float(w) == 0:
                continue
            tag = f"{ds}-lgrid-w{w}-s{args.scale}-f{f}"
            s = load(tag)
            if s is None:
                print(f"missing: {tag}", file=sys.stderr)
                continue
            points.append((w, f, s))

    base = next((s for w, f, s in points if float(w) == 0), None)

    print(f"\n{ds}: small-variant genotype F1, and structural-variant F1 with per-type recall")
    print("(SV per-type is recall, not F1 -- sv_by_class has no false-positive counts)\n")
    hdr = (f"{'w':>4s} {'f':>5s} | {'smallALL':>8s} {'SNV':>7s} {'smINS':>7s} {'smDEL':>7s} "
           f"| {'SV F1':>7s} {'SVrec':>7s} {'SVprec':>7s} {'INSrec':>7s} {'DELrec':>7s}")
    print(hdr)
    print("-" * len(hdr))
    for w, f, s in points:
        sv = s.get("sv") or {}
        ins, n_ins = sv_recall(s, "INS")
        dele, n_del = sv_recall(s, "DEL")
        lab_f = "-" if float(w) == 0 else f
        print(f"{w:>4s} {lab_f:>5s} | {fmt(gt_f1(s, 'ALL'), 8)} {fmt(gt_f1(s, 'Snv'))} "
              f"{fmt(gt_f1(s, 'Insertion'))} {fmt(gt_f1(s, 'Deletion'))} "
              f"| {fmt(sv.get('f1'))} {fmt(sv.get('recall'))} {fmt(sv.get('precision'))} "
              f"{fmt(ins)} {fmt(dele)}")

    # The frequency prior against its own control, at matched weight. Against the w = 0 baseline
    # every cell is dominated by what the transition did, and `f`'s contribution -- an order of
    # magnitude smaller -- is invisible inside it. Holding `w` fixed is the only way to see it.
    print("\nfreq-prior effect: delta against f = 0 at the same weight")
    print(hdr)
    print("-" * len(hdr))
    by_w = {}
    for w, f, s in points:
        if float(w) != 0:
            by_w.setdefault(w, {})[f] = s
    for w in args.weights:
        if float(w) == 0 or w not in by_w:
            continue
        ref = by_w[w].get("0")
        if ref is None:
            continue
        r_sv = ref.get("sv") or {}
        r_ins, _ = sv_recall(ref, "INS")
        r_del, _ = sv_recall(ref, "DEL")
        for f in args.freq_priors:
            if float(f) == 0 or f not in by_w[w]:
                continue
            s = by_w[w][f]
            sv = s.get("sv") or {}
            ins, _ = sv_recall(s, "INS")
            dele, _ = sv_recall(s, "DEL")

            def d(a, b):
                return a - b if isinstance(a, float) and isinstance(b, float) else None

            print(f"{w:>4s} {f:>5s} | "
                  f"{fmt(d(gt_f1(s, 'ALL'), gt_f1(ref, 'ALL')), 8)} "
                  f"{fmt(d(gt_f1(s, 'Snv'), gt_f1(ref, 'Snv')))} "
                  f"{fmt(d(gt_f1(s, 'Insertion'), gt_f1(ref, 'Insertion')))} "
                  f"{fmt(d(gt_f1(s, 'Deletion'), gt_f1(ref, 'Deletion')))} "
                  f"| {fmt(d(sv.get('f1'), r_sv.get('f1')))} "
                  f"{fmt(d(sv.get('recall'), r_sv.get('recall')))} "
                  f"{fmt(d(sv.get('precision'), r_sv.get('precision')))} "
                  f"{fmt(d(ins, r_ins))} {fmt(d(dele, r_del))}")

    if base:
        print(f"\ndeltas against w = 0 (n = {n_ins} INS, {n_del} DEL in the SV truth set)")
        print(hdr)
        print("-" * len(hdr))
        b_sv = base.get("sv") or {}
        b_ins, _ = sv_recall(base, "INS")
        b_del, _ = sv_recall(base, "DEL")
        for w, f, s in points:
            if float(w) == 0:
                continue
            sv = s.get("sv") or {}
            ins, _ = sv_recall(s, "INS")
            dele, _ = sv_recall(s, "DEL")

            def d(a, b):
                return a - b if isinstance(a, float) and isinstance(b, float) else None

            print(f"{w:>4s} {f:>5s} | "
                  f"{fmt(d(gt_f1(s, 'ALL'), gt_f1(base, 'ALL')), 8)} "
                  f"{fmt(d(gt_f1(s, 'Snv'), gt_f1(base, 'Snv')))} "
                  f"{fmt(d(gt_f1(s, 'Insertion'), gt_f1(base, 'Insertion')))} "
                  f"{fmt(d(gt_f1(s, 'Deletion'), gt_f1(base, 'Deletion')))} "
                  f"| {fmt(d(sv.get('f1'), b_sv.get('f1')))} "
                  f"{fmt(d(sv.get('recall'), b_sv.get('recall')))} "
                  f"{fmt(d(sv.get('precision'), b_sv.get('precision')))} "
                  f"{fmt(d(ins, b_ins))} {fmt(d(dele, b_del))}")


if __name__ == "__main__":
    main()
