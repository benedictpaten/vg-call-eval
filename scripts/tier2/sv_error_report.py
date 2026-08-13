#!/usr/bin/env python3
"""Generate every table in docs/tier2-sv-errors.md from the atlas and the sensitivity run.

Inputs, both produced by other scripts in this directory:
  work/sv-atlas/{truth,calls}.tsv        sv_error_atlas.py
  work/sv-atlas/metric_sensitivity.json  sv_metric_sensitivity.py

Prints markdown. The page is written by hand around these tables rather than emitted
whole, because the argument is prose; but no number in it is transcribed.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
ATLAS = REPO / "work/sv-atlas"
DS = ["chr6-4hap", "chr6-34hap", "chr20-4hap", "chr20-34hap"]
BINS = ["50-99", "100-299", "300-999", "1k+"]


def load():
    truth = list(csv.DictReader(open(ATLAS / "truth.tsv"), delimiter="\t"))
    calls = list(csv.DictReader(open(ATLAS / "calls.tsv"), delimiter="\t"))
    sens = {}
    p = ATLAS / "metric_sensitivity.json"
    if p.exists():
        sens = json.loads(p.read_text())
    return truth, calls, sens


def is_hom(gt: str) -> bool:
    return gt.replace("|", "/") == "1/1"


def h(*cols):
    print("| " + " | ".join(str(c) for c in cols) + " |")


def rule(n):
    print("|" + "|".join(["---"] * n) + "|")


def recall(rows, **f):
    sel = [r for r in rows if all(r[k] == v for k, v in f.items())]
    tp = sum(1 for r in sel if r["outcome"] == "TP")
    return tp, len(sel)


def sec(title):
    print(f"\n### {title}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.parse_args()
    truth, calls, sens = load()

    sec("Recall by type, size and zygosity")
    for arm in ["poisson-z", "readlik"]:
        print(f"\n**`{arm}`**\n")
        h("dataset", "type", "zyg", *BINS)
        rule(3 + len(BINS))
        for ds in DS:
            for t in ["DEL", "INS"]:
                for z, lab in ((False, "het"), (True, "hom")):
                    cells = []
                    for b in BINS:
                        sel = [r for r in truth if r["dataset"] == ds and r["arm"] == arm
                               and r["svtype"] == t and r["sizebin"] == b
                               and is_hom(r["gt_truth"]) == z]
                        tp = sum(1 for r in sel if r["outcome"] == "TP")
                        cells.append(f"{tp/len(sel):.3f} ({len(sel)})" if sel else "-")
                    h(ds, t, lab, *cells)

    sec("Where the read model's SV recall deficit lives")
    h("dataset", "net deficit", "het DEL >=300", "het DEL <300", "hom DEL", "INS")
    rule(6)
    for ds in DS:
        P = {(r["chrom"], r["pos"], r["svlen"]): r for r in truth
             if r["dataset"] == ds and r["arm"] == "poisson-z"}
        R = {(r["chrom"], r["pos"], r["svlen"]): r for r in truth
             if r["dataset"] == ds and r["arm"] == "readlik"}
        b = collections.Counter()
        for k in P:
            if k not in R:
                continue
            d = (P[k]["outcome"] == "TP") - (R[k]["outcome"] == "TP")
            if not d:
                continue
            r = R[k]
            if r["svtype"] == "DEL" and not is_hom(r["gt_truth"]):
                key = "hetDEL>=300" if r["sizebin"] in ("300-999", "1k+") else "hetDEL<300"
            elif r["svtype"] == "DEL":
                key = "homDEL"
            elif r["svtype"] == "INS":
                key = "INS"
            else:
                key = "other"
            b[key] += d
        h(ds, sum(b.values()), b["hetDEL>=300"], b["hetDEL<300"], b["homDEL"], b["INS"])

    sec("What the caller did at each missed truth SV")
    CL = ["called-large", "split-merge", "fragmented", "small-only", "no-call"]
    h("dataset", "arm", "FN", *CL)
    rule(3 + len(CL))
    fn = [r for r in truth if r["outcome"] == "FN"]
    for ds in DS:
        for arm in ["poisson-z", "readlik"]:
            sel = [r for r in fn if r["dataset"] == ds and r["arm"] == arm]
            c = collections.Counter(r["call_class"] for r in sel)
            h(ds, f"`{arm}`", len(sel),
              *[f"{c[k]} ({c[k]/len(sel)*100:.0f}%)" for k in CL])

    sec("Why truvari rejected each false positive")
    MC = ["placement", "consumed", "dissimilar", "none"]
    h("dataset", "arm", "FP", *MC)
    rule(3 + len(MC))
    for ds in DS:
        for arm in ["poisson-z", "readlik"]:
            sel = [r for r in calls if r["dataset"] == ds and r["arm"] == arm
                   and r["outcome"] == "FP"]
            c = collections.Counter(r["match_class"] for r in sel)
            h(ds, f"`{arm}`", len(sel),
              *[f"{c[k]} ({c[k]/len(sel)*100:.0f}%)" for k in MC])

    sec("Metric sensitivity: refdist sweep and phab refinement")
    if not sens:
        print("_(no metric_sensitivity.json)_")
    else:
        h("dataset", "arm", "setting", "TP-base", "FP", "recall", "precision", "F1")
        rule(8)
        for key in sorted(sens):
            ds, arm, setting = key.split("|")
            s = sens[key]
            h(ds, f"`{arm}`", setting, s["TP-base"], s["FP"],
              f"{s['recall']:.4f}", f"{s['precision']:.4f}", f"{s['f1']:.4f}")

    sec("34-haplotype false positives: carried over or new")
    def idx(ds, arm, out):
        return {(r["chrom"], int(r["pos"]), int(r["svlen"])): r for r in calls
                if r["dataset"] == ds and r["arm"] == arm and r["outcome"] == out}

    def near(k, S, tol=100, ltol=0.2):
        c, p, l = k
        return any(c2 == c and abs(p2 - p) <= tol and l * l2 > 0
                   and abs(abs(l2) - abs(l)) <= max(20, ltol * abs(l))
                   for (c2, p2, l2) in S)

    h("contig", "arm", "4-hap FP", "34-hap FP", "carried over", "new",
      "4-hap FP that vanish", "new with no truth candidate")
    rule(8)
    for contig, a, b_ in [("chr6", "chr6-4hap", "chr6-34hap"),
                          ("chr20", "chr20-4hap", "chr20-34hap")]:
        for arm in ["poisson-z", "readlik"]:
            f4, f34 = idx(a, arm, "FP"), idx(b_, arm, "FP")
            s4, s34 = set(f4), set(f34)
            carried = [k for k in f34 if near(k, s4)]
            new = [k for k in f34 if k not in set(carried)]
            gone = [k for k in f4 if not near(k, s34)]
            newnone = [k for k in new if f34[k]["match_class"] == "none"]
            h(contig, f"`{arm}`", len(f4), len(f34), len(carried), len(new),
              len(gone), len(newnone))

    sec("Support profile of SV calls: true against false")
    h("dataset", "outcome", "n", "median share", "median GQ", "median DP")
    rule(6)
    import statistics as st
    for ds in DS:
        for out in ["TP", "FP"]:
            sel = [r for r in calls if r["dataset"] == ds and r["arm"] == "readlik"
                   and r["outcome"] == out]
            sh = [float(r["share"]) for r in sel if r["share"]]
            gq = [float(r["gq"]) for r in sel if r["gq"]]
            dp = [float(r["dp"]) for r in sel if r["dp"]]
            h(ds, out, len(sel),
              f"{st.median(sh):.3f}" if sh else "n/a",
              f"{st.median(gq):.0f}" if gq else "n/a",
              f"{st.median(dp):.0f}" if dp else "n/a")


if __name__ == "__main__":
    main()
