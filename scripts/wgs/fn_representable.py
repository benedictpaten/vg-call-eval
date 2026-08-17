#!/usr/bin/env python3
"""Of the SNVs vg misses, how many are not in the graph at all?

small_swallowed.py shows 43.2% of vg's SNV false negatives sit inside a large allele vg emitted --
recall lost to representation. That leaves the complementary question, which bounds how much recall
is reachable by *any* change to the caller: how many of these variants the graph cannot express, so
that no traversal, scope or scoring change could ever produce them.

The instrument is the panel VCF PanGenie was run from. It is a deconstruction of the same graph into
biallelic records, so it is that graph's own allele inventory: a truth variant with no matching
record there is not carried by the graph, and is uncallable by construction for either tool.

Worth noting what that comparison also demonstrates. `vg deconstruct` produces the nested small
variants as separate records from the same graph that `vg call` collapses into whole-snarl alleles.
The information is present; the default traversal scope is what loses it.

SNVs only, deliberately. A single-base substitution has an unambiguous position, so REF/ALT matching
between two independently normalised VCFs is safe. Indels would need left-alignment reconciliation
and the answer would be about normalisation as much as about the graph.

Three-way split of every vg SNV false negative:

  * **not in the panel** -- the graph does not carry the allele. Irreducible.
  * **in the panel, swallowed** -- carried, and vg emitted an allele spanning it, with no record of
    its own. Representation.
  * **in the panel, not swallowed** -- carried, offered as its own site, and still missed. Scoring,
    or a traversal that was never enumerated.
"""

from __future__ import annotations

import argparse
import bisect
import subprocess
from collections import defaultdict
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
LONG = 50


def query(path: str, fmt: str, region: str | None = None) -> list[str]:
    cmd = ["bcftools", "query", "-f", fmt] + (["-r", region] if region else []) + [path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", nargs="*", default=AUTOSOMES)
    ap.add_argument("--vg-score", default="work/wgs/score")
    ap.add_argument("--pg-score", default="work/pangenie/score")
    ap.add_argument("--vg-vcf", default="work/wgs/HG002.vcf.gz")
    ap.add_argument("--panel", required=True, help="the biallelic panel VCF, indexed")
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    counts = defaultdict(int)
    by_pg = defaultdict(lambda: defaultdict(int))

    for c in args.contigs:
        # Truth SNVs vg missed, and what PanGenie made of each.
        fns = []
        for ln in query(f"{args.vg_score}/{c}.aardvark/truth.vcf.gz",
                        "%CHROM\t%POS\t%REF\t%ALT\t[%BD]\n"):
            f = ln.split("\t")
            if len(f) < 5 or f[4] != "FN":
                continue
            ref, alt = f[2], f[3].split(",")[0]
            if len(ref) != 1 or len(alt) != 1:
                continue
            fns.append((int(f[1]), ref, alt))
        if not fns:
            continue
        pgbd = {}
        for ln in query(f"{args.pg_score}/{c}.aardvark/truth.vcf.gz",
                        "%POS\t%REF\t%ALT\t[%BD]\n"):
            f = ln.split("\t")
            if len(f) >= 4:
                pgbd[(int(f[0]), f[1], f[2].split(",")[0])] = f[3]

        # The graph's allele inventory on this contig. Two forms, because the panel does not always
        # express a variant as its own SNV record: an exact (pos, ref, alt) triple, and the spans of
        # longer panel records, since a base substitution can sit inside a multi-base panel allele.
        # Testing only for the exact triple marked 86 chr20 variants as absent from the panel that
        # PanGenie had nonetheless called -- a contradiction, and the reason the span test is here.
        panel = set()
        panel_spans = []
        for ln in query(args.panel, "%POS\t%REF\t%ALT\n", region=c):
            f = ln.split("\t")
            p, ref = int(f[0]), f[1]
            if len(ref) == 1:
                for a in f[2].split(","):
                    if len(a) == 1:
                        panel.add((p, ref, a))
            if len(ref) > 1:
                panel_spans.append((p, p + len(ref) - 1))
        panel_spans.sort()
        pstarts = [s[0] for s in panel_spans]
        plongest = max((e - s + 1 for s, e in panel_spans), default=1)

        def in_panel_span(pos: int) -> bool:
            lo = bisect.bisect_left(pstarts, pos - plongest)
            return any(st <= pos <= en
                       for st, en in panel_spans[lo:bisect.bisect_right(pstarts, pos)])

        # Records vg emitted, for the swallowed test.
        spans = []
        for ln in query(args.vg_vcf, "%POS\t%REF\n", region=c):
            f = ln.split("\t")
            p, L = int(f[0]), len(f[1])
            spans.append((p, p + L - 1, L))
        spans.sort()
        starts = [s[0] for s in spans]
        longest = max((s[2] for s in spans), default=1)

        for pos, ref, alt in fns:
            in_panel = (pos, ref, alt) in panel or in_panel_span(pos)
            lo = bisect.bisect_left(starts, pos - longest)
            swallowed = any(st < pos <= en and L >= LONG
                            for st, en, L in spans[lo:bisect.bisect_right(starts, pos)])
            if not in_panel:
                state = "not in the panel (graph cannot express it)"
            elif swallowed:
                state = "in the panel, swallowed by a large vg allele"
            else:
                state = "in the panel, its own site, still missed"
            counts[state] += 1
            by_pg[state][pgbd.get((pos, ref, alt), "no verdict")] += 1

    total = sum(counts.values()) or 1
    STATES = ["not in the panel (graph cannot express it)",
              "in the panel, swallowed by a large vg allele",
              "in the panel, its own site, still missed"]
    L = ["# How much of vg's SNV recall loss is reachable at all", "",
         f"Every SNV aardvark marks FN for vg, over {len(args.contigs)} contig(s): {total:,} variants.",
         "",
         "The panel VCF is a deconstruction of the same graph vg calls on, so a truth variant absent",
         "from it is not carried by the graph and no caller change could produce it.",
         "",
         "| | n | share |", "|---|---|---|"]
    for s in STATES:
        L.append(f"| {s} | {counts[s]:,} | **{100*counts[s]/total:.1f}%** |")
    L += ["", "## Split by what PanGenie made of the same variant", "",
          "| | PanGenie TP | PanGenie FN | no verdict |", "|---|---|---|---|"]
    for s in STATES:
        d = by_pg[s]
        L.append(f"| {s} | {d.get('TP',0):,} | {d.get('FN',0):,} | {d.get('no verdict',0):,} |")
    L += ["",
          "A variant absent from the panel that PanGenie nonetheless called TP would be a",
          "contradiction, so that cell is a check on the method rather than a result."]
    Path(args.out).write_text("\n".join(L) + "\n")
    for s in STATES:
        print(f"  {s:48s} {counts[s]:>8,}  {100*counts[s]/total:5.1f}%")


if __name__ == "__main__":
    main()
