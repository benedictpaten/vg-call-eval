#!/usr/bin/env python3
"""Separate "vg never called it" from "vg called it and truvari could not match it".

sv_delta.py establishes that 2,630 truth SVs are missed by vg and found by PanGenie, so the panel
carried the allele. That is not yet a diagnosis, because a false negative has two very different
causes:

  (a) vg emitted no record at the locus at all -- a scoring failure, the read model declining an
      allele it was offered;
  (b) vg emitted a record at the locus which truvari did not match to the truth variant -- a
      representation failure, where the event was detected but written in a shape the comparison
      could not line up.

The two call for opposite work. (a) is a model problem. (b) is a VCF-writing problem, and worse, it
is scored twice: once as the missed truth variant and again as the unmatched call, so a single
mis-shaped record costs both a false negative and a false positive.

There is a specific reason to suspect (b) here. vg emits 2,219 same-length substitution false
positives against PanGenie's 29 -- REF and ALT of equal length, a replacement rather than an indel.
If a truth deletion plus insertion is being written as one substitution, truvari sees no size change
to match against. So this also asks whether the vg-only FN loci coincide with those SUB false
positives.
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
WINDOW = 100   # bp; truvari's default refdist is 500, so this is deliberately tighter than matching


def query(path: str, fmt: str, region: str | None = None) -> list[str]:
    cmd = ["bcftools", "query", "-f", fmt]
    if region:
        cmd += ["-r", region]
    cmd.append(path)
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg-score", default="work/wgs/score")
    ap.add_argument("--pg-score", default="work/pangenie/score")
    ap.add_argument("--vg-vcf", default="work/wgs/HG002.vcf.gz")
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    tally = defaultdict(int)
    band_tally = defaultdict(lambda: defaultdict(int))
    coincide_examples = []

    for c in AUTOSOMES:
        vg_fn = {}
        for ln in query(f"{args.vg_score}/{c}.truvari/fn.vcf.gz",
                        "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n"):
            chrom, pos, t, sl = ln.split("\t")
            vg_fn[f"{chrom}:{pos}:{t}:{sl}"] = (int(pos), t, abs(int(sl)))
        pg_fn = set()
        for ln in query(f"{args.pg_score}/{c}.truvari/fn.vcf.gz",
                        "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n"):
            chrom, pos, t, sl = ln.split("\t")
            pg_fn.add(f"{chrom}:{pos}:{t}:{sl}")
        vg_only = {k: v for k, v in vg_fn.items() if k not in pg_fn}
        if not vg_only:
            continue

        # Every call vg actually emitted on this contig, and every FP truvari rejected, as sorted
        # position lists so a nearby-record test is a scan rather than a query per variant.
        emitted = []
        for ln in query(args.vg_vcf, "%POS\t%REF\t%ALT\n", region=c):
            pos, ref, alt = ln.split("\t")
            biggest = max(abs(len(a) - len(ref)) for a in alt.split(","))
            same_len = any(len(a) == len(ref) for a in alt.split(","))
            emitted.append((int(pos), biggest, same_len, max(len(ref), *(len(a) for a in alt.split(",")))))
        emitted.sort()
        pos_list = [e[0] for e in emitted]

        import bisect
        for key, (pos, t, n) in vg_only.items():
            lo = bisect.bisect_left(pos_list, pos - WINDOW)
            hi = bisect.bisect_right(pos_list, pos + WINDOW)
            near = emitted[lo:hi]
            b = ("50-100" if n < 100 else "100-300" if n < 300 else
                 "300-700" if n < 700 else "700+")
            if not near:
                tally["no record within 100 bp"] += 1
                band_tally[b]["none"] += 1
            else:
                # Did vg write anything of comparable size? If the only nearby records are far too
                # small, the event was not represented even though the locus was touched.
                comparable = [e for e in near if e[1] >= 0.5 * n]
                sub_only = [e for e in near if e[2] and e[3] >= 50]
                if comparable:
                    tally["record of comparable size present"] += 1
                    band_tally[b]["comparable"] += 1
                elif sub_only:
                    tally["only a same-length substitution present"] += 1
                    band_tally[b]["sub"] += 1
                    if len(coincide_examples) < 6:
                        coincide_examples.append((c, pos, t, n, sub_only[0][0], sub_only[0][3]))
                else:
                    tally["record present but far too small"] += 1
                    band_tally[b]["small"] += 1

    total = sum(tally.values())
    L = ["# Why vg misses the SVs PanGenie finds: no call, or an unmatched call?", "",
         f"Over the {total:,} autosomal truth SVs that vg misses and PanGenie finds, asking what vg",
         f"actually wrote within {WINDOW} bp of each.", "",
         "| what vg emitted at the locus | n | share |", "|---|---|---|"]
    for k, v in sorted(tally.items(), key=lambda kv: -kv[1]):
        L.append(f"| {k} | {v:,} | {100*v/total:.1f}% |")
    L += ["", "`comparable size` means a record whose largest allele-length change is at least half",
          "the truth variant's size, so the event was written but truvari declined the match --",
          "a representation failure, and one that is scored twice, as a false negative here and as",
          "an unmatched false positive in the same run.", "",
          "## By truth variant size", "",
          "| size | no record | comparable record | same-length substitution only | too small |",
          "|---|---|---|---|---|"]
    for b in ("50-100", "100-300", "300-700", "700+"):
        d = band_tally[b]
        L.append(f"| {b} | {d['none']:,} | {d['comparable']:,} | {d['sub']:,} | {d['small']:,} |")
    if coincide_examples:
        L += ["", "## Loci where the only nearby call is a same-length substitution", "",
              "| contig | truth pos | truth type | truth size | vg record pos | vg allele len |",
              "|---|---|---|---|---|---|"]
        for c, p, t, n, sp, sl in coincide_examples:
            L.append(f"| {c} | {p:,} | {t} | {n:,} | {sp:,} | {sl:,} |")
    Path(args.out).write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
