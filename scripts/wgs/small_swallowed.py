#!/usr/bin/env python3
"""How many small variants does vg lose inside large snarl alleles?

A `SnarlTraversal` spans a snarl from end to end through every interior node, so nested variation is
already baked into each top-level allele. When the caller genotypes a top-level snarl successfully it
emits one record whose ALT is a whole-allele replacement, and `RecurseOnFail` then never descends --
a successful parent call is a dead end regardless of what it decided. A SNP inside that snarl is
therefore either one base of difference inside a long ALT string, or gone.

That prediction is already visible on the structural-variant side: 90.6% of vg's same-length SV false
positives are long REF/ALT pairs differing at ten bases or fewer, some at two or three in an allele
thousands of bases long. Those are nested small variants sized as structural ones.

This asks the converse and more consequential question: of the small variants vg misses, how many sit
inside a long record vg *did* emit? Those are recall lost to representation rather than to evidence --
the allele was called, the variant is inside it, and no record says so.

aardvark writes FORMAT/BD (TP/FP/FN) per truth variant, so the sets come straight out of the existing
scoring. Every vg false negative is counted, and split by what PanGenie made of the same variant,
because that separates two different claims: that swallowing explains the *competitive* gap (what
PanGenie got and vg did not) and that it is a general ceiling on vg's recall (all of them, however
PanGenie fared). Those turn out to be very different rates.

The control matters as much as the measurement. Long records are common in a 32-haplotype graph, so
"most missed variants sit under a long record" means nothing unless correctly-called variants sit
under them less often. Both rates are reported.
"""

from __future__ import annotations

import argparse
import bisect
import subprocess
from collections import defaultdict
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
LONG = 50          # a record at least this long in REF is "a large allele"


def query(path: str, fmt: str, region: str | None = None) -> list[str]:
    # Region-restricted, or a 22-contig run re-reads the whole 4.7M-record genome VCF each time.
    cmd = ["bcftools", "query", "-f", fmt] + (["-r", region] if region else []) + [path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def vtype(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    return "indel"


def load_bd(path: str) -> dict:
    """key -> (BD, type) for every truth record aardvark labelled."""
    out = {}
    for ln in query(path, "%CHROM\t%POS\t%REF\t%ALT\t[%BD]\n"):
        f = ln.split("\t")
        if len(f) < 5:
            continue
        chrom, pos, ref, alt, bd = f[:5]
        first = alt.split(",")[0]
        out[f"{chrom}:{pos}:{ref}:{alt}"] = (bd, vtype(ref, first), int(pos))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", nargs="*", default=AUTOSOMES)
    ap.add_argument("--vg-score", default="work/wgs/score")
    ap.add_argument("--pg-score", default="work/pangenie/score")
    ap.add_argument("--vg-vcf", default="work/wgs/HG002.vcf.gz")
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    # counts[(population, type)][bucket] -> n
    counts = defaultdict(lambda: defaultdict(int))
    spanning_len = defaultdict(list)

    for c in args.contigs:
        vg = load_bd(f"{args.vg_score}/{c}.aardvark/truth.vcf.gz")
        pg = load_bd(f"{args.pg_score}/{c}.aardvark/truth.vcf.gz")
        if not vg:
            continue

        # Every record vg emitted on this contig, as (start, end, reflen). A record's REF spans
        # [POS, POS+len(REF)-1]; a nested variant is "swallowed" if it falls strictly inside one.
        spans = []
        for ln in query(args.vg_vcf, "%POS\t%REF\n", region=c):
            f = ln.split("\t")
            p, L = int(f[0]), len(f[1])
            spans.append((p, p + L - 1, L))
        spans.sort()
        starts = [s[0] for s in spans]
        longest = max((s[2] for s in spans), default=1)

        def swallowing_record(pos: int):
            """The longest vg record strictly containing pos, if any is at least LONG."""
            lo = bisect.bisect_left(starts, pos - longest)
            best = 0
            for st, en, L in spans[lo:bisect.bisect_right(starts, pos)]:
                if st < pos <= en and L >= LONG:
                    best = max(best, L)
            return best

        # Every vg false negative is counted once in "all", and again in whichever sub-population
        # PanGenie's own decision puts it in. Splitting on that separates two different claims: that
        # swallowing explains the *competitive* gap (the variants PanGenie got and vg did not), and
        # that it is a general ceiling on vg's recall (all of them, however PanGenie fared).
        for key, (bd, vt, pos) in vg.items():
            pbd = pg.get(key, (None, None, None))[0]
            if bd == "FN":
                pops = ["vg FN (all)"]
                if pbd == "TP":
                    pops.append("  of which PanGenie called it")
                elif pbd == "FN":
                    pops.append("  of which PanGenie missed it too")
                else:
                    pops.append("  of which PanGenie has no verdict")
            elif bd == "TP":
                pops = ["vg TP (control)"]
            else:
                continue
            L = swallowing_record(pos)
            for pop in pops:
                counts[(pop, vt)]["total"] += 1
                if L:
                    counts[(pop, vt)]["swallowed"] += 1
                    spanning_len[(pop, vt)].append(L)

    L = ["# Small variants vg loses inside large snarl alleles", "",
         f"Generated by `scripts/wgs/small_swallowed.py` over {len(args.contigs)} contig(s).",
         "",
         "A truth variant counts as *swallowed* when a record vg emitted strictly contains its",
         f"position and that record's REF is at least {LONG} bp -- the variant is inside an allele vg",
         "called, with no record of its own.",
         "",
         "`vg FN (all)` is every truth variant aardvark marked FN for vg, split by what PanGenie made of",
         "the same variant. `vg TP` is the control: long records are common on a 32-haplotype graph, so a",
         "swallowed rate only means something against the rate among variants vg got right.",
         "",
         "| population | type | n | inside a large vg allele | rate |",
         "|---|---|---|---|---|"]
    for pop in ("vg FN (all)", "  of which PanGenie called it", "  of which PanGenie missed it too", "  of which PanGenie has no verdict", "vg TP (control)"):
        for vt in ("SNV", "indel"):
            d = counts[(pop, vt)]
            t, s = d["total"], d["swallowed"]
            if not t:
                continue
            L.append(f"| {pop} | {vt} | {t:,} | {s:,} | **{100*s/t:.1f}%** |")
    L += ["", "## Size of the swallowing allele", "",
          "| population | type | median REF bp | 90th percentile | max |", "|---|---|---|---|---|"]
    for pop in ("vg FN (all)", "  of which PanGenie called it", "  of which PanGenie missed it too", "  of which PanGenie has no verdict", "vg TP (control)"):
        for vt in ("SNV", "indel"):
            v = sorted(spanning_len[(pop, vt)])
            if not v:
                continue
            L.append(f"| {pop} | {vt} | {v[len(v)//2]:,} | {v[int(len(v)*0.9)]:,} | {v[-1]:,} |")
    Path(args.out).write_text("\n".join(L) + "\n")
    for pop in ("vg FN (all)", "  of which PanGenie called it", "  of which PanGenie missed it too", "  of which PanGenie has no verdict", "vg TP (control)"):
        for vt in ("SNV", "indel"):
            d = counts[(pop, vt)]
            if d["total"]:
                print(f"  {pop:18s} {vt:6s} {d['swallowed']:>8,} / {d['total']:>9,} "
                      f"= {100*d['swallowed']/d['total']:.1f}%")


if __name__ == "__main__":
    main()
