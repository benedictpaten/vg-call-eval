#!/usr/bin/env python3
"""Localise the structural-variant gap between vg call and PanGenie.

The headline is that PanGenie leads SV F1 0.5739 to 0.5152 on the autosomes with both more true
calls and fewer false ones. That is a summary, not a diagnosis, and a difference spread evenly over
24k variants would call for different work than one concentrated in a size band or a repeat class.

Both call sets were scored by the same truvari invocations against byte-identical truth, so the
per-contig fn/tp-base/fp VCFs are directly pairable on the truth side: an `fn` record in one and a
`tp-base` record in the other is the *same* truth variant, judged differently. That pairing is what
makes this more than two histograms side by side, and it supports one inference that a single call
set cannot:

    A truth SV that PanGenie called and vg missed proves the panel carried the allele.

So the vg-only FN set is not a panel limitation. It is a set of variants whose alleles were
available and which the read model declined to call, which makes it the actionable population. The
converse does not follow for PanGenie-only FNs from these files alone.

Three things are measured:

  1. FN decomposition -- missed by both, vg only, PanGenie only -- profiled by type, size band and
     truth zygosity. Answers where the recall gap lives.
  2. Genotype-aware scoring. truvari matches on locus, size and sequence, *not* genotype, and it
     annotates the difference as GTMatch. Recomputing with GTMatch==0 required tests whether the
     advantage survives having to get the genotype right -- which matters because the read model's
     known SV defect was a mis-genotyping one that this metric is blind to.
  3. Repeat context, from the truth VCF's own TRF and LCR annotations, for the FN sets and the FP
     sets. Answers whether the gap is a tandem-repeat phenomenon.
"""

from __future__ import annotations

import argparse
import subprocess
from collections import defaultdict
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

# Bands chosen to straddle the ~700 bp break-even the interior-vs-junction analysis identified
# (docs/tier2-sv-errors.md): if the residual defect is still that mechanism, the gap should widen
# above it.
BANDS = [(50, 100), (100, 300), (300, 700), (700, 2000), (2000, 10000), (10000, 10**9)]


def band(n: int) -> str:
    for lo, hi in BANDS:
        if lo <= n < hi:
            return f"{lo}-{hi}" if hi < 10**9 else f"{lo}+"
    return "?"


def f1(tp: int, fp: int, fn: int) -> float:
    return 2 * tp / (2 * tp + fp + fn) if tp else float("nan")


def query(path: Path, fmt: str) -> list[list[str]]:
    """bcftools query, or an empty list if the file is missing or empty."""
    if not path.exists():
        return []
    out = subprocess.run(["bcftools", "query", "-f", fmt, str(path)],
                         capture_output=True, text=True)
    if out.returncode != 0:
        return []
    return [ln.split("\t") for ln in out.stdout.splitlines() if ln]


# Truth-side records carry the annotations we profile on; the call side carries only what the caller
# emitted, so FP profiling uses size and type alone.
TRUTH_FMT = ("%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t%INFO/GTMatch"
             "\t%INFO/TRF\t%INFO/LCR\n")
# Neither call set annotates SVTYPE/SVLEN -- both are sequence-resolved -- so the call side is
# classified from allele lengths, the same way truvari infers them.
CALL_FMT = "%CHROM\t%POS\t%REF\t%ALT\n"


def load_truth(score_dir: Path, contig: str, which: str) -> tuple[dict, int]:
    """(key -> record, raw row count) for one of fn / tp-base.

    The row count is returned separately because a truth variant can appear on more than one row --
    truvari emits a row per match, so a multi-matched variant is repeated -- and the dict collapses
    those. Set arithmetic needs the dict; the headline totals need the rows, and reporting the dict
    size as a TP count would silently undercount against the published figures.
    """
    rows = query(score_dir / f"{contig}.truvari" / f"{which}.vcf.gz", TRUTH_FMT)
    recs = {}
    for chrom, pos, svtype, svlen, gt, gtmatch, trf, lcr in rows:
        try:
            n = abs(int(svlen))
        except ValueError:
            continue
        try:
            lcr_hi = float(lcr) >= 0.9
        except ValueError:
            lcr_hi = False
        recs[f"{chrom}:{pos}:{svtype}:{svlen}"] = {
            "chrom": chrom, "type": svtype, "len": n, "band": band(n),
            # A truth genotype with two distinct alleles is heterozygous. The separator varies.
            "het": len({a for a in gt.replace("|", "/").split("/") if a != "."}) > 1,
            # GTMatch is an allele-count difference: 0 means the genotypes agree. Absent on fn.
            "gt_ok": gtmatch == "0",
            "trf": trf not in (".", ""),
            "lcr": lcr_hi,
        }
    return recs, len(rows)


def classify(ref: str, alt: str) -> tuple[str, int]:
    """Type and size of a sequence-resolved allele.

    A same-length REF/ALT pair is a substitution, not an indel. Its length *change* is zero, so a
    size filter on the change would discard it -- but truvari sizes such a record by its allele
    length and scores it as a structural variant, which is why they show up in fp.vcf.gz at all.
    They are called SUB here and kept, both to reproduce the published FP totals and because they
    are a known non-variant population in their own right (see tier2-sv-errors.md): a third of the
    4-haplotype false calls were placement or bookkeeping artefacts rather than biology.
    """
    d = len(alt) - len(ref)
    if d == 0:
        return "SUB", max(len(ref), len(alt))
    return ("INS" if d > 0 else "DEL"), abs(d)


def load_calls(score_dir: Path, contig: str) -> list[dict]:
    """Every record in fp.vcf.gz, classified. No size filter: truvari already applied its own, and
    re-deriving it here only risks disagreeing with the totals being explained."""
    out = []
    for chrom, pos, ref, alt in query(score_dir / f"{contig}.truvari" / "fp.vcf.gz", CALL_FMT):
        # A multiallelic record's ALT field can hold several alleles; take the largest, which is the
        # one that made the record structural.
        t, n = max((classify(ref, a) for a in alt.split(",")), key=lambda ta: ta[1])
        out.append({"type": t, "len": n, "band": band(n)})
    return out


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f}%" if b else "--"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg", default="work/wgs/score")
    ap.add_argument("--pg", default="work/pangenie/score")
    ap.add_argument("--out", default="docs/sv-delta.md")
    args = ap.parse_args()
    vgd, pgd = Path(args.vg), Path(args.pg)

    vg_fn, pg_fn, vg_tp, pg_tp = {}, {}, {}, {}
    vg_fp, pg_fp = [], []
    rows = defaultdict(int)
    for c in AUTOSOMES:
        for tag, dest, d in (("vg_fn", vg_fn, vgd), ("pg_fn", pg_fn, pgd)):
            r, n = load_truth(d, c, "fn"); dest.update(r); rows[tag] += n
        for tag, dest, d in (("vg_tp", vg_tp, vgd), ("pg_tp", pg_tp, pgd)):
            r, n = load_truth(d, c, "tp-base"); dest.update(r); rows[tag] += n
        vg_fp += load_calls(vgd, c)
        pg_fp += load_calls(pgd, c)

    both = set(vg_fn) & set(pg_fn)
    vg_only = set(vg_fn) - set(pg_fn)      # panel carried it, vg declined: the actionable set
    pg_only = set(pg_fn) - set(vg_fn)

    L = []
    add = L.append
    add("# Where the structural-variant gap against PanGenie actually is")
    add("")
    add("Generated by `scripts/wgs/sv_delta.py`. Autosomes only, SVs >=50 bp, truth T2T-Q100.")
    add("")
    add("## The gap as truvari scores it")
    add("")
    add("| | TP | FP | FN | F1 |")
    add("|---|---|---|---|---|")
    for name, tp, fp, fn in (("vg call", rows["vg_tp"], len(vg_fp), rows["vg_fn"]),
                             ("PanGenie", rows["pg_tp"], len(pg_fp), rows["pg_fn"])):
        add(f"| {name} | {tp:,} | {fp:,} | {fn:,} | **{f1(tp, fp, fn):.4f}** |")
    add("")
    add(f"Distinct truth variants behind those TP rows: {len(vg_tp):,} for vg and {len(pg_tp):,} for")
    add("PanGenie -- truvari emits a row per match, so a multi-matched variant repeats. The set")
    add("arithmetic below is over distinct variants; the totals above are rows, to stay comparable")
    add("with the published figures.")
    add("")

    # --- 1. where the recall gap lives -------------------------------------------------------
    add("## 1. The recall gap, and the part of it that is actionable")
    add("")
    add(f"- missed by **both**: {len(both):,} ({pct(len(both), len(set(vg_fn) | set(pg_fn)))} of all missed)")
    add(f"- **vg only** (PanGenie found it, so the panel carried the allele): {len(vg_only):,}")
    add(f"- **PanGenie only**: {len(pg_only):,}")
    add("")
    add("`vg only` is the population worth working on: the allele was demonstrably available and the")
    add("read model declined it. The same inference does not run the other way from these files.")
    add("")
    add("### vg-only misses by type, size and zygosity")
    add("")
    add("| type | size | truth zygosity | vg-only FN | PanGenie-only FN | net to vg |")
    add("|---|---|---|---|---|---|")
    keys = defaultdict(lambda: [0, 0])
    for k in vg_only:
        r = vg_fn[k]
        keys[(r["type"], r["band"], "het" if r["het"] else "hom")][0] += 1
    for k in pg_only:
        r = pg_fn[k]
        keys[(r["type"], r["band"], "het" if r["het"] else "hom")][1] += 1
    for (t, b, z), (v, p) in sorted(keys.items(), key=lambda kv: kv[1][0] - kv[1][1], reverse=True):
        if v + p < 20:
            continue
        add(f"| {t} | {b} | {z} | {v:,} | {p:,} | {p - v:+,} |")
    add("")
    add("Rows with fewer than 20 variants on both sides are omitted; `net to vg` is what vg would")
    add("gain by matching PanGenie in that cell, so negative means vg is already ahead there.")
    add("")

    # --- 2. genotype-aware -------------------------------------------------------------------
    add("## 2. Does the advantage survive requiring the right genotype?")
    add("")
    add("truvari's default match is locus, size and sequence -- **not** genotype. It records the")
    add("genotype difference as `GTMatch`, so the same comparison can be rescored with a correct")
    add("genotype required. This matters here specifically: the read model's known SV defect was a")
    add("mis-genotyping one, invisible to a genotype-blind metric.")
    add("")
    add("| | TP (locus match) | of those, GT correct | TP (GT required) | F1 (GT required) |")
    add("|---|---|---|---|---|")
    for name, tp, fp, fn in (("vg call", vg_tp, len(vg_fp), rows["vg_fn"]),
                             ("PanGenie", pg_tp, len(pg_fp), rows["pg_fn"])):
        ok = sum(1 for r in tp.values() if r["gt_ok"])
        wrong = len(tp) - ok
        add(f"| {name} | {len(tp):,} | {ok:,} ({pct(ok, len(tp))}) | {ok:,} | "
            f"**{f1(ok, fp + wrong, fn + wrong):.4f}** |")
    add("")
    add("A locus-matched call with the wrong genotype is counted as both a false positive and a")
    add("false negative under the stricter rule, which is what a genotyper getting the copy number")
    add("wrong has actually done.")
    add("")

    # --- 3. repeat context -------------------------------------------------------------------
    add("## 3. Is it a tandem-repeat phenomenon?")
    add("")
    add("From the truth VCF's own annotations, so identical for both tools.")
    add("")
    add("| set | n | in a tandem repeat (TRF) | low-complexity (LCR>=0.9) |")
    add("|---|---|---|---|")
    for label, keyset, src in (("missed by both", both, vg_fn),
                              ("vg-only FN", vg_only, vg_fn),
                              ("PanGenie-only FN", pg_only, pg_fn),
                              ("vg TP", set(vg_tp), vg_tp),
                              ("PanGenie TP", set(pg_tp), pg_tp)):
        n = len(keyset)
        t = sum(1 for k in keyset if src[k]["trf"])
        c = sum(1 for k in keyset if src[k]["lcr"])
        add(f"| {label} | {n:,} | {t:,} ({pct(t, n)}) | {c:,} ({pct(c, n)}) |")
    add("")

    # --- 4. false positives ------------------------------------------------------------------
    add("## 4. The false-positive excess")
    add("")
    add(f"vg emits {len(vg_fp):,} SV false positives against PanGenie's {len(pg_fp):,}, "
        f"a difference of {len(vg_fp) - len(pg_fp):+,}.")
    add("")
    add("| type | size | vg FP | PanGenie FP | excess to vg |")
    add("|---|---|---|---|---|")
    fpk = defaultdict(lambda: [0, 0])
    for r in vg_fp:
        fpk[(r["type"], r["band"])][0] += 1
    for r in pg_fp:
        fpk[(r["type"], r["band"])][1] += 1
    for (t, b), (v, p) in sorted(fpk.items(), key=lambda kv: kv[1][0] - kv[1][1], reverse=True):
        if v + p < 20:
            continue
        add(f"| {t} | {b} | {v:,} | {p:,} | {v - p:+,} |")
    add("")

    # The excess is dominated by one type, so split it out and rescore without it. This is the
    # single most load-bearing number in the comparison, because it decides whether vg's
    # false-positive deficit is about evidence or about how records are written.
    vsub = sum(1 for r in vg_fp if r["type"] == "SUB")
    psub = sum(1 for r in pg_fp if r["type"] == "SUB")
    vind, pind = len(vg_fp) - vsub, len(pg_fp) - psub
    add("### Splitting the excess by type")
    add("")
    add("| | same-length substitutions | genuine INS/DEL | total |")
    add("|---|---|---|---|")
    add(f"| vg call | {vsub:,} | {vind:,} | {len(vg_fp):,} |")
    add(f"| PanGenie | {psub:,} | {pind:,} | {len(pg_fp):,} |")
    add(f"| difference | **{vsub - psub:+,}** | **{vind - pind:+,}** | {len(vg_fp) - len(pg_fp):+,} |")
    add("")
    add(f"**On genuine insertions and deletions vg emits {pind - vind:,} *fewer* false positives than")
    add("PanGenie.** The entire false-positive deficit, and more, is same-length substitutions: REF and")
    add("ALT of equal length, which truvari sizes by allele length and therefore scores as structural")
    add("variants. vg's output is multiallelic and carries these; PanGenie's biallelic-split output")
    add("essentially does not. It is a difference in what gets written, not in what the evidence")
    add("supports.")
    add("")
    add("Rescoring with substitutions excluded from both sides:")
    add("")
    add("| | TP | FP | FN | F1 |")
    add("|---|---|---|---|---|")
    vgf = f1(rows["vg_tp"], vind, rows["vg_fn"])
    pgf = f1(rows["pg_tp"], pind, rows["pg_fn"])
    add(f"| vg call | {rows['vg_tp']:,} | {vind:,} | {rows['vg_fn']:,} | **{vgf:.4f}** |")
    add(f"| PanGenie | {rows['pg_tp']:,} | {pind:,} | {rows['pg_fn']:,} | **{pgf:.4f}** |")
    add("")
    full_gap = f1(rows["pg_tp"], len(pg_fp), rows["pg_fn"]) - f1(rows["vg_tp"], len(vg_fp), rows["vg_fn"])
    add(f"The gap falls from {full_gap:.4f} to {pgf - vgf:.4f}, so **{100 * (1 - (pgf - vgf) / full_gap):.0f}% "
        "of the headline SV F1 gap is this one representation artefact** and the rest is recall.")
    add("")

    Path(args.out).write_text("\n".join(L) + "\n")
    print(f"wrote {args.out}")
    print(f"  missed by both {len(both):,}   vg-only {len(vg_only):,}   pangenie-only {len(pg_only):,}")


if __name__ == "__main__":
    main()
