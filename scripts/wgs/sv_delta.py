#!/usr/bin/env python3
"""Localise the structural-variant gap between vg call and PanGenie.

PanGenie leads SV F1 on the autosomes. That is a summary, not a diagnosis, and a difference spread evenly over 24k variants would call for different
work than one concentrated in a size band or a repeat class.

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
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import bench_metrics as bm  # noqa: E402

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
TRUTH_FMT = ("%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t%INFO/GTMatch"
             "\t%INFO/TRF\t%INFO/LCR\n")
# Neither call set annotates SVTYPE/SVLEN -- both are sequence-resolved -- so the call side is
# classified from allele lengths, the same way truvari infers them.
CALL_FMT = "%CHROM\t%POS\t%REF\t%ALT\n"


def load_truth(score_dir: Path, contig: str, which: str) -> tuple[dict, int]:
    """(key -> record, raw row count) for one of fn / tp-base.

    A truth variant is keyed on (CHROM, POS, REF, ALT). POS, type and length are not an identity: the
    two alleles of a compound het can share all three, and keying on them merged 146 of vg's TP rows
    and 28 of its FN rows into their neighbours. The row count is still returned separately, so a
    genuinely repeated record would show up as rows exceeding the dict rather than vanish.
    """
    rows = query(score_dir / f"{contig}.truvari" / f"{which}.vcf.gz", TRUTH_FMT)
    recs = {}
    for chrom, pos, ref, alt, svtype, svlen, gt, gtmatch, trf, lcr in rows:
        try:
            n = abs(int(svlen))
        except ValueError:
            continue
        try:
            lcr_hi = float(lcr) >= 0.9
        except ValueError:
            lcr_hi = False
        recs[(chrom, pos, ref, alt)] = {
            "chrom": chrom, "type": svtype, "len": n, "band": band(n),
            # A truth genotype with two distinct alleles is heterozygous. The separator varies.
            "het": len({a for a in gt.replace("|", "/").split("/") if a != "."}) > 1,
            # GTMatch is an allele-count difference: 0 means the genotypes agree. Absent on fn.
            "gt_ok": gtmatch == "0",
            "trf": trf not in (".", ""),
            "lcr": lcr_hi,
        }
    return recs, len(rows)


def load_comp(score_dir: Path, contig: str) -> tuple[int, int]:
    """(tp-comp rows, of which genotype-matched). Precision is over calls, so it needs the call
    side's own TP count; the truth side's tp-base count is a different number."""
    rows = query(score_dir / f"{contig}.truvari" / "tp-comp.vcf.gz", "%INFO/GTMatch\n")
    return len(rows), sum(1 for (g,) in rows if g == "0")


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
    ap.add_argument("--vg", default="work/wgs-mm095/score")
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
        for tag, d in (("vg", vgd), ("pg", pgd)):
            n, ok = load_comp(d, c)
            rows[f"{tag}_tpc"] += n
            rows[f"{tag}_tpc_gt"] += ok
            rows[f"{tag}_tp_gt"] += sum(1 for r in query(d / f"{c}.truvari" / "tp-base.vcf.gz",
                                                         "%INFO/GTMatch\n") if r == ["0"])
        vg_fp += load_calls(vgd, c)
        pg_fp += load_calls(pgd, c)

    both = set(vg_fn) & set(pg_fn)
    vg_only = set(vg_fn) - set(pg_fn)      # panel carried it, vg declined: the actionable set
    pg_only = set(pg_fn) - set(vg_fn)

    L = []
    add = L.append
    add("# Where the structural-variant gap against PanGenie actually is")
    add("")
    add(f"Generated by `scripts/wgs/sv_delta.py` from `{args.vg}` and `{args.pg}`. Autosomes only, SVs")
    add("of 50 bp up to truvari's default 50 kb size cap, truth T2T-Q100 v1.1 (GIAB defrabb V0.019")
    add("draft benchmark). The `10000+` size band therefore ends below 50 kb.")
    add("")
    add("## The gap as truvari scores it")
    add("")
    add("| | TP-base | FN | TP-comp | FP | F1 |")
    add("|---|---|---|---|---|---|")
    for name, t in (("vg call", "vg"), ("PanGenie", "pg")):
        fp = len(vg_fp if t == "vg" else pg_fp)
        c = bm.Counts(rows[f"{t}_tp"], rows[f"{t}_fn"], rows[f"{t}_tpc"], fp)
        add(f"| {name} | {c.truth_tp:,} | {c.truth_fn:,} | {c.query_tp:,} | {fp:,} | "
            f"**{c.f1:.4f}** |")
    add("")
    add("Recall is over truth records (TP-base, FN) and precision over calls (TP-comp, FP), as")
    add("truvari's own F1 is; one call can match several truth records and the reverse.")
    add("")
    if len(vg_tp) == rows["vg_tp"] and len(pg_tp) == rows["pg_tp"] and \
            len(vg_fn) == rows["vg_fn"] and len(pg_fn) == rows["pg_fn"]:
        add("Every truth row is a distinct variant keyed on (CHROM, POS, REF, ALT), so the set")
        add("arithmetic below counts the same things as the totals above.")
    else:
        add(f"Distinct truth variants behind those TP rows: {len(vg_tp):,} for vg and {len(pg_tp):,} for")
        add("PanGenie, keyed on (CHROM, POS, REF, ALT). The set arithmetic below is over distinct")
        add("variants; the totals above are rows, to stay comparable with the published figures.")
    add("")

    # --- 1. where the recall gap lives -------------------------------------------------------
    add("## 1. The misses, and the part of them that is actionable")
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
    add("| type | size | truth zygosity | vg-only FN | PanGenie-only FN | vg behind by |")
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
        add(f"| {t} | {b} | {z} | {v:,} | {p:,} | {v - p:+,} |")
    add("")
    add("Rows whose two sides sum to fewer than 20 are omitted. `vg behind by` is vg-only minus")
    add("PanGenie-only misses: the net truth variants vg would gain by matching PanGenie in that cell.")
    add("Negative means vg is already ahead there. Rows run from where vg is furthest behind.")
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
    for name, t, tp, fp in (("vg call", "vg", vg_tp, len(vg_fp)),
                            ("PanGenie", "pg", pg_tp, len(pg_fp))):
        ok = sum(1 for r in tp.values() if r["gt_ok"])
        # Rows, not distinct variants, so the ungated rates are the table's above. A matched pair
        # whose genotypes disagree moves to FN on the truth side and to FP on the call side.
        tb, tc = rows[f"{t}_tp"], rows[f"{t}_tpc"]
        tb_ok, tc_ok = rows[f"{t}_tp_gt"], rows[f"{t}_tpc_gt"]
        c = bm.Counts(tb_ok, rows[f"{t}_fn"] + tb - tb_ok, tc_ok, fp + tc - tc_ok)
        add(f"| {name} | {len(tp):,} | {ok:,} ({pct(ok, len(tp))}) | {ok:,} | **{c.f1:.4f}** |")
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
    add("Rows whose two sides sum to fewer than 20 are omitted, so the rows do not add up to the totals.")
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
    # The sign here has flipped once already, when nested calling removed most of the substitutions,
    # so the sentence is derived rather than asserted. A template that only reads correctly one way
    # round is how a stale conclusion survives a rerun.
    add("Same-length substitutions -- REF and ALT of equal length -- are a representation artefact "
        "rather than an evidence one: truvari sizes such a record by its allele length and so scores "
        "it as structural. vg's output carries them and PanGenie's essentially does not.")
    add("")
    if vind <= pind:
        add(f"**On genuine insertions and deletions vg emits {pind - vind:,} fewer false positives "
            f"than PanGenie**, so the whole of its false-positive deficit, and more, is those "
            f"substitutions.")
    else:
        add(f"**vg emits {vind - pind:,} more genuine insertion and deletion false positives than "
            f"PanGenie**, so the substitutions account for only {vsub - psub:,} of the "
            f"{len(vg_fp) - len(pg_fp):,} excess and the rest is real disagreement about what is "
            f"there.")
    add("")
    add("Rescoring with substitutions excluded from both sides:")
    add("")
    add("| | TP-base | FN | TP-comp | FP | F1 |")
    add("|---|---|---|---|---|---|")

    def score(t, fp):
        return bm.Counts(rows[f"{t}_tp"], rows[f"{t}_fn"], rows[f"{t}_tpc"], fp).f1

    vgf, pgf = score("vg", vind), score("pg", pind)
    for name, t, fp, f in (("vg call", "vg", vind, vgf), ("PanGenie", "pg", pind, pgf)):
        add(f"| {name} | {rows[t + '_tp']:,} | {rows[t + '_fn']:,} | {rows[t + '_tpc']:,} | "
            f"{fp:,} | **{f:.4f}** |")
    add("")
    full_gap = score("pg", len(pg_fp)) - score("vg", len(vg_fp))
    add(f"The gap falls from {full_gap:.4f} to {pgf - vgf:.4f}, so **{100 * (1 - (pgf - vgf) / full_gap):.0f}% "
        "of the headline SV F1 gap is this one representation artefact**.")
    add("")
    # Which side the remainder is on has to be derived: it was recall once, and it is not now.
    rec = {t: rows[f"{t}_tp"] / (rows[f"{t}_tp"] + rows[f"{t}_fn"]) for t in ("vg", "pg")}
    prec = {"vg": rows["vg_tpc"] / (rows["vg_tpc"] + vind), "pg": rows["pg_tpc"] / (rows["pg_tpc"] + pind)}
    behind = [m for m, d in (("recall", rec), ("precision", prec)) if d["vg"] < d["pg"]]
    what = {0: "vg is ahead on both", 1: f"the rest is {behind[0]}", 2: "vg is behind on both"}[len(behind)]
    add(f"With the substitutions excluded, {what}: vg's recall is {rec['vg']:.4f} against PanGenie's "
        f"{rec['pg']:.4f}, and its precision {prec['vg']:.4f} against {prec['pg']:.4f}.")
    add("")

    Path(args.out).write_text("\n".join(L) + "\n")
    print(f"wrote {args.out}")
    print(f"  missed by both {len(both):,}   vg-only {len(vg_only):,}   pangenie-only {len(pg_only):,}")


if __name__ == "__main__":
    main()
