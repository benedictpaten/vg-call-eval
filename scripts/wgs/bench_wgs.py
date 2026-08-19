#!/usr/bin/env python3
"""Benchmark the whole-genome call set, per contig, then aggregate.

Scored contig by contig rather than genome at once. aardvark's whole-genome memory profile is
unknown here and the point of this exercise is that it runs on a laptop, so per-contig scoring is
bounded by construction. The cost is that a genome-wide F1 has to be recomputed from summed
TP/FP/FN rather than read off a tool's output -- exact for counts, which is why the aggregation
sums counts and never averages per-contig F1s.

Reuses the harness's own `aardvark.compare` and truvari invocation rather than reconstructing the
command lines. Those have accumulated specifics -- `--pick ac`, matched `--sizemin/--sizefilt`,
multiallelic splitting before truvari, the sample-name mapping -- and a whole-genome run that
quietly used different ones would not be comparable with any tier-2 number.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))

from vgcalleval.engines import aardvark  # noqa: E402

sys.path.insert(0, str(HERE.parent / "tier2"))
from truvari_sv import split_multiallelic  # noqa: E402

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
ALL_CONTIGS = AUTOSOMES + ["chrX", "chrY"]


def run(cmd, **kw):
    proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        print(" ".join(str(c) for c in cmd), file=sys.stderr)
        print(proc.stderr[-2500:], file=sys.stderr)
        return False
    return True


def score_contig(work: Path, contig: str, sample: str, threads: int, truvari: str) -> dict:
    d = work / contig
    score = work / "score"
    score.mkdir(parents=True, exist_ok=True)
    out = {"contig": contig}

    calls = d / f"{contig}.vcf.gz"
    if not calls.exists():
        out["error"] = "no calls"
        return out

    # aardvark pairs by sample name, and the caller writes whatever -s was given.
    # Freshness, not existence, for the same reason assemble_wgs.sh now checks it: a cached
    # rename or an old aardvark summary will happily answer for a VCF that has since been
    # recalled, and the only trace is a file timestamp.
    def stale(path):
        return not path.exists() or path.stat().st_mtime < calls.stat().st_mtime

    renamed = score / f"{contig}.renamed.vcf.gz"
    if stale(renamed):
        names = score / f"{contig}.sample.txt"
        names.write_text(f"{sample}\n")
        run(["bcftools", "reheader", "-s", names, "-o", renamed, calls])
        run(["bcftools", "index", "-f", "-t", renamed])

    adir = score / f"{contig}.aardvark"
    if stale(adir / "summary.tsv"):
        try:
            aardvark.compare(
                aardvark="aardvark",
                reference=d / f"{contig}.fa",
                truth_vcf=d / f"truth.{contig}.smvar.vcf.gz",
                query_vcf=renamed,
                regions_bed=d / f"truth.{contig}.smvar.bed",
                out_dir=adir,
                truth_sample=sample,
                query_sample=sample,
                label=contig,
                options=aardvark.AardvarkOptions(threads=threads),
            )
        except Exception as exc:   # noqa: BLE001
            out["aardvark_error"] = str(exc)
    if (adir / "summary.tsv").exists():
        out["aardvark"] = aardvark.read_summary(adir)

    # truvari, matching truvari_sv.py: split multiallelics first, matched size bounds, --pick ac.
    tdir = score / f"{contig}.truvari"
    if stale(tdir / "summary.json"):
        ref = d / f"{contig}.fa"
        norm = score / f"{contig}.norm.vcf.gz"
        truth_norm = score / f"{contig}.truth.norm.vcf.gz"
        # The harness's own splitter, not a local bcftools norm. It sorts after splitting, which
        # is not optional: left-aligning a split allele can move it upstream of the record that
        # followed, and tabix then refuses to index. Reimplementing it here reproduced exactly
        # that failure, which the tier-2 version carries a comment about having already hit.
        try:
            for src, dst in ((renamed, norm), (d / f"truth.{contig}.stvar.vcf.gz", truth_norm)):
                split_multiallelic(src, dst, ref)
        except (Exception, SystemExit) as exc:   # noqa: BLE001
            # SystemExit as well as Exception: tier2's run() reports a failed command with
            # sys.exit rather than by raising, so `except Exception` does not catch it and the
            # first draft of this guard let chrY kill the run exactly as before.
            #
            # One contig's failure must not end the run. chrY's did, after 23 contigs had already
            # been scored, and the script died before writing any summary at all.
            out["truvari_error"] = f"normalisation failed: {exc}"
            return out
        if tdir.exists():
            run(["rm", "-rf", tdir])
        ok = run([truvari, "bench", "-b", truth_norm, "-c", norm, "-f", ref,
                  "-o", tdir, "--includebed", d / f"truth.{contig}.stvar.bed",
                  "--sizemin", 50, "--sizefilt", 50,
                  "--bSample", sample, "--cSample", sample, "--pick", "ac"])
        if not ok:
            out["truvari_error"] = "bench failed"
    if (tdir / "summary.json").exists():
        out["truvari"] = json.loads((tdir / "summary.json").read_text())
    return out


def pick(rows, comparison, vtype):
    for r in rows or []:
        if r.get("comparison", "").upper() == comparison and r.get("variant_type") == vtype:
            return r
    return None


def f1(tp, fp, fn):
    if tp + fp == 0 or tp + fn == 0:
        return float("nan")
    prec = tp / (tp + fp)
    rec = tp / (tp + fn)
    return 0.0 if prec + rec == 0 else 2 * prec * rec / (prec + rec)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", default="work/wgs")
    p.add_argument("--out", default="work/wgs/score/wgs-summary.md")
    p.add_argument("--sample", default="HG002")
    p.add_argument("--threads", type=int, default=6)
    # The harness keeps truvari in its own venv; truvari_sv.py defaults to the same path, and a
    # different build would not be comparable with the tier-2 SV numbers.
    p.add_argument("--truvari", default=str(REPO / "work/truvari-venv/bin/truvari"))
    p.add_argument("--contigs", nargs="*", default=ALL_CONTIGS)
    # chrY is called but not scoreable against this truth. The graph's CHM13 chrY path is
    # 57,686,750 bp where the truth's chrY runs past 62,111,784, and the sequences do not
    # correspond at any constant offset -- REF alleles match the graph FASTA at chance level
    # (15/60) whether shifted by a PAR length or not at all. CHM13v2.0's chrY is HG002-derived at
    # 62.46 Mb, and this graph's matches neither it nor GRCh38's 57.23 Mb. So chrY calls are on a
    # different coordinate system than the truth, and scoring them measures the mismatch: aardvark
    # returns recall 0.09 at precision 0.000 against 0.9361 for chrX.
    #
    # Excluded by name rather than dropped by a filter, so the exclusion is visible in the output
    # instead of being a silent hole in a genome-wide number.
    p.add_argument("--unscoreable", nargs="*", default=["chrY"])
    args = p.parse_args()

    work = Path(args.work)
    results = []
    for c in args.contigs:
        print(f"[score] {c}", flush=True)
        r = score_contig(work, c, args.sample, args.threads, args.truvari)
        r["scoreable"] = c not in args.unscoreable
        results.append(r)

    (work / "score" / "per-contig.json").write_text(json.dumps(results, indent=2))

    # Genome-wide totals from summed counts. Averaging per-contig F1s would weight chr21 like
    # chr1; summing the counts is the only aggregation that means anything.
    lines = ["# Whole-genome results: HG002 against T2T-Q100", "",
             "Called per contig on the 34-haplotype HPRC graph, `--read-likelihood` with panel",
             "enumeration, phasing and mosaic on. chrY haploid; chrX haploid outside the",
             "pseudoautosomal regions and diploid inside them, in one run via --ploidy-bed.", "",
             "**chrY is called but excluded from every total below.** The graph's CHM13 chrY path",
             "is 57,686,750 bp where the truth's chrY runs past 62,111,784, and the two do not",
             "correspond at any constant offset -- REF alleles match the graph's own FASTA at",
             "chance level whether shifted by a PAR length or not at all. CHM13v2.0's chrY is",
             "HG002-derived at 62.46 Mb and this graph's matches neither it nor GRCh38's 57.23 Mb.",
             "Scored anyway it returns recall 0.09 at precision 0.000, which measures the",
             "coordinate mismatch and not the caller. The calls remain in the VCF and the mosaic.",
             "",
             # These cross-links live here rather than in the .md because this script rewrites the
             # whole file. Added by hand once, they were silently deleted by the next rescore --
             # which is the failure mode of hand-editing a generated file, and it does not announce
             # itself: the numbers were byte-identical, so the diff looked like a pure deletion of
             # prose nobody had asked about. Anything meant to survive a rescore belongs in here.
             "**How to run this, and how long it takes**: [wgs-performance.md](wgs-performance.md).",
             "**Behaviour across coverage and ploidy**: [coverage.md](coverage.md).",
             "",
             "**Compared against PanGenie on the same graph and reads**: see",
             "[pangenie-comparison.md](pangenie-comparison.md). Briefly, on the autosomes vg is ahead on small",
             "variants (ALL F1 0.9630 against 0.9505, driven by indels) and PanGenie is ahead on structural",
             "variants (0.5739 against 0.5152) -- though 42% of that gap is a representation artefact",
             "rather than an evidence one, and the recall remainder sits at 50-300 bp; see",
             "[sv-delta.md](sv-delta.md).",
             "",
             "**The mosaic** this run also emits: 182,328 segments over 5,041,066 sites, 14 MB.",
             "See wgs-performance.md for why assembling it is not `cat`.",
             "",
             "**Nested calling and phasing are the defaults** as of this run, which is why these",
             "numbers moved: SNV F1 0.9752 -> 0.9833, ALL F1 0.9626 -> 0.9699, SV F1 0.5134 -> 0.5467,",
             "with 59,413 SNV false negatives recovered, at no runtime or memory cost. `--no-nested`",
             "and `--no-phased` restore the old behaviour. See",
             "[nested-calling-design.md](nested-calling-design.md).",
             "",
             "**Two caveats that belong with these numbers.** The gain is a rich-panel effect: on the",
             "4-haplotype tier-2 graphs nested calling is flat to 0.0005 *down* on ALL F1, because a",
             "small panel enumerates few of the long collapsing ALTs it exists to break up while the",
             "extra-records cost still applies. And 0.15% of records carry a ploidy-coherence FILTER",
             "(`nested_diploid` 2,458, `nested_unreachable` 5,000, `nested_haploid` 0), meaning the",
             "child's ploidy and its parent's final genotype disagree; those calls are flagged rather",
             "than corrected.",
             "", "## Small variants (aardvark, GT)", ""]

    # JointIndel, not Indel: aardvark's plain Indel row is query-only (truth_total 0), so summing
    # it reports FPs against no truth at all and an F1 of nan. The tier-2 pages use the joint row
    # for the same reason.
    for vtype, label in (("ALL", "ALL"), ("Snv", "SNV"), ("JointIndel", "Indel")):
        tp = fp = fn = 0
        for r in results:
            if not r.get("scoreable", True):
                continue
            row = pick(r.get("aardvark"), "GT", vtype)
            if row:
                tp += int(row.get("truth_tp", 0) or 0)
                fp += int(row.get("query_fp", 0) or 0)
                fn += int(row.get("truth_fn", 0) or 0)
        lines.append(f"- **{label}**: TP {tp:,}  FP {fp:,}  FN {fn:,}  "
                     f"recall {tp/(tp+fn) if tp+fn else float('nan'):.4f}  "
                     f"precision {tp/(tp+fp) if tp+fp else float('nan'):.4f}  "
                     f"**F1 {f1(tp, fp, fn):.4f}**")

    lines += ["", "## Structural variants (truvari, >=50 bp)", ""]
    tp = fp = fn = 0
    for r in results:
        if not r.get("scoreable", True):
            continue
        s = r.get("truvari") or {}
        tp += int(s.get("TP-base", 0) or 0)
        fp += int(s.get("FP", 0) or 0)
        fn += int(s.get("FN", 0) or 0)
    lines.append(f"- TP {tp:,}  FP {fp:,}  FN {fn:,}  **F1 {f1(tp, fp, fn):.4f}**")

    lines += ["", "## Per contig", "",
              "| contig | small F1 | SV F1 | notes |", "|---|---|---|---|"]
    for r in results:
        a = pick(r.get("aardvark"), "GT", "ALL")
        sm = "-"
        if a:
            sm = f"{f1(int(a.get('truth_tp',0) or 0), int(a.get('query_fp',0) or 0), int(a.get('truth_fn',0) or 0)):.4f}"
        t = r.get("truvari") or {}
        sv = f"{t['f1']:.4f}" if isinstance(t.get("f1"), (int, float)) else "-"
        notes = "; ".join(k for k in ("error", "aardvark_error", "truvari_error") if k in r)
        if not r.get("scoreable", True):
            notes = ("excluded: reference mismatch with truth" + ("; " + notes if notes else ""))
        lines.append(f"| {r['contig']} | {sm} | {sv} | {notes} |")

    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
