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
    renamed = score / f"{contig}.renamed.vcf.gz"
    if not renamed.exists():
        names = score / f"{contig}.sample.txt"
        names.write_text(f"{sample}\n")
        run(["bcftools", "reheader", "-s", names, "-o", renamed, calls])
        run(["bcftools", "index", "-f", "-t", renamed])

    adir = score / f"{contig}.aardvark"
    if not (adir / "summary.tsv").exists():
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
    if not (tdir / "summary.json").exists():
        ref = d / f"{contig}.fa"
        norm = score / f"{contig}.norm.vcf.gz"
        truth_norm = score / f"{contig}.truth.norm.vcf.gz"
        for src, dst in ((renamed, norm), (d / f"truth.{contig}.stvar.vcf.gz", truth_norm)):
            if not dst.exists():
                p1 = subprocess.run(["bcftools", "norm", "-m", "-any", "-f", str(ref),
                                     "-Oz", "-o", str(dst), str(src)],
                                    capture_output=True, text=True)
                if p1.returncode != 0:
                    out["truvari_error"] = p1.stderr[-500:]
                    return out
                run(["bcftools", "index", "-f", "-t", dst])
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
    p.add_argument("--truvari", default="truvari")
    p.add_argument("--contigs", nargs="*", default=ALL_CONTIGS)
    args = p.parse_args()

    work = Path(args.work)
    results = []
    for c in args.contigs:
        print(f"[score] {c}", flush=True)
        results.append(score_contig(work, c, args.sample, args.threads, args.truvari))

    (work / "score" / "per-contig.json").write_text(json.dumps(results, indent=2))

    # Genome-wide totals from summed counts. Averaging per-contig F1s would weight chr21 like
    # chr1; summing the counts is the only aggregation that means anything.
    lines = ["# Whole-genome results: HG002 against T2T-Q100", "",
             "Called per contig on the 34-haplotype HPRC graph, `--read-likelihood` with panel",
             "enumeration, phasing and mosaic on. chrY haploid; chrX haploid outside the",
             "pseudoautosomal regions and diploid inside them, spliced from two runs.", "",
             "## Small variants (aardvark, GT)", ""]

    for vtype, label in (("ALL", "ALL"), ("Snv", "SNV"), ("Indel", "Indel")):
        tp = fp = fn = 0
        for r in results:
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
        lines.append(f"| {r['contig']} | {sm} | {sv} | {notes} |")

    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
