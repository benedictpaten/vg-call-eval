#!/usr/bin/env python3
"""Score an arbitrary call set against both tier-2 benchmarks.

`run_arms.py` scores the five fixed arms. This scores one VCF produced by any
invocation -- a new flag, a swept parameter, a patched build -- against *both* the
small-variant benchmark (aardvark) and the structural one (truvari), and reports the
heterozygous-deletion breakdown that the SV investigation turns on.

Both benchmarks, always, and that is the point. A change aimed at structural variants
can wreck small-variant genotyping without touching the SV numbers, and a change aimed
at heterozygotes can move the het/hom balance genome-wide while every summary F1 looks
respectable. Scoring only the class a change was designed for is how a regression ships.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))
WORK = REPO / "work"
TRUVARI = WORK / "truvari-venv/bin/truvari"

DATASETS = {
    "chr6-4hap": ("tier2-chr6", "chr6"),
    "chr6-34hap": ("tier2-chr6-hap32", "chr6"),
    "chr20-4hap": ("tier2-chr20", "chr20"),
    "chr20-34hap": ("tier2-chr20-hap32", "chr20"),
}


def run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"failed: {' '.join(str(c) for c in cmd)}\n{r.stderr[-2000:]}")
    return r


def normalise(src: Path, dst: Path, ref: Path) -> None:
    """Split multi-allelics and sort, as truvari_sv.py does.

    Truvari matches one ALT at a time, so a multi-allelic record would otherwise be
    judged on its first ALT alone. The sort is required, not tidiness: left-aligning a
    split allele can move it upstream of the record that followed it.
    """
    if dst.exists():
        return
    tmp = dst.with_suffix(".unsorted.vcf.gz")
    run(["bcftools", "norm", "-m-any", "-f", ref, "-Oz", "-o", tmp, src])
    run(["bcftools", "sort", "-Oz", "-o", dst, tmp])
    run(["bcftools", "index", "-t", "-f", dst])
    tmp.unlink(missing_ok=True)


def sv_breakdown(bench_dir: Path) -> dict:
    """Recall by type, size class and zygosity, from truvari's base-side files."""
    c = collections.Counter()
    for side, out in (("tp-base", "TP"), ("fn", "FN")):
        p = bench_dir / f"{side}.vcf.gz"
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                t = info.get("SVTYPE")
                if t not in ("INS", "DEL"):
                    continue
                sl = abs(int(info.get("SVLEN", 0) or 0))
                b = ("1k+" if sl >= 1000 else "300-999" if sl >= 300
                     else "100-299" if sl >= 100 else "50-99")
                z = "hom" if f[9].split(":")[0].replace("|", "/") == "1/1" else "het"
                c[(t, b, z, out)] += 1
    return c


def gt_mix(vcf: Path) -> dict:
    """Het/hom balance over the whole call set.

    The max-allele variant is expected to shift this hard, and it is invisible in any
    per-class F1, so it is reported for every scored VCF.
    """
    n = collections.Counter()
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            gt = line.split("\t")[9].split(":")[0].replace("|", "/")
            a = gt.split("/")
            if len(a) != 2:
                continue
            n["het" if a[0] != a[1] else "hom"] += 1
    total = sum(n.values()) or 1
    return {"het": n["het"], "hom": n["hom"], "het_frac": round(n["het"] / total, 4)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--dataset", default="chr6-4hap")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--skip-aardvark", action="store_true")
    args = ap.parse_args()

    sub, contig = DATASETS[args.dataset]
    w = WORK / sub
    res = w / "results"
    ref = w / f"{contig}.fa"
    vcf = Path(args.vcf)

    out = {"label": args.label, "dataset": args.dataset, "vcf": str(vcf)}
    out["genotype_mix"] = gt_mix(vcf)
    print(f"genotype mix: {out['genotype_mix']}", flush=True)

    # ---- structural, truvari ----
    norm = res / "truvari-norm" / f"{args.label}.norm.vcf.gz"
    normalise(vcf, norm, ref)
    bench = res / f"truvari-{args.label}"
    if bench.exists():
        shutil.rmtree(bench)
    run([TRUVARI, "bench", "-b", res / "truvari-norm/truth.stvar.norm.vcf.gz",
         "-c", norm, "-f", ref, "-o", bench,
         "--includebed", w / f"truth.{contig}.stvar.bed",
         "--sizemin", 50, "--sizefilt", 50,
         "--bSample", "HG002", "--cSample", "SAMPLE", "--pick", "ac"])
    s = json.loads((bench / "summary.json").read_text())
    out["sv"] = {k: s.get(k) for k in
                 ("TP-base", "TP-comp", "FP", "FN", "recall", "precision", "f1")}
    print(f"SV: TP-base {s['TP-base']} FP {s['FP']} FN {s['FN']} "
          f"R {s['recall']:.4f} P {s['precision']:.4f} F1 {s['f1']:.4f}", flush=True)

    c = sv_breakdown(bench)
    out["sv_by_class"] = {}
    print("\n  recall by type/size/zygosity:")
    for t in ("DEL", "INS"):
        for b in ("50-99", "100-299", "300-999", "1k+"):
            for z in ("het", "hom"):
                tp, fn = c[(t, b, z, "TP")], c[(t, b, z, "FN")]
                if tp + fn == 0:
                    continue
                out["sv_by_class"][f"{t} {b} {z}"] = [tp, tp + fn]
                print(f"    {t} {b:8s} {z}  {tp:4d}/{tp+fn:4d} = {tp/(tp+fn):.3f}")

    # ---- small variants, aardvark ----
    if not args.skip_aardvark:
        from vgcalleval.engines import aardvark  # noqa: E402
        adir = res / f"aardvark-{args.label}"
        try:
            aardvark.compare(
                aardvark="aardvark", reference=ref,
                truth_vcf=w / f"truth.{contig}.smvar.vcf.gz", query_vcf=vcf,
                regions_bed=w / f"truth.{contig}.smvar.bed", out_dir=adir,
                truth_sample="HG002", query_sample="SAMPLE", label=args.label,
                options=aardvark.AardvarkOptions(threads=args.threads))
            out["smallvar"] = aardvark.read_summary(adir)
            print("\nsmall variants (aardvark) written to", adir, flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"aardvark failed: {exc}", file=sys.stderr)

    dest = WORK / "sv-atlas" / f"score-{args.label}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
