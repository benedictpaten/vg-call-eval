#!/usr/bin/env python3
"""H2 of plan §9.22: an SV metric that does not care how a large event is represented.

Why this exists. Aardvark compares by reconstructing haplotypes, and the SV precision
reported since §9.18 is recomputed from its per-variant BD labels over >=50 bp query
records -- necessary, because aardvark zeroes its own query columns for the Sv*
categories. That metric charges a >=50 bp record as a false positive whenever it is not
*individually* credited, which is exactly what happens when a caller emits one large
event as several smaller records. On the 34-haplotype graph aardvark credits 74 more
truth SVs while the count of individually-credited >=50 bp records does not move at all,
so that is not a hypothetical.

Truvari matches by locus proximity, size similarity and sequence similarity instead, so a
correct event split across records, or merged, or shifted within `--refdist`, still
matches. If truvari and aardvark disagree about the direction of the 4-hap to 34-hap
change, the disagreement is the result: it means the deficit is in the metric.

Records are split with `bcftools norm -m-any` first because truvari matches one ALT at a
time and a multi-allelic record would otherwise be judged on its first ALT alone.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
VENV_TRUVARI = REPO / "work/truvari-venv/bin/truvari"


def run(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"failed: {' '.join(str(c) for c in cmd)}\n{r.stderr[-2000:]}")
    return r


def split_multiallelic(src: Path, dst: Path, reference: Path) -> None:
    """Split multi-allelic records, then sort.

    The sort is not optional. Left-aligning a split allele can move it upstream of the
    record that followed it, and bcftools then refuses to index -- "Unsorted positions on
    sequence #1". Splitting without sorting produced exactly that on the 34-haplotype
    graph while happening to succeed on the 4-haplotype one, which is the kind of
    difference that silently becomes a comparison between two different pipelines.
    """
    # Regenerate whenever the source is newer. Skipping purely on existence silently
    # compares a freshly called VCF against a *stale* normalised copy: on the run that
    # made the length-weighted mixture the default, this reported chr20 SV numbers
    # byte-identical to the previous model, because truvari was handed a two-day-old
    # file. Nothing errored, and the numbers were entirely plausible.
    if dst.exists() and dst.stat().st_mtime >= src.stat().st_mtime:
        return
    tmp = dst.with_suffix(".unsorted.vcf.gz")
    run(["bcftools", "norm", "-m-any", "-f", str(reference), "-Oz", "-o", str(tmp), str(src)])
    run(["bcftools", "sort", "-Oz", "-o", str(dst), str(tmp)])
    run(["bcftools", "index", "-t", "-f", str(dst)])
    tmp.unlink(missing_ok=True)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work", required=True)
    p.add_argument("--arms", nargs="*", default=["poisson-z", "readlik-z"])
    p.add_argument("--label", required=True)
    p.add_argument("--truvari", default=str(VENV_TRUVARI))
    p.add_argument("--sizemin", type=int, default=50)
    p.add_argument("--contig", default="chr20")
    args = p.parse_args()

    W = Path(args.work)
    res = W / "results"
    ref = W / f"{args.contig}.fa"
    truth = W / f"truth.{args.contig}.stvar.vcf.gz"
    bed = W / f"truth.{args.contig}.stvar.bed"

    norm_dir = res / "truvari-norm"
    norm_dir.mkdir(exist_ok=True)
    truth_norm = norm_dir / "truth.stvar.norm.vcf.gz"
    split_multiallelic(truth, truth_norm, ref)

    out = {}
    for arm in args.arms:
        q = res / f"{arm}.vcf.gz"
        if not q.exists():
            print(f"skip {arm}: {q} missing", flush=True)
            continue
        qn = norm_dir / f"{arm}.norm.vcf.gz"
        split_multiallelic(q, qn, ref)

        odir = res / f"truvari-{arm}"
        if odir.exists():
            run(["rm", "-rf", str(odir)])
        run([args.truvari, "bench", "-b", str(truth_norm), "-c", str(qn),
             "-f", str(ref), "-o", str(odir), "--includebed", str(bed),
             "--sizemin", str(args.sizemin), "--sizefilt", str(args.sizemin),
             "--bSample", "HG002", "--cSample", "SAMPLE", "--pick", "ac"])
        summary = json.loads((odir / "summary.json").read_text())
        out[arm] = summary
        print(f"{args.label} {arm}: recall {summary.get('recall'):.4f} "
              f"precision {summary.get('precision'):.4f} f1 {summary.get('f1'):.4f} "
              f"(TP-base {summary.get('TP-base')}, FP {summary.get('FP')}, FN {summary.get('FN')})",
              flush=True)

    (res / "truvari-summary.json").write_text(json.dumps(out, indent=2))
    print(f"wrote {res / 'truvari-summary.json'}")


if __name__ == "__main__":
    main()
