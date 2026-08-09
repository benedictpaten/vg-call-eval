#!/usr/bin/env python3
"""How much of the SV error budget is the ruler rather than the caller.

Two knobs, run over the arms and graphs that matter.

**`--refdist` sweep.** Truvari will not match a query record to a base record more
than `--refdist` bp away, default 500 here. In tandem repeats the same event can be
written anywhere in the array, and the atlas finds ~11% of false positives carrying a
near-perfect near-miss (>=0.7 sequence and size similarity) at a distance beyond that.
Widening the window puts a number on how much of the FP/FN budget is placement. This
changes the *metric*, not the caller, so it is reported as a sensitivity band and never
as a headline accuracy.

**`truvari refine`.** The principled version of the same question: re-align the
candidate regions truvari itself flags (`candidate.refine.bed`) with MAFFT so that the
benchmark and the call set are written in a common representation, then re-compare.
Where the sweep bounds the placement component, refine settles it -- a pair that
survives harmonisation is a real disagreement.

Both are needed because they fail differently. A wider `--refdist` can manufacture
matches between genuinely different events that happen to be nearby and similar;
refine cannot, but it only looks at regions truvari already suspected.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"
TRUVARI = REPO / "work/truvari-venv/bin/truvari"

DATASETS = {
    "chr6-4hap": ("tier2-chr6", "chr6"),
    "chr6-34hap": ("tier2-chr6-hap32", "chr6"),
    "chr20-4hap": ("tier2-chr20", "chr20"),
    "chr20-34hap": ("tier2-chr20-hap32", "chr20"),
}


def run(cmd, **kw):
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, **kw)
    return r


def bench(truth: Path, query: Path, ref: Path, bed: Path, odir: Path,
          refdist: int, threads: int = 4) -> dict | None:
    if odir.exists():
        shutil.rmtree(odir)
    # truvari refuses --chunksize below --refdist: records are gathered into chunks
    # first, so a match window wider than the chunk could never be found. The default
    # chunksize is 1000, which silently caps the sweep at refdist 1000.
    r = run([TRUVARI, "bench", "-b", truth, "-c", query, "-f", ref, "-o", odir,
             "--includebed", bed, "--sizemin", 50, "--sizefilt", 50,
             "--bSample", "HG002", "--cSample", "SAMPLE", "--pick", "ac",
             "--refdist", refdist, "--chunksize", max(1000, refdist)])
    s = odir / "summary.json"
    if not s.exists():
        print(f"  bench failed (refdist={refdist}): {r.stderr[-400:]}", flush=True)
        return None
    return json.loads(s.read_text())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["poisson-z", "readlik-z"])
    ap.add_argument("--datasets", nargs="*", default=list(DATASETS))
    ap.add_argument("--refdists", nargs="*", type=int, default=[500, 1000, 2000])
    ap.add_argument("--refine", action="store_true", help="also run truvari refine")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    out = {}
    for ds in args.datasets:
        sub, contig = DATASETS[ds]
        w = WORK / sub
        res = w / "results"
        ref = w / f"{contig}.fa"
        bed = w / f"truth.{contig}.stvar.bed"
        truth = res / "truvari-norm" / "truth.stvar.norm.vcf.gz"
        for arm in args.arms:
            q = res / "truvari-norm" / f"{arm}.norm.vcf.gz"
            if not q.exists():
                print(f"skip {ds} {arm}: no normalised VCF", file=sys.stderr)
                continue
            for rd in args.refdists:
                odir = res / f"truvari-{arm}-rd{rd}"
                # refdist 500 is the shipped setting; reuse its existing directory
                # rather than recomputing an identical answer.
                if rd == 500:
                    odir = res / f"truvari-{arm}"
                    s = json.loads((odir / "summary.json").read_text())
                else:
                    s = bench(truth, q, ref, bed, odir, rd, args.threads)
                if s is None:
                    continue
                out[f"{ds}|{arm}|rd{rd}"] = {
                    k: s.get(k) for k in
                    ("TP-base", "TP-comp", "FP", "FN", "precision", "recall", "f1",
                     "base cnt", "comp cnt")}
                print(f"{ds:13s} {arm:11s} refdist {rd:5d}: "
                      f"TP-base {s['TP-base']:5d} FP {s['FP']:5d} FN {s['FN']:5d} "
                      f"P {s['precision']:.4f} R {s['recall']:.4f} F1 {s['f1']:.4f}",
                      flush=True)

            if args.refine:
                base = res / f"truvari-{arm}"
                cand = base / "candidate.refine.bed"
                if not cand.exists() or cand.stat().st_size == 0:
                    print(f"  {ds} {arm}: no refine candidates", flush=True)
                    continue
                # refine writes in place; copy so the original bench stays intact.
                rdir = res / f"truvari-{arm}-refine"
                if rdir.exists():
                    shutil.rmtree(rdir)
                shutil.copytree(base, rdir)
                r = run([TRUVARI, "refine", "-f", ref, "-t", args.threads,
                         "--align", "mafft", "--use-original-vcfs", rdir])
                rs = rdir / "refine.variant_summary.json"
                if rs.exists():
                    s = json.loads(rs.read_text())
                    out[f"{ds}|{arm}|refine"] = {
                        k: s.get(k) for k in
                        ("TP-base", "TP-comp", "FP", "FN", "precision", "recall", "f1")}
                    print(f"{ds:13s} {arm:11s} refine     : "
                          f"TP-base {s.get('TP-base')} FP {s.get('FP')} FN {s.get('FN')} "
                          f"P {s.get('precision'):.4f} R {s.get('recall'):.4f} "
                          f"F1 {s.get('f1'):.4f}", flush=True)
                else:
                    print(f"  {ds} {arm}: refine produced no summary\n"
                          f"    {r.stderr[-500:]}", flush=True)

    dest = WORK / "sv-atlas" / "metric_sensitivity.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
