#!/usr/bin/env python3
"""Score a parameter sweep and say, per coverage and ploidy, where the optimum sits.

The question this exists to answer is whether a default fitted at ~30x diploid is still right at
5x, and at ploidy 1. A parameter whose optimum moves with coverage cannot have one shipped value;
one whose optimum does not move can keep the value it has and be documented as coverage-invariant,
which is the more useful outcome because it is one less thing for a user to tune.

Reports the best value per arm **and the size of the gain**, which matters more than the argmax.
An optimum that is 0.0003 better than the default is not a reason to change a default -- it is
noise plus a different arm, and changing on that basis would invalidate every published number for
nothing.
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

# arm -> (truth dir, contig, aardvark dir for the truth-total lookup)
ARMS = {
    "chr20.5":   ("work/tier2-chr20-hap32", "chr20"),
    "chr20.30":  ("work/tier2-chr20-hap32", "chr20"),
    "chrX.2.5":  ("work/coverage/chrX/truthdir", "chrX"),
    "chrX.14.6": ("work/coverage/chrX/truthdir", "chrX"),
}


def f1(tp: int, fp: int, fn: int) -> float:
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def score(vcf: Path, truth_dir: Path, contig: str, out_dir: Path, sample: str,
          threads: int) -> dict | None:
    if not vcf.exists():
        return None
    renamed = out_dir.parent / (vcf.stem.replace(".vcf", "") + ".renamed.vcf.gz")
    if not renamed.exists():
        names = out_dir.parent / "sample.txt"
        names.write_text(f"{sample}\n")
        subprocess.run(["bcftools", "reheader", "-s", str(names), "-o", str(renamed), str(vcf)],
                       check=True)
        subprocess.run(["bcftools", "index", "-f", "-t", str(renamed)], check=True)
    if not (out_dir / "summary.tsv").exists():
        try:
            aardvark.compare(
                aardvark="aardvark", reference=truth_dir / f"{contig}.fa",
                truth_vcf=truth_dir / f"truth.{contig}.smvar.vcf.gz", query_vcf=renamed,
                regions_bed=truth_dir / f"truth.{contig}.smvar.bed", out_dir=out_dir,
                truth_sample=sample, query_sample=sample, label=out_dir.name,
                options=aardvark.AardvarkOptions(threads=threads))
        except (Exception, SystemExit) as exc:   # noqa: BLE001
            print(f"    aardvark failed for {vcf.name}: {exc}")
            return None
    if not (out_dir / "summary.tsv").exists():
        return None
    rows = aardvark.read_summary(out_dir)
    for r in rows or []:
        if r.get("comparison", "").upper() == "GT" and r.get("variant_type") == "ALL":
            try:
                return {"f1": float(r["metric_f1"]), "recall": float(r["metric_recall"]),
                        "precision": float(r["metric_precision"])}
            except (KeyError, TypeError, ValueError):
                return None
    return None


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tag", required=True, help="sweep tag, e.g. lw")
    p.add_argument("--param", required=True, help="the flag swept, for the report header")
    p.add_argument("--values", nargs="+", required=True)
    p.add_argument("--default", required=True, help="the shipped value, marked in the table")
    p.add_argument("--work", default=None)
    p.add_argument("--sample", default="HG002")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    work = Path(args.work) if args.work else Path(f"work/coverage/sweep/{args.tag}")
    results = {}
    for arm, (truth_dir, contig) in ARMS.items():
        results[arm] = {}
        for v in args.values:
            vcf = work / f"{arm}.{args.tag}{v}.vcf.gz"
            print(f"[score] {arm} {args.param} {v}", flush=True)
            results[arm][v] = score(vcf, Path(truth_dir), contig,
                                    work / f"aardvark.{arm}.{v}", args.sample, args.threads)

    print(f"\n== F1 by {args.param} ==\n")
    print(f"{'arm':12s}" + "".join(f"{v:>10}" for v in args.values)
          + f"{'best':>8}{'gain vs default':>17}")
    summary = {}
    for arm, per in results.items():
        cells = []
        for v in args.values:
            r = per.get(v)
            cells.append(f"{r['f1']:10.4f}" if r else f"{'-':>10}")
        got = {v: per[v]["f1"] for v in args.values if per.get(v)}
        if not got:
            print(f"{arm:12s}" + "".join(cells) + f"{'-':>8}{'-':>17}")
            continue
        best = max(got, key=lambda v: got[v])
        base = got.get(args.default)
        gain = (got[best] - base) if base is not None else float("nan")
        summary[arm] = {"best": best, "best_f1": got[best], "default_f1": base, "gain": gain}
        print(f"{arm:12s}" + "".join(cells) + f"{best:>8}{gain:>+17.4f}")

    print(f"\nShipped default: {args.param} {args.default}")
    print("Read the gain column, not the argmax. A different best value worth 0.0003 is noise and")
    print("a different arm; changing a default on that basis invalidates every published number")
    print("for nothing.")
    if summary:
        moves = {a: s["best"] for a, s in summary.items()}
        if len(set(moves.values())) == 1:
            print(f"\nThe optimum is {next(iter(moves.values()))} on every arm: this parameter is")
            print("coverage- and ploidy-invariant over the range measured, and needs no table.")
        else:
            print(f"\nThe optimum differs by arm: {moves}")
            print("Check whether the gains are large enough to be worth a coverage-dependent value.")

    if args.out:
        Path(args.out).write_text(json.dumps({"param": args.param, "default": args.default,
                                              "results": results, "summary": summary}, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
