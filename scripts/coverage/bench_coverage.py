#!/usr/bin/env python3
"""Score the coverage titration and produce the three curves Stage 1 is judged against.

  1. F1 vs coverage. The headline, and the least interesting of the three -- it will go down as
     coverage falls and that is not news.

  2. The GQ-vs-depth slope, per zygosity. This is the thing that has to be fixed. Measured on the
     whole-genome run, chr7 het median GQ runs 24 at DP 8-14 to 151 at DP 34-44 -- about 3.5x DP,
     linear. If that slope holds across the titration then GQ is a pure depth multiple and no fixed
     threshold can mean the same thing at two coverages.

  3. **The calibration curve: observed precision against claimed GQ, one line per coverage.** This
     is the artifact that decides Stage 1. If the lines lie on top of each other, GQ is already
     calibrated and only the *distribution* shifts with coverage -- so the fix is a gate, not a
     rescale. If they fan out, GQ means different things at different depths and needs the
     normalisation. Those are different fixes and the curve tells us which.

Scoring reuses the harness's own `aardvark.compare` rather than a fresh command line, for the same
reason bench_wgs.py does: the tier-2 invocation has accumulated specifics, and a run that quietly
used different ones would not be comparable with any existing number.

The full-coverage arm doubles as a control. It uses the same whole-genome reads database the tier-2
run used, so its F1 must reproduce the published chr20 figure. If it does not, the harness is
wrong and nothing else here can be trusted.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))

from vgcalleval.engines import aardvark  # noqa: E402


def query(vcf: Path, fmt: str) -> list[list[str]]:
    out = subprocess.run(["bcftools", "query", "-f", fmt, str(vcf)],
                         capture_output=True, text=True).stdout
    return [l.split("\t") for l in out.splitlines()]


def f1(tp: int, fp: int, fn: int) -> float:
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def score(work: Path, tag: str, truth_vcf: Path, truth_bed: Path, ref: Path,
          sample: str, threads: int) -> dict | None:
    """aardvark for one coverage level, cached on summary.tsv."""
    calls = work / f"call.{tag}x.vcf.gz"
    if not calls.exists():
        print(f"  {tag}x: no calls, skipping")
        return None
    renamed = work / f"call.{tag}x.renamed.vcf.gz"
    if not renamed.exists():
        names = work / "sample.txt"
        names.write_text(f"{sample}\n")
        subprocess.run(["bcftools", "reheader", "-s", str(names), "-o", str(renamed), str(calls)],
                       check=True)
        subprocess.run(["bcftools", "index", "-f", "-t", str(renamed)], check=True)
    adir = work / f"aardvark.{tag}x"
    if not (adir / "summary.tsv").exists():
        try:
            aardvark.compare(
                aardvark="aardvark", reference=ref, truth_vcf=truth_vcf, query_vcf=renamed,
                regions_bed=truth_bed, out_dir=adir, truth_sample=sample, query_sample=sample,
                label=f"cov{tag}", options=aardvark.AardvarkOptions(threads=threads),
            )
        except (Exception, SystemExit) as exc:   # noqa: BLE001
            # SystemExit as well as Exception: the harness reports a failed command by exiting
            # rather than raising, so a bare `except Exception` misses it -- the same trap
            # bench_wgs.py hit, where one contig's failure ended a run after 23 had been scored.
            print(f"  {tag}x: aardvark failed: {exc}")
            return None
    if not (adir / "summary.tsv").exists():
        print(f"  {tag}x: no summary produced")
        return None
    return {"tag": tag, "aardvark": aardvark.read_summary(adir), "adir": adir,
            "calls": calls}


def pick(rows, comparison, vtype):
    for r in rows or []:
        if r.get("comparison", "").upper() == comparison and r.get("variant_type") == vtype:
            return r
    return None


def per_call(adir: Path, calls: Path) -> list[dict]:
    """BD decisions joined to GQ/DP/GT, keyed by position."""
    bd = {int(r[0]): r[1] for r in query(adir / "query.vcf.gz", "%POS\t[%BD]\n")}
    rows = []
    for r in query(calls, "%POS\t[%GT\t%GQ\t%DP]\n"):
        pos = int(r[0])
        try:
            gq, dp = float(r[2]), float(r[3])
        except ValueError:
            continue
        rows.append({"pos": pos, "gt": r[1], "gq": gq, "dp": dp, "bd": bd.get(pos)})
    return rows


def med(v):
    v = [x for x in v if x == x]
    return statistics.median(v) if v else float("nan")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", default="work/coverage/chr20", type=Path)
    p.add_argument("--truth-dir", default="work/tier2-chr20-hap32", type=Path)
    p.add_argument("--contig", default="chr20")
    p.add_argument("--tags", nargs="*", default=["5", "10", "15", "20", "25", "30"])
    p.add_argument("--sample", default="HG002")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--out", default=None, help="write a JSON summary here")
    args = p.parse_args()

    truth_vcf = args.truth_dir / f"truth.{args.contig}.smvar.vcf.gz"
    truth_bed = args.truth_dir / f"truth.{args.contig}.smvar.bed"
    ref = args.truth_dir / f"{args.contig}.fa"

    results = []
    for tag in args.tags:
        print(f"[score] {tag}x", flush=True)
        r = score(args.work, tag, truth_vcf, truth_bed, ref, args.sample, args.threads)
        if r:
            r["rows"] = per_call(r["adir"], r["calls"])
            results.append(r)
    if not results:
        sys.exit("nothing scored")

    # --- 1. F1 vs coverage ---------------------------------------------------
    print("\n== 1. accuracy vs coverage (aardvark GT) ==")
    print(f"{'cov':>5} {'calls':>9} {'medDP':>6} {'ALL F1':>8} {'SNV F1':>8} {'Indel F1':>9} "
          f"{'recall':>7} {'prec':>7}")
    table = []
    for r in results:
        a = pick(r["aardvark"], "GT", "ALL")
        s = pick(r["aardvark"], "GT", "Snv")
        i = pick(r["aardvark"], "GT", "JointIndel")
        def g(row, k):
            try:
                return float(row.get(k, "nan"))
            except (TypeError, ValueError):
                return float("nan")
        row = {"coverage": r["tag"], "n_calls": len(r["rows"]),
               "median_dp": med([x["dp"] for x in r["rows"]]),
               "f1_all": g(a, "metric_f1"), "f1_snv": g(s, "metric_f1"),
               "f1_indel": g(i, "metric_f1"),
               "recall": g(a, "metric_recall"), "precision": g(a, "metric_precision")}
        table.append(row)
        print(f"{r['tag']+'x':>5} {row['n_calls']:9,d} {row['median_dp']:6.0f} "
              f"{row['f1_all']:8.4f} {row['f1_snv']:8.4f} {row['f1_indel']:9.4f} "
              f"{row['recall']:7.4f} {row['precision']:7.4f}")

    # --- 2. GQ vs depth ------------------------------------------------------
    print("\n== 2. median GQ by zygosity: is GQ a pure multiple of depth? ==")
    HET = {"0/1", "1/0", "0|1", "1|0", "1/2", "1|2", "2/1", "2|1"}
    HOM = {"1/1", "1|1", "1", "2/2", "2|2"}
    # A haploid contig has essentially no het calls -- the only ones are in the PAR, which the
    # T2T-Q100 confident regions exclude -- so the het column is empty there and the hom column
    # carries the signal. Reported as "-" rather than nan so an empty cell reads as "not
    # applicable at this ploidy" instead of as a failed computation.
    print(f"{'cov':>5} {'medDP':>6} {'het GQ':>7} {'hom GQ':>7} {'GQ/DP het':>10} "
          f"{'GQ/DP hom':>10} {'%GQ<10':>7} {'%at cap':>8}")
    for r in results:
        rows = r["rows"]
        dp = med([x["dp"] for x in rows])
        h = med([x["gq"] for x in rows if x["gt"] in HET])
        m = med([x["gq"] for x in rows if x["gt"] in HOM])
        low = 100 * sum(1 for x in rows if x["gq"] < 10) / len(rows)
        cap = 100 * sum(1 for x in rows if x["gq"] >= 256) / len(rows)
        def cell(v):
            return "      -" if v != v else f"{v:7.0f}"
        def ratio(v):
            return "         -" if (v != v or not dp) else f"{v/dp:10.2f}"
        print(f"{r['tag']+'x':>5} {dp:6.0f} {cell(h)} {cell(m)} {ratio(h)} {ratio(m)} "
              f"{low:6.1f}% {cap:7.1f}%")
        for t in table:
            if t["coverage"] == r["tag"]:
                t.update(gq_het=h, gq_hom=m,
                         gq_per_dp_het=h / dp if (dp and h == h) else None,
                         gq_per_dp_hom=m / dp if (dp and m == m) else None,
                         pct_gq_lt10=low, pct_at_cap=cap)
    print("Read the GQ/DP columns with care. Taken across arms they compare *different populations*")
    print("-- which sites get called het changes with depth -- and that alone made the diploid")
    print("column look strongly superlinear (2.00 to 3.79) when a comparison paired on identical")
    print("sites, which the nested subsampling makes possible, gives a nearly flat 3.2 to 3.8.")
    print("'%at cap' is the share censored by the 256 clamp: 23% of haploid calls at full depth,")
    print("which is why a normalised field has to be computed inside the caller, pre-clamp.")

    # --- 3. calibration ------------------------------------------------------
    print("\n== 3. calibration: observed precision at a claimed GQ, per coverage ==")
    buckets = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 160), (160, 10**9)]
    labels = [f"{lo}-{hi if hi < 10**9 else ''}" for lo, hi in buckets]
    print(f"{'cov':>5} " + " ".join(f"{l:>11}" for l in labels))
    calib = {}
    for r in results:
        cells, series = [], []
        for lo, hi in buckets:
            k = [x for x in r["rows"] if lo <= x["gq"] < hi and x["bd"] in ("TP", "FP")]
            if not k:
                cells.append(f"{'-':>11}")
                series.append(None)
                continue
            prec = sum(1 for x in k if x["bd"] == "TP") / len(k)
            cells.append(f"{prec:.3f}({len(k)//1000:>3d}k)")
            series.append({"precision": prec, "n": len(k)})
        calib[r["tag"]] = dict(zip(labels, series))
        print(f"{r['tag']+'x':>5} " + " ".join(f"{c:>11}" for c in cells))
    print("Lines on top of each other  -> GQ is already calibrated; only its distribution moves,")
    print("                               so the fix is an emission gate, not a rescale.")
    print("Lines fanning out           -> GQ means different things at different depths, and the")
    print("                               normalised companion field is needed.")

    if args.out:
        Path(args.out).write_text(json.dumps(
            {"per_coverage": table, "calibration": calib}, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
