#!/usr/bin/env python3
"""A two-condition guard for implausible read pile-ups, scored as a rule rather than a fit.

Falls out of the coverage analysis (`coverage_model.py`). On the 34-haplotype graph, false
structural calls carry a median depth ~30% *above* true ones (53.5 vs 41 on chr6, 51 vs 40
on chr20) while only ~63% of their reads fit any called allele, against 100% for true calls.
On the 4-haplotype graph neither difference exists. So the signature is specific to the rich
graph, and it is not "too few reads" -- the intuition the investigation started from -- but
too many, fitting nothing.

That is §5.3.3's depth-*plausibility* guard rather than a ranking term: it does not ask which
genotype is likelier, it asks whether this pile-up could have been produced by any genotype
at this site. Expressed as a rule with two thresholds rather than a fitted model, because a
guard has to be explainable and has to transfer, and because a logistic fit on ~400 SV calls
would not be worth the standard error.

Scored on both benchmarks: a guard that improves SVs by damaging small variants is not a
guard, it is a trade, and small variants outnumber SVs three hundred to one.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
from pathlib import Path


def labels(work: Path, kind: str) -> dict:
    bd = {}
    if kind == "truvari":
        for fn, lab in (("tp-comp.vcf.gz", "TP"), ("fp.vcf.gz", "FP")):
            with gzip.open(work / "results/truvari-readlik" / fn, "rt") as fh:
                for line in fh:
                    if not line.startswith("#"):
                        bd[int(line.split("\t", 2)[1])] = lab
    else:
        with gzip.open(work / "results/aardvark-readlik/query.vcf.gz", "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                v = dict(zip(f[8].split(":"), f[9].split(":"))).get("BD")
                if v in ("TP", "FP"):
                    bd[int(f[1])] = v
    return bd


def rolling_median(vals: list[float], window: int = 201) -> list[float]:
    half = window // 2
    out = []
    for i in range(len(vals)):
        w = sorted(vals[max(0, i - half):min(len(vals), i + half + 1)])
        out.append(w[len(w) // 2] if w else 0.0)
    return out


def collect(work: Path, kind: str):
    bd = labels(work, kind)
    q = subprocess.run(
        ["bcftools", "query", "-f", "%POS[\t%DP\t%AD]\n",
         str(work / "results/readlik.vcf.gz")], capture_output=True, text=True)
    # The local-depth baseline must come from *all* calls, not only labelled ones: with
    # truvari labels only a few hundred calls are labelled and a median over those would
    # be a median over megabases.
    pos, dps, ads = [], [], []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        try:
            dp = float(f[1])
            ad = [int(x) for x in f[2].split(",")]
        except (ValueError, IndexError):
            continue
        if dp <= 0:
            continue
        pos.append(int(f[0]))
        dps.append(dp)
        ads.append(sum(ad) / dp)
    med = rolling_median(dps)
    return [(p, bd[p], d / m if m > 0 else 1.0, s)
            for p, d, s, m in zip(pos, dps, ads, med) if p in bd]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work", nargs="+", required=True)
    p.add_argument("--depth-ratio", type=float, default=1.3,
                   help="flag when DP exceeds this multiple of the local median")
    p.add_argument("--min-share", type=float, default=0.8,
                   help="flag when the fraction of reads fitting a called allele is below this")
    args = p.parse_args()

    print(f"guard: DP > {args.depth_ratio} x local median  AND  sum(AD)/DP < {args.min_share}\n")
    hdr = (f"  {'dataset':<24}{'benchmark':<10}{'flagged':>9}{'of which FP':>13}"
           f"{'FP removed':>12}{'TP lost':>10}{'prec before':>13}{'prec after':>12}")
    print(hdr)
    for w in args.work:
        W = Path(w)
        for kind in ("aardvark", "truvari"):
            try:
                rows = collect(W, kind)
            except FileNotFoundError:
                continue
            if not rows:
                continue
            flagged = [r for r in rows if r[2] > args.depth_ratio and r[3] < args.min_share]
            tp = sum(1 for r in rows if r[1] == "TP")
            fp = len(rows) - tp
            f_fp = sum(1 for r in flagged if r[1] == "FP")
            f_tp = len(flagged) - f_fp
            before = tp / (tp + fp) if tp + fp else 0.0
            after = ((tp - f_tp) / (tp - f_tp + fp - f_fp)) if (tp - f_tp + fp - f_fp) else 0.0
            frac = f"{f_fp / len(flagged):.0%}" if flagged else "-"
            print(f"  {W.name:<24}{kind:<10}{len(flagged):>9,}{frac:>13}"
                  f"{f_fp:>12,}{f_tp:>10,}{before:>13.4f}{after:>12.4f}")


if __name__ == "__main__":
    main()
