#!/usr/bin/env python3
"""Score candidate GQ normalisers against the titration, before any of them is built in C++.

The question Stage 1 has to answer is: what quantity can a consumer threshold and get the same
meaning at any coverage and any ploidy? This scores candidates by how much the precision-at-a-
claimed-score curve moves between arms. Lower spread is better calibrated.

**Always evaluate pooled across ploidies, never on one contig.** `GQ/DP` looks like a clear win on
the diploid series alone -- it halves the spread, 0.101 to 0.050 -- and it is a clear loss once the
haploid series is included:

    series                    raw GQ     GQ/DP
    chr20 (diploid)            0.101     0.050
    chrX  (haploid)            0.150     0.161
    POOLED                     0.348     0.496

At a matched GQ/DP the two ploidies sit 0.6 apart in precision (chr20 0.75-0.95, chrX 0.14-0.25),
so dividing by depth removes the depth axis and leaves the larger ploidy axis untouched -- and by
compressing the score range it actually makes the pooled picture worse. Validating on chr20 alone
would have shipped a field that degrades the thing it exists to fix.

Why ploidy dominates: at ploidy 1 the runner-up genotype is a different allele outright, so each
read discriminates fully. At ploidy 2 a het's runner-up differs on one strand, so a read
discriminates about half as much. The per-read gap scale is therefore a function of ploidy, which
is exactly what a normaliser has to divide out and what `1/DP` cannot see.

That points at the remaining candidate: divide the observed gap by the gap achievable at this
site -- the gap a noise-free pileup would produce given lambda, the ploidy, and the allele set.
That quantity needs the per-read likelihood matrix, so it cannot be evaluated from a VCF; use
`vg call --dump-likelihoods` to get the matrix and extend this script rather than guessing.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

RAW_BUCKETS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 80), (80, 1e9)]
NORM_BUCKETS = [(0, .25), (.25, .75), (.75, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 1e9)]


def load(work: Path, tag: str) -> list[tuple[float, float, str]]:
    """(GQ, DP, TP/FP) for every scored call in one arm."""
    bd = subprocess.run(["bcftools", "query", "-f", "%POS\t[%BD]\n",
                         str(work / f"aardvark.{tag}x" / "query.vcf.gz")],
                        capture_output=True, text=True).stdout
    B = {int(l.split("\t")[0]): l.split("\t")[1] for l in bd.splitlines()}
    out = subprocess.run(["bcftools", "query", "-f", "%POS\t[%GQ\t%DP]\n",
                          str(work / f"call.{tag}x.vcf.gz")],
                         capture_output=True, text=True).stdout
    rows = []
    for l in out.splitlines():
        f = l.split("\t")
        try:
            gq, dp = float(f[1]), float(f[2])
        except ValueError:
            continue
        b = B.get(int(f[0]))
        if b in ("TP", "FP"):
            rows.append((gq, dp, b))
    return rows


def spread(score, buckets, arms: dict) -> float:
    """Mean over buckets of (max - min) precision across arms. 0 = perfectly comparable."""
    per = [[] for _ in buckets]
    for rows in arms.values():
        for i, (lo, hi) in enumerate(buckets):
            k = [r for r in rows if lo <= score(r) < hi]
            if len(k) < 200:      # a bucket with too few calls is noise, not a data point
                continue
            per[i].append(sum(1 for r in k if r[2] == "TP") / len(k))
    s = [max(p) - min(p) for p in per if len(p) >= 4]
    return sum(s) / len(s) if s else float("nan")


CANDIDATES = {
    "raw GQ (status quo)": (lambda r: r[0], RAW_BUCKETS),
    "GQ / DP": (lambda r: r[0] / r[1] if r[1] else 0.0, NORM_BUCKETS),
    "GQ / sqrt(DP)": (lambda r: r[0] / (r[1] ** 0.5) if r[1] else 0.0,
                      [(0, 1), (1, 2), (2, 4), (4, 8), (8, 14), (14, 1e9)]),
}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--series", nargs="+",
                   default=["chr20 (diploid)=work/coverage/chr20=5,10,15,20,25,30",
                            "chrX (haploid)=work/coverage/chrX=2.5,5,7.5,10,12.5,14.6"],
                   help="NAME=WORKDIR=tag,tag,...")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    data = {}
    for spec in args.series:
        name, work, tags = spec.split("=")
        data[name] = {t: load(Path(work), t) for t in tags.split(",")}

    names = list(CANDIDATES)
    print("Mean precision spread across arms (lower = better calibrated)\n")
    print(f"{'series':26s} " + " ".join(f"{n:>20s}" for n in names))
    result = {}
    for s, arms in data.items():
        row = {n: spread(f, b, arms) for n, (f, b) in CANDIDATES.items()}
        result[s] = row
        print(f"{s:26s} " + " ".join(f"{row[n]:20.3f}" for n in names))

    # The test that matters. A field comparable across coverage but not across ploidy passes
    # every row above and fails this one.
    pooled = {f"{s}|{t}": rows for s, arms in data.items() for t, rows in arms.items()}
    row = {n: spread(f, b, pooled) for n, (f, b) in CANDIDATES.items()}
    result["POOLED across ploidies"] = row
    print(f"\n{'POOLED across ploidies':26s} " + " ".join(f"{row[n]:20.3f}" for n in names))
    print("\nThe pooled row is the acceptance test for Stage 1. A candidate that does not beat")
    print("raw GQ there is not a fix, however good it looks on one contig.")

    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2))
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
