#!/usr/bin/env python3
"""What fraction of a het call's reads actually goes to the non-reference allele, as a
function of how long that allele is?

This is the measurement that has to come before any length-aware binomial prior. The
obvious geometric model -- a 151 bp read can only span an insertion of length L if
151 > L + 2*flank, so expected ALT support falls to nothing by L ~ 130 -- turns out to
predict a skew so extreme that scoring true het SV calls against it makes them look
impossible. That model is wrong because the caller does not require a read to span an
allele: a read covering one breakpoint still aligns better to one allele than the other.

Rather than guess a better geometry, measure it. Restricted to calls the benchmark says
are true, so the curve describes what a *correct* het looks like, which is the only thing
a null hypothesis can usefully be built from.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
from pathlib import Path


def labels_aardvark(work: Path) -> dict:
    bd = {}
    with gzip.open(work / "results/aardvark-readlik-z/query.vcf.gz", "rt") as fh:
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


def labels_truvari(work: Path) -> dict:
    bd = {}
    d = work / "results/truvari-readlik-z"
    for fname, lab in (("tp-comp.vcf.gz", "TP"), ("fp.vcf.gz", "FP")):
        with gzip.open(d / fname, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    bd[int(line.split("\t", 2)[1])] = lab
    return bd


BUCKETS = [(0, 0), (1, 1), (2, 3), (4, 7), (8, 15), (16, 31), (32, 63),
           (64, 127), (128, 255), (256, 1000), (1001, 10 ** 9)]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--labels", choices=["aardvark", "truvari"], default="aardvark")
    args = p.parse_args()

    W = Path(args.work)
    bd = labels_aardvark(W) if args.labels == "aardvark" else labels_truvari(W)

    q = subprocess.run(
        ["bcftools", "query", "-f", "%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD]\n",
         str(W / "results/readlik-z.vcf.gz")], capture_output=True, text=True)

    # bucket -> [sum of alt fraction, n, sum of AD total / DP]
    acc = {b: [0.0, 0, 0.0] for b in BUCKETS}
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6 or bd.get(int(f[0])) != "TP":
            continue
        gt = [int(x) for x in f[3].replace("|", "/").split("/") if x.isdigit()]
        if len(gt) != 2 or gt[0] == gt[1] or 0 not in gt:
            continue  # only ref/alt hets: with two alt alleles there is no "the" alt
        try:
            dp = float(f[4])
            ad = [int(x) for x in f[5].split(",")]
        except ValueError:
            continue
        alt_i = gt[0] or gt[1]
        alts = f[2].split(",")
        if alt_i - 1 >= len(alts) or alt_i >= len(ad) or dp <= 0:
            continue
        a = alts[alt_i - 1]
        if a[0] in "<*":
            continue
        d = abs(len(a) - len(f[1]))
        tot = ad[0] + ad[alt_i]
        if tot < 4:
            continue
        for b in BUCKETS:
            if b[0] <= d <= b[1]:
                acc[b][0] += ad[alt_i] / tot
                acc[b][1] += 1
                acc[b][2] += sum(ad) / dp
                break

    print(f"=== {args.label} ({args.labels} labels): true ref/het calls ===")
    print(f"  {'|len ALT-REF|':<16}{'n':>8}{'mean ALT share':>16}{'mean AD/DP':>13}")
    for b in BUCKETS:
        s, n, dpf = acc[b]
        if n == 0:
            continue
        rng = f"{b[0]}" if b[0] == b[1] else f"{b[0]}-{b[1] if b[1] < 10**8 else '+'}"
        print(f"  {rng:<16}{n:>8,}{s / n:>16.3f}{dpf / n:>13.3f}")


if __name__ == "__main__":
    main()
