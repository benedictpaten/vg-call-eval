#!/usr/bin/env python3
"""Reduce `vg call --dump-likelihoods` to per-site, per-allele read support.

The genotype likelihood uses only *which* alleles a read fits, through a mixture that
assumes each haplotype in the genotype contributed exactly 1/|G| of the reads. It has no
sampling model for that split and no length adjustment, so two signals the reads carry are
currently unused:

  - **allele balance**: how the reads actually divide between the called alleles, against
    what a binomial with a length-aware expectation would predict;
  - **total depth**: whether this site has as many reads as its neighbourhood implies.

This table is the raw material for asking whether either separates false calls from true
ones. For each site it emits the number of reads whose best-fitting allele is `a`, for
every `a`, plus the total and the mean absolute fit (`best_ln`, the row divisor, which is
the only surviving record of how well a read matched anything at all).

Ties are counted fractionally rather than assigned to the first allele: at exactly the
sites of interest many reads fit several alleles equally, and awarding those to index 0
would manufacture the very skew we are trying to measure.
"""

from __future__ import annotations

import argparse
import sys


def flush(out, site, n_alleles, counts, n_reads, best_ln_sum):
    if site is None:
        return
    support = ",".join(f"{c:.2f}" for c in counts[:n_alleles])
    out.write(f"{site}\t{n_alleles}\t{n_reads}\t"
              f"{best_ln_sum / n_reads if n_reads else 0:.4f}\t{support}\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    args = p.parse_args()

    site = None
    n_alleles = 0
    counts: list[float] = []
    n_reads = 0
    best_ln_sum = 0.0
    written = 0

    with open(args.out, "w") as out:
        out.write("site\tn_alleles\tn_reads\tmean_best_ln\tsupport\n")
        for line in sys.stdin:
            if line.startswith("#site"):
                if site is not None:
                    flush(out, site, n_alleles, counts, n_reads, best_ln_sum)
                    written += 1
                n_alleles = len(line.rstrip("\n").split("\t")) - 4
                site, counts, n_reads, best_ln_sum = None, [0.0] * n_alleles, 0, 0.0
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 + n_alleles or not f[0]:
                continue
            if site is None:
                site = f[0]
            n_reads += 1
            try:
                best_ln_sum += float(f[3])
                rel = [float(x) for x in f[4:4 + n_alleles]]
            except ValueError:
                continue
            top = max(rel)
            if top <= 0:
                continue
            winners = [i for i, v in enumerate(rel) if v >= top - 1e-12]
            share = 1.0 / len(winners)
            for i in winners:
                counts[i] += share
        flush(out, site, n_alleles, counts, n_reads, best_ln_sum)
        written += 1

    print(f"wrote {written:,} sites to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
