#!/usr/bin/env python3
"""Reduce a `vg call --dump-likelihoods` stream to one row per site.

H5 of plan §9.22 needs the number of alleles the genotyper actually *enumerated* at each
site, which is not recoverable from the VCF: `AT` lists only the alleles that survived
de-duplication and made it into the record, while the argmax ran over `A(A+1)/2` genotypes
built from everything enumerated.

The dump is ~32 M lines for one chromosome, so it is consumed from a FIFO and reduced on
the fly rather than written out. Each site emits a `#site` header whose trailing columns
are `allele_0 .. allele_{n-1}` -- that is where the allele count comes from -- followed by
one row per read, which is where the site name is. The writer holds an OpenMP critical
section for a whole site, so blocks never interleave.
"""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    args = p.parse_args()

    n_alleles = 0
    site = None
    n_reads = 0
    written = 0

    with open(args.out, "w") as out:
        out.write("site\tn_alleles\tn_reads\n")
        for line in sys.stdin:
            if line.startswith("#site"):
                if site is not None:
                    out.write(f"{site}\t{n_alleles}\t{n_reads}\n")
                    written += 1
                # trailing columns after site/read/mismap_prob/best_ln are the alleles
                n_alleles = len(line.rstrip("\n").split("\t")) - 4
                site, n_reads = None, 0
                continue
            f = line.split("\t", 1)
            if not f or not f[0]:
                continue
            if site is None:
                site = f[0]
            n_reads += 1
        if site is not None:
            out.write(f"{site}\t{n_alleles}\t{n_reads}\n")
            written += 1

    print(f"wrote {written:,} sites to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
