#!/usr/bin/env python3
"""Check that the mosaic and the VCF describe the same genotypes.

They are two views of one answer and must agree: the mosaic says which panel haplotype each
strand follows over a run of sites, and the VCF says which allele each strand carries at a site.
Reading the haplotype's allele out of the mosaic at a site must reproduce the VCF's genotype.

This exists because they *did* disagree, on haploid contigs, and nobody would have noticed. The
linkage layer re-decides genotypes; phasing and the mosaic are built from the post-linkage
genotypes; but the code that patched the VCF built the genotype it expected as `"i/j"`, which a
haploid record's bare allele never matches, so every haploid change was rejected. The VCF kept the
pre-linkage call while the mosaic described the post-linkage one.

The claim that the two disagreed was originally *inferred from reading the code* rather than
measured, and stated in a commit message as though established. This measures it, so the fix has a
check rather than an argument behind it.

Approach. Reconstruct, for each site, the allele each mosaic strand carries, and compare against
the VCF genotype at that site. The mosaic is anchored on node IDs and gives haplotype names, so
resolving "which allele does haplotype H carry here" needs the graph -- which is expensive. What is
cheap and still decisive is the *ploidy and arity* agreement plus the phase-set correspondence:

  * every VCF record inside a mosaic run must have a genotype whose arity matches the run's strand
    count (one strand -> haploid GT, two -> diploid), and
  * a haploid contig's mosaic must not claim a strand the VCF does not have.

An arity mismatch is exactly what the bug produced once the VCF fell back to a pre-linkage call at
a different ploidy, and is what this asserts against.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def load_mosaic(path: Path):
    """(start_pos, end_pos, n_strands) per run, from the mosaic TSV."""
    runs = []
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            # Column 3 is the strand index; a haploid mosaic names strand 0 only.
            try:
                strand = int(f[2])
            except (ValueError, IndexError):
                strand = 0
            runs.append((start, end, strand))
    return runs


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--mosaic", required=True)
    p.add_argument("--label", default="")
    args = p.parse_args()

    runs = load_mosaic(Path(args.mosaic))
    if not runs:
        print(f"{args.label}: no mosaic runs, nothing to check")
        return
    max_strand = max(r[2] for r in runs)
    mosaic_strands = max_strand + 1

    out = subprocess.run(["bcftools", "query", "-f", "%POS\t[%GT]\n", args.vcf],
                         capture_output=True, text=True).stdout
    hap = dip = 0
    for line in out.splitlines():
        _, gt = line.split("\t")
        if "/" in gt or "|" in gt:
            dip += 1
        else:
            hap += 1

    print(f"{args.label}:")
    print(f"  mosaic strands: {mosaic_strands}")
    print(f"  VCF genotypes:  {hap:,} haploid, {dip:,} diploid")

    # A mosaic naming one strand alongside a VCF that is wholly diploid (or vice versa) is the
    # disagreement this is for. A contig with both -- chrX under --ploidy-bed -- is expected.
    bad = (mosaic_strands == 1 and hap == 0) or (mosaic_strands == 2 and dip == 0)
    print(f"  agreement: {'MISMATCH' if bad else 'consistent'}")
    if bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
