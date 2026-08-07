#!/usr/bin/env python3
"""H6 of plan §9.22: re-genotype offline under a prior, from the VCF's own GL field.

`ReadLikelihoodSnarlCaller::genotype` takes the argmax of the raw likelihood with no
genotype prior at all (ties break to reference, nothing else discriminates). §9.22's H5
found the false-SV rate rising steeply with the number of enumerated alleles -- 0.436 at
two alleles to 0.768 at ten or more -- which is what an unprioritised argmax over
`A(A+1)/2` hypotheses should do: the maximum of many noisy scores is biased upward, and
the bias grows with the number of hypotheses.

A Hardy-Weinberg prior is the principled correction, and it is a multiplicity correction
by construction. With reference frequency `1-q` and the non-reference mass `q` spread over
the `A-1` alternates, each *specific* non-reference genotype gets less prior mass as `A`
grows, so a site offering many alternates demands more evidence to leave the reference.
A flat "penalise heterozygotes" term would not do this -- it does not scale with `A`.

Testable with no calling run because `GL` is in the VCF. **Its reach is limited and the
limit matters:** `GL` covers the genotypes over the alleles that *survived into the
record*, not the full enumerated set, so this bounds the achievable effect from below. A
site where 34 alleles were scored but two were emitted is corrected here as a two-allele
site. Read a positive result as "a prior helps, and would help more in the caller"; read a
null result as inconclusive rather than as refuting the prior.
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys


def genotype_index(i: int, j: int) -> int:
    """VCF GL ordering for diploid: index of unordered (i, j) with i <= j."""
    return j * (j + 1) // 2 + i


def log10_prior(i: int, j: int, n_alleles: int, q: float) -> float:
    """Hardy-Weinberg log10 P(genotype), reference frequency 1-q, alternates sharing q."""
    n_alt = max(1, n_alleles - 1)
    f_ref = 1.0 - q
    f_alt = q / n_alt
    fi = f_ref if i == 0 else f_alt
    fj = f_ref if j == 0 else f_alt
    p = fi * fj * (1.0 if i == j else 2.0)
    return math.log10(max(p, 1e-300))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--q", type=float, default=0.01,
                   help="prior non-reference allele frequency; smaller favours reference")
    args = p.parse_args()

    opener = gzip.open if args.vcf.endswith(".gz") else open
    changed = total = no_gl = 0

    with opener(args.vcf, "rt") as fh, open(args.out, "w") as out:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            if line.startswith("#CHROM"):
                out.write(f'##genotypePrior=HWE,q={args.q}\n')
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                out.write(line)
                continue
            fmt = f[8].split(":")
            val = f[9].split(":")
            d = dict(zip(fmt, val))
            gl = d.get("GL")
            n_alleles = 1 + len([a for a in f[4].split(",") if a not in (".",)])
            if not gl or gl == "." or n_alleles < 2:
                no_gl += 1
                out.write(line)
                continue
            try:
                lls = [float(x) for x in gl.split(",")]
            except ValueError:
                no_gl += 1
                out.write(line)
                continue

            best = None
            second = -math.inf
            for j in range(n_alleles):
                for i in range(j + 1):
                    k = genotype_index(i, j)
                    if k >= len(lls):
                        continue
                    s = lls[k] + log10_prior(i, j, n_alleles, args.q)
                    if best is None or s > best[0]:
                        if best is not None:
                            second = best[0]
                        best = (s, i, j)
                    elif s > second:
                        second = s
            if best is None:
                no_gl += 1
                out.write(line)
                continue

            total += 1
            _, i, j = best
            old_gt = d.get("GT", "./.")
            new_gt = f"{i}/{j}"
            if set(old_gt.replace("|", "/").split("/")) != {str(i), str(j)}:
                changed += 1
            d["GT"] = new_gt
            # GQ is the phred gap between best and second best, as in the caller.
            if math.isfinite(second):
                d["GQ"] = str(max(0, min(256, int(round(10.0 * (best[0] - second))))))
            f[9] = ":".join(d.get(k, ".") for k in fmt)
            out.write("\t".join(f) + "\n")

    print(f"q={args.q}: {total:,} records re-genotyped, {changed:,} changed "
          f"({100*changed/total if total else 0:.2f}%), {no_gl:,} untouched (no usable GL)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
