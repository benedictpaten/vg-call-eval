#!/usr/bin/env python3
"""Why the read-likelihood model loses large heterozygous deletions.

The claim. At a heterozygous deletion the reads that lie *inside* the deleted interval
come from the intact haplotype and fit a long allele; the only reads supporting the
deletion are the few that span its junction. Under

    ln P(reads | G) = sum_r ln[ (1 - e_r) * mean_{h in G} rel(r,h) + e_r ]

an interior read scores `1` under a homozygous long genotype and `0.5` under the correct
heterozygote, so it favours the wrong answer by `ln 2 = 0.693` nats. A junction read
favours the heterozygote by `ln(0.5 / e_r)`, capped at `ln(0.5 / mismap-min) = 3.22` nats
at the shipped floor of 0.02 -- and by *zero* for a read at the `--mismap-max` cap.

Interior reads grow with deletion length; junction reads do not. Break-even is at
`3.22 / 0.693 = 4.64` junction-equivalents, so with 151 bp reads the model should start
losing heterozygous deletions around `4.64 * 151 = 701 bp` and lose them badly above
1 kb. That is exactly the observed shape: het deletion recall 0.586 -> 0.395 in the
300-999 bin and 0.670 -> 0.064 above 1 kb, against the Poisson caller.

This script tests the claim quantitatively rather than by eye. For each site it splits
the reads into interior and junction, predicts the margin between the homozygous long
genotype and the best genotype containing the deletion, and compares that against the
margin the model actually produces.

Identifying the deletion allele. `--dump-likelihoods` records no allele identity -- only
`rel(read, allele)` columns -- and the allele order does *not* agree with the order
`-T/--traversals` reports, so the column cannot be looked up. It is identified from the
data instead: at these sites exactly one allele is fit well by under 15% of reads while
the others are fit by 45% or more, which is the signature of a deletion and is
unambiguous. `-T` independently confirms a deletion-scale traversal is enumerated at
every one of these snarls (for example `[296, 2943, 2944, 2945]` bp where the benchmark
records a 2648 bp deletion), so this is a scoring failure and not a missing allele.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import math
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
WORK = REPO / "work"

LN2 = math.log(2)
WELL_FIT = 0.9      # rel above this counts as "this allele explains the read"
POOR_FIT = 0.1      # rel below this counts as "this allele does not"
DEL_FRAC = 0.15     # an allele fit by fewer than this fraction of reads is the deletion


def load_dump(path: Path, wanted: set[str]) -> dict[str, list[list[str]]]:
    sites: dict[str, list] = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if not wanted or f[0] in wanted:
                sites[f[0]].append(f)
    return sites


def discordant_snarls(contig: str, dataset: str) -> tuple[set[str], dict]:
    """Snarl IDs of large heterozygous deletions poisson-z recovers and readlik misses."""
    truth = list(csv.DictReader(open(WORK / "sv-atlas/truth.tsv"), delimiter="\t"))

    def key(r):
        return (r["chrom"], r["pos"], r["svlen"])

    P = {key(r): r for r in truth if r["dataset"] == dataset and r["arm"] == "poisson-z"
         and r["svtype"] == "DEL" and r["sizebin"] == "1k+"}
    R = {key(r): r for r in truth if r["dataset"] == dataset and r["arm"] == "readlik"
         and r["svtype"] == "DEL" and r["sizebin"] == "1k+"}
    disc = {int(k[1]) for k in P if P[k]["outcome"] == "TP" and R[k]["outcome"] == "FN"}

    sub = {"chr6-4hap": "tier2-chr6", "chr20-4hap": "tier2-chr20"}[dataset]
    ids, meta = set(), {}
    with gzip.open(WORK / sub / "results/poisson-z.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if int(f[1]) in disc and f[2].startswith(">") and len(f[4]) - len(f[3]) <= -1000:
                n = f[2].strip(">").split(">")
                sid = f"{n[0]}+_{n[-1]}+"
                ids.add(sid)
                meta[sid] = (int(f[1]), len(f[4]) - len(f[3]))
    return ids, meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", default=str(WORK / "sv-atlas/chr6-large.dump.tsv"))
    ap.add_argument("--dataset", default="chr6-4hap")
    ap.add_argument("--contig", default="chr6")
    ap.add_argument("--mismap-min", type=float, default=0.02,
                    help="only used for the printed break-even length")
    args = ap.parse_args()

    ids, meta = discordant_snarls(args.contig, args.dataset)
    sites = load_dump(Path(args.dump), ids)
    if not sites:
        print("no discordant snarls found in the dump; nothing to test")
        return

    cap = math.log(0.5 / args.mismap_min)
    print(f"interior read favours the wrong genotype by  {LN2:.3f} nats")
    print(f"junction read favours the right one by up to {cap:.3f} nats "
          f"(--mismap-min {args.mismap_min})")
    print(f"break-even at {cap/LN2:.2f} interior per junction read, "
          f"so ~{cap/LN2*151:.0f} bp with 151 bp reads\n")

    print(f"{'snarl':24s} {'svlen':>7s} {'reads':>6s} {'interior':>9s} {'junction':>9s} "
          f"{'observed':>9s} {'predicted':>10s} {'ratio':>6s}")
    ratios = []
    for s in sorted(sites):
        rows = sites[s]
        na = len(rows[0]) - 4
        rel = [[float(r[4 + a]) for a in range(na)] for r in rows]
        e = [float(r[2]) for r in rows]

        def LL(G):
            return sum(math.log((1 - e[i]) * sum(rel[i][h] for h in G) / len(G) + e[i])
                       for i in range(len(rel)))

        fit = [sum(1 for r in rel if r[a] > WELL_FIT) for a in range(na)]
        d = min(range(na), key=lambda a: fit[a])
        if fit[d] > DEL_FRAC * len(rows):
            print(f"{s:24s}  no allele looks like a deletion; skipped")
            continue
        long_a = max((a for a in range(na) if a != d), key=lambda a: fit[a])

        observed = LL([long_a, long_a]) - LL([d, long_a])
        interior = [i for i, r in enumerate(rel) if r[long_a] > WELL_FIT and r[d] < POOR_FIT]
        junction = [i for i, r in enumerate(rel) if r[d] > WELL_FIT]
        predicted = len(interior) * LN2 - sum(math.log(0.5 / e[i]) for i in junction)
        ratio = observed / predicted if predicted else float("nan")
        ratios.append(ratio)
        svlen = meta.get(s, (0, 0))[1]
        print(f"{s:24s} {svlen:7d} {len(rows):6d} {len(interior):9d} {len(junction):9d} "
              f"{observed:9.1f} {predicted:10.1f} {ratio:6.2f}")

    good = [r for r in ratios if 0.85 <= r <= 1.15]
    print(f"\n{len(good)}/{len(ratios)} sites within 15% of the two-term prediction")
    print("A margin *smaller* than predicted means something partly offsets the interior")
    print("reads -- reads that fit both alleles in part, so they neither vote cleanly nor")
    print("cancel. Those sites are the ones the model comes closest to getting right.")


if __name__ == "__main__":
    main()
