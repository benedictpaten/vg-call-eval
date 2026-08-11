#!/usr/bin/env python3
"""How often does the call set imply a recombination the haplotype panel has never seen?

`vg call` genotypes each snarl independently. The emitted call set is the concatenation of
per-site argmaxes, which corresponds to a pair of haplotypes free to switch panel haplotype
at every site, at no cost. Nothing in the objective notices when consecutive called alleles
are carried by no single panel haplotype.

This measures how often that happens, before any HMM exists to fix it. The prize is bounded
by two things and both are measured here: how often an apparent recombination occurs at all,
and whether it occurs where the reads were undecided -- because a linkage prior can only pay
where the emission is flat. If apparent recombinations cluster at high-`GQ` sites, the reads
already disagree with the panel and a prior that overrides them would make things worse.

Method. `vg deconstruct` on the same graph gives one column per panel haplotype and uses the
same snarl IDs as `vg call`, so the panel matrix joins to the call set on the ID column with
no coordinate matching. For adjacent called sites, take the distinct called alleles at each,
and ask whether any single panel haplotype carries one from each. Zero means no haplotype
explains any combination we called.

Reading the numbers. The panel in a haplotype-sampled GBZ is synthetic: haplotypes are
recombinations of real assemblies chosen in blocks, so linkage within a block is real and
block boundaries carry switches that are an artefact of sampling rather than biology. With
~10 kb blocks and sites ~1 kb apart, most adjacent pairs sit inside a block; the artefact
should show as a rise in apparent recombination with inter-site distance, saturating around
the block length. That stratification is reported, because it is what would calibrate a
position-dependent recombination rate later -- and because a *flat* profile would mean the
signal is not linkage at all.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_panel(vcf: Path) -> dict:
    """snarl ID -> {allele index: frozenset of haplotype (column, phase) pairs}.

    A haplotype is identified by its column and phase, not by name: a sampled GBZ puts every
    recombinant under one sample name, so the phase index is the only thing distinguishing
    them. Missing genotypes are simply absent -- a haplotype that does not traverse the
    snarl carries no allele there, which is different from carrying the reference.
    """
    panel = {}
    with open_maybe_gz(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or not f[2].startswith(">"):
                continue
            by_allele = collections.defaultdict(set)
            for col, cell in enumerate(f[9:]):
                gt = cell.split(":")[0]
                for phase, tok in enumerate(gt.replace("/", "|").split("|")):
                    if tok.isdigit():
                        by_allele[int(tok)].add((col, phase))
            panel[f[2]] = {a: frozenset(h) for a, h in by_allele.items()}
    return panel


def load_calls(vcf: Path) -> list:
    """Called records in reference order: (snarl, pos, distinct called allele indices, GQ)."""
    out = []
    with open_maybe_gz(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or not f[2].startswith(">"):
                continue
            sample = dict(zip(f[8].split(":"), f[9].split(":")))
            gt = sample.get("GT", "")
            alleles = {int(t) for t in gt.replace("/", "|").split("|") if t.isdigit()}
            if not alleles:
                continue
            try:
                gq = float(sample.get("GQ", "nan"))
            except ValueError:
                gq = float("nan")
            out.append((f[2], int(f[1]), alleles, gq))
    return out


DIST_BINS = [(0, 200, "<200"), (200, 1000, "200-1k"), (1000, 5000, "1k-5k"),
             (5000, 20000, "5k-20k"), (20000, 10 ** 12, ">20k")]


def dist_bin(d: int) -> str:
    for lo, hi, name in DIST_BINS:
        if lo <= d < hi:
            return name
    return ">20k"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="a vg call VCF")
    ap.add_argument("--panel", required=True, help="vg deconstruct VCF for the same graph")
    ap.add_argument("--label", required=True)
    args = ap.parse_args()

    panel = load_panel(Path(args.panel))
    calls = load_calls(Path(args.calls))
    joined = [c for c in calls if c[0] in panel]
    print(f"\n=== {args.label} ===")
    print(f"  panel sites {len(panel):,}, called records {len(calls):,}, "
          f"joined on snarl ID {len(joined):,} ({len(joined) / max(len(calls), 1) * 100:.1f}%)")
    n_hap = len({h for site in panel.values() for hs in site.values() for h in hs})
    print(f"  panel haplotypes seen: {n_hap}")

    by_dist = collections.defaultdict(lambda: [0, 0])
    by_gq = collections.defaultdict(lambda: [0, 0])
    shared_counts = []
    lifts = []
    total = broken = 0
    for (s1, p1, a1, gq1), (s2, p2, a2, gq2) in zip(joined, joined[1:]):
        if p2 < p1:
            continue
        best = 0
        best_lift = None
        for x in a1:
            hx = panel[s1].get(x, frozenset())
            if not hx:
                continue
            for y in a2:
                hy = panel[s2].get(y, frozenset())
                inter = len(hx & hy)
                if inter > best:
                    best = inter
                # Enrichment over independence, which is the quantity that decides whether
                # linkage carries information the per-site model is throwing away. A raw
                # count is misleading: if both alleles are carried by half the panel and
                # they are in perfect linkage, the intersection is still half the panel.
                # lift = P(allele at site 2 | allele at site 1) / P(allele at site 2).
                if hx and hy:
                    lift = (inter / len(hx)) / (len(hy) / n_hap)
                    if best_lift is None or lift > best_lift:
                        best_lift = lift
        total += 1
        shared_counts.append(best)
        if best_lift is not None:
            lifts.append(best_lift)
        d = dist_bin(p2 - p1)
        by_dist[d][0] += 1
        by_dist[d][1] += (best == 0)
        # The prize only exists where the reads were undecided.
        g = min(gq1, gq2)
        key = "GQ<10" if g < 10 else ("GQ 10-40" if g < 40 else "GQ>=40")
        by_gq[key][0] += 1
        by_gq[key][1] += (best == 0)
        broken += (best == 0)

    if not total:
        print("  no adjacent pairs joined")
        return
    print(f"\n  adjacent called pairs: {total:,}")
    print(f"  no panel haplotype carries any called combination: {broken:,} "
          f"({broken / total * 100:.1f}%)")
    # How *constraining* the linkage is, which decides whether a prior is worth having.
    # Two of three panel haplotypes co-carrying tells you almost nothing; two of
    # thirty-four is a strong constraint. Reported as a fraction of the panel, so the two
    # graphs are comparable.
    nz = [c for c in shared_counts if c]
    if nz:
        print(f"  co-carrying haplotypes when >0: median {statistics.median(nz):.0f} "
              f"of {n_hap} ({statistics.median(nz) / n_hap * 100:.0f}% of the panel)")
        dist = collections.Counter(min(c, 5) for c in shared_counts)
        print("  " + "  ".join(f"{k if k < 5 else '5+'}:{dist[k] / total * 100:.1f}%"
                               for k in sorted(dist)))
    if lifts:
        lifts.sort()
        print(f"  linkage lift P(b|a)/P(b): median {statistics.median(lifts):.2f}x, "
              f"p25 {lifts[len(lifts) // 4]:.2f}x, p75 {lifts[3 * len(lifts) // 4]:.2f}x")
        print(f"    pairs where the called combination is *more* likely than independence: "
              f"{sum(1 for x in lifts if x > 1.05) / len(lifts) * 100:.1f}%; "
              f"at or below independence: {sum(1 for x in lifts if x <= 1.05) / len(lifts) * 100:.1f}%")

    print(f"\n  {'distance':>10s} {'pairs':>9s} {'apparent recomb':>17s}")
    for _, _, name in DIST_BINS:
        n, b = by_dist[name]
        if n:
            print(f"  {name:>10s} {n:>9,} {b / n * 100:>16.1f}%")
    print(f"\n  {'min GQ':>10s} {'pairs':>9s} {'apparent recomb':>17s}")
    for name in ("GQ<10", "GQ 10-40", "GQ>=40"):
        n, b = by_gq[name]
        if n:
            print(f"  {name:>10s} {n:>9,} {b / n * 100:>16.1f}%")

    dest = WORK / "sv-atlas" / f"apparent-recombination-{args.label}.json"
    dest.write_text(json.dumps(
        {"pairs": total, "broken": broken,
         "by_distance": {k: v for k, v in by_dist.items()},
         "by_gq": {k: v for k, v in by_gq.items()}}, indent=2))
    print(f"\n  wrote {dest}")


if __name__ == "__main__":
    main()
