#!/usr/bin/env python3
"""Stage 0 of the depth term: predict what it would do, before writing any of it.

The model scores P(reads | G) conditioned on the reads it was handed, and never asks
whether that many reads should be there. A complete generative model factorises as

    P(data | G) = P(N | G) * P(reads | N, G)

so depth enters as an *additive* term, not a filter:

    ln P(data | G) = ln Poisson(N ; lambda_G) + sum_r ln[ (1-e_r) * sum_h w_h rel(r,h) + e_r ]

    lambda_G = c * sum_{h in G} (L_h + R - 1)

**The footprint is used two different ways, and they are not the same quantity.** The
mixture weight `w_h` wants sequence *unique* to an allele, because only reads over unique
sequence can separate two genotypes -- reads in shared sequence fit everything and cancel.
`lambda_G` wants the *whole* traversal length, because every base of the haplotype
generates reads whether or not those reads discriminate. Confusing the two would be an
easy and invisible mistake, so they are computed separately here.

`c` is calibrated as the median of N / sum_h(L_h + R - 1) over the called genotype across
all dumped sites. Calibrating against the same N the model actually sees matters: N is
"rows in the likelihood matrix", which depends on the read-fetch window and the placement
filter, and is not the same as coverage. A `c` derived from a pack file would be in
different units and would silently mis-scale every site.

What this script decides. Two questions, and the second is the kill criterion:

  1. Do the large heterozygous deletions the caller currently misses flip to correct?
  2. Do the collapsed-repeat pile-ups get *worse*? lambda_G grows with allele length, so
     an anomalously large N mechanically favours whichever genotype presents the most
     sequence. That is a preference for long alleles, not a rejection of the site, and it
     is the way this term could do harm.

Approximation, stated. The read-side margins here use whole-traversal mixture weights
rather than the shipped unique-content weights, because the dump records no allele
identity and reconstructing unique node content needs per-node lengths the dump does not
carry. The difference between the two weightings is single-digit nats; the depth term is
tens to hundreds, so it does not change any conclusion below. It would need fixing before
this became a measurement rather than a prediction.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import math
import re
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"


def load_traversals(gaf: Path) -> dict[str, dict[int, int]]:
    """snarl id -> {allele index: traversal length in bp}."""
    out: dict[str, dict[int, int]] = collections.defaultdict(dict)
    with open(gaf) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.split("\t")
            parts = f[0].split("#")
            snarl = next((p for p in parts if p.startswith(">")), None)
            if snarl is None:
                continue
            out[snarl][int(parts[-1])] = int(f[1])
    return out


def load_dump(path: Path) -> dict[str, list[list[str]]]:
    sites: dict[str, list] = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            sites[f[0]].append(f)
    return sites


def dump_key_to_snarl(key: str) -> str:
    """`170518348+_170519407+` -> `>170518348>170519407`."""
    n = key.replace("+", "").split("_")
    return ">" + n[0] + ">" + n[-1]


def ln_poisson(n: int, lam: float) -> float:
    if lam <= 0:
        return -math.inf
    # Stirling for the factorial; n is in the hundreds here so it is exact enough.
    if n == 0:
        return -lam
    return n * math.log(lam) - lam - (n * math.log(n) - n + 0.5 * math.log(2 * math.pi * n))


def genotypes(n_alleles: int):
    return [(a, b) for b in range(n_alleles) for a in range(b + 1)]


def read_ll(rel, e, G, weights):
    total = 0.0
    for i in range(len(rel)):
        mix = sum(weights[k] * rel[i][h] for k, h in enumerate(G))
        total += math.log((1 - e[i]) * mix + e[i])
    return total


def footprint_weights(G, lengths, R):
    eff = [max(lengths.get(h, 0) + R - 1, 1.0) for h in G]
    s = sum(eff)
    return [x / s for x in eff], s


def analyse(dump: Path, trav: Path, R: float, label: str, depth_weights):
    sites = load_dump(dump)
    travs = load_traversals(trav)

    # ---- calibrate c against the called genotype, in the same units as N ----
    ratios = []
    per_site = {}
    for key, rows in sites.items():
        snarl = dump_key_to_snarl(key)
        L = travs.get(snarl)
        if not L:
            continue
        na = len(rows[0]) - 4
        if na < 2 or any(a not in L for a in range(na)):
            continue
        rel = [[float(r[4 + a]) for a in range(na)] for r in rows]
        e = [float(r[2]) for r in rows]
        N = len(rows)
        best, best_ll = None, -math.inf
        for G in genotypes(na):
            w, _ = footprint_weights(G, L, R)
            ll = read_ll(rel, e, G, w)
            if ll > best_ll:
                best, best_ll = G, ll
        _, sum_called = footprint_weights(best, L, R)
        ratios.append(N / sum_called)
        per_site[key] = (snarl, L, rel, e, N, na, best, best_ll)

    if not ratios:
        print(f"{label}: no usable sites")
        return {}, 0.0
    c = statistics.median(ratios)
    print(f"{label}: {len(per_site)} sites, calibrated c = {c:.4f} reads per position "
          f"(IQR {statistics.quantiles(ratios, n=4)[0]:.3f}-{statistics.quantiles(ratios, n=4)[2]:.3f})")
    return per_site, c


def combined_best(L, rel, e, N, na, c, R, dw):
    best, best_score, parts = None, -math.inf, None
    for G in genotypes(na):
        w, s = footprint_weights(G, L, R)
        rl = read_ll(rel, e, G, w)
        dl = ln_poisson(N, c * s)
        score = rl + dw * dl
        if score > best_score:
            best, best_score, parts = G, score, (rl, dl, s)
    return best, parts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--read-length", type=float, default=151.0)
    ap.add_argument("--depth-weights", nargs="+", type=float, default=[0.0, 0.25, 0.5, 1.0])
    args = ap.parse_args()

    # ---------- chr6: the heterozygous deletions the caller misses ----------
    per_site, c = analyse(WORK / "sv-atlas/chr6-large.v2.dump.tsv",
                          WORK / "sv-atlas/chr6-trav.gaf",
                          args.read_length, "chr6-4hap", args.depth_weights)

    # Snarls of large heterozygous deletions poisson-z recovers and readlik-z misses.
    truth = list(csv.DictReader(open(WORK / "sv-atlas/truth.tsv"), delimiter="\t"))
    def key(r):
        return (r["chrom"], r["pos"], r["svlen"])
    P = {key(r): r for r in truth if r["dataset"] == "chr6-4hap" and r["arm"] == "poisson-z"
         and r["svtype"] == "DEL" and r["sizebin"] == "1k+"}
    Rz = {key(r): r for r in truth if r["dataset"] == "chr6-4hap" and r["arm"] == "readlik-z"
          and r["svtype"] == "DEL" and r["sizebin"] == "1k+"}
    disc = {int(k[1]) for k in P if P[k]["outcome"] == "TP" and Rz[k]["outcome"] == "FN"}
    want = set()
    with gzip.open(WORK / "tier2-chr6/results/poisson-z.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if int(f[1]) in disc and f[2].startswith(">") and len(f[4]) - len(f[3]) <= -1000:
                n = f[2].strip(">").split(">")
                want.add(f"{n[0]}+_{n[-1]}+")

    rows = [(k, v) for k, v in per_site.items() if k in want]
    print(f"\n=== Q1: {len(rows)} large heterozygous deletions the caller currently misses ===")
    print(f"{'snarl':24s} {'N':>5s} {'called':>8s} " +
          " ".join(f"{'w=' + str(w):>10s}" for w in args.depth_weights))
    flips = {w: 0 for w in args.depth_weights}
    for k, (snarl, L, rel, e, N, na, called, _) in sorted(rows):
        # The deletion allele is the shortest traversal; the correct genotype pairs it
        # with the best-supported long one.
        short = min(L, key=lambda a: L[a])
        cells = []
        for w in args.depth_weights:
            G, _ = combined_best(L, rel, e, N, na, c, args.read_length, w)
            has_del = short in G and len(set(G)) > 1
            flips[w] += has_del
            cells.append(("HET-DEL" if has_del else str(G)))
        print(f"{k:24s} {N:5d} {str(called):>8s} " + " ".join(f"{x:>10s}" for x in cells))
    print("\nsites calling the heterozygous deletion, by depth weight:")
    for w in args.depth_weights:
        print(f"  w={w:<5} {flips[w]:2d}/{len(rows)}")

    # ---------- chr20: the collapsed-repeat pile-ups ----------
    per20, c20 = analyse(WORK / "sv-atlas/chr20-large.v2.dump.tsv",
                         WORK / "sv-atlas/chr20-trav.gaf",
                         args.read_length, "chr20-4hap", args.depth_weights)
    if per20:
        depths = sorted(v[4] for v in per20.values())
        med = depths[len(depths) // 2]
        pileups = [(k, v) for k, v in per20.items() if v[4] > 3 * med]
        print(f"\n=== Q2: {len(pileups)} pile-up sites (N > 3x the median of {med}) ===")
        print("Does the depth term push these toward *more* sequence, i.e. make them worse?")
        print(f"{'snarl':24s} {'N':>6s} {'called':>8s} " +
              " ".join(f"{'w=' + str(w):>12s}" for w in args.depth_weights))
        longer = {w: 0 for w in args.depth_weights}
        for k, (snarl, L, rel, e, N, na, called, _) in sorted(pileups)[:15]:
            _, called_sum = footprint_weights(called, L, args.read_length)
            cells = []
            for w in args.depth_weights:
                G, parts = combined_best(L, rel, e, N, na, c20, args.read_length, w)
                longer[w] += parts[2] > called_sum
                cells.append(f"{str(G)}")
            print(f"{k:24s} {N:6d} {str(called):>8s} " + " ".join(f"{x:>12s}" for x in cells))
        print("\npile-up sites moved to a LONGER-footprint genotype, by depth weight:")
        for w in args.depth_weights:
            print(f"  w={w:<5} {longer[w]:2d}/{len(pileups)}")


if __name__ == "__main__":
    main()
