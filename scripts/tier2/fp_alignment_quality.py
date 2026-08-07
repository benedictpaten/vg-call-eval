#!/usr/bin/env python3
"""H4 of plan §9.19: is the alignment worse at the sites that go wrong?

Globally the two read sets are near-identical -- mean MAPQ 55.31 vs 54.65, and the
32-haplotype alignments are actually *better* on identity and divergence. That rules out
a broad mapping regression but says nothing about whether the damage is concentrated at
the sites that produce spurious calls, which is the question that matters.

Method. Each VCF record carries an `AT` field listing the node IDs of every traversal, so
a site's node set is available without touching the graph. Aardvark's annotated query VCF
says which records are TP and which FP. So: collect the node sets of FP SNV sites and of
a matched sample of TP SNV sites, stream the chromosome's GAF once, and accumulate the
mapping quality of every read touching each set.

The comparison that matters is FP-vs-TP *within* one graph, then how that gap differs
between graphs. Comparing FP sites across graphs directly would confound "these sites are
hard" with "this graph is harder", because the two FP sets are largely different sites.
"""

from __future__ import annotations

import argparse
import gzip
import random
import re
import subprocess
import sys
from pathlib import Path

NODE_RE = re.compile(r"\d+")


def decisions(aardvark_dir: Path) -> dict[int, str]:
    """POS -> BD for single-base substitutions only."""
    out: dict[int, str] = {}
    with gzip.open(aardvark_dir / "query.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            alt = f[4].split(",")[0]
            if alt.startswith("<") or alt == "*" or len(f[3]) != 1 or len(alt) != 1:
                continue
            bd = dict(zip(f[8].split(":"), f[9].split(":"))).get("BD")
            if bd in ("TP", "FP"):
                out[int(f[1])] = bd
    return out


def site_nodes(vcf: Path, wanted: set[int]) -> dict[int, set[int]]:
    """POS -> set of node IDs, read from the AT field."""
    q = subprocess.run(["bcftools", "query", "-f", "%POS\t%INFO/AT\n", str(vcf)],
                       capture_output=True, text=True)
    out: dict[int, set[int]] = {}
    for line in q.stdout.splitlines():
        pos, _, at = line.partition("\t")
        try:
            p = int(pos)
        except ValueError:
            continue
        if p not in wanted or at in (".", ""):
            continue
        out[p] = {int(n) for n in NODE_RE.findall(at)}
    return out


def scan_gaf(gaf: Path, groups: dict[str, set[int]]) -> dict[str, dict]:
    """One pass over the GAF, accumulating MAPQ stats for reads touching each node set."""
    # The MAPQ<10 bucket is called out because that is exactly where --mismap-max binds:
    # phred 10 is p=0.1, the default cap. Below it the cap is what decides how much a
    # read is discounted, not its MAPQ.
    stats = {k: {"n": 0, "sum": 0, "q60": 0, "q0": 0, "lt10": 0, "lt30": 0,
                 "id_sum": 0.0, "id_n": 0} for k in groups}
    opener = gzip.open if str(gaf).endswith(".gz") else open
    with opener(gaf, "rt") as fh:
        for line in fh:
            f = line.split("\t", 12)
            if len(f) < 12:
                continue
            try:
                mapq = int(f[11])
                matches, blocklen = int(f[9]), int(f[10])
            except ValueError:
                continue
            nodes = {int(n) for n in NODE_RE.findall(f[5])}
            for key, wanted in groups.items():
                if nodes & wanted:
                    s = stats[key]
                    s["n"] += 1
                    s["sum"] += mapq
                    if mapq == 60:
                        s["q60"] += 1
                    if mapq == 0:
                        s["q0"] += 1
                    if mapq < 10:
                        s["lt10"] += 1
                    if mapq < 30:
                        s["lt30"] += 1
                    if blocklen > 0:
                        s["id_sum"] += matches / blocklen
                        s["id_n"] += 1
    return stats


def report(label: str, stats: dict[str, dict]) -> None:
    print(f"\n=== {label} ===")
    print(f"{'site set':<10}{'reads':>12}{'mean MAPQ':>11}{'%MAPQ60':>10}{'%MAPQ<30':>10}"
          f"{'%MAPQ<10':>10}{'%MAPQ0':>9}{'identity':>10}")
    for key, s in stats.items():
        if not s["n"]:
            print(f"{key:<10}{'no reads':>12}")
            continue
        print(f"{key:<10}{s['n']:>12,}{s['sum']/s['n']:>11.2f}{100*s['q60']/s['n']:>10.2f}"
              f"{100*s['lt30']/s['n']:>10.2f}{100*s['lt10']/s['n']:>10.2f}"
              f"{100*s['q0']/s['n']:>9.2f}"
              f"{s['id_sum']/s['id_n'] if s['id_n'] else 0:>10.5f}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--work", required=True, help="tier2 work dir")
    p.add_argument("--arm", default="readlik-z")
    p.add_argument("--gaf", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--tp-sample", type=int, default=2000)
    p.add_argument("--seed", type=int, default=11)
    args = p.parse_args()

    work = Path(args.work)
    dec = decisions(work / "results" / f"aardvark-{args.arm}")
    fps = [p for p, b in dec.items() if b == "FP"]
    tps = [p for p, b in dec.items() if b == "TP"]
    random.Random(args.seed).shuffle(tps)
    # Sampled rather than exhaustive: TPs outnumber FPs ~45:1, and a read touching any
    # node of any TP site would otherwise sweep in most of the chromosome, which would
    # measure the chromosome rather than the sites.
    tps = tps[:args.tp_sample]
    print(f"{args.label}: {len(fps):,} FP SNV sites, {len(tps):,} sampled TP SNV sites",
          file=sys.stderr, flush=True)

    nodes = site_nodes(work / "results" / f"{args.arm}.vcf.gz", set(fps) | set(tps))
    groups = {
        "FP": set().union(*(nodes[p] for p in fps if p in nodes)) if fps else set(),
        "TP": set().union(*(nodes[p] for p in tps if p in nodes)) if tps else set(),
    }
    # A node shared by both sets would be counted for both and blur the contrast.
    both = groups["FP"] & groups["TP"]
    groups["FP"] -= both
    groups["TP"] -= both
    print(f"  nodes: FP {len(groups['FP']):,}  TP {len(groups['TP']):,}  "
          f"(dropped {len(both):,} shared)", file=sys.stderr, flush=True)

    report(args.label, scan_gaf(Path(args.gaf), groups))


if __name__ == "__main__":
    main()
