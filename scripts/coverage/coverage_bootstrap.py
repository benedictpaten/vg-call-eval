#!/usr/bin/python3
"""Confidence intervals for docs/coverage.md: paired 1 Mb block bootstrap of F1, recall and precision
deltas between two points of the coverage sweep, pooled over chr20 and chr6.

A point is (technology, depth tag), e.g. ("sr", "5x") or ("ont", "full"). Per-variant tags are binned
into 1 Mb blocks keyed (contig, pos // 1 Mb): aardvark's FORMAT/BD on truth.vcf.gz (TP/FN) and
query.vcf.gz (TP/FP) for small variants, truvari's four VCFs for SVs. Every run's block totals are
CHECKED against its summary (aardvark's Snv and JointIndel rows, never the query-only Indel row; the
SV counts in score.<tag>.json) before use, so a block count can never disagree with the tables.

Recall comes from truth TP, precision from query TP (four-count). Blocks are resampled jointly for
both points, so the CI is on the paired difference. Comparing technologies pairs different graphs
(32 vs 16 haplotypes); the interval is on the callsets, not on the reads alone.

Needs numpy (the system /usr/bin/python3 has it).  Usage: coverage_bootstrap.py [--nboot N]
"""
import collections
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path.home() / "PycharmProjects/vg-call-eval"
BLOCK, SEED = 1_000_000, 20260923
CONTIGS = ("chr20", "chr6")
CLASSES = ("ALL", "SNV", "INDEL", "SV")

# (name, point a, point b, class): each reports b - a
COMPARISONS = [
    ("short small variants 5x->full", ("sr", "5x"), ("sr", "full"), "ALL"),
    ("short SNV 5x->full", ("sr", "5x"), ("sr", "full"), "SNV"),
    ("short indel 5x->full", ("sr", "5x"), ("sr", "full"), "INDEL"),
    ("short SV 5x->full", ("sr", "5x"), ("sr", "full"), "SV"),
    ("ONT small variants 25x->full", ("ont", "25x"), ("ont", "full"), "ALL"),
    ("ONT SNV 5x->full", ("ont", "5x"), ("ont", "full"), "SNV"),
    ("ONT indel 15x->full", ("ont", "15x"), ("ont", "full"), "INDEL"),
    ("ONT indel 20x->full", ("ont", "20x"), ("ont", "full"), "INDEL"),
    ("ONT indel 25x->full", ("ont", "25x"), ("ont", "full"), "INDEL"),
    ("ONT SV 5x->full", ("ont", "5x"), ("ont", "full"), "SV"),
    ("ONT SV 15x->full", ("ont", "15x"), ("ont", "full"), "SV"),
    ("short->ONT SNV 10x", ("sr", "10x"), ("ont", "10x"), "SNV"),
    ("short->ONT SNV 15x", ("sr", "15x"), ("ont", "15x"), "SNV"),
    ("short->ONT SNV 20x", ("sr", "20x"), ("ont", "20x"), "SNV"),
    ("short->ONT SNV full", ("sr", "full"), ("ont", "full"), "SNV"),
    ("short->ONT SV 10x", ("sr", "10x"), ("ont", "10x"), "SV"),
    ("short->ONT SV full", ("sr", "full"), ("ont", "full"), "SV"),
    # appended, so the seeds (SEED + index) of the rows above do not move
    ("ONT SV 5x->15x", ("ont", "5x"), ("ont", "15x"), "SV"),
    ("ONT SV 10x->15x", ("ont", "10x"), ("ont", "15x"), "SV"),
    ("short SV 10x->full", ("sr", "10x"), ("sr", "full"), "SV"),
]
# The same test on one contig at a time, for claims made per contig.
PER_CONTIG = [
    ("ONT indel 20x->full", ("ont", "20x"), ("ont", "full"), "INDEL"),
]


def bq(path, fmt):
    r = subprocess.run(["bcftools", "query", "-f", fmt, str(path)], capture_output=True, text=True, check=True)
    return r.stdout.splitlines()


def is_snv(ref, alt):
    return len(ref) == 1 and all(len(a) == 1 for a in alt.split(",") if a not in (".", "*"))


def run_blocks(tech, contig, tag):
    """(contig, block) -> class -> [tpb, fn, tpc, fp], checked against the run's summaries."""
    label = f"cov-{tech}-{contig}-{tag}"
    res = REPO / f"work/tier2-{contig}-hap32/results"
    B = collections.defaultdict(lambda: {c: np.zeros(4, dtype=np.int64) for c in CLASSES})
    for side, fname, tags in ((0, "truth.vcf.gz", ("TP", "FN")), (2, "query.vcf.gz", ("TP", "FP"))):
        for line in bq(res / f"aardvark-{label}" / fname, "%POS\t%REF\t%ALT[\t%BD]\n"):
            f = line.split("\t")
            if len(f) < 4 or f[3] not in tags:
                continue
            k, i = (contig, int(f[0]) // BLOCK), side + tags.index(f[3])
            B[k]["ALL"][i] += 1
            B[k]["SNV" if is_snv(f[1], f[2]) else "INDEL"][i] += 1
    for i, name in enumerate(("tp-base", "fn", "tp-comp", "fp")):
        for line in bq(res / f"truvari-{label}" / f"{name}.vcf.gz", "%POS\n"):
            if line:
                B[(contig, int(line) // BLOCK)]["SV"][i] += 1

    want = {}
    for r in csv.DictReader(open(res / f"aardvark-{label}" / "summary.tsv"), delimiter="\t"):
        if r["comparison"] == "GT" and r["region_label"] == "ALL" and r["filter"] == "ALL":
            want[r["variant_type"]] = tuple(int(r[k]) for k in ("truth_tp", "truth_fn", "query_tp", "query_fp"))
    sv = json.loads((REPO / f"work/cov/{tech}-{contig}/score.{tag}.json").read_text())["sv"]
    want = {"ALL": want["ALL"], "SNV": want["Snv"], "INDEL": want["JointIndel"],
            "SV": (sv["TP-base"], sv["FN"], sv["TP-comp"], sv["FP"])}
    for c in CLASSES:
        got = tuple(int(x) for x in sum(v[c] for v in B.values()))
        if got != want[c]:
            sys.exit(f"{label} {c}: block totals {got} != summary {want[c]}")
    return B


def metrics(c):
    """c: (..., 4) counts -> recall, precision, F1 (four-count)."""
    tpb, fn, tpc, fp = (c[..., i].astype(float) for i in range(4))
    with np.errstate(invalid="ignore", divide="ignore"):
        r, p = tpb / (tpb + fn), tpc / (tpc + fp)
        return {"F1": np.nan_to_num(2 * r * p / (r + p)), "R": r, "P": p}


def compare(cache, a, b, cls, seed, nboot, contigs=CONTIGS):
    A, B = {}, {}
    for contig in contigs:
        for pt, dst in ((a, A), (b, B)):
            key = (pt[0], contig, pt[1])
            if key not in cache:
                cache[key] = run_blocks(*key)
            dst.update({k: v[cls] for k, v in cache[key].items()})
    keys = sorted(set(A) | set(B))
    xa = np.array([A.get(k, np.zeros(4)) for k in keys])
    xb = np.array([B.get(k, np.zeros(4)) for k in keys])
    idx = np.random.default_rng(seed).integers(0, len(keys), size=(nboot, len(keys)))
    pa, pb = metrics(xa.sum(0)), metrics(xb.sum(0))
    ba, bb = metrics(xa[idx].sum(1)), metrics(xb[idx].sum(1))
    out = {}
    for m in ("F1", "R", "P"):
        lo, hi = np.percentile(bb[m] - ba[m], [2.5, 97.5])
        out[m] = (float(pb[m] - pa[m]), float(lo), float(hi))
    return out, len(keys)


def main():
    nboot = int(sys.argv[sys.argv.index("--nboot") + 1]) if "--nboot" in sys.argv else 10000
    cache = {}
    print(f"# paired 1 Mb block bootstrap, {nboot} replicates, pooled over {'+'.join(CONTIGS)}; b - a")
    for i, (name, a, b, cls) in enumerate(COMPARISONS):
        res, n = compare(cache, a, b, cls, SEED + i, nboot)
        cells = "  ".join(f"d{m} {pt:+.4f} [{lo:+.4f}, {hi:+.4f}] {'SIG ' if lo > 0 or hi < 0 else 'n.s.'}"
                          for m, (pt, lo, hi) in res.items())
        print(f"{name:32} {n:3} blocks  {cells}", flush=True)
    for i, (name, a, b, cls) in enumerate(PER_CONTIG):
        for j, contig in enumerate(CONTIGS):
            res, n = compare(cache, a, b, cls, SEED + 1000 + 10 * i + j, nboot, (contig,))
            cells = "  ".join(f"d{m} {pt:+.4f} [{lo:+.4f}, {hi:+.4f}] {'SIG ' if lo > 0 or hi < 0 else 'n.s.'}"
                              for m, (pt, lo, hi) in res.items())
            print(f"{name + ' ' + contig:32} {n:3} blocks  {cells}", flush=True)


if __name__ == "__main__":
    main()
