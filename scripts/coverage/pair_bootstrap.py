#!/usr/bin/python3
"""Paired 1 Mb block bootstrap between aardvark results, pooled over the contigs given.

  boot.py <name> <contig>:<aardA>:<aardB> [<contig>:<aardA>:<aardB> ...] [--nboot N]
Reports B - A for ALL, SNV and INDEL (JointIndel counts, reconstructed per record from FORMAT/BD):
F1, recall and precision, each with a 95% interval. Block totals are checked against summary.tsv.
"""
import collections
import csv
import subprocess
import sys

import numpy as np

BLOCK, SEED = 1_000_000, 20260924
CLASSES = ("ALL", "SNV", "INDEL")


def bq(path, fmt):
    return subprocess.run(["bcftools", "query", "-f", fmt, path], capture_output=True, text=True,
                          check=True).stdout.splitlines()


def is_snv(ref, alt):
    return len(ref) == 1 and all(len(a) == 1 for a in alt.split(",") if a not in (".", "*"))


def blocks(contig, aard):
    B = collections.defaultdict(lambda: {c: np.zeros(4, dtype=np.int64) for c in CLASSES})
    for side, fname, tags in ((0, "truth.vcf.gz", ("TP", "FN")), (2, "query.vcf.gz", ("TP", "FP"))):
        for line in bq(f"{aard}/{fname}", "%POS\t%REF\t%ALT[\t%BD]\n"):
            f = line.split("\t")
            if len(f) < 4 or f[3] not in tags:
                continue
            k, i = (contig, int(f[0]) // BLOCK), side + tags.index(f[3])
            B[k]["ALL"][i] += 1
            B[k]["SNV" if is_snv(f[1], f[2]) else "INDEL"][i] += 1
    want = {}
    for r in csv.DictReader(open(f"{aard}/summary.tsv"), delimiter="\t"):
        if r["comparison"] == "GT" and r["region_label"] == "ALL" and r["filter"] == "ALL":
            want[r["variant_type"]] = tuple(int(r[k]) for k in ("truth_tp", "truth_fn", "query_tp", "query_fp"))
    for c, vt in (("ALL", "ALL"), ("SNV", "Snv"), ("INDEL", "JointIndel")):
        got = tuple(int(x) for x in sum(v[c] for v in B.values()))
        if got != want[vt]:
            sys.exit(f"{aard} {c}: block totals {got} != summary {want[vt]}")
    return B


def metrics(c):
    tpb, fn, tpc, fp = (c[..., i].astype(float) for i in range(4))
    with np.errstate(invalid="ignore", divide="ignore"):
        r, p = tpb / (tpb + fn), tpc / (tpc + fp)
        return {"F1": np.nan_to_num(2 * r * p / (r + p)), "R": r, "P": p}


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    nboot = int(sys.argv[sys.argv.index("--nboot") + 1]) if "--nboot" in sys.argv else 10000
    name, specs = args[0], [a for a in args[1:] if ":" in a]
    A, Bb = {}, {}
    for spec in specs:
        contig, a, b = spec.split(":")
        for src, dst in ((a, A), (b, Bb)):
            for k, v in blocks(contig, src).items():
                dst[k] = v
    keys = sorted(set(A) | set(Bb))
    z = {c: np.zeros(4) for c in CLASSES}
    out = [name]
    for c in ("INDEL", "SNV", "ALL"):
        xa = np.array([A.get(k, z)[c] for k in keys]); xb = np.array([Bb.get(k, z)[c] for k in keys])
        idx = np.random.default_rng(SEED).integers(0, len(keys), size=(nboot, len(keys)))
        pa, pb = metrics(xa.sum(0)), metrics(xb.sum(0))
        ba, bb = metrics(xa[idx].sum(1)), metrics(xb[idx].sum(1))
        cells = []
        for m in ("F1", "R", "P"):
            lo, hi = np.percentile(bb[m] - ba[m], [2.5, 97.5])
            pt = float(pb[m] - pa[m])
            cells.append(f"d{m} {pt:+.4f} [{lo:+.4f},{hi:+.4f}]{'*' if lo > 0 or hi < 0 else ''}")
        out.append(f"{c}: A {float(pa['F1']):.4f} B {float(pb['F1']):.4f}  " + "  ".join(cells))
    print("\n  ".join(out), flush=True)


if __name__ == "__main__":
    main()
