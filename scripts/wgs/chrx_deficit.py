#!/usr/bin/env python3
"""Why does chrX score lower than the autosomes? Reproduces the analysis in docs/wgs-results.md.

The headline is that chrX's deficit is a precision deficit made of low-GQ calls, and that a GQ
filter *raises* chrX while *lowering* the autosomes. That asymmetry is the load-bearing measurement:
it says the caller ranks chrX calls correctly and the problem is how much mass sits in the bad bin,
not that the bin is misidentified.

The mechanism is ploidy meeting paralogy. Under ploidy 1 a 50/50 pileup has no genotype that
explains it -- the model must pick one allele and the likelihood gap collapses to nothing -- so
every locus where stray reads from a paralogous copy manufacture a balanced pileup becomes a
coin-flip hom-alt call at GQ 0. On a diploid contig the same pileup is just a heterozygote.

Two caveats worth carrying, because both cut against the easy reading:

  * The matched-depth table is not evidence that chrX is easier than chr7. On chr7 a DP-12 call is
    abnormal -- coverage failed there for a reason -- so chr7's low-depth calls are a self-selected
    hard population. The claim it supports is only the weaker one: chrX is not anomalous once depth
    is controlled for.

  * The k-mer copy-number test must be run on *distinct canonical* k-mers with self-hits excluded by
    position. A first version sampled every 7th base when building the sets and every 7th base when
    streaming, so the offsets did not align and the self-hit counts were nonsense; a second version
    reported the mean external count, which a handful of satellite k-mers dominate. The fraction of
    k-mers with any external copy is the statistic that survives both problems.

Usage:  chrx_deficit.py [--work work/wgs] [--contig chrX] [--against chr7 chr20]
"""

from __future__ import annotations

import argparse
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

K = 31
COMP = str.maketrans("ACGT", "TGCA")

# The two paralogous windows that carry 29% of chrX's false positives. Found by binning FP rate
# along the chromosome at 5 Mb and then 100 kb, not assumed.
HOTSPOTS = [(47_600_000, 47_800_000), (48_700_000, 48_900_000)]
CONTROLS = [(80_000_000, 80_200_000), (140_000_000, 140_200_000)]


def query(vcf: Path, fmt: str, region: str | None = None) -> list[list[str]]:
    cmd = ["bcftools", "query", "-f", fmt]
    if region:
        cmd += ["-r", region]
    cmd.append(str(vcf))
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    return [l.split("\t") for l in out.splitlines()]


def f1(tp: int, fp: int, fn: int) -> float:
    """Symmetric F1 from counts.

    aardvark's own F1 uses truth_tp for recall and query_tp for precision, which differ slightly;
    this one is used for every contig here so the comparisons are internally consistent.
    """
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def load(work: Path, contig: str) -> dict:
    """BD decisions joined to the caller's own FORMAT fields, keyed by position."""
    adir = work / "score" / f"{contig}.aardvark"
    bd = {int(r[0]): r[1] for r in query(adir / "query.vcf.gz", "%POS\t[%BD]\n")}
    truth_bd = [(int(r[0]), r[1]) for r in query(adir / "truth.vcf.gz", "%POS\t[%BD]\n")]
    rows = []
    for r in query(work / contig / f"{contig}.vcf.gz", "%POS\t[%GT\t%DP\t%AD\t%DR\t%GQ]\n"):
        pos = int(r[0])
        def num(i, default=float("nan")):
            try:
                return float(r[i])
            except (ValueError, IndexError):
                return default
        ad = [int(x) for x in r[3].split(",") if x.isdigit()]
        rows.append({
            "pos": pos, "gt": r[1], "dp": num(2), "dr": num(4), "gq": num(5, 0.0),
            # Minor-allele fraction: 0 for a unanimous pileup, ->0.5 for a balanced one.
            "bal": min(ad) / sum(ad) if len(ad) >= 2 and sum(ad) else float("nan"),
            "bd": bd.get(pos),
        })
    return {"rows": rows, "truth_bd": truth_bd}


def med(v):
    v = [x for x in v if x == x]
    return statistics.median(v) if v else float("nan")


def report_gq(data: dict) -> None:
    print("\n== F1 vs GQ threshold (symmetric F1, so columns are comparable) ==")
    contigs = list(data)
    print(f"{'GQ>=':>6} " + " ".join(f"{c:>9}" for c in contigs))
    for thr in (0, 10, 20, 30):
        cells = []
        for c in contigs:
            d = data[c]
            tp = sum(1 for r in d["rows"] if r["bd"] == "TP" and r["gq"] >= thr)
            fp = sum(1 for r in d["rows"] if r["bd"] == "FP" and r["gq"] >= thr)
            truth_total = sum(1 for _, b in d["truth_bd"] if b in ("TP", "FN"))
            cells.append(f"{f1(tp, fp, truth_total - tp):9.4f}")
        print(f"{thr:6d} " + " ".join(cells))
    print("\nThe point is the asymmetry: a GQ filter should raise chrX and lower the autosomes.")

    print("\n== GQ / DP distribution ==")
    for c, d in data.items():
        gq = [r["gq"] for r in d["rows"]]
        dp = [r["dp"] for r in d["rows"] if r["dp"] == r["dp"]]
        low = 100 * sum(1 for g in gq if g < 10) / len(gq)
        print(f"  {c:6s} n={len(gq):7d}  median GQ {med(gq):6.0f}  median DP {med(dp):5.0f}  "
              f"GQ<10 {low:5.1f}%")


def report_balance(d: dict) -> None:
    print("\n== allele balance: a haploid call on a split pileup has no right answer ==")
    def hot(p):
        return any(lo <= p < hi for lo, hi in HOTSPOTS)
    groups = defaultdict(list)
    for r in d["rows"]:
        if r["bd"] not in ("TP", "FP"):
            continue
        groups[("hotspot" if hot(r["pos"]) else "background", r["bd"])].append(r)
    print(f"{'region':12s} {'cls':4s} {'n':>7s} {'medDP':>7s} {'medDR':>7s} "
          f"{'med minor-AF':>13s} {'%>0.3':>7s}")
    for k in sorted(groups):
        v = groups[k]
        bals = [r["bal"] for r in v if r["bal"] == r["bal"]]
        pb = 100 * sum(1 for b in bals if b > 0.3) / len(bals) if bals else float("nan")
        print(f"{k[0]:12s} {k[1]:4s} {len(v):7d} {med([r['dp'] for r in v]):7.1f} "
              f"{med([r['dr'] for r in v]):7.2f} {med([r['bal'] for r in v]):13.3f} {pb:6.1f}%")


def report_hotspots(d: dict) -> None:
    """Note the TP here is aardvark's *truth-side* count, so the whole-chrX F1 printed below
    (0.9364) matches aardvark's summary, while report_gq's GQ>=0 row (0.9334) uses the query-side
    count -- it has to, since GQ is a property of the query record. The two differ by the usual
    truth_tp/query_tp gap and neither is wrong; do not read the difference as drift."""
    print("\n== hotspot contribution ==")
    def hot(p):
        return any(lo <= p < hi for lo, hi in HOTSPOTS)
    tp = sum(1 for p, b in d["truth_bd"] if b == "TP")
    fn = sum(1 for p, b in d["truth_bd"] if b == "FN")
    fp = sum(1 for r in d["rows"] if r["bd"] == "FP")
    htp = sum(1 for p, b in d["truth_bd"] if b == "TP" and hot(p))
    hfn = sum(1 for p, b in d["truth_bd"] if b == "FN" and hot(p))
    hfp = sum(1 for r in d["rows"] if r["bd"] == "FP" and hot(r["pos"]))
    span = sum(hi - lo for lo, hi in HOTSPOTS)
    print(f"  whole chrX      TP {tp:6d} FP {fp:5d} FN {fn:5d}  F1 {f1(tp, fp, fn):.4f}")
    print(f"  two hotspots    TP {htp:6d} FP {hfp:5d} FN {hfn:5d}  "
          f"({span/1e6:.1f} Mb, {100*hfp/fp:.0f}% of all FPs)")
    print(f"  chrX minus them TP {tp-htp:6d} FP {fp-hfp:5d} FN {fn-hfn:5d}  "
          f"F1 {f1(tp-htp, fp-hfp, fn-hfn):.4f}")


def report_kmers(work: Path, contig: str) -> None:
    """Copy number of the hotspots' k-mers across chrX and chrY.

    Distinct canonical k-mers, self-hits excluded by position, and reported as the *fraction* with
    any external copy -- the mean external count is dominated by a few satellite k-mers.
    """
    def seq_of(p):
        return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()

    def canon(km):
        r = km.translate(COMP)[::-1]
        return km if km <= r else r

    x = seq_of(work / contig / f"{contig}.fa")
    others = [(f.stem, seq_of(f)) for f in [work / "chrY" / "chrY.fa"] if f.exists()]
    regions = {f"hot {lo/1e6:.1f}-{hi/1e6:.1f}Mb": (lo, hi) for lo, hi in HOTSPOTS}
    regions.update({f"ctl {lo/1e6:.1f}-{hi/1e6:.1f}Mb": (lo, hi) for lo, hi in CONTROLS})

    sets = {}
    for n, (lo, hi) in regions.items():
        s = x[lo:hi]
        sets[n] = {canon(s[i:i+K]) for i in range(len(s)-K+1) if "N" not in s[i:i+K]}
    ext = {n: dict.fromkeys(sets[n], 0) for n in regions}
    for src, seqn in [(contig, x)] + others:
        for i in range(len(seqn)-K+1):
            km = seqn[i:i+K]
            if "N" in km:
                continue
            c = canon(km)
            for n, (lo, hi) in regions.items():
                if src == contig and lo <= i < hi:
                    continue                      # self-hit, excluded by position not by count
                if c in ext[n]:
                    ext[n][c] += 1
    print(f"\n== 31-mer copy number vs {contig}"
          + (" + " + " + ".join(n for n, _ in others) if others else "") + " ==")
    print(f"{'region':22s} {'distinct/200kb':>14s} {'% with external copy':>21s}")
    for n in regions:
        v = list(ext[n].values())
        print(f"{n:22s} {len(v):14d} {100*sum(1 for q in v if q>0)/len(v):20.1f}%")
    print("Hotspots should run 3-4x the control paralogy rate; a low distinct count is a tandem "
          "array.")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--work", default="work/wgs", type=Path)
    p.add_argument("--contig", default="chrX")
    p.add_argument("--against", nargs="*", default=["chr7", "chr20"])
    p.add_argument("--skip-kmers", action="store_true",
                   help="the k-mer pass streams both chromosomes and takes a few minutes")
    args = p.parse_args()

    data = {c: load(args.work, c) for c in [args.contig] + args.against}
    report_gq(data)
    report_balance(data[args.contig])
    report_hotspots(data[args.contig])
    if not args.skip_kmers:
        report_kmers(args.work, args.contig)


if __name__ == "__main__":
    main()
