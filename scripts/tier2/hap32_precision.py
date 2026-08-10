#!/usr/bin/env python3
"""Why the 34-haplotype graph emits more false structural variants at the same true yield.

The framing matters, because the summary number is misleading. Going 4-hap -> 34-hap,
`readlik-z` loses 0.034 to 0.056 structural-variant F1 -- but recall is flat on chr6
(-0.005) and *better* on chr20 (+0.005). True positives barely move (852 -> 844,
376 -> 380) while false positives rise 35-52% (625 -> 949, 351 -> 476). So the question is
not "is it worse at finding structural variants" but "why does it emit more false ones".

And it is not the genotyper: every arm loses roughly the same amount, and on chr20
`readlik-z` loses the least of the five. Whatever this is lives in the graph-and-alignment
substrate. This script therefore measures the *population* of false calls rather than
comparing callers.

Phase 0 -- where the extra false calls sit: by type, by size, and how many
structural-variant-sized records each graph emitted at all. That last one is a control for
a mechanical explanation: truvari matches with `--pick ac`, so several query records can
compete for one base record and the losers become false positives. If the richer graph
simply emits more records per event, precision falls without the caller being any less
accurate per event.

Phase 1 -- what the largest unexplained bucket actually contains. Truvari annotates each
false positive with the best candidate it considered and rejected, so `PctSeqSimilarity`,
`PctSizeSimilarity` and `StartDistance` can be read straight off the record with no
re-alignment. The four buckets:

  no candidate  no annotation at all -- nothing comparable within the search window
  placement     sequence and size both pass, but the match sits beyond --refdist
  consumed      everything passes and it is still false, because the base record was
                matched to a different query record under --pick ac
  dissimilar    a candidate exists and the sequence or size test fails

`dissimilar` is half of the 34-haplotype false calls and nobody has looked inside it. If
its similarity scores cluster just under the 0.7 threshold, these are the same events
spelled differently and the remedy is representation harmonisation. If they are spread
low, the caller is emitting sequence that is genuinely not there.

Every evaluated call is inside the GIAB confident region, because truvari runs with
`--includebed`. "No truth candidate" therefore does not mean "the benchmark is silent
here" -- it means the benchmark claims to characterise this region and has no record.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
WORK = REPO / "work"

DATASETS = [("chr20-4hap", "tier2-chr20"), ("chr20-34hap", "tier2-chr20-hap32"),
            ("chr6-4hap", "tier2-chr6"), ("chr6-34hap", "tier2-chr6-hap32")]

SIZE_BINS = [(50, 99, "50-99"), (100, 299, "100-299"), (300, 999, "300-999"),
             (1000, 9999, "1k-10k"), (10000, 10 ** 12, "10k+")]

SEQ_PASS = 0.7      # truvari's default --pctseq
SIZE_PASS = 0.7     # --pctsize
REFDIST = 500       # --refdist


def size_bin(n: int) -> str:
    for lo, hi, name in SIZE_BINS:
        if lo <= n <= hi:
            return name
    return "<50"


def info_map(field: str) -> dict:
    out = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        elif kv:
            out[kv] = ""
    return out


def records(path: Path) -> list[dict]:
    out = []
    if not path.exists():
        return out
    with gzip.open(path, "rt") as fh:
        keys = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            ref, alt = f[3], f[4].split(",")[0]
            svlen = len(alt) - len(ref)
            info = info_map(f[7])
            fmt = {}
            if len(f) >= 10:
                fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            out.append({"chrom": f[0], "pos": int(f[1]), "svlen": svlen,
                        "svtype": "DEL" if svlen < 0 else ("INS" if svlen > 0 else "OTHER"),
                        "bin": size_bin(abs(svlen)), "info": info, "fmt": fmt})
    return out


def bucket(r: dict) -> str:
    i = r["info"]
    if "PctSeqSimilarity" not in i and "PctSizeSimilarity" not in i:
        return "no candidate"
    def num(k):
        try:
            return float(i[k])
        except (KeyError, ValueError):
            return None
    seq, size = num("PctSeqSimilarity"), num("PctSizeSimilarity")
    dist = num("StartDistance")
    seq_ok = seq is not None and seq >= SEQ_PASS
    size_ok = size is not None and size >= SIZE_PASS
    if seq_ok and size_ok:
        if dist is not None and abs(dist) > REFDIST:
            return "placement"
        return "consumed"
    return "dissimilar"


def fnum(d: dict, key: str):
    v = d.get(key)
    if v in (None, ".", ""):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def share(fmt: dict):
    dp, ad = fnum(fmt, "DP"), fmt.get("AD")
    if not dp or not ad:
        return None
    try:
        return min(1.0, sum(int(x) for x in ad.split(",")) / dp)
    except ValueError:
        return None


def med(xs):
    xs = [x for x in xs if x is not None]
    return statistics.median(xs) if xs else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="readlik-z")
    args = ap.parse_args()

    data = {}
    for ds, sub in DATASETS:
        d = WORK / sub / "results" / f"truvari-{args.arm}"
        data[ds] = {"fp": records(d / "fp.vcf.gz"),
                    "tp": records(d / "tp-comp.vcf.gz"),
                    "fn": records(d / "fn.vcf.gz")}

    # ---------------- Phase 0a: record inflation, the mechanical control ----------------
    print(f"=== Phase 0a: how many SV-sized query records each graph emitted ({args.arm}) ===")
    print(f"{'dataset':12s} {'TP-comp':>8s} {'FP':>6s} {'records':>8s} "
          f"{'precision':>10s} {'FN':>5s}")
    for ds, _ in DATASETS:
        n_tp, n_fp = len(data[ds]["tp"]), len(data[ds]["fp"])
        print(f"{ds:12s} {n_tp:8d} {n_fp:6d} {n_tp + n_fp:8d} "
              f"{n_tp / (n_tp + n_fp) if n_tp + n_fp else 0:10.4f} {len(data[ds]['fn']):5d}")

    # ---------------- Phase 0b: where the extra false calls sit ----------------
    for key, label in (("svtype", "type"), ("bin", "size")):
        print(f"\n=== Phase 0b: false positives by {label} ===")
        cats = sorted({r[key] for ds, _ in DATASETS for r in data[ds]["fp"]})
        print(f"{'dataset':12s} " + " ".join(f"{c:>10s}" for c in cats))
        for ds, _ in DATASETS:
            c = collections.Counter(r[key] for r in data[ds]["fp"])
            print(f"{ds:12s} " + " ".join(f"{c.get(x, 0):>10d}" for x in cats))
        print("  per-record precision within stratum (TP / (TP+FP)):")
        for ds, _ in DATASETS:
            cf = collections.Counter(r[key] for r in data[ds]["fp"])
            ct = collections.Counter(r[key] for r in data[ds]["tp"])
            cells = []
            for x in cats:
                n = ct.get(x, 0) + cf.get(x, 0)
                cells.append(f"{ct.get(x, 0) / n:10.3f}" if n else f"{'-':>10s}")
            print(f"  {ds:10s} " + " ".join(cells))

    # ---------------- Phase 1: inside the buckets ----------------
    print(f"\n=== Phase 1: false-positive buckets, from truvari's rejected candidate ===")
    order = ["no candidate", "placement", "consumed", "dissimilar"]
    print(f"{'dataset':12s} {'FP':>6s} " + " ".join(f"{b:>13s}" for b in order))
    for ds, _ in DATASETS:
        c = collections.Counter(bucket(r) for r in data[ds]["fp"])
        n = len(data[ds]["fp"])
        print(f"{ds:12s} {n:6d} " +
              " ".join(f"{c.get(b, 0):5d} ({c.get(b, 0) / n * 100 if n else 0:4.1f}%)"
                       for b in order))

    print("\n=== Phase 1b: inside `dissimilar` -- near-threshold, or genuinely not there? ===")
    print("PctSeqSimilarity of the rejected candidate, for FPs in the dissimilar bucket")
    print(f"{'dataset':12s} {'n':>5s} {'median':>7s} " +
          " ".join(f"{lab:>12s}" for lab in ("seq>=0.5", "0.2-0.5", "<0.2", "no seq")))
    for ds, _ in DATASETS:
        dis = [r for r in data[ds]["fp"] if bucket(r) == "dissimilar"]
        seqs = [fnum(r["info"], "PctSeqSimilarity") for r in dis]
        have = [s for s in seqs if s is not None]
        hi = sum(1 for s in have if s >= 0.5)
        mid = sum(1 for s in have if 0.2 <= s < 0.5)
        lo = sum(1 for s in have if s < 0.2)
        none = sum(1 for s in seqs if s is None)
        m = med(have)
        print(f"{ds:12s} {len(dis):5d} {m if m is not None else float('nan'):7.3f} "
              f"{hi:12d} {mid:12d} {lo:12d} {none:12d}")

    print("\n=== Phase 1c: do the shipped quality signals already know? ===")
    print("medians over SV records, true vs false, per bucket")
    for ds, _ in DATASETS:
        print(f"\n  {ds}")
        print(f"    {'population':<16s} {'n':>5s} {'GQ':>7s} {'share':>7s} {'DR':>7s} {'BL':>9s}")
        groups = [("TP", data[ds]["tp"])]
        by = collections.defaultdict(list)
        for r in data[ds]["fp"]:
            by[bucket(r)].append(r)
        groups += [(f"FP {b}", by[b]) for b in order if by[b]]
        for label, rs in groups:
            print(f"    {label:<16s} {len(rs):5d} "
                  f"{med([fnum(r['fmt'], 'GQ') for r in rs]) or float('nan'):7.1f} "
                  f"{med([share(r['fmt']) for r in rs]) or float('nan'):7.3f} "
                  f"{med([fnum(r['fmt'], 'DR') for r in rs]) or float('nan'):7.3f} "
                  f"{med([fnum(r['fmt'], 'BL') for r in rs]) or float('nan'):9.1f}")

    dest = WORK / "sv-atlas" / f"hap32-precision-{args.arm}.json"
    dest.write_text(json.dumps(
        {ds: {"n_fp": len(data[ds]["fp"]), "n_tp": len(data[ds]["tp"]),
              "n_fn": len(data[ds]["fn"]),
              "buckets": dict(collections.Counter(bucket(r) for r in data[ds]["fp"])),
              "fp_by_type": dict(collections.Counter(r["svtype"] for r in data[ds]["fp"])),
              "fp_by_size": dict(collections.Counter(r["bin"] for r in data[ds]["fp"]))}
         for ds, _ in DATASETS}, indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
