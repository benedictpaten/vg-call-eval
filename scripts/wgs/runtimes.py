#!/usr/bin/env python3
"""Per-contig wall, CPU and peak RSS for whole-genome runs, from their `/usr/bin/time -l` blocks.

Wall clock is NOT comparable between two runs at different concurrency -- the short-read run packs
contigs under schedule_wgs.py's memory budget (up to 3 at -t 5), the ONT run is `xargs -P 2` at
-t 5 -- so the comparison this prints is CPU seconds, and only over contigs BOTH runs have. An
earlier version of this divided one run's 6 contigs by the other's 24 and reported 1.13x where the
like-for-like figure is 2.96x.
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path

CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

def parse(log: Path):
    try:
        t = log.read_text(errors="ignore")
    except FileNotFoundError:
        return None
    m = re.search(r"([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys", t)
    if not m:
        return None
    r = re.search(r"(\d+)\s+maximum resident set size", t)
    return (float(m.group(1)), float(m.group(2)) + float(m.group(3)),
            int(r.group(1)) / 2**30 if r else 0.0)

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="work/wgs-0912", help="first run's work dir")
    ap.add_argument("--b", default="work/ont-wgs", help="second run's work dir")
    ap.add_argument("--a-label", default="short")
    ap.add_argument("--b-label", default="ONT")
    args = ap.parse_args()

    rows = [(c, parse(Path(args.a) / c / f"{c}.log"), parse(Path(args.b) / c / f"{c}.log"))
            for c in CONTIGS]
    print(f"{'contig':7} | {args.a_label+' wall':>11} {'cpu':>8} {'RSS':>6} | "
          f"{args.b_label+' wall':>11} {'cpu':>8} {'RSS':>6} | {'cpu x':>6}")
    print("-" * 76)
    for c, a, b in rows:
        fa = f"{a[0]:>11.0f} {a[1]:>8.0f} {a[2]:>6.1f}" if a else f"{'-':>11} {'-':>8} {'-':>6}"
        fb = f"{b[0]:>11.0f} {b[1]:>8.0f} {b[2]:>6.1f}" if b else f"{'-':>11} {'-':>8} {'-':>6}"
        rt = f"{b[1]/a[1]:>6.2f}" if a and b and a[1] else f"{'-':>6}"
        print(f"{c:7} | {fa} | {fb} | {rt}")

    both = [(c, a, b) for c, a, b in rows if a and b]
    if not both:
        print("\n  no contig has both runs yet")
        return 0
    ca = sum(a[1] for _, a, _ in both); cb = sum(b[1] for _, _, b in both)
    print("-" * 76)
    print(f"\n  over the {len(both)} contigs both runs have:")
    print(f"    CPU  {args.a_label} {ca/3600:.2f} h   {args.b_label} {cb/3600:.2f} h   "
          f"ratio {cb/ca:.2f}x")
    print(f"    peak RSS  {args.a_label} {max(a[2] for _, a, _ in both):.1f} GB   "
          f"{args.b_label} {max(b[2] for _, _, b in both):.1f} GB")
    print(f"    {args.a_label}: {sum(1 for _, a, _ in rows if a)}/24 contigs, "
          f"{args.b_label}: {sum(1 for _, _, b in rows if b)}/24")
    return 0

if __name__ == "__main__":
    sys.exit(main())
