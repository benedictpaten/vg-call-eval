#!/usr/bin/env python3
"""Subsample a contig's GAF to several target coverages in one pass.

Two properties this gets deliberately, both of which matter for reading the titration curves:

**The levels are nested.** A read is kept at level c if `hash(name)/2^32 < c/source_coverage`, so
the 5x set is a subset of the 10x set is a subset of the 20x set. The alternative -- an independent
draw per level -- would put sampling noise into every pairwise comparison, and the whole point of
the exercise is to attribute an F1 difference to coverage rather than to which reads got drawn.
Nesting makes the comparison paired.

**Mates stay together.** The key is the read name, and paired mates share a name in this GAF (a
fact `extract_reads_from_db.py` documents having been bitten by: 20,000 records carry 10,000
distinct names). Hashing the name keeps or drops both mates as a unit. The read-likelihood model
scores reads independently so this changes little for the caller, but it keeps the subsets honest
as read sets and keeps `vg pack` counts sane if an arm ever needs them.

Coverage is measured from the data rather than assumed: the source's coverage is the summed query
alignment span (GAF columns 3 and 4, 0-based half-open) divided by the contig length. Assuming
151 bp per read would have been close here but wrong in general, and the ratio is what sets every
sampling fraction.
"""

from __future__ import annotations

import argparse
import zlib
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gaf", required=True, type=Path, help="source GAF for one contig")
    p.add_argument("--contig-length", required=True, type=int)
    p.add_argument("--out-prefix", required=True, type=Path,
                   help="writes <prefix>.<cov>x.gaf per level")
    p.add_argument("--levels", type=float, nargs="+", default=[5, 10, 15, 20],
                   help="target coverages; levels at or above source are skipped, since the "
                        "source itself is the full-coverage arm and should be used directly")
    p.add_argument("--measure-only", action="store_true",
                   help="report the source coverage and exit, writing nothing")
    args = p.parse_args()

    # Pass 1: measure. Cheap relative to writing, and the fractions cannot be chosen without it.
    total_span = 0
    n_reads = 0
    with open(args.gaf) as fh:
        for line in fh:
            f = line.split("\t", 5)
            if len(f) < 5:
                continue
            try:
                total_span += int(f[3]) - int(f[2])
            except ValueError:
                continue
            n_reads += 1
    source_cov = total_span / args.contig_length
    print(f"source: {n_reads:,} reads, {total_span/1e9:.3f} Gbp aligned span, "
          f"{args.contig_length/1e6:.1f} Mb contig -> {source_cov:.2f}x")
    if args.measure_only:
        return

    levels = sorted(c for c in args.levels if c < source_cov)
    skipped = [c for c in args.levels if c >= source_cov]
    if skipped:
        print(f"skipping {skipped} at or above source coverage {source_cov:.2f}x -- "
              f"use the source GAF directly as the full-coverage arm")
    if not levels:
        print("nothing to write")
        return

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fracs = [(c, c / source_cov) for c in levels]
    handles = []
    for c, fr in fracs:
        path = Path(f"{args.out_prefix}.{c:g}x.gaf")
        handles.append((fr, open(path, "w"), path, [0]))
        print(f"  {c:g}x -> fraction {fr:.4f} -> {path}")

    with open(args.gaf) as fh:
        for line in fh:
            name = line.split("\t", 1)[0]
            h = zlib.crc32(name.encode()) / 4294967296.0
            for fr, out, _, count in handles:
                if h < fr:
                    out.write(line)
                    count[0] += 1

    print()
    for (c, _), (_, out, path, count) in zip(fracs, handles):
        out.close()
        print(f"  {c:g}x: {count[0]:,} reads -> {path}")


if __name__ == "__main__":
    main()
