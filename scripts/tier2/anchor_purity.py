#!/usr/bin/env python3
"""Purity and calibration of `vg call --anchors-out`, against reads of known haplotype origin.

Simulate separately from two haplotype paths, prefix the read names, merge, map, call. Then this
answers the question the whole boundary-pin design rests on: an anchor's reads are assigned by a
model rather than by exact match, so how often is that assignment wrong, and does the per-read score
say so?

  purity       fraction of anchor members whose true haplotype is the anchor's majority, over
               anchors at sites that HAVE two slots
  yield        anchors, reads per anchor, sites whose two slots resolve to DIFFERENT truth
               haplotypes (a site whose slots agree has not partitioned anything)
  calibration  observed error rate against the claimed per-read phred, in buckets. A score of 10
               should mean about 10% wrong. This is the one that matters: an uncalibrated score is
               worse than none, because a downstream filter would trust it.

**Single-slot anchors are counted but excluded from purity, and that is not a convenience.** A
homozygous or haploid site produces one anchor holding every read, from both haplotypes; its
"majority purity" is therefore about 0.5 by construction and measures nothing. Since `vg` emits those
by default, folding them in would halve the headline number while nothing had actually got worse.
"""

import argparse
import math
import re
import sys
from collections import Counter, defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--pattern", default=r"^(h[12])_",
                    help="regex whose first group is the read's true haplotype [^(h[12])_]")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--by-score", action="store_true", help="print the calibration table")
    args = ap.parse_args()

    truth_re = re.compile(args.pattern)

    # The interning table: R rows carry an id into it.
    id_name = {}
    with open(args.anchors) as handle:
        for line in handle:
            if line.startswith("#read\t"):
                f = line.rstrip("\n").split("\t")
                if len(f) == 3:
                    id_name[f[1]] = f[2]
            elif not line.startswith("#"):
                break
    if not id_name:
        sys.exit(f"{args.anchors} has no #read table; expected an anchors-version 2 file")

    anchors = []           # (node, snarl, slot, gqn, [(hap, score), ...])
    current = None
    members = []
    unknown = 0

    def flush():
        if current is not None and len(members) >= args.min_reads:
            anchors.append((*current, members[:]))

    with open(args.anchors) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "A":
                flush()
                members = []
                gqn = None if f[4] == "." else float(f[4])
                current = (int(f[1]), f[2], int(f[3]), gqn)
            elif f[0] == "R" and current is not None:
                m = truth_re.match(id_name.get(f[1], ""))
                if not m:
                    unknown += 1
                    continue
                members.append((m.group(1), float(f[4])))
    flush()

    if not anchors:
        sys.exit("no anchors with a recognisable truth haplotype; check --pattern")

    # Which sites actually offer a partition? A site with one slot holds every read in one anchor,
    # so there is nothing for its members to be pure or impure about.
    slot_count = defaultdict(set)
    for node, snarl, slot, gqn, reads in anchors:
        slot_count[(node, snarl)].add(slot)

    total_reads = 0
    correct_reads = 0
    per_anchor_purity = []
    buckets = defaultdict(lambda: [0, 0])     # score bucket -> [n, wrong]
    slots_by_site = defaultdict(dict)         # (node, snarl) -> slot -> majority hap
    single_anchors = 0
    single_reads = 0

    for node, snarl, slot, gqn, reads in anchors:
        if len(slot_count[(node, snarl)]) < 2:
            single_anchors += 1
            single_reads += len(reads)
            continue
        counts = Counter(hap for hap, _ in reads)
        majority, majority_n = counts.most_common(1)[0]
        per_anchor_purity.append(majority_n / len(reads))
        slots_by_site[(node, snarl)][slot] = majority
        for hap, score in reads:
            total_reads += 1
            wrong = hap != majority
            correct_reads += not wrong
            b = min(int(score), 30)
            buckets[b][0] += 1
            buckets[b][1] += wrong

    if not per_anchor_purity:
        sys.exit("no anchor sits at a site with two slots, so there is no partition to score")

    partitioned = 0
    collapsed = 0
    for site, slots in slots_by_site.items():
        if len(slots) < 2:
            continue
        if len(set(slots.values())) == 2:
            partitioned += 1
        else:
            collapsed += 1

    # Per-site partition accuracy: the two slots are assigned to the two truth haplotypes the way
    # that agrees best, then accuracy is the share of reads sitting in the slot matching their own
    # origin. Unlike per-anchor purity this is a statement about the PARTITION rather than about
    # within-anchor agreement, so a site whose two slots both fill up with h1 scores about 0.5 here
    # where purity would call both anchors pure. Chance is 0.5 for both.
    site_reads = defaultdict(lambda: defaultdict(Counter))
    for node, snarl, slot, gqn, reads in anchors:
        if len(slot_count[(node, snarl)]) < 2:
            continue
        for hap, score in reads:
            site_reads[(node, snarl)][slot][hap] += 1
    acc_right = 0
    acc_total = 0
    for site, slots in site_reads.items():
        keys = sorted(slots)
        if len(keys) != 2:
            continue
        a, b = slots[keys[0]], slots[keys[1]]
        direct = a["h1"] + b["h2"]
        swapped = a["h2"] + b["h1"]
        acc_right += max(direct, swapped)
        acc_total += sum(a.values()) + sum(b.values())

    scored_anchors = len(anchors) - single_anchors
    print(f"anchors                 {len(anchors):,}"
          + (f"   ({single_anchors:,} single-slot, excluded from purity)" if single_anchors else ""))
    print(f"read placements         {total_reads + single_reads:,}"
          + (f"   ({unknown:,} with no truth label)" if unknown else ""))
    print(f"reads per anchor        {(total_reads + single_reads) / len(anchors):.1f}")
    print()
    print(f"-- over the {scored_anchors:,} anchors at sites that offer a partition --")
    print(f"read placements scored  {total_reads:,}")
    print(f"purity, read-weighted   {correct_reads / total_reads:.4f}")
    print(f"purity, anchor-mean     {sum(per_anchor_purity) / len(per_anchor_purity):.4f}")
    if acc_total:
        print(f"partition accuracy      {acc_right / acc_total:.4f}   "
              f"(chance 0.5; slots assigned to haplotypes the way that agrees best)")
    both = partitioned + collapsed
    if both:
        print(f"sites with two slots    {both:,}")
        print(f"  slots name 2 haps     {partitioned:,} ({100 * partitioned / both:.1f}%)")
        print(f"  slots collapse to 1   {collapsed:,} ({100 * collapsed / both:.1f}%)")

    if args.by_score:
        print()
        print("score  n          wrong    observed   claimed")
        for b in sorted(buckets):
            n, wrong = buckets[b]
            observed = wrong / n
            obs_phred = 99.0 if observed == 0 else -10 * math.log10(observed)
            print(f"{b:>5}  {n:<10,} {wrong:<8,} {observed:.4f}     {obs_phred:.1f}")


if __name__ == "__main__":
    sys.exit(main())
