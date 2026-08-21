#!/usr/bin/env python3
"""Check which haplotype a nested haploid call was placed on, against a phased truth.

A nested site is called at ploidy 1 because exactly one parent allele crosses the child chain, so
which haplotype its allele sits on is *determined* by the parent's phase rather than estimated. That
makes it checkable, and it is not checkable any other way: whatshap refuses a VCF of mixed ploidy, so
every phasing number published for `--nested` so far was computed on the diploid records alone and is
blind to exactly the sites this measures.

The strand is not in the VCF -- a haploid `GT` of `1` with a `PS` says which block the site is in and
nothing about which strand -- so it is recovered from the mosaic. A nested site is a wildcard on the
strand it is *not* on, and a wildcard breaks a segment, so the mosaic isolates it: the site falls
inside a `*` segment on one strand and a named-haplotype segment on the other.

Then the orientation against the truth comes from the nearest matched diploid het site, read as a
local frame. Not one orientation per block: at 2.4% switch error per adjacent het pair a chromosome-length
block has no single orientation -- blockwise hamming is 49%, a coin flip -- and a first version of this
script took a block-wide majority and measured that noise. Only sites where the call and the truth agree
on POS, REF and ALT are counted, and only where the truth is heterozygous, since a homozygous truth site
carries the allele on both haplotypes and cannot tell a right strand from a wrong one.

**Read the controls, not the percentage.** This script cannot compare two strand conventions, and
finding that out is the main thing it has produced. On chr20 it scores the traversal-order slot the
caller records at 73.9% and the derivation from the parent's phased allele pair at 47.1% -- but "put
every nested site on strand 1" scores **74.5%** on the same sites. The recorded slot is 1 for about 82%
of nested sites, so its apparent accuracy is that skew and not information, and a balanced convention
loses to a constant against a one-sided truth subset. Both controls are printed for exactly this reason.

What can compare them is relative phase, which a constant cannot game: run whatshap over the call set
including the nested records (they are readable now that the caller writes `a|.` and `.|a`, and
`phasing_benchmark.py --half-missing ref` will feed them in). That says the two conventions are
indistinguishable -- 1,655 switches against 1,661 on ~58,900 pairs -- while the nested sites themselves
switch at about 21% against a 2.77% baseline. See the Stage 7 notes in docs/nested-calling-design.md.
"""

from __future__ import annotations

import argparse
import bisect
import statistics
import subprocess
import sys
from collections import defaultdict


def read_records(path, contig):
    """(pos, ref, alt, gt) for every biallelic record on one contig."""
    out = []
    cmd = ["bcftools", "view", "-H", "-r", contig, path]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stderr[-2000:], file=sys.stderr)
        sys.exit(f"bcftools failed on {path}")
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10 or "," in f[4]:
            continue
        gt = f[9].split(":")[0]
        out.append((int(f[1]), f[3], f[4], gt))
    return out


def strands_from_vcf(calls_records):
    """position -> strand, read straight off a half-missing phased genotype.

    `vg call --nested --phased` writes a nested haploid record as `a|.` or `.|a`, so the strand is in
    the record and needs no recovery. This supersedes the mosaic route below, which loses a quarter of
    the sites: wildcard segments merge across neighbouring sites, so a position can fall inside a
    wildcard run on both strands or on neither.
    """
    out = {}
    for pos, ref, alt, gt in calls_records:
        if gt.endswith("|."):
            out[pos] = 0
        elif gt.startswith(".|"):
            out[pos] = 1
    return out


def wildcard_intervals(mosaic, contig):
    """Per strand, the [start, end] ranges the panel names no haplotype over.

    Requires mosaic version 3 or later. Version 2 wrote * for two different things -- a strand the
    panel cannot name a haplotype for, and a strand carrying no sequence at all -- so on a v2 file
    this function silently returned the union of the two and over-reported.
    """
    spans = {0: [], 1: []}
    seen_version = None
    with open(mosaic) as fh:
        for line in fh:
            if line.startswith("#mosaic-version"):
                seen_version = line.rstrip("\n").split("\t")[1]
                if seen_version != "3":
                    raise SystemExit(
                        f"{mosaic}: mosaic-version {seen_version}, need 3 -- in version 2 the "
                        "haplotype column conflated 'panel cannot explain' with 'no sequence here', "
                        "so this measurement would over-report"
                    )
            if not line.startswith("H\t"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or f[1] != contig or f[7] != "*":
                continue
            spans[int(f[2])].append((int(f[3]), int(f[4])))
    for s in spans:
        spans[s].sort()
    return spans


def covered(spans, pos):
    """Whether pos falls in one of the sorted, non-overlapping ranges. Linear scan is fine: there
    are a few hundred of them and the callers walk positions in order anyway."""
    for start, end in spans:
        if start <= pos <= end:
            return True
        if start > pos:
            return False
    return False


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls", required=True)
    p.add_argument("--mosaic", required=True)
    p.add_argument("--truth", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--other-mosaic",
                   help="a second run's mosaic; report how many nested sites the two disagree on, "
                        "which is how the strand recovery itself gets checked -- vg reports how many "
                        "strands it moved, and a sound recovery has to find that many")
    args = p.parse_args()

    truth = {}
    for pos, ref, alt, gt in read_records(args.truth, args.contig):
        if "|" in gt:
            a, b = gt.split("|")
            if a.isdigit() and b.isdigit():
                truth[(pos, ref, alt)] = (int(a), int(b))

    calls = read_records(args.calls, args.contig)
    spans = wildcard_intervals(args.mosaic, args.contig)
    vcf_strand = strands_from_vcf(calls)
    if vcf_strand:
        print(f"strand read from the VCF for {len(vcf_strand)} nested records (half-missing "
              f"genotypes); the mosaic is a fallback only where the VCF has none")

    # Anchors: diploid het calls whose truth record matches, and whether our first strand is the
    # truth's first haplotype *there*.
    #
    # Local, not one orientation per block, because a block does not have one. chr20's phase carries
    # a 2.77% switch rate over 58,807 adjacent het pairs, so within a 248 Mb block the orientation
    # re-randomises every forty sites or so and the blockwise hamming rate is 49.3% -- a coin flip.
    # A block-wide majority is therefore noise, and an earlier version of this script measured
    # exactly that: 41.4% on 157 sites, indistinguishable from chance in either direction.
    anchors = []
    same = flip = 0
    for pos, ref, alt, gt in calls:
        if "|" not in gt:
            continue
        a, b = gt.split("|")
        if not (a.isdigit() and b.isdigit()) or a == b:
            continue
        t = truth.get((pos, ref, alt))
        if t is None or t[0] == t[1]:
            continue
        if (int(a), int(b)) == t:
            anchors.append((pos, 0))
            same += 1
        elif (int(b), int(a)) == t:
            anchors.append((pos, 1))
            flip += 1
    anchors.sort()
    print(f"anchors: {len(anchors)} matched diploid het sites "
          f"({same} unflipped, {flip} flipped against the truth) -- used as local frames, since "
          f"the block has no single orientation")

    anchor_pos = [a[0] for a in anchors]

    def local_frame(pos):
        """The nearest anchor's frame, and how far away it was."""
        if not anchors:
            return None, None
        k = bisect.bisect_left(anchor_pos, pos)
        best = None
        for cand in (k - 1, k):
            if 0 <= cand < len(anchors):
                d = abs(anchors[cand][0] - pos)
                if best is None or d < best[1]:
                    best = (anchors[cand][1], d)
        return best

    if args.other_mosaic:
        other = wildcard_intervals(args.other_mosaic, args.contig)

        def strand_of(sp, pos):
            w0, w1 = covered(sp[0], pos), covered(sp[1], pos)
            if w0 == w1:
                return None      # both or neither: no single strand to read off
            return 1 if w0 else 0

        agree = differ = unusable = 0
        for pos, ref, alt, gt in calls:
            if "|" in gt or "/" in gt or not gt.isdigit():
                continue
            a, b = strand_of(spans, pos), strand_of(other, pos)
            if a is None or b is None:
                unusable += 1
            elif a == b:
                agree += 1
            else:
                differ += 1
        print(f"against {args.other_mosaic}: {differ} nested sites on a different strand, "
              f"{agree} the same, {unusable} not comparable")

        # Paired against the truth on the sites where both runs have a strand and the truth can
        # tell them apart. Unpaired totals compare different subsets, since the two runs do not
        # make the same sites recoverable.
        both = other_right = this_right = 0
        for pos, ref, alt, gt in calls:
            if "|" in gt or "/" in gt or not gt.isdigit():
                continue
            a, b = strand_of(spans, pos), strand_of(other, pos)
            if a is None or b is None:
                continue
            t = truth.get((pos, ref, alt))
            if t is None or t[0] == t[1]:
                continue
            frame, _ = local_frame(pos)
            if frame is None:
                continue
            both += 1
            this_right += t[frame if a == 0 else 1 - frame] == int(gt)
            other_right += t[frame if b == 0 else 1 - frame] == int(gt)
        if both:
            print(f"  paired on {both} sites decisive for both: this mosaic right {this_right} "
                  f"({100.0 * this_right / both:.1f}%), the other {other_right} "
                  f"({100.0 * other_right / both:.1f}%)")

    stats = defaultdict(int)
    distances = []
    for pos, ref, alt, gt in calls:
        if pos in vcf_strand:
            # The record says which strand, so no recovery and no losses.
            strand = vcf_strand[pos]
            gt = gt.replace("|.", "").replace(".|", "")
            stats["haploid_records"] += 1
        elif "|" in gt or "/" in gt or not gt.isdigit():
            continue          # diploid, or no call
        else:
            stats["haploid_records"] += 1
            wild0, wild1 = covered(spans[0], pos), covered(spans[1], pos)
            if wild0 and wild1:
                stats["no_strand"] += 1
                continue
            if not wild0 and not wild1:
                stats["strand_not_recoverable"] += 1
                continue
            strand = 1 if wild0 else 0
        t = truth.get((pos, ref, alt))
        if t is None:
            stats["no_matching_truth_record"] += 1
            continue
        if t[0] == t[1]:
            stats["truth_homozygous"] += 1
            continue
        frame, distance = local_frame(pos)
        if frame is None:
            stats["no_local_frame"] += 1
            continue
        # Our strand 0 is the truth's haplotype `frame` in this neighbourhood.
        claimed = frame if strand == 0 else 1 - frame
        if t[claimed] == int(gt):
            stats["correct"] += 1
        else:
            stats["wrong"] += 1
        # Controls. The recorded slot is 1 for about 82% of nested sites on chr20, so "it agrees with
        # the truth 74% of the time" is only a result if a constant does worse. Scored on exactly the
        # same sites and frames.
        for name, const in (("always_strand_0", 0), ("always_strand_1", 1)):
            c = const if frame == 0 else 1 - const
            if t[c] == int(gt):
                stats[name] += 1
        distances.append(distance)

    decisive = stats["correct"] + stats["wrong"]
    for k in ("haploid_records", "no_strand", "strand_not_recoverable",
              "no_matching_truth_record", "truth_homozygous", "correct", "wrong"):
        print(f"  {k:28s} {stats[k]}")
    if decisive:
        print(f"\nstrand correct on {stats['correct']} of {decisive} decisive sites "
              f"({100.0 * stats['correct'] / decisive:.1f}%)")
        for name in ("always_strand_0", "always_strand_1"):
            print(f"  control {name:16s} {stats[name]:4d} "
                  f"({100.0 * stats[name] / decisive:.1f}%)")
        if distances:
            print(f"  nearest anchor: median {int(statistics.median(distances))} bp, "
                  f"max {max(distances)} bp")
    else:
        print("\nno decisive sites")


if __name__ == "__main__":
    main()
