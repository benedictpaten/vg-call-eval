#!/usr/bin/env python3
"""Is the hom-split gate selective? Simulated on het sites, where the answer is known.

A homozygous site is split on an inference it cannot check. A heterozygous site makes the same
inference and also knows the answer -- its own alleles partition the reads -- so running the hom
gate on het sites asks exactly what the gate does blind, against held-out ground truth, at 77k
sites instead of 32 switch events.

Input is the per (het site, read) dump: the leave-one-out strand log-odds and whether its sign
agrees with the slot the site's OWN alleles put the read in. Requires --no-anchors-phase-hets,
or the slot partly follows the strand and the check measures itself.
"""
import argparse, collections, math, statistics, sys


def gate(los, min_q, frac, floor=2):
    s0 = s1 = 0
    for v in los:
        if v == 0.0 or abs(v) < min_q:
            continue
        if v > 0.0:
            s0 += 1
        else:
            s1 += 1
    need = max(floor, math.ceil(frac * (s0 + s1)))
    return s0 >= need and s1 >= need


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hetcheck", required=True)
    ap.add_argument("--settings", default="0.5:0.0,1:0.0,2:0.0,4:0.0,8:0.0,"
                                          "0.5:0.10,0.5:0.20,0.5:0.33,0.5:0.40,"
                                          "2:0.20,2:0.33,4:0.33,8:0.33")
    a = ap.parse_args()

    los = collections.defaultdict(list)
    agr = collections.defaultdict(list)
    for line in open(a.hetcheck):
        if line[0] == "#":
            continue
        s, v, g = line.rstrip("\n").split("\t")
        los[s].append(float(v))
        agr[s].append(g == "1")
    sites = list(los)
    print(f"{len(sites)} het sites, {sum(len(v) for v in los.values())} read observations",
          file=sys.stderr)

    tot_r = sum(len(v) for v in agr.values())
    tot_ok = sum(sum(v) for v in agr.values())
    print(f"baseline agreement over all reads at all het sites: "
          f"{tot_ok}/{tot_r} = {100*tot_ok/tot_r:.4f}%\n")

    print("min_q\tfrac\tsites kept\tkept %\t| reads kept\tagreement kept\t| reads dropped\t"
          "agreement dropped\tdisagreements removed")
    base_bad = tot_r - tot_ok
    for tok in a.settings.split(","):
        q, f = tok.split(":")
        q, f = float(q), float(f)
        kr = ko = dr = do = 0
        nk = 0
        for s in sites:
            if gate(los[s], q, f):
                nk += 1
                kr += len(agr[s]); ko += sum(agr[s])
            else:
                dr += len(agr[s]); do += sum(agr[s])
        rem = (dr - do)
        print(f"{q}\t{f}\t{nk}\t{100*nk/len(sites):.1f}%\t| {kr}\t"
              f"{100*ko/kr if kr else 0:.4f}%\t| {dr}\t"
              f"{100*do/dr if dr else 0:.4f}%\t{rem} ({100*rem/base_bad:.1f}% of {base_bad})")

    # Where does the disagreement actually live? Per-site agreement, sites with >= 10 reads.
    print("\nper-site agreement distribution (sites with >= 10 reads):")
    per = sorted((sum(agr[s]) / len(agr[s]), s) for s in sites if len(agr[s]) >= 10)
    n = len(per)
    for label, lo_i, hi_i in [("worst 1%", 0, n // 100), ("worst 5%", 0, n // 20),
                              ("worst 10%", 0, n // 10), ("all", 0, n)]:
        sl = per[lo_i:hi_i]
        if not sl:
            continue
        reads = sum(len(agr[s]) for _, s in sl)
        bad = sum(len(agr[s]) - sum(agr[s]) for _, s in sl)
        print(f"  {label}: {len(sl)} sites, mean agreement {100*statistics.mean(x for x,_ in sl):.2f}%, "
              f"{bad} of {reads} reads disagree ({100*bad/base_bad:.1f}% of all disagreements)")


if __name__ == "__main__":
    main()
