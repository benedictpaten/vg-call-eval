#!/usr/bin/env python3
"""Refusing a weak relink is the only thing that actually shortens the chain.

Suppressing hom splits around a probable switch does not: the flanking het sites are linked
directly by reads at every switch junction measured (28/28 chr20, 26/26 chr6), so the chain
stays connected through them and only loses some homozygous labels. What shortens a chain that
contains a switch is a PHASE-SET BOUNDARY, and stage 2 is where to put one -- it is the only
place the phasing admits it had no cascade evidence, and 79% of chr20's switches are there.

Today stage 2 always relinks (unless the vote is exactly zero, where the panel's frame stands).
This prices the alternative: refuse the weakest relinks, and report the resulting block N50
against the switches left inside blocks.
"""
import argparse, bisect, collections, pickle, sys


def load_relinks(path):
    rows, best = [], -1
    for line in open(path):
        if line[0] == "#" or not line.startswith("relink"):
            continue
        c = line.rstrip("\n").split("\t")
        best = max(best, int(c[1]))
        rows.append((int(c[1]), int(c[2]), int(c[3]), float(c[4]), int(c[5]), int(c[6]),
                     int(c[7]), float(c[8]), float(c[9])))
    return [r for r in rows if r[0] == best]


def switches(path):
    rows = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(path)
                  if not l.startswith(("track", "#")))
    runs, cur = [], []
    for r in rows:
        if cur and r[0] == cur[-1][1]:
            cur.append(r)
        else:
            if cur:
                runs.append(cur)
            cur = [r]
    if cur:
        runs.append(cur)
    return ([(r[0][0], r[0][1]) for r in runs if len(r) == 1],
            [(r[0][0], r[-1][1]) for r in runs if len(r) == 2],
            [(r[0][0], r[-1][1]) for r in runs if len(r) > 2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--switch-bed", required=True)
    ap.add_argument("--contig-span", type=int, required=True)
    a = ap.parse_args()

    rel = load_relinks(a.junctions)
    iso, paired, longer = switches(a.switch_bed)
    targets = iso + longer
    print(f"{len(rel)} stage-2 relinks; {len(targets)} true switches, {len(paired)} flips\n",
          file=sys.stderr)

    rows = []
    for (_, pa, pb, tot, n, npos, nneg, mp, mn) in rel:
        lo, hi = min(pa, pb), max(pa, pb)
        rows.append({"lo": lo, "hi": hi, "mid": (lo + hi) // 2, "tot": abs(tot), "n": n,
                     "minf": min(npos, nneg) / n if n else 0.5,
                     "gap": hi - lo})

    def evaluate(refused, optimistic=True):
        """optimistic: a switch is separated if its interval OVERLAPS a refused gap at all --
        an upper bound, since a site inside the gap is hung by stage 3 and may land either side.
        pessimistic: only if the gap's midpoint falls inside the switch interval."""
        cuts = sorted(r["mid"] for r in refused)
        gaps = sorted((r["lo"], r["hi"]) for r in refused)
        gstart = [g[0] for g in gaps]
        inside = 0
        for t in targets:
            sep = False
            if optimistic:
                i = bisect.bisect_right(gstart, t[1]) - 1
                for q in (i, i - 1):
                    if 0 <= q < len(gaps) and gaps[q][1] >= t[0]:
                        sep = True
                        break
            else:
                i = bisect.bisect_left(cuts, t[0])
                sep = i < len(cuts) and cuts[i] <= t[1]
            if not sep:
                inside += 1
        edges = [0] + cuts + [a.contig_span]
        spans = sorted((edges[i + 1] - edges[i] for i in range(len(edges) - 1)), reverse=True)
        tot = sum(spans)
        acc = 0
        n50 = 0
        for x in spans:
            acc += x
            if acc >= tot / 2:
                n50 = x
                break
        return len(spans), n50, inside

    print("rule\trefused\tblocks\tblock N50\tswitches left inside (optimistic / strict)")
    b, n50, ins = evaluate([])
    print(f"today (always relink)\t0\t{b}\t{n50/1e6:.2f} Mb\t{ins}/{len(targets)}")
    # the bound: an oracle that cuts at exactly the true switches
    ocuts = [{"lo": t[0], "hi": t[1], "mid": (t[0]+t[1])//2} for t in targets]
    b, n50, ins = evaluate(ocuts)
    print(f"ORACLE (cut at every true switch)\t{len(targets)}\t{b}\t{n50/1e6:.2f} Mb\t"
          f"{ins}/{len(targets)}")
    for name, key, desc in [("weakest |vote|", lambda r: r["tot"], False),
                            ("most split vote", lambda r: r["minf"], True),
                            ("fewest votes", lambda r: r["n"], False),
                            ("widest gap", lambda r: r["gap"], True)]:
        order = sorted(rows, key=key, reverse=desc)
        for f in (0.002, 0.01, 0.05, 0.25, 1.0):
            k = max(1, int(f * len(rows)))
            b, n50, ins = evaluate(order[:k])
            _, _, strict = evaluate(order[:k], optimistic=False)
            print(f"{name}\t{k} ({f*100:g}%)\t{b}\t{n50/1e6:.2f} Mb\t"
                  f"{ins} / {strict} of {len(targets)}")


if __name__ == "__main__":
    main()
