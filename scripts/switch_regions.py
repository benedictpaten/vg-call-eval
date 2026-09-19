#!/usr/bin/env python3
"""Switches cluster. Can the regions be found without truth?

On chr20, 82% of true switches fall in 5.4% of the contig -- pericentromeric (26.5 Mb) and
subtelomeric (65.3-65.5 Mb). Picking 28 junctions out of 120,372 is hopeless at a 4x lift;
picking 3.6 Mb out of 66 Mb is a different and much easier problem, and it is the one that
matches "avoid splits AROUND them".

Every statistic here is available at run time. None uses truth.
"""
import argparse, bisect, collections, pickle, statistics, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--switch-bed", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--positions", required=True)
    ap.add_argument("--lambda-table", required=True)
    ap.add_argument("--window", type=int, default=100000)
    a = ap.parse_args()

    sys.path.insert(0, "scripts")
    from relink_pareto import switches
    iso, paired, longer = switches(a.switch_bed)
    targets = iso + longer

    W = a.window
    # --- per-junction statistics from the phasing ---
    jd = collections.defaultdict(list)     # window -> [|d|]
    jn = collections.defaultdict(list)     # window -> [n shared reads]
    jmin = collections.defaultdict(list)   # window -> [minority pair fraction]
    jbrk = collections.defaultdict(lambda: [0, 0])
    jcoh = collections.defaultdict(list)
    lastp = 0
    rows = []
    for line in open(a.junctions):
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        rows.append(c)
        lastp = max(lastp, int(c[1]) if c[0] == "link" else 0)
    for c in rows:
        kind, pas, pa = c[0], int(c[1]), int(c[2])
        w = pa // W
        if kind == "link" and pas == lastp:
            d, n, nc, nt = abs(float(c[4])), int(c[5]), int(c[6]), int(c[7])
            jd[w].append(d); jn[w].append(n)
            jmin[w].append(min(nc, nt) / n if n else 0.5)
            jbrk[w][0] += 1
            jbrk[w][1] += 1 if c[10] == "break" else 0
        elif kind == "coh":
            jcoh[w].append(float(c[4]))

    # --- per-site statistics from the anchors ---
    idx = pickle.load(open(a.index, "rb"))
    posn = pickle.load(open(a.positions, "rb"))
    los = {}
    for line in open(a.lambda_table):
        if line[0] == "#":
            continue
        n_, v = line.rstrip("\n").split("\t")
        los[n_] = abs(float(v))
    names = idx["names"]
    sdepth = collections.defaultdict(list)
    slam = collections.defaultdict(list)
    splits_at = collections.defaultdict(int)
    allsplits = []
    for s, k in idx["kind"].items():
        p = posn.get(s)
        if p is None:
            continue
        w = p // W
        rd = idx["reads_of"].get(s, ())
        sdepth[w].append(len(rd))
        if k == "homsplit":
            splits_at[w] += 1
            allsplits.append(p)
            for r in rd:
                slam[w].append(los.get(names[r], 0.0))

    wins = sorted(set(jd) | set(sdepth))
    nsw = collections.Counter()
    for t in targets:
        nsw[t[0] // W] += 1

    def agg(d, w, fn=statistics.mean, default=None):
        v = d.get(w)
        return fn(v) if v else default

    feats = {}
    for w in wins:
        feats[w] = {
            "break rate": (jbrk[w][1] / jbrk[w][0]) if jbrk[w][0] else 0.0,
            "median |d|": agg(jd, w, statistics.median, 1e9),
            "median shared reads": agg(jn, w, statistics.median, 1e9),
            "mean minority frac": agg(jmin, w, statistics.mean, 0.0),
            "mean coherence": agg(jcoh, w, statistics.mean, 1.0),
            "median site depth": agg(sdepth, w, statistics.median, 1e9),
            "mean |lambda| at splits": agg(slam, w, statistics.mean, 1e9),
            "site density": (len(sdepth.get(w, ())) ),
        }
    total_splits = len(allsplits)
    total_sw = sum(nsw.values())
    print(f"{len(wins)} windows of {W//1000} kb; {total_sw} switches; {total_splits} splits\n",
          file=sys.stderr)

    print("signal\tflagged windows\tswitches covered\tsplits suppressed\tlift")
    order_specs = [("break rate", True), ("median |d|", False), ("median shared reads", False),
                   ("mean minority frac", True), ("mean coherence", False),
                   ("median site depth", False), ("mean |lambda| at splits", False),
                   ("site density", False)]
    for name, desc in order_specs:
        order = sorted(wins, key=lambda w: feats[w][name], reverse=desc)
        for frac in (0.02, 0.05, 0.10, 0.20):
            k = max(1, int(frac * len(wins)))
            sel = set(order[:k])
            cov = sum(nsw[w] for w in sel)
            sp = sum(splits_at[w] for w in sel)
            c = 100 * cov / total_sw
            s = 100 * sp / total_splits
            print(f"{name}\t{frac*100:g}% ({k})\t{cov}/{total_sw} ({c:.0f}%)\t{s:.1f}%\t"
                  f"{c/s if s else 0:.1f}x")


if __name__ == "__main__":
    main()
