#!/usr/bin/env python3
"""What long-range linkage do the split homozygotes actually carry?

The companion to homsplit_switch.py. That one asks whether splits bridge a junction where the
phasing is WRONG; this one asks whether they bridge anything at all, by comparing the
haplotype-carrying graph built from het sites alone against the same graph with split
homozygotes admitted, at each gate setting.

Connectivity is built from the reads directly: a read links the consecutive admitted carriers
along its own placement list, and a link counts once it has --min-link reads. That is exact for
a chain -- any connection a read makes between two admitted carriers passes through its own
consecutive admitted pairs -- and costs one pass over the placements instead of the ~10M set
intersections a windowed pairwise scan would need.
"""
import argparse, collections, math, pickle, statistics, sys


def gate(los, min_q, frac, floor=2):
    s0 = s1 = 0
    for v in los:
        if v == 0.0 or math.isnan(v) or abs(v) < min_q:
            continue
        if v > 0.0:
            s0 += 1
        else:
            s1 += 1
    need = max(floor, math.ceil(frac * (s0 + s1)))
    return s0 >= need and s1 >= need


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def components(idx, posn, admitted, min_link):
    order = sorted(admitted, key=lambda s: (posn[s], s))
    rank = {s: i for i, s in enumerate(order)}
    pair = collections.Counter()
    for r, sites in idx["read_sites"].items():
        mine = sorted((rank[s] for s in sites if s in rank))
        for a, b in zip(mine, mine[1:]):
            if a != b:
                pair[(a, b)] += 1
    dsu = DSU(len(order))
    kept = 0
    for (a, b), n in pair.items():
        if n >= min_link:
            dsu.union(a, b)
            kept += 1
    spans = collections.defaultdict(lambda: [10**12, -1, 0])
    for s in order:
        root = dsu.find(rank[s])
        e = spans[root]
        p = posn[s]
        e[0] = min(e[0], p)
        e[1] = max(e[1], p)
        e[2] += 1
    sizes = sorted(((e[1] - e[0]) for e in spans.values()), reverse=True)
    total = sum(sizes)
    acc, n50 = 0, 0
    for x in sizes:
        acc += x
        if acc >= total / 2:
            n50 = x
            break
    return len(order), kept, len(spans), n50, (sizes[0] if sizes else 0), total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--lambda-table", required=True)
    ap.add_argument("--positions", required=True, help="pickle of snarl -> position")
    ap.add_argument("--min-link", type=int, default=2)
    ap.add_argument("--settings", default="none,0.5:0.0,2:0.0,4:0.0,8:0.0,0.5:0.33,4:0.33,8:0.33")
    a = ap.parse_args()

    idx = pickle.load(open(a.index, "rb"))
    posn = pickle.load(open(a.positions, "rb"))
    los = {}
    for line in open(a.lambda_table):
        if line[0] == "#":
            continue
        n, v = line.rstrip("\n").split("\t")
        los[n] = float(v)
    names = idx["names"]
    kind = idx["kind"]
    hets = {s for s, k in kind.items() if k == "het" and s in posn}
    homsplits = {s for s, k in kind.items() if k == "homsplit" and s in posn}
    site_los = {s: [los.get(names[r], 0.0) for r in idx["reads_of"].get(s, ())] for s in homsplits}

    print("setting\tcarriers\tlinks\tcomponents\tN50 (bp)\tlargest (bp)\tspanned (bp)")
    for tok in a.settings.split(","):
        if tok == "none":
            admitted, label = hets, "hets only (no hom split)"
        else:
            q, f = tok.split(":")
            live = {s for s in homsplits if gate(site_los[s], float(q), float(f))}
            admitted = hets | live
            label = f"q={q} frac={f} ({len(live)} split)"
        n, links, comp, n50, big, total = components(idx, posn, admitted, a.min_link)
        print(f"{label}\t{n}\t{links}\t{comp}\t{n50:,}\t{big:,}\t{total:,}")


if __name__ == "__main__":
    main()
