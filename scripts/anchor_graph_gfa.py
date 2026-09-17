#!/usr/bin/env python3
"""Condense an anchor file into a GFA small enough to open in Bandage.

THE GRAPH.  Every pin -- one `A` row, a (boundary node, snarl, slot) triple -- is one bidirected
node: a read passes through it, entering one side and leaving the other, and the row's `strand`
says which way round.  Not every site-slot owns two pins (some are dropped), so the pin, not the
site, is the unit that always exists.  Walking one read's pins in `offset` order gives its path;
consecutive pins give the links, weighted by how many reads carry them.

WHAT THE PICTURE SHOWS.  A phased site carries two slots, so it contributes two parallel strands,
one per haplotype.  A site that failed to split carries ONE slot that every read passes through, so
the two strands PINCH into a single node there.  Those pinches are the run breaks, and in Bandage
they are the visible feature: long parallel ladders knotted at each collapsed site.  Colour by
`Class` from the CSV to pick them out.

CONDENSATION.  Unitig compaction over the bidirected graph: a link is unique when the node side at
BOTH of its ends has degree one, and unique links chain nodes into paths that collapse to a single
segment.  A read that merely ENDS contributes no link, so termination never blocks a merge -- only
a real branch does.  --min-edge-reads prunes thin links first; without it a few chimeric reads
shatter the condensation.

ANNOTATION.  LN:i: is the reference span in bases of the chain a segment represents (needs --vcf;
falls back to the pin count).  DP:f: is mean read depth over its pins.  The --csv carries pin
count, span, depth, class and locus so Bandage can label and colour by any of them.
"""
from __future__ import annotations

import argparse, array, collections, sys

SEG_BITS, OFF_SHIFT = 21, 22          # packed = off<<22 | seg<<1 | strand
SEG_MASK = (1 << SEG_BITS) - 1


def parse_region(s):
    if not s:
        return None
    contig, _, rng = s.partition(":")
    if not rng:
        return (contig, 0, 1 << 62)
    lo, _, hi = rng.partition("-")
    return (contig, int(lo.replace(",", "")), int(hi.replace(",", "")))


def load_vcf(path):
    pos, gt = {}, {}
    if path:
        for line in open(path):
            if line[0] == "#":
                continue
            c = line.split("\t", 10)
            pos[c[2]] = (c[0], int(c[1]))
            gt[c[2]] = c[9].split(":", 1)[0]
    return pos, gt


def read_anchors(path, pos, region, max_reads=None):
    """One segment per pin. Returns segments, per-read packed placements, and per-segment depth."""
    seg_id, seg_key = {}, []
    slots_of_site = collections.defaultdict(set)
    per_read = collections.defaultdict(lambda: array.array("q"))
    depth = []
    cur, keep_cur = -1, False
    rc, rlo, rhi = region if region else (None, 0, 0)
    for line in open(path):
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] == "A":
            snarl, slot = c[2], int(c[3])
            slots_of_site[snarl].add(slot)
            if region:
                p = pos.get(snarl)
                keep_cur = bool(p and p[0] == rc and rlo <= p[1] <= rhi)
            else:
                keep_cur = True
            if not keep_cur:
                cur = -1
                continue
            key = (int(c[1]), snarl, slot)
            cur = seg_id.get(key)
            if cur is None:
                cur = seg_id[key] = len(seg_key)
                seg_key.append(key)
                depth.append(0)
                if len(seg_key) > SEG_MASK:
                    sys.exit(f"more than {SEG_MASK} pins; widen SEG_BITS")
        elif c[0] == "R" and keep_cur:
            depth[cur] += 1
            rid = int(c[1])
            if max_reads is not None and rid >= max_reads:
                continue
            per_read[rid].append((int(c[3]) << OFF_SHIFT) | (cur << 1) | (int(c[2]) & 1))
    return seg_key, seg_id, slots_of_site, per_read, depth


def build_links(per_read):
    """(seg,side) pairs joined by reads. side 1 = the node's right end, side 0 = its left."""
    links = collections.Counter()
    for packed in per_read.values():
        if len(packed) < 2:
            continue
        v = sorted(packed)
        prev = None
        for x in v:
            seg, strand = (x >> 1) & SEG_MASK, x & 1
            if prev is not None and prev[0] != seg:
                # leaving prev: exit its right side if it was read forward, else its left
                u = (prev[0], 1 if prev[1] == 0 else 0)
                w = (seg, 0 if strand == 0 else 1)     # entering seg's left if forward
                links[(u, w) if u <= w else (w, u)] += 1
            prev = (seg, strand)
    return links


def prune(links, min_w, min_frac):
    """Keep a link unless it is a minority continuation at BOTH of the sides it joins.

    The rule we want is "merge when every read either continues to the same place or terminates".
    Taken literally a single chimeric read is a branch, so a link is dropped when it carries less
    than `min_frac` of the busiest link at each of its two sides -- relative, so it adapts to local
    depth instead of guessing an absolute count. A link that dominates either side survives, which
    stops a thin-but-only continuation from being cut loose.
    """
    best = collections.defaultdict(int)
    for (u, w), n in links.items():
        if u[0] == w[0]:
            continue
        best[u] = max(best[u], n)
        best[w] = max(best[w], n)
    kept = {}
    for (u, w), n in links.items():
        if u[0] == w[0] or n < min_w:
            continue
        if n < min_frac * best[u] and n < min_frac * best[w]:
            continue
        kept[(u, w)] = n
    return kept


def compact(nseg, links, min_w, min_frac=0.0):
    """Chain segments across links whose node side has degree one at BOTH ends."""
    kept = prune(links, min_w, min_frac)
    adj = collections.defaultdict(set)
    for (u, w) in kept:
        adj[u].add(w)
        adj[w].add(u)
    nxt = {}
    for (u, w) in kept:
        if len(adj[u]) == 1 and len(adj[w]) == 1:
            nxt[u], nxt[w] = w, u

    seen = bytearray(nseg)
    unitigs = []
    for s in range(nseg):
        if seen[s]:
            continue
        seen[s] = 1
        chain = collections.deque([(s, 0, 1)])          # (seg, in_side, out_side)
        cur, side = s, 1
        while (cur, side) in nxt:
            v, vs = nxt[(cur, side)]
            if seen[v]:
                break
            seen[v] = 1
            chain.append((v, vs, 1 - vs))
            cur, side = v, 1 - vs
        cur, side = s, 0
        while (cur, side) in nxt:
            v, vs = nxt[(cur, side)]
            if seen[v]:
                break
            seen[v] = 1
            chain.appendleft((v, 1 - vs, vs))
            cur, side = v, 1 - vs
        unitigs.append(list(chain))
    return unitigs, kept, nxt


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--vcf", help="call VCF; supplies positions (join on ID) so LN is real bases")
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--csv", help="Bandage annotation CSV")
    ap.add_argument("--region", help="contig[:start-end]; restrict to positioned sites inside it")
    ap.add_argument("--min-edge-reads", type=int, default=2,
                    help="absolute floor: drop links carried by fewer reads than this [2]")
    ap.add_argument("--min-edge-frac", type=float, default=0.15,
                    help="drop a link carrying less than this fraction of the busiest link at "
                         "BOTH of its sides; this is the 'all reads continue or terminate' rule [0.15]")
    ap.add_argument("--min-unitig-pins", type=int, default=1,
                    help="omit unitigs holding fewer pins than this [1]")
    a = ap.parse_args()

    region = parse_region(a.region)
    pos, gt = load_vcf(a.vcf)
    seg_key, seg_id, slots_of_site, per_read, depth = read_anchors(a.anchors, pos, region)
    print(f"pins {len(seg_key)}   reads {len(per_read)}   placements {sum(depth)}", file=sys.stderr)

    links = build_links(per_read)
    print(f"raw links {len(links)}", file=sys.stderr)
    unitigs, kept, nxt = compact(len(seg_key), links, a.min_edge_reads, a.min_edge_frac)
    print(f"links kept {len(kept)}   unitigs {len(unitigs)}",
          file=sys.stderr)

    where = {}                                   # (seg, side) -> (unitig, end) for terminal sides
    uinfo = []
    for ui, chain in enumerate(unitigs):
        where[(chain[0][0], chain[0][1])] = (ui, 0)
        where[(chain[-1][0], chain[-1][2])] = (ui, 1)
        segs = [c[0] for c in chain]
        ps = [pos[seg_key[s][1]] for s in segs if seg_key[s][1] in pos]
        span = (max(p[1] for p in ps) - min(p[1] for p in ps)) if len(ps) > 1 else 0
        d = sum(depth[s] for s in segs) / len(segs)
        cls = collections.Counter()
        for s in segs:
            _, snarl, slot = seg_key[s]
            n = len(slots_of_site[snarl])
            g = gt.get(snarl)
            if n == 2:
                cls[f"hap{slot}"] += 1
            elif g and g.replace("|", "/").split("/").count(".") == 1:
                cls["haploid"] += 1
            else:
                cls["collapsed"] += 1
        locus = (f"{ps[0][0]}:{min(p[1] for p in ps)}-{max(p[1] for p in ps)}") if ps else "off-ref"
        uinfo.append((len(segs), span, d, cls.most_common(1)[0][0], locus, ps))

    sz = sorted((u[0] for u in uinfo), reverse=True)
    tot = sum(sz)
    acc = 0
    n50 = 0
    for L in sz:
        acc += L
        if acc >= tot / 2:
            n50 = L
            break
    spans = sorted((u[1] for u in uinfo if u[1] > 0), reverse=True)
    print(f"unitig pins: largest {sz[0]}  N50 {n50}  median {sz[len(sz)//2]}  "
          f"singletons {sum(1 for x in sz if x == 1)}", file=sys.stderr)
    if spans:
        print(f"unitig span: largest {spans[0]/1000:.1f} kb  "
              f"median {spans[len(spans)//2]/1000:.2f} kb", file=sys.stderr)
    cls = collections.Counter(u[3] for u in uinfo)
    print(f"classes: {dict(cls.most_common())}", file=sys.stderr)

    big = [i for i, u in enumerate(uinfo) if u[0] >= a.min_unitig_pins]
    bigset = set(big)
    with open(a.gfa, "w") as f:
        f.write("H\tVN:Z:1.0\n")
        for i in big:
            n, span, d, cls, locus, _ = uinfo[i]
            f.write(f"S\tu{i}\t*\tLN:i:{max(1, span if span else n)}\tDP:f:{d:.1f}\tRC:i:{n}\n")
        nl = 0
        for (u, w), cnt in kept.items():
            if u in nxt and nxt[u] == w:
                continue                          # consumed inside a unitig
            if u not in where or w not in where:
                continue
            (ua, ue), (wa, we) = where[u], where[w]
            if ua == wa or ua not in bigset or wa not in bigset:
                continue
            f.write(f"L\tu{ua}\t{'+' if ue == 1 else '-'}\tu{wa}\t{'+' if we == 0 else '-'}"
                    f"\t0M\tRC:i:{cnt}\n")
            nl += 1
    print(f"wrote {len(big)} segments and {nl} links to {a.gfa}", file=sys.stderr)

    if a.csv:
        with open(a.csv, "w") as f:
            f.write("Name,Class,Pins,SpanBp,Depth,Locus\n")
            for i in big:
                n, span, d, cls, locus, _ = uinfo[i]
                f.write(f"u{i},{cls},{n},{span},{d:.1f},{locus}\n")
        print(f"wrote annotations to {a.csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
