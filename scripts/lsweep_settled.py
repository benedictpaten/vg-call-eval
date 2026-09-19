#!/usr/bin/env python3
"""The L-sweep, compared against the orientation the algorithm ACTUALLY settled on.

The span rows are emitted during stage 1 and carry the parity implied by the signs of `d`. That
is the wrong baseline wherever stage 2 relinked: there the orientation came from the K x K vote,
not from `d`. So the first pass of this analysis silently excluded every relink event -- which is
exactly the population the containment idea targets. With the settled `o[]` dumped after stage 3,
the parity between any two sites is `o[a] ^ o[b]` and every span can be judged against what the
algorithm decided.
"""
import argparse, bisect, collections, gzip, pickle, sys


def second_batch(rows, keyfn):
    """read_phase_flips runs twice -- once inside re-genotyping, once at render. Take the render
    traversal: the second time the chain index resets to zero."""
    starts = [i for i, r in enumerate(rows) if keyfn(r) == 0]
    return rows[starts[-1]:] if len(starts) > 1 else rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--keymap", required=True)
    ap.add_argument("--switch-bed", required=True)
    ap.add_argument("--correct-vcf", required=True)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()

    link, span, settled = [], [], []
    for line in open(a.junctions):
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        if c[0] == "link" and int(c[1]) == 2:
            link.append((int(c[2]), int(c[3]), int(c[4]), int(c[5]), float(c[6]), c[8] == "1"))
        elif c[0] == "span" and int(c[1]) == 2:
            span.append((int(c[2]), int(c[3]), int(c[6]), float(c[7]), int(c[8]), int(c[10])))
        elif c[0] == "settled":
            settled.append((int(c[2]), int(c[4]), int(c[6]), int(c[7])))
    d = [i for i in range(1, len(link)) if link[i][2] < link[i - 1][2]]
    link = link[d[0]:] if d else link
    settled = second_batch(settled, lambda r: r[2])
    o = {r[0]: r[3] for r in settled}            # record_key -> settled orientation
    chain = {r[0]: r[2] for r in settled}        # record_key -> chain index
    print(f"{a.tag}: {len(link):,} junctions, {sum(1 for r in link if r[5]):,} breaks, "
          f"{len(settled):,} settled sites", file=sys.stderr)

    order = [link[0][0]] + [r[1] for r in link]
    idx = {k: i for i, k in enumerate(order)}
    brk = [r[5] for r in link]

    km = {}
    for line in open(a.keymap):
        if line[0] == "#":
            continue
        k, s = line.rstrip("\n").split("\t")
        km[int(k)] = s
    inv = {v: k for k, v in km.items()}
    ass = []
    for line in gzip.open(a.correct_vcf, "rt"):
        if line[0] == "#":
            continue
        c = line.split("\t", 4)
        k = inv.get(c[2])
        if k is not None and k in idx:
            ass.append((int(c[1]), idx[k]))
    ass.sort()
    apos = [x[0] for x in ass]

    raw = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(a.switch_bed)
                 if not l.startswith(("track", "#")))
    runs, cur = [], []
    for r in raw:
        if cur and r[0] == cur[-1][1]:
            cur.append(r)
        else:
            if cur:
                runs.append(cur)
            cur = [r]
    if cur:
        runs.append(cur)
    sw = [(r[0][0], r[-1][1]) for r in runs if len(r) % 2 == 1]

    events = []
    for t0, t1 in sw:
        i = min(max(bisect.bisect_left(apos, t0 + 1), 0), len(ass) - 1)
        j = min(max(bisect.bisect_left(apos, t1 + 1), 0), len(ass) - 1)
        ci, cj = sorted((ass[i][1], ass[j][1]))
        js = list(range(ci, cj))
        events.append((js, any(brk[x] for x in js)))
    allbad = set(x for js, _ in events for x in js)

    # spans, keyed by junction chain index
    kp = {(r[0], r[1]): i for i, r in enumerate(link)}
    sp = collections.defaultdict(dict)
    for ka, kb, L, skip, shared, bp in span:
        i = kp.get((ka, kb))
        if i is None:
            continue
        aidx, bidx = i + 1 - L, i + L
        if aidx < 0 or bidx >= len(order):
            continue
        ka_, kb_ = order[aidx], order[bidx]
        if ka_ not in o or kb_ not in o:
            continue
        parity = o[ka_] ^ o[kb_]
        if skip == 0.0:
            verdict = 2
        else:
            verdict = 1 if ((skip < 0.0) == (parity == 1)) else 0
        sp[i][L] = (verdict, shared, bp, abs(skip))

    Ls = [2, 5, 10, 20, 50]
    print(f"\n== {a.tag}: spans judged against the SETTLED orientation ==")
    print("  L    tests     disagree (background)   at switch junctions (informative only)")
    for L in Ls:
        have = [(i, v[L]) for i, v in sp.items() if L in v and v[L][0] != 2]
        bg = sum(1 for _, x in have if x[0] == 0)
        n = dis = 0
        for i in allbad:
            x = sp.get(i, {}).get(L)
            if x is None or x[0] == 2:
                continue
            lo, hi = i + 1 - L, i + L - 1
            if sum(1 for j in allbad if lo <= j <= hi) % 2 == 0:
                continue
            n += 1
            dis += 1 if x[0] == 0 else 0
        print(f"  {L:<3}  {len(have):>7,}   {bg:>5,} ({100*bg/len(have):.2f}%)"
              f"{'':>10}{dis}/{n}" + (f" ({100*dis/n:.0f}%)" if n else ""))

    nrel = sum(1 for _, r in events if r)
    hit = hitrel = hitcas = 0
    for js, isrel in events:
        got = any(sp.get(i, {}).get(L, (1,))[0] == 0 for i in js for L in Ls)
        hit += got
        if isrel:
            hitrel += got
        else:
            hitcas += got
    flagged = len({i for i, v in sp.items() if any(v.get(L, (1,))[0] == 0 for L in Ls)})
    print(f"\n  events with >=1 disagreeing span: {hit}/{len(events)}")
    print(f"    of the {nrel} RELINK events:  {hitrel}")
    print(f"    of the {len(events)-nrel} CASCADE events: {hitcas}")
    print(f"  junctions flagged: {flagged:,}/{len(link):,} ({100*flagged/len(link):.2f}%)"
          f"  -> lift {(hit/len(events))/(flagged/len(link)):.0f}x")


if __name__ == "__main__":
    main()
