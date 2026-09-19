#!/usr/bin/env python3
"""Where are the switches, and can splits be suppressed around them affordably?

Stage 1 cascades over reliable sites and BREAKS where |d| < --phase-break; stage 2 then relinks
the blocks across each break with a K x K pairwise vote. A switch is a wrong orientation, so it
is either a wrong sign in the cascade or a wrong relink -- and it turns out to be overwhelmingly
the relink. That localises the problem to 5.7% of junctions before any new signal is computed.

This prices suppression in splits, not base pairs: for a flagged junction, suppress every split
inside its interval plus the `--adjacent` nearest splits on each side, which is what "directly
between them or adjacent to them" means for an anchor file.
"""
import argparse, bisect, collections, pickle, statistics, sys


def load(path):
    link, relink, coh = [], [], {}
    passes = collections.defaultdict(set)
    for line in open(path):
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        kind, pas = c[0], int(c[1])
        passes[kind].add(pas)
        rec = (pas, int(c[2]), int(c[3]), float(c[4]), int(c[5]), int(c[6]), int(c[7]),
               float(c[8]), float(c[9]), c[10])
        if kind == "link":
            link.append(rec)
        elif kind == "relink":
            relink.append(rec)
        else:
            coh.setdefault(pas, {})[int(c[2])] = float(c[4])
    lp = max(passes["link"]) if passes["link"] else 0
    rp = max(passes["relink"]) if passes["relink"] else 0
    cp = max(coh) if coh else None
    return ([r for r in link if r[0] == lp], [r for r in relink if r[0] == rp],
            coh.get(cp, {}) if cp is not None else {})


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
    ap.add_argument("--split-positions", required=True)
    ap.add_argument("--split-index", required=True)
    ap.add_argument("--adjacent", type=int, default=2,
                    help="splits suppressed on EACH side of a flagged junction, beyond those inside")
    ap.add_argument("--include-longer", action="store_true",
                    help="count runs of >2 BED intervals as switches too")
    a = ap.parse_args()

    link, relink, coh = load(a.junctions)
    iso, paired, longer = switches(a.switch_bed)
    targets = iso + (longer if a.include_longer else [])
    posn = pickle.load(open(a.split_positions, "rb"))
    idx = pickle.load(open(a.split_index, "rb"))
    splits = sorted(posn[s] for s, k in idx["kind"].items() if k == "homsplit" and s in posn)
    print(f"{len(link)} stage-1 junctions, {len(relink)} stage-2 relinks, "
          f"{len(coh)} coherence rows, {len(splits)} splits", file=sys.stderr)
    print(f"{len(targets)} true switches ({len(iso)} isolated, {len(longer)} longer runs), "
          f"{len(paired)} flips\n", file=sys.stderr)

    def mk(rows, is_relink):
        out = []
        for (_, pa, pb, d, n, npos, nneg, mp, mn, flag) in rows:
            lo, hi = min(pa, pb), max(pa, pb)
            out.append({"lo": lo, "hi": hi, "d": abs(d), "n": n,
                        "minf": min(npos, nneg) / n if n else 0.5,
                        "minm": min(mp, mn) / (mp + mn) if (mp + mn) > 0 else 0.5,
                        "gap": hi - lo,
                        "coh": min(coh.get(pa, 1.0), coh.get(pb, 1.0)),
                        "brk": flag == "break", "relink": is_relink})
        return out

    L, R = mk(link, False), mk(relink, True)

    def price(sel):
        """splits suppressed, and switches covered, for a selection."""
        iv = []
        for r in sel:
            i = bisect.bisect_left(splits, r["lo"]) - a.adjacent
            j = bisect.bisect_right(splits, r["hi"]) + a.adjacent
            i, j = max(0, i), min(len(splits), j)
            if j > i:
                iv.append((i, j, splits[i], splits[j - 1]))
        iv.sort()
        merged, nsp = [], 0
        for i, j, s, e in iv:
            if merged and i <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], j)
                merged[-1][3] = max(merged[-1][3], e)
            else:
                merged.append([i, j, s, e])
        nsp = sum(m[1] - m[0] for m in merged)
        # a switch is covered if it falls inside a flagged junction interval (not the padding:
        # the padding is cost, not reach)
        ivals = sorted((r["lo"], r["hi"]) for r in sel)
        st = [x[0] for x in ivals]
        cov = 0
        for t in targets:
            k = bisect.bisect_right(st, t[1]) - 1
            for q in (k, k - 1):
                if 0 <= q < len(ivals) and ivals[q][1] >= t[0]:
                    cov += 1
                    break
        return nsp, cov

    def report(label, rows, sigs, fracs):
        print(f"\n== {label} ({len(rows)} junctions) ==")
        print("signal\tflagged\tswitches covered\tsplits suppressed\tlift")
        nsp_all, cov_all = price(rows)
        print(f"ALL\t{len(rows)} (100%)\t{cov_all}/{len(targets)}\t"
              f"{100*nsp_all/len(splits):.1f}%\t"
              f"{(100*cov_all/len(targets))/(100*nsp_all/len(splits)):.2f}x")
        for name, key, desc in sigs:
            order = sorted(rows, key=key, reverse=desc)
            for f in fracs:
                k = max(1, int(f * len(rows)))
                nsp, cov = price(order[:k])
                c = 100 * cov / len(targets)
                s = 100 * nsp / len(splits)
                print(f"{name}\t{f*100:g}% ({k})\t{cov}/{len(targets)} ({c:.0f}%)\t{s:.1f}%\t"
                      f"{c/s if s else 0:.2f}x")

    sigs = [("|d| asc", lambda r: r["d"], False),
            ("minority pair frac", lambda r: r["minf"], True),
            ("minority mass frac", lambda r: r["minm"], True),
            ("n votes asc", lambda r: r["n"], False),
            ("gap bp desc", lambda r: r["gap"], True),
            ("coherence asc", lambda r: r["coh"], False)]
    report("stage-2 relinks -- the 5.7% of junctions where stage 1 declined to link",
           R, sigs, (0.01, 0.05, 0.10, 0.25, 0.50))
    report("all stage-1 junctions", L, sigs, (0.01, 0.05, 0.25))
    print("\nlift 1.00x is random: covering X% of switches by suppressing X% of the splits.")


if __name__ == "__main__":
    main()
