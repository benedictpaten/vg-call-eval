#!/usr/bin/env python3
"""Can the phasing locate its own switch errors well enough to steer hom splitting?

The harm from splitting homozygous anchors is not that the split bridges a junction nothing
else bridges -- measured, it never is -- but that it lengthens a haplotype-labelled chain that
CONTAINS a switch. Contiguity bought against correctness. The remedy would be to suppress splits
around probable switches, which needs the switches located from the run's own evidence.

Evaluates each candidate signal as a ranking of stage-1 junctions, and prices it in the currency
that matters: what fraction of homozygous splits would have to be suppressed to cover what
fraction of the true switches.

`phase_link` returns a SUM of per-read log10 odds, so a junction can look confident while its
reads are split -- two confident readers disagreeing cancel. The vote decomposition is therefore
carried alongside |d|, and is the reason to expect it to beat |d|.
"""
import argparse, bisect, collections, statistics, sys


def load_junctions(path, want_pass=None):
    links, coh = [], {}
    passes = set()
    for line in open(path):
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        kind, pas = c[0], int(c[1])
        passes.add(pas)
        if kind == "link":
            links.append((pas, int(c[2]), int(c[3]), float(c[4]), int(c[5]), int(c[6]),
                          int(c[7]), float(c[8]), float(c[9]), c[10]))
        elif kind == "coh":
            coh.setdefault(pas, {})[int(c[2])] = (float(c[4]), int(c[5]))
    last = max(passes) if passes else 0
    links = [l for l in links if l[0] == last]
    # the coherence pass runs on passes BEFORE the last, so take the latest one that has rows
    cpass = max(coh) if coh else None
    return links, (coh.get(cpass, {}) if cpass is not None else {}), last, cpass


def load_switch_bed(path):
    rows = []
    for line in open(path):
        if line.startswith(("track", "#")):
            continue
        c = line.split()
        rows.append((int(c[1]), int(c[2])))
    rows.sort()
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
    isolated = [(r[0][0], r[0][1]) for r in runs if len(r) == 1]
    paired = [(r[0][0], r[-1][1]) for r in runs if len(r) == 2]
    longer = [(r[0][0], r[-1][1]) for r in runs if len(r) > 2]
    return isolated, paired, longer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--switch-bed", required=True)
    ap.add_argument("--split-positions", help="pickle of snarl -> pos, for pricing in splits")
    ap.add_argument("--split-index", help="pickle from homsplit_index.py, for the kind map")
    ap.add_argument("--windows", default="0,1000,5000",
                    help="bp added each side of a flagged junction: splits are suppressed there")
    a = ap.parse_args()

    links, coh, lastp, cpass = load_junctions(a.junctions)
    iso, paired, longer = load_switch_bed(a.switch_bed)
    print(f"stage-1 junctions (pass {lastp}): {len(links)}; coherence rows (pass {cpass}): {len(coh)}",
          file=sys.stderr)
    print(f"switch bed: {len(iso)} isolated (true switches), {len(paired)} paired (flips), "
          f"{len(longer)} longer runs", file=sys.stderr)

    targets = iso                       # long-range harm: a flip is local and self-correcting
    tstart = [t[0] for t in targets]

    def contains_switch(lo, hi):
        i = bisect.bisect_left(tstart, lo)
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(targets) and targets[j][0] < hi and lo <= targets[j][1]:
                return True
        return False

    rows = []
    for (_, pa, pb, d, n, nc, nt, mc, mt, flag) in links:
        lo, hi = min(pa, pb), max(pa, pb)
        minority_n = min(nc, nt) / n if n else 0.5
        minority_m = min(mc, mt) / (mc + mt) if (mc + mt) > 0 else 0.5
        c = min(coh.get(pa, (1.0, 0))[0], coh.get(pb, (1.0, 0))[0])
        rows.append({"lo": lo, "hi": hi, "d": abs(d), "n": n,
                     "minority_n": minority_n, "minority_m": minority_m,
                     "coh": c if c >= 0 else 1.0,
                     "hit": contains_switch(lo, hi), "brk": flag == "break"})
    covered_any = sum(1 for t in targets if any(r["lo"] < t[1] and t[0] <= r["hi"] for r in rows))
    print(f"{covered_any} of {len(targets)} true switches fall inside some stage-1 junction "
          f"interval ({sum(r['hit'] for r in rows)} junctions are hit)\n", file=sys.stderr)

    span = max(r["hi"] for r in rows) - min(r["lo"] for r in rows)
    splits = None
    if a.split_positions and a.split_index:
        import pickle
        posn = pickle.load(open(a.split_positions, "rb"))
        idx = pickle.load(open(a.split_index, "rb"))
        splits = sorted(posn[s] for s, k in idx["kind"].items()
                        if k == "homsplit" and s in posn)
        print(f"pricing against {len(splits)} split homozygotes", file=sys.stderr)

    def merge(sel, win):
        iv = sorted((r["lo"] - win, r["hi"] + win) for r in sel)
        merged = []
        for s, e in iv:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        return merged

    def price(sel, win):
        """bp flagged, splits suppressed, and SWITCHES covered, for a selection + window."""
        merged = merge(sel, win)
        bp = sum(e - s for s, e in merged)
        nsp = 0
        if splits is not None:
            for s, e in merged:
                nsp += bisect.bisect_right(splits, e) - bisect.bisect_left(splits, s)
        starts = [m[0] for m in merged]
        cov = 0
        for t in targets:
            i = bisect.bisect_right(starts, t[1]) - 1
            if i >= 0 and merged[i][1] >= t[0]:
                cov += 1
        return bp, nsp, cov

    signals = [("|d| (ascending)", lambda r: r["d"], False),
               ("minority read fraction", lambda r: r["minority_n"], True),
               ("minority evidence mass", lambda r: r["minority_m"], True),
               ("shared reads (ascending)", lambda r: r["n"], False),
               ("coherence of flanks", lambda r: r["coh"], False)]

    wins = [int(w) for w in a.windows.split(",")]
    print(f"\n{len(targets)} true switches, {covered_any} inside a junction interval; "
          f"{len(rows)} junctions over {span/1e6:.1f} Mb (one per {span//len(rows)} bp); "
          f"{len(splits) if splits else 0} splits")
    print("\nsignal\twindow\tflagged\tswitches covered\tbp flagged\tsplits suppressed")
    for name, key, desc in signals:
        order = sorted(rows, key=key, reverse=desc)
        for win in wins:
            for frac in (0.001, 0.01, 0.05, 0.25):
                k = max(1, int(frac * len(rows)))
                bp, nsp, cov = price(order[:k], win)
                sp = f"{100*nsp/len(splits):.1f}%" if splits else "-"
                print(f"{name}\t+/-{win}\t{frac*100:g}% ({k})\t{cov}/{len(targets)} "
                      f"({100*cov/len(targets):.0f}%)\t{100*bp/span:.1f}%\t{sp}")

    for win in wins:
        sel = [r for r in rows if r["brk"]]
        bp, nsp, cov = price(sel, win)
        sp = f"{100*nsp/len(splits):.1f}%" if splits else "-"
        print(f"ALREADY BROKEN (|d| < --phase-break)\t+/-{win}\t{len(sel)} "
              f"({100*len(sel)/len(rows):.1f}%)\t{cov}/{len(targets)} "
              f"({100*cov/len(targets):.0f}%)\t{100*bp/span:.1f}%\t{sp}")


if __name__ == "__main__":
    main()
