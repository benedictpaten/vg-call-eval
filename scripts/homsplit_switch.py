#!/usr/bin/env python3
"""Does hom splitting bridge the phasing's switch events?

A hom split labels its slots by sign(lambda), the read's cross-site strand against the settled
chain. It therefore always FOLLOWS the chain and can never invent a switch. What it can do is
carry an existing switch across a junction the anchor graph would otherwise have refused to
link: where het density is low -- often why the phasing switched there in the first place --
the splits supply the only haplotype-carrying path, and the graph then states a confident,
wrong pairing instead of breaking.

So the question is binary per junction:

    direct         some read pins at BOTH bracketing het sites -- hom splits are irrelevant
    split-bridged  no read does, but a chain through split homozygotes connects them
    unbridged      neither

and the cost of raising a gate is how many junctions move from split-bridged to unbridged,
measured at switch junctions against every other junction as the background.

`lo` is one number per read: read_strand_log_odds' leave-one-out never fires at a homozygous
site, because phase_sites skips trav_first == trav_second, so the value is the read's global
calibrated lambda and is the same at every homozygous site it crosses. That is what makes the
whole (--split-min-q, --split-min-side) sweep computable from a single run.
"""
import argparse, collections, gzip, math, pickle, random, statistics, sys


def load_vcf(path):
    """snarl id -> (pos, is_phased_het). The ID column is the snarl, which is the anchor key."""
    pos, phased_het = {}, set()
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            c = line.split("\t", 10)
            sid = c[2]
            if sid == ".":
                continue
            pos[sid] = int(c[1])
            gt = c[9].split(":", 1)[0]
            if "|" in gt:
                a, b = gt.split("|")[:2]
                if a != b and a != "." and b != ".":
                    phased_het.add(sid)
    return pos, phased_het


def load_lambda(path):
    lo = {}
    with open(path) as f:
        for line in f:
            if line[0] == "#":
                continue
            n, v = line.rstrip("\n").split("\t")
            lo[n] = float(v)
    return lo


def load_switch_bed(path):
    """whatshap --switch-error-bed. A switch is one interval; a FLIP is two contiguous ones."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith(("track", "#")):
                continue
            c = line.split()
            rows.append((c[0], int(c[1]), int(c[2])))
    rows.sort()
    runs, cur = [], []
    for r in rows:
        if cur and r[0] == cur[-1][0] and r[1] == cur[-1][2]:
            cur.append(r)
        else:
            if cur:
                runs.append(cur)
            cur = [r]
    if cur:
        runs.append(cur)
    switches = [(r[0][0], r[0][1], r[0][2]) for r in runs if len(r) == 1]
    flips = [(r[0][0], r[0][1], r[-1][2]) for r in runs if len(r) == 2]
    other = [r for r in runs if len(r) > 2]
    return switches, flips, other


def site_positions(idx, vcf_pos):
    """Every site gets a coordinate. Positioned sites take the VCF's; the rest -- 60.6% of
    splits have no VCF line -- take the MEDIAN position of the positioned sites their own
    reads also touch. Read-mediated rather than node-id interpolation: it needs no assumption
    that node ids run monotone with the reference, and it is bounded by a read length by
    construction."""
    posn = {}
    for s in idx["kind"]:
        p = vcf_pos.get(s)
        if p is not None:
            posn[s] = p
    unplaced = [s for s in idx["kind"] if s not in posn]
    for s in unplaced:
        seen = []
        for r in idx["reads_of"].get(s, ()):
            for t in idx["read_sites"].get(r, ()):
                p = posn.get(t)
                if p is not None:
                    seen.append(p)
        if seen:
            posn[s] = int(statistics.median(seen))
    return posn, len(unplaced)


def carrier_order(posn, carriers):
    order = sorted(carriers, key=lambda s: (posn[s], s))
    return order, {s: i for i, s in enumerate(order)}


def local_links(idx, order, rl, rr):
    """Shared-read counts among the carriers in [rl, rr], built on demand.

    Local rather than global because the windows are tiny: consecutive assessed het sites sit
    ~1,090 bp apart on chr20 while carriers run one per ~390 bp, so a junction spans a handful
    of carriers. Globally the same counts would be ~55M pairs -- every pair a 15 kb read makes
    among the ~38 carriers it crosses -- for no extra information.
    """
    ids = [set(idx["reads_of"].get(order[i], ())) for i in range(rl, rr + 1)]
    links = {}
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            n = len(ids[a] & ids[b])
            if n:
                links[(rl + a, rl + b)] = n
    return links


def reachable(rl, rr, links, allowed, min_link):
    """A chain of allowed carriers from rl to rr. Endpoints are always allowed: they are the
    het sites the junction is defined by, and the question is what lies between them."""
    if rl == rr:
        return True
    front, seen = [rl], {rl}
    while front:
        nxt = []
        for i in front:
            for j in range(i + 1, rr + 1):
                if j in seen or (j != rr and j not in allowed):
                    continue
                if links.get((i, j), 0) >= min_link:
                    if j == rr:
                        return True
                    seen.add(j)
                    nxt.append(j)
        front = nxt
    return False


def gate(los, min_q, frac, floor):
    """The hom-split gate, recomputed. side counts over reads clearing --split-min-q, then
    each side must reach max(floor, ceil(frac * confident))."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", required=True)
    ap.add_argument("--lambda-table", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--switch-bed", required=True)
    ap.add_argument("--contig", default="chr20")
    ap.add_argument("--min-link", type=int, default=2)
    ap.add_argument("--background", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--settings", default="0.5:0.0,2:0.0,4:0.0,8:0.0,"
                                          "0.5:0.10,0.5:0.20,0.5:0.33,"
                                          "2:0.20,4:0.20,4:0.33,8:0.33")
    a = ap.parse_args()

    idx = pickle.load(open(a.index, "rb"))
    vcf_pos, phased_het = load_vcf(a.vcf)
    los = load_lambda(a.lambda_table)
    switches, flips, other = load_switch_bed(a.switch_bed)
    print(f"switch bed: {len(switches)} switches, {len(flips)} flips, {len(other)} longer runs",
          file=sys.stderr)

    posn, n_unplaced = site_positions(idx, vcf_pos)
    placed = sum(1 for s in idx["kind"] if s in posn)
    print(f"sites {len(idx['kind'])}: {placed} placed "
          f"({n_unplaced} had no VCF line, read-mediated)", file=sys.stderr)

    kind = idx["kind"]
    names = idx["names"]
    homsplits = [s for s, k in kind.items() if k == "homsplit" and s in posn]
    hets = [s for s, k in kind.items() if k == "het" and s in posn]
    print(f"{len(homsplits)} split homozygotes, {len(hets)} het sites", file=sys.stderr)

    # per split site, its reads' lambdas -- the gate's whole input
    site_los = {}
    for s in homsplits:
        site_los[s] = [los.get(names[r], 0.0) for r in idx["reads_of"].get(s, ())]

    carriers = set(homsplits) | set(hets)
    order, rank = carrier_order(posn, carriers)
    print(f"{len(order)} haplotype-carrying sites in coordinate order", file=sys.stderr)

    # the junctions: consecutive ASSESSED het sites bracketing each switch, plus a background
    assessed = sorted((posn[s], s) for s in hets if s in phased_het)
    apos = [p for p, _ in assessed]
    import bisect

    def junction_for(lo_bp, hi_bp):
        i = bisect.bisect_left(apos, lo_bp)
        if i == 0 or i >= len(apos):
            return None
        return assessed[i - 1][1], assessed[i][1]

    switch_j = set()
    for c, lo_bp, hi_bp in switches:
        j = junction_for(lo_bp + 1, hi_bp + 1)   # whatshap BED is 0-based
        if j:
            switch_j.add(j)
    jset = []
    for i in range(len(assessed) - 1):
        j = (assessed[i][1], assessed[i + 1][1])
        jset.append(("switch" if j in switch_j else "background", j))
    print(f"{sum(1 for t, _ in jset if t == 'switch')} switch junctions resolved of "
          f"{len(switches)}, {sum(1 for t, _ in jset if t == 'background')} background",
          file=sys.stderr)

    # link counts per junction, computed once and reused across every setting
    windows = []
    for tag, (sl, sr) in jset:
        rl, rr = sorted((rank[sl], rank[sr]))
        gap = abs(posn[sr] - posn[sl])
        windows.append((tag, rl, rr, gap, local_links(idx, order, rl, rr)))
    gaps = sorted(w[3] for w in windows)
    print(f"{len(windows)} junction windows, median {statistics.median([w[2]-w[1] for w in windows]):.0f}"
          f" carriers and {statistics.median(gaps)} bp apart; "
          f"90th pct {gaps[int(0.9*len(gaps))]} bp, max {gaps[-1]} bp", file=sys.stderr)

    settings = []
    for tok in a.settings.split(","):
        q, f = tok.split(":")
        settings.append((float(q), float(f)))

    print("\nmin_q\tfrac\tsplit\tcollapsed\t| switch: direct  split-bridged  unbridged"
          "\t| background: direct  split-bridged  unbridged")
    for min_q, frac in settings:
        live = {s for s in homsplits if gate(site_los[s], min_q, frac, 2)}
        allowed = {rank[s] for s in hets} | {rank[s] for s in live}
        tally = {"switch": collections.Counter(), "background": collections.Counter()}
        for tag, rl, rr, gap, links in windows:
            if links.get((rl, rr), 0) >= a.min_link:
                tally[tag]["direct"] += 1
            elif reachable(rl, rr, links, allowed, a.min_link):
                tally[tag]["split"] += 1
            else:
                tally[tag]["unbridged"] += 1
        sw, bg = tally["switch"], tally["background"]
        print(f"{min_q}\t{frac}\t{len(live)}\t{len(homsplits)-len(live)}\t| "
              f"{sw['direct']}\t{sw['split']}\t{sw['unbridged']}\t| "
              f"{bg['direct']}\t{bg['split']}\t{bg['unbridged']}")

    # Stratified by how far apart the bracketing hets are. The aggregate is dominated by the
    # dense majority, where a 15 kb read spans both hets and nothing else can matter; if hom
    # splits ever carry a junction it is where that fails.
    buckets = [(0, 2000), (2000, 5000), (5000, 15000), (15000, 40000), (40000, 10**9)]
    for min_q, frac in [settings[0], settings[-1]]:
        live = {s for s in homsplits if gate(site_los[s], min_q, frac, 2)}
        allowed = {rank[s] for s in hets} | {rank[s] for s in live}
        print(f"\n-- by het gap, --split-min-q {min_q} --split-min-side max(2, {frac}*n) --")
        print("gap (bp)\tswitch: n  direct  split  unbr\tbackground: n  direct  split  unbr")
        for lo_b, hi_b in buckets:
            t = {"switch": collections.Counter(), "background": collections.Counter()}
            for tag, rl, rr, gap, links in windows:
                if not (lo_b <= gap < hi_b):
                    continue
                t[tag]["n"] += 1
                if links.get((rl, rr), 0) >= a.min_link:
                    t[tag]["direct"] += 1
                elif reachable(rl, rr, links, allowed, a.min_link):
                    t[tag]["split"] += 1
                else:
                    t[tag]["unbridged"] += 1
            sw, bg = t["switch"], t["background"]
            label = f"{lo_b}-{hi_b}" if hi_b < 10**9 else f">{lo_b}"
            print(f"{label}\t{sw['n']}\t{sw['direct']}\t{sw['split']}\t{sw['unbridged']}"
                  f"\t\t{bg['n']}\t{bg['direct']}\t{bg['split']}\t{bg['unbridged']}")


if __name__ == "__main__":
    main()
