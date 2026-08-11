#!/usr/bin/env python3
"""Does panel linkage actually break at haplotype-sampling subchain boundaries?

This was the gate on `--linkage-block-switch`. It did not pass, and the flag is gone; this is
kept as the record of why, and as the measurement to redo if a graph with a genuinely different
sampling structure ever makes the question live again.

Why it is needed. As implemented, `beta` spreads its switch mass as `crossings = gap /
block_length`, and that makes it *algebraically identical* to shortening the distance scale:

    1 - rho' = (1 - rho_bio) * (1 - beta)^(gap/L)
             = (1 - rho_min) * exp(-gap/s) * exp(-gap * -ln(1-beta)/L)
             = (1 - rho_min) * exp(-gap/s_eff),   1/s_eff = 1/s + -ln(1-beta)/L

So `beta = 0.57` at `L = 10 kb` is exactly `--linkage-scale 5423`, and the whole beta axis is a
duplicate of an axis the model already has. It cannot be tuned to anything new. Measured, that
axis is close to flat -- 10 kb to 40 kb moves four-dataset F1 by about 0.001 -- so no value of
smeared beta, in either direction, is worth anything.

Beta becomes a distinct mechanism only when the crossing count is an integer drawn from known
boundary positions, so that it is 0 for the ~90% of adjacent site pairs that sit inside one
block and >= 1 for the rest. That is a different shape, not a reparameterisation. But it costs a
new required input (a `.hapl`), so it should be paid for only if the effect is there.

What is measured. Linkage in the *panel matrix itself*, not in the calls. The claim under test is
about the panel: within a subchain a sampled haplotype is a contiguous piece of one assembly, so
allele sharing between haplotypes is stable; at a boundary the sampler continues the same
assembly only some of the time, so the correspondence reshuffles. That is a property of the
haplotype-by-site matrix and needs no genotyping to see. Measuring it in the calls instead would
add caller noise to a question the calls have no bearing on.

The statistic is normalised mutual information between the two sites' allele assignments over
the haplotypes typed at both. NMI rather than an r^2-style measure because panel sites here are
routinely multi-allelic, and rather than the "is any single haplotype consistent" indicator used
in `apparent_recombination.py` because that saturates: with 32 haplotypes almost every allele
pair is carried by someone, so it cannot see a *weakening* of linkage, only its total absence.

**The comparison must be gap-matched, and coarse bins are not enough.** Boundary-crossing pairs
are farther apart than non-crossing pairs -- that is what a boundary is -- so a pooled comparison
rediscovers the distance decay and calls it a boundary effect. Wide bins do not fix this: the
crossing probability rises with gap *within* a bin too, so crossing pairs pile up at the top of
each bin while NMI falls across it. Measured, that residual confound is larger than the effect
being looked for: with six decade-ish bins, permuted boundaries produced a bigger apparent drop
below 200 bp (0.090) than the real ones did (0.059).

So the strata here are narrow -- every distinct gap up to 64 bp, then multiplicative steps of 6%
-- and the estimate is the n-weighted mean of the within-stratum difference. The permutation
control is the check that this worked: with matched gaps it should sit at zero, and if it does
not, the design is still leaking and the real number cannot be read.

Reading a null. A drop at boundaries is conclusive: beta is real and the `.hapl` plumbing is
worth writing. *No* drop is ambiguous, because the boundaries used here are recomputed on the
post-sampling graph rather than the pre-sampling pangenome the sampler actually partitioned
(that graph is not on disk). Same 10 kb target and same snarl-walking rule, but a different
topology to walk, so the partition can be offset. A misphased partition would wash out a real
effect. So a null sends the question to the full pangenome, it does not close it.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import gzip
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_boundaries(tsv: Path, contig: str | None) -> list[int]:
    """Subchain boundary positions along the reference, from `vg haplotypes --statistics`.

    Rows are `type, id, start, end, length, kmers, sequences` after a `C` line naming the chain;
    `start`/`end` are a semiopen reference range for the subchain *interior*, so consecutive
    subchains do not abut -- they are separated by the boundary node and any unary path, a few
    hundred bp here. One junction is therefore one switch opportunity, not two, and it is placed
    at the midpoint of that gap. Taking both endpoints instead would leave a pair spanning a
    single junction with a crossing count of 2, and `(1-beta)^k` squares the penalty.
    """
    intervals, current, keep = [], None, True
    with open_maybe_gz(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or f[0] == "H":
                continue
            if f[0] == "C":
                current = f[2] if len(f) > 2 else None
                keep = contig is None or current == contig
                continue
            if not keep or len(f) < 4 or f[0] not in ("N", "P", "S", "F"):
                continue
            try:
                intervals.append((int(f[2]), int(f[3])))
            except ValueError:
                continue
    intervals.sort()
    return [(a[1] + b[0]) // 2 for a, b in zip(intervals, intervals[1:]) if b[0] >= a[1]]


def load_panel(vcf: Path) -> list[tuple[int, dict]]:
    """[(position, {haplotype -> allele})] in reference order, polymorphic sites only.

    A haplotype is (column, phase): a sampled GBZ files every recombinant under one sample name,
    so the phase index is the only thing separating them. An untyped haplotype is omitted rather
    than called reference -- a haplotype that does not traverse the site carries no allele, and
    counting that as reference would manufacture the very linkage this is trying to measure.
    """
    sites = []
    with open_maybe_gz(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            alleles = {}
            for col, cell in enumerate(f[9:]):
                gt = cell.split(":")[0]
                for phase, tok in enumerate(gt.replace("/", "|").split("|")):
                    if tok.isdigit():
                        alleles[(col, phase)] = int(tok)
            if len(set(alleles.values())) < 2:
                continue          # monomorphic: entropy zero, NMI undefined
            sites.append((int(f[1]), alleles))
    sites.sort(key=lambda s: s[0])
    return sites


def nmi(a: dict, b: dict) -> float | None:
    """Normalised mutual information over the haplotypes typed at both sites.

    Normalised by min(H(A), H(B)) so that 1.0 means one site determines the other even when the
    two differ in allele count -- the asymmetric case is common here and dividing by the mean or
    the joint entropy would report perfect linkage as partial.
    """
    shared = [h for h in a if h in b]
    n = len(shared)
    if n < 6:
        return None
    ca, cb, cab = collections.Counter(), collections.Counter(), collections.Counter()
    for h in shared:
        ca[a[h]] += 1
        cb[b[h]] += 1
        cab[(a[h], b[h])] += 1
    if len(ca) < 2 or len(cb) < 2:
        return None

    def ent(counter):
        return -sum((c / n) * math.log(c / n) for c in counter.values() if c)

    ha, hb = ent(ca), ent(cb)
    if ha <= 0 or hb <= 0:
        return None
    mi = 0.0
    for (x, y), c in cab.items():
        pxy = c / n
        mi += pxy * math.log(pxy / ((ca[x] / n) * (cb[y] / n)))
    return max(0.0, min(1.0, mi / min(ha, hb)))


def gap_stratum(d: int) -> int:
    """Narrow gap strata: exact below 64 bp, then ~6% multiplicative steps.

    Fine enough that NMI is near-constant inside a stratum, which is the whole point -- the
    difference between crossing and non-crossing pairs must not be able to hide a gap difference.
    """
    if d < 64:
        return d
    return 64 + int(math.log(d / 64.0) / math.log(1.06))


def matched_drop(pairs, bounds) -> tuple[float, int]:
    """n-weighted mean of (NMI inside - NMI crossing) computed within each narrow gap stratum.

    Strata with no crossing pair or no non-crossing pair contribute nothing: there is no matched
    comparison to make, and including them as zero would dilute toward the null in proportion to
    how rare crossings are.
    """
    acc = collections.defaultdict(lambda: {"in": [], "cross": []})
    for gap, pa, pb, v in pairs:
        k = bisect.bisect_right(bounds, pb) - bisect.bisect_right(bounds, pa)
        acc[gap_stratum(gap)]["cross" if k > 0 else "in"].append(v)
    num = den = 0.0
    for d in acc.values():
        if not d["in"] or not d["cross"]:
            continue
        w = len(d["cross"])
        num += w * (sum(d["in"]) / len(d["in"]) - sum(d["cross"]) / len(d["cross"]))
        den += w
    return (num / den if den else float("nan")), int(den)


REPORT_BINS = [(0, 200, "<200"), (200, 500, "200-500"), (500, 1000, "500-1k"),
               (1000, 2000, "1k-2k"), (2000, 5000, "2k-5k"), (5000, 20000, "5k-20k")]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, type=Path)
    ap.add_argument("--subchains", required=True, type=Path)
    ap.add_argument("--contig", default=None)
    ap.add_argument("--permutations", type=int, default=10)
    ap.add_argument("--out", type=Path, default=WORK / "sv-atlas" / "subchain-linkage.json")
    args = ap.parse_args()

    bounds = load_boundaries(args.subchains, args.contig)
    sites = load_panel(args.panel)
    if not bounds or not sites:
        sys.exit(f"nothing to do: {len(bounds)} boundaries, {len(sites)} panel sites")

    lo, hi = sites[0][0], sites[-1][0]
    covered = sum(1 for b in bounds if lo <= b <= hi)
    print(f"boundaries: {len(bounds)} ({covered} inside the panel's span {lo}-{hi})")
    print(f"panel sites: {len(sites)} polymorphic")
    if covered < 10:
        sys.exit("boundary positions do not overlap the panel: wrong contig or wrong coordinates")
    print(f"median block: {(hi - lo) / max(covered, 1):.0f} bp")

    # NMI is the expensive part and does not depend on the boundaries, so it is computed once and
    # reused across the real run and every permutation.
    pairs = []
    for (pa, aa), (pb, ab) in zip(sites, sites[1:]):
        gap = pb - pa
        if gap <= 0 or gap >= 20000:
            continue
        v = nmi(aa, ab)
        if v is not None:
            pairs.append((gap, pa, pb, v))
    print(f"adjacent pairs scored: {len(pairs)}")

    rng = random.Random(20260811)
    real, n_real = matched_drop(pairs, bounds)
    ctrl = []
    for _ in range(args.permutations):
        fake = sorted(rng.randrange(lo, hi) for _ in range(covered))
        ctrl.append(matched_drop(pairs, fake)[0])
    cmean = sum(ctrl) / len(ctrl)
    csd = (sum((c - cmean) ** 2 for c in ctrl) / max(len(ctrl) - 1, 1)) ** 0.5

    print(f"\ngap-matched drop (real):    {real:+.4f}   over {n_real} matched crossing pairs")
    print(f"gap-matched drop (control): {cmean:+.4f} +/- {csd:.4f}  ({len(ctrl)} permutations)")
    z = (real - cmean) / csd if csd > 0 else float("nan")
    print(f"z against the permutation null: {z:+.1f}")

    # Reported per range as well, because a break that is real should be concentrated in the gap
    # range where sites actually sit -- if it appears only in ranges holding a handful of pairs,
    # it is noise wearing the right sign.
    print(f"\n{'gap':>9s} {'matched drop':>13s} {'control':>9s} {'n cross':>9s}")
    print("-" * 44)
    rows = []
    for lo_g, hi_g, name in REPORT_BINS:
        sub = [p for p in pairs if lo_g <= p[0] < hi_g]
        if not sub:
            continue
        r, n = matched_drop(sub, bounds)
        c = sum(matched_drop(sub, sorted(rng.randrange(lo, hi) for _ in range(covered)))[0]
                for _ in range(3)) / 3.0
        print(f"{name:>9s} {r:13.4f} {c:9.4f} {n:9d}")
        rows.append({"bin": name, "matched_drop": r, "control": c, "n_cross": n})

    verdict = ("boundaries carry a linkage break the gap-matched control does not reproduce"
               if z > 3 and real > 0 else
               "no boundary effect separable from gap composition -- see the docstring on why a "
               "null here is ambiguous")
    print(f"\nverdict: {verdict}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"boundaries": len(bounds), "sites": len(sites), "pairs": len(pairs),
         "matched_drop_real": real, "matched_drop_control_mean": cmean,
         "matched_drop_control_sd": csd, "z": z, "bins": rows, "verdict": verdict}, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
