#!/usr/bin/env python3
"""Stage 0 of the mosaic output: how many switches does the Viterbi path actually make?

The mosaic file records switch *points*, not sites, so its size is set by how piecewise the
inferred mosaic is -- not by the 105k sites on chr20. If the path switches rarely the file is
tens of KB; if it switches at a large fraction of sites the run-length encoding buys nothing and
the format should be per-site instead. Nothing measured so far answers that, because the shipped
model does marginal posterior decoding and never forms a path at all.

Two decodings are measured here, and the difference between them is the point:

  unconstrained   argmax over state paths of the joint probability. The most likely mosaic.
  constrained     the same, restricted at each site to states spelling the **called** genotype.

The constrained one is what the feature needs, because the emitted genome has to agree with the
emitted VCF. Its *path probability* can only be lower, since it maximises over a subset of the
states. Its *switch count* is not monotone in the same way and had to be measured: constraining
removes states, which takes away opportunities to switch as readily as it forces them.

Measured on chr20-34hap, 105,251 sites: 2,115 unconstrained switches against 2,064 constrained,
so consistency with the VCF costs essentially nothing (0.98x) rather than the inflation this was
written to look for. About 2% of sites are switch points, so run-length encoding buys a factor of
50 and Format A is viable: ~141 KB for chr20, ~6.6 MB genome-wide.

Reuses linkage_hmm_offline for loading, so the emission, the panel matrix and the chunking are
literally the same code the offline forward-backward was validated with.
"""

from __future__ import annotations

import argparse
import collections
import math
from pathlib import Path

import numpy as np

import linkage_hmm_offline as off

NEG = -np.inf


def _top2(x, axis):
    """(best, second best, argbest) along an axis. Second best is -inf if the axis has one entry."""
    arg = np.argmax(x, axis=axis)
    best = np.max(x, axis=axis)
    masked = np.copy(x)
    if axis == 1:
        masked[np.arange(x.shape[0]), arg] = NEG
    else:
        masked[arg, np.arange(x.shape[1])] = NEG
    second = np.max(masked, axis=axis)
    return best, second, arg


def _excl(best, second, arg, index):
    """max over the axis excluding position `index`, from precomputed top-2."""
    return np.where(arg == index, second, best)


def viterbi_path(recs, n_hap, w_t, scale, rho_min, escape, constrain):
    """Max-product over ordered pairs of panel haplotypes. Returns one (a, b) state per record.

    The transition factorises per strand -- `T = (1-rho)I + (rho/m)1`, so T(x->y) takes only two
    values -- but the *maximisation* does not separate, because delta(a,b) couples the strands.
    The standard reduction is by cases on which strands stayed:

        delta'(a',b') = e(a',b') + max of
            delta(a',b')                              + 2S      both stayed
            max_{b != b'}   delta(a',b)               + S + J   strand 1 stayed
            max_{a != a'}   delta(a,b')               + J + S   strand 2 stayed
            max_{a != a', b != b'} delta(a,b)         + 2J      both jumped

    with S = ln(1 - rho + rho/m) and J = ln(rho/m). Each leave-one-out maximum comes from a
    precomputed top-2 along the relevant axis, so the step is O(m^2) rather than the O(m^4) a
    literal pairwise maximisation would cost. That is what makes m = 35 over 100k sites tractable.

    In logs throughout: sum-product needs per-site rescaling to avoid underflow, max-product does
    not, and logs turn "stay versus jump" into a comparison of sums.
    """
    m = n_hap + 1
    n = len(recs)

    ems = []
    for r in recs:
        e = off.emission_matrix(r, n_hap, escape)
        with np.errstate(divide="ignore"):
            le = np.log(e)
        if constrain:
            le = np.where(_constraint_mask(r, n_hap, m), le, NEG)
        ems.append(le)

    delta = np.full((m, m), -math.log(m * m)) + ems[0]
    back_a = np.zeros((n, m, m), dtype=np.int8)
    back_b = np.zeros((n, m, m), dtype=np.int8)
    idx = np.arange(m)

    for t in range(1, n):
        gap = max(recs[t].pos - recs[t - 1].pos, 1)
        rho = off.weighted_switch(gap, w_t, scale, rho_min, 0.0, 0.0)
        rho = min(max(rho, 1e-12), 1.0 - 1e-12)
        S = math.log(1.0 - rho + rho / m)
        J = math.log(rho / m)

        D = delta
        # Row-wise top-2: for case "strand 1 stayed", exclude the arriving b'.
        rbest, rsecond, rarg = _top2(D, axis=1)
        # rowExcl[a, b'] = max_{b != b'} D[a, b]
        rowExcl = np.where(rarg[:, None] == idx[None, :], rsecond[:, None], rbest[:, None])
        # Column-wise top-2: for case "strand 2 stayed", exclude the arriving a'.
        cbest, csecond, carg = _top2(D, axis=0)
        colExcl = np.where(carg[None, :] == idx[:, None], csecond[None, :], cbest[None, :])
        # Both jumped: max over a != a' of rowExcl[a, b'].
        bbest, bsecond, barg = _top2(rowExcl, axis=0)
        bothExcl = np.where(barg[None, :] == idx[:, None], bsecond[None, :], bbest[None, :])

        c1 = D + 2.0 * S
        c2 = rowExcl + S + J
        c3 = colExcl + J + S
        c4 = bothExcl + 2.0 * J

        stacked = np.stack([c1, c2, c3, c4])
        which = np.argmax(stacked, axis=0)
        delta = np.take_along_axis(stacked, which[None], axis=0)[0] + ems[t]

        # Backpointers, resolved per case.
        rargExcl = np.where(rarg[:, None] == idx[None, :],
                            _argsecond(D, rarg, axis=1)[:, None], rarg[:, None])
        cargExcl = np.where(carg[None, :] == idx[:, None],
                            _argsecond(D, carg, axis=0)[None, :], carg[None, :])
        bargExcl = np.where(barg[None, :] == idx[:, None],
                            _argsecond(rowExcl, barg, axis=0)[None, :], barg[None, :])

        A = np.empty((m, m), dtype=np.int64)
        B = np.empty((m, m), dtype=np.int64)
        aa = np.broadcast_to(idx[:, None], (m, m))
        bb = np.broadcast_to(idx[None, :], (m, m))
        # case 1: (a', b')
        A = np.where(which == 0, aa, 0)
        B = np.where(which == 0, bb, 0)
        # case 2: (a', argmax_{b != b'})
        A = np.where(which == 1, aa, A)
        B = np.where(which == 1, rargExcl, B)
        # case 3: (argmax_{a != a'}, b')
        A = np.where(which == 2, cargExcl, A)
        B = np.where(which == 2, bb, B)
        # case 4: a from bargExcl; b is that row's own best excluding b'.
        A = np.where(which == 3, bargExcl, A)
        B = np.where(which == 3, rargExcl[bargExcl, bb], B)

        back_a[t] = A.astype(np.int8)
        back_b[t] = B.astype(np.int8)

    finite = np.where(np.isfinite(delta), delta, NEG)
    cur = np.unravel_index(int(np.argmax(finite)), delta.shape)
    path = [None] * n
    path[n - 1] = (int(cur[0]), int(cur[1]))
    for t in range(n - 1, 0, -1):
        a, b = path[t]
        path[t - 1] = (int(back_a[t, a, b]), int(back_b[t, a, b]))
    return path


def _argsecond(x, arg, axis):
    """Index of the second largest along an axis, given the argmax."""
    masked = np.copy(x)
    if axis == 1:
        masked[np.arange(x.shape[0]), arg] = NEG
        return np.argmax(masked, axis=1)
    masked[arg, np.arange(x.shape[1])] = NEG
    return np.argmax(masked, axis=0)


def _parse_gt(gt):
    """VCF GT string -> (i, j), or None where it is missing or not diploid."""
    if not gt or gt.startswith("."):
        return None
    toks = gt.replace("|", "/").split("/")
    if len(toks) != 2 or not all(t.isdigit() for t in toks):
        return None
    return int(toks[0]), int(toks[1])


def _constraint_mask(r, n_hap, m):
    """States whose two alleles spell the called genotype.

    A haplotype absent from the site, and the wildcard, may carry anything and so satisfy any
    constraint -- the same convention the emission uses for absence, and what keeps the
    constrained problem feasible where the panel cannot spell a call.
    """
    # Record uses __slots__, so the parse happens here rather than being cached on the record.
    gt = _parse_gt(r.gt)
    if gt is None:
        return np.full((m, m), True)
    alleles = np.full(m, -1, dtype=np.int64)
    if r.hap_allele:
        for h, a in r.hap_allele.items():
            alleles[h] = a
    free = alleles < 0
    free[n_hap] = True
    gi, gj = gt
    mask = np.outer((alleles == gi) | free, (alleles == gj) | free)
    if gi != gj:
        mask |= np.outer((alleles == gj) | free, (alleles == gi) | free)
    return mask


def count_switches(path):
    sw = [0, 0]
    for t in range(1, len(path)):
        for s in (0, 1):
            if path[t][s] != path[t - 1][s]:
                sw[s] += 1
    return sw


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls", required=True, help="emitted VCF with GL and AT")
    p.add_argument("--panel", required=True, help="vg deconstruct VCF for the panel matrix")
    p.add_argument("--max-sites", type=int, default=0,
                   help="stop after this many records (0 = all)")
    p.add_argument("--weight", type=float, default=2.0)
    p.add_argument("--scale", type=float, default=10000.0)
    p.add_argument("--rho-min", type=float, default=1e-3)
    p.add_argument("--escape", type=float, default=1e-2)
    p.add_argument("--max-gap", type=int, default=1_000_000)
    p.add_argument("--chunk", type=int, default=2000)
    p.add_argument("--margin", type=int, default=250)
    args = p.parse_args()

    stats = collections.Counter()
    panel, n_hap = off.load_panel(Path(args.panel))
    recs, _header = off.load_calls(Path(args.calls), panel, stats)
    if args.max_sites:
        recs = recs[: args.max_sites]
    mapped = sum(1 for r in recs if r.hap_allele)
    print(f"{len(recs)} records, {n_hap} panel haplotypes, "
          f"{mapped} with a panel mapping ({mapped / max(len(recs), 1) * 100:.1f}%)", flush=True)
    if stats:
        print("  loader notes:", dict(stats), flush=True)

    totals = {"unconstrained": [0, 0], "constrained": [0, 0]}
    sites = 0
    n_seg = 0
    for seg in off.segments(recs, args.max_gap):
        n_seg += 1
        for sub, lo, hi in off.chunks(seg, args.chunk, args.margin):
            sites += hi - lo
            for name, constrain in (("unconstrained", False), ("constrained", True)):
                path = viterbi_path(sub, n_hap, args.weight, args.scale, args.rho_min,
                                    args.escape, constrain)
                # Count only inside the interior, so chunk seams are not counted as switches --
                # the artifact Stage 2 of the plan exists to prevent in the real implementation.
                sw = count_switches(path[lo:hi])
                totals[name][0] += sw[0]
                totals[name][1] += sw[1]
        print(f"  segment {n_seg}: {sites} interior sites so far", flush=True)

    print(f"\ninterior sites scored: {sites} over {n_seg} chain segment(s)")
    for name in ("unconstrained", "constrained"):
        a, b = totals[name]
        tot = a + b
        print(f"  {name:14s} switches: strand0={a} strand1={b} total={tot} "
              f"({tot / max(sites, 1) * 100:.3f}% of sites)")
    u = sum(totals["unconstrained"]) or 1
    c = sum(totals["constrained"])
    print(f"\nconstraining inflates switches by {c / u:.2f}x")
    seg_lines = c + 2 * n_seg
    print(f"mosaic segments implied: {seg_lines}  ->  ~{seg_lines * 70 / 1024:.0f} KB at 70 B/line")


if __name__ == "__main__":
    main()
