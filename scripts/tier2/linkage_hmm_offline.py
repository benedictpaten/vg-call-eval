#!/usr/bin/env python3
"""Stage 0 of the linkage HMM: run the whole model offline, with no changes to vg.

`vg call` already emits everything the HMM consumes, which is why this stage costs a day
rather than a week:

    per-site ln P(reads | G) for every genotype   GL, log10-scaled, per record
    haplotype -> allele matrix                    vg deconstruct, joined on snarl ID
    site order and gaps                           POS
    allele identity across the two files          AT (allele traversal)

So the emission is fixed and cached, and `w_t`, `beta` and `L` can be swept freely without a
single run of `vg call`. If this does not improve genotype concordance where the reads were
undecided, no C++ gets written. See planning/vg-call-linkage-hmm.md sections 7 and 9.

The model, as specified in that note. Hidden states are ordered pairs of panel haplotypes plus
a wildcard for alleles the panel does not carry; the emission for a state is the existing
per-site likelihood of the genotype that state implies; transitions are Li-Stephens with a
distance term and a boost at haplotype-sampling block boundaries. Inference is forward-backward,
and the call is `argmax` over **genotypes** of the summed state posterior -- not the genotype of
the most likely state, which is a different and wrong answer when many states imply one genotype.

Ordered rather than unordered pairs, deliberately: it doubles the state count to N^2, which is
nothing at N = 34, and it makes the transition factorise exactly, so the forward step is O(N^2)
per site instead of O(N^4). Genotypes are symmetrised when the posterior is summed.

Approximations, both of which only *understate* the model:

  * `GL` covers the alleles that reached the record, not every allele enumerated at the site. A
    site where many alleles were enumerated and few emitted gives the HMM fewer states than the
    real implementation would have.
  * Chain structure is approximated by reference order with a gap cutoff, because the snarl tree
    is not available here. Fragment boundaries are not modelled, so linkage is allowed to cross
    places where a real implementation would cut it -- again conservative for the *gain*, though
    it means a positive result needs re-confirming in Stage 1.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def info_field(info: str, key: str):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv[len(key) + 1:]
    return None


# --------------------------------------------------------------------------- panel

def load_panel(vcf: Path):
    """snarl -> (allele-traversal list, {haplotype index: allele index}), plus the haplotype list.

    A haplotype is (column, phase): a sampled GBZ puts every recombinant under one sample name,
    so phase is the only thing that distinguishes them.
    """
    per_snarl = {}
    haps = set()
    with open_maybe_gz(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or not f[2].startswith(">"):
                continue
            at = info_field(f[7], "AT")
            if at is None:
                continue
            assign = {}
            for col, cell in enumerate(f[9:]):
                for phase, tok in enumerate(cell.split(":")[0].replace("/", "|").split("|")):
                    if tok.isdigit():
                        assign[(col, phase)] = int(tok)
                        haps.add((col, phase))
            if assign:
                per_snarl[f[2]] = (at.split(","), assign)
    hap_list = sorted(haps)
    hap_index = {h: i for i, h in enumerate(hap_list)}
    out = {}
    for snarl, (at, assign) in per_snarl.items():
        out[snarl] = (at, {hap_index[h]: a for h, a in assign.items()})
    return out, len(hap_list)


# --------------------------------------------------------------------------- calls

def gt_index(i: int, j: int) -> int:
    """VCF diploid genotype ordering: index of (i,j) with i <= j."""
    if i > j:
        i, j = j, i
    return j * (j + 1) // 2 + i


class Record:
    __slots__ = ("chrom", "pos", "snarl", "n_alleles", "gl", "gt", "gqi", "line",
                 "hap_allele", "fields")

    def __init__(self, fields, chrom, pos, snarl, n_alleles, gl, gt, gqi):
        self.fields = fields
        self.chrom, self.pos, self.snarl = chrom, pos, snarl
        self.n_alleles, self.gl, self.gt, self.gqi = n_alleles, gl, gt, gqi
        self.hap_allele = None


def load_calls(vcf: Path, panel: dict, stats: collections.Counter):
    """Records in file order, each carrying its per-genotype likelihood and panel mapping."""
    out = []
    with open_maybe_gz(vcf) as fh:
        header = []
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            keys = f[8].split(":")
            vals = f[9].split(":")
            sample = dict(zip(keys, vals))
            gl_raw = sample.get("GL")
            if gl_raw is None or gl_raw == ".":
                stats["no GL"] += 1
                continue
            try:
                gl = [float(x) if x not in (".", "") else -np.inf for x in gl_raw.split(",")]
            except ValueError:
                stats["bad GL"] += 1
                continue
            n_alleles = 1 + len(f[4].split(","))
            if len(gl) != n_alleles * (n_alleles + 1) // 2:
                stats["GL/allele count mismatch"] += 1
                continue
            try:
                gqi = float(sample.get("GQI", sample.get("GQ", "nan")))
            except ValueError:
                gqi = float("nan")
            rec = Record(f, f[0], int(f[1]), f[2], n_alleles, np.array(gl, dtype=np.float64),
                         sample.get("GT", ""), gqi)

            # Map panel alleles onto this record's alleles through the traversal strings. Both
            # files are vg output over the same graph, so the strings are directly comparable;
            # anything that fails to line up is skipped and counted rather than guessed at.
            entry = panel.get(rec.snarl)
            at = info_field(f[7], "AT")
            if entry is not None and at is not None:
                call_at = at.split(",")
                index_of = {a: i for i, a in enumerate(call_at)}
                panel_at, assign = entry
                mapped = {}
                ok = True
                for hap, pa in assign.items():
                    if pa >= len(panel_at):
                        ok = False
                        break
                    ci = index_of.get(panel_at[pa])
                    if ci is None or ci >= n_alleles:
                        ok = False
                        break
                    mapped[hap] = ci
                if ok and mapped:
                    rec.hap_allele = mapped
                    stats["panel-linked"] += 1
                else:
                    stats["allele mapping failed"] += 1
            else:
                stats["no panel entry"] += 1
            out.append(rec)
    return out, header


# --------------------------------------------------------------------------- HMM

def emission_matrix(rec: Record, n_hap: int, escape: float) -> np.ndarray:
    """(n_hap+1) x (n_hap+1) matrix of relative P(reads | genotype implied by the state).

    Index n_hap is the wildcard, whose allele is unknown; its emission marginalises uniformly
    over the record's alleles. Without it, a genotype the panel cannot spell is unreachable, and
    the model would silently suppress novel alleles.
    """
    gl = rec.gl - np.max(rec.gl[np.isfinite(rec.gl)], initial=0.0)
    p = np.power(10.0, gl)              # relative, per genotype, in VCF ordering
    p[~np.isfinite(p)] = 0.0
    n = rec.n_alleles

    table = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            table[i, j] = p[gt_index(i, j)]

    alleles = np.full(n_hap + 1, -1, dtype=np.int64)
    if rec.hap_allele:
        for h, a in rec.hap_allele.items():
            alleles[h] = a

    e = np.zeros((n_hap + 1, n_hap + 1), dtype=np.float64)
    known = alleles >= 0
    if known.any():
        idx = alleles[known]
        sub = table[np.ix_(idx, idx)]
        e[np.ix_(known, known)] = sub
    # Wildcard rows/columns: marginalise the unknown strand uniformly over alleles.
    marg_rows = table.mean(axis=1)      # by the *known* strand's allele
    if known.any():
        e[known, n_hap] = marg_rows[alleles[known]] * escape
        e[n_hap, known] = marg_rows[alleles[known]] * escape
    e[n_hap, n_hap] = table.mean() * escape * escape
    # A haplotype that does not traverse this snarl carries no allele: treat it as the wildcard
    # rather than as evidence, which is what "absent" has to mean.
    absent = (~known)
    absent[n_hap] = False
    if absent.any() and known.any():
        e[np.ix_(absent, known)] = e[n_hap, known] * np.ones((absent.sum(), 1))
        e[np.ix_(known, absent)] = e[known, n_hap][:, None] * np.ones((1, absent.sum()))
    if absent.any():
        e[np.ix_(absent, absent)] = e[n_hap, n_hap]
    return e


def switch_prob(gap: int, scale: float, rho_min: float, beta: float, block: float) -> float:
    """Per-strand switch probability: distance term, then a boost at sampling block boundaries."""
    rho = rho_min + (1.0 - rho_min) * (1.0 - math.exp(-gap / scale))
    if beta > 0.0 and block > 0.0:
        # Expected number of block boundaries crossed. The positions are not available offline,
        # so the mass is spread rather than placed -- only the product of beta and block count is
        # identifiable from aggregate data anyway (see the note, section 4).
        crossings = gap / block
        rho = 1.0 - (1.0 - rho) * ((1.0 - beta) ** crossings)
    return min(max(rho, 0.0), 1.0)


def weighted_switch(gap, w_t, scale, rho_min, beta, block) -> float:
    """The switch probability, tempered by the transition weight.

    `rho ** w_t` has the endpoints we need: at `w_t = 0` it is 1, so every transition is uniform,
    the chain is memoryless and the posterior collapses to the per-site emission -- the current
    caller, recovered exactly. At `w_t = 1` it is the modelled rate. Above 1 linkage tightens.
    Tempering the *switch* probability rather than the log-transition keeps it a probability at
    every setting, which a weight on the log would not.
    """
    rho = switch_prob(max(gap, 1), scale, rho_min, beta, block)
    if w_t == 1.0:
        return rho
    return min(max(rho ** w_t, 1e-12), 1.0)


def forward_backward(recs, n_hap, w_t, scale, rho_min, beta, block, escape, w_freq=1.0):
    """Returns per-record genotype posteriors, as dicts {(i,j): probability}.

    `w_freq` separates two effects the state space silently bundles together, which the
    inertness check caught. Summing the state posterior over the states implying each genotype
    weights that genotype by its **multiplicity** -- how many haplotype pairs spell it -- which
    is a panel allele-frequency prior. So `w_t = 0` does *not* recover the current caller: it
    recovers the current caller plus that prior, and at `GQI < 10` it moves 40% of genotypes on
    its own, with no linkage involved.

    Dividing by multiplicity^(1 - w_freq) makes the two separable. `w_freq = 1` keeps the
    frequency prior, `w_freq = 0` removes it, and `(w_freq = 0, w_t = 0)` reproduces the per-site
    argmax exactly -- which is the inertness property the plan asked for and now has a test.

    Worth separating for its own sake: a frequency prior is a per-site quantity needing no HMM,
    no chains and no transitions, so if it turns out to do most of the work it is a far cheaper
    thing to ship.
    """
    m = n_hap + 1
    ems = [emission_matrix(r, n_hap, escape) for r in recs]

    def transition_apply(a, rho):
        """Sum over previous ordered pairs of a * T x T, exploiting T = (1-rho)I + rho/m."""
        stay = (1.0 - rho)
        jump = rho / m
        row = a.sum(axis=1, keepdims=True)
        col = a.sum(axis=0, keepdims=True)
        tot = a.sum()
        return (stay * stay * a
                + stay * jump * (row + col)
                + jump * jump * tot)

    n = len(recs)
    alphas = np.zeros((n, m, m))
    scales = np.zeros(n)
    a = np.full((m, m), 1.0 / (m * m)) * ems[0]
    s = a.sum() or 1.0
    a /= s
    alphas[0], scales[0] = a, s
    for t in range(1, n):
        rho = weighted_switch(recs[t].pos - recs[t - 1].pos, w_t, scale, rho_min, beta, block)
        a = transition_apply(alphas[t - 1], rho) * ems[t]
        s = a.sum() or 1.0
        a /= s
        alphas[t], scales[t] = a, s

    betas = np.zeros((n, m, m))
    betas[n - 1] = 1.0
    for t in range(n - 2, -1, -1):
        rho = weighted_switch(recs[t + 1].pos - recs[t].pos, w_t, scale, rho_min, beta, block)
        b = transition_apply(betas[t + 1] * ems[t + 1], rho)
        s = b.sum() or 1.0
        betas[t] = b / s

    out = []
    for t, rec in enumerate(recs):
        g = alphas[t] * betas[t]
        tot = g.sum()
        if tot <= 0 or not np.isfinite(tot):
            out.append(None)
            continue
        g /= tot
        alleles = np.full(m, -1, dtype=np.int64)
        if rec.hap_allele:
            for h, al in rec.hap_allele.items():
                alleles[h] = al
        post = collections.defaultdict(float)
        mult = collections.defaultdict(int)
        for i in range(m):
            ai = alleles[i]
            for j in range(m):
                aj = alleles[j]
                if ai < 0 or aj < 0:
                    continue          # wildcard/absent: no specific genotype implied
                key = (min(ai, aj), max(ai, aj))
                post[key] += g[i, j]
                mult[key] += 1
        if post and w_freq != 1.0:
            for key in post:
                if mult[key] > 0:
                    post[key] /= mult[key] ** (1.0 - w_freq)
        out.append(dict(post) if post else None)
    return out


# --------------------------------------------------------------------------- driver

def segments(recs, max_gap: int):
    """Split into chains on contig change or a gap wide enough that linkage is spent."""
    seg = []
    for r in recs:
        if seg and (r.chrom != seg[-1].chrom or r.pos - seg[-1].pos > max_gap):
            yield seg
            seg = []
        seg.append(r)
    if seg:
        yield seg


def chunks(seg, size: int, margin: int):
    """Overlapping chunks, yielding (records, lo, hi) where [lo,hi) is the interior to keep.

    Exact inference over a whole chain would need ~1 GB of betas at 100k sites, and the plan
    calls for a windowed implementation regardless (section 5): lift is ~1.50 by 10-30 kb, so a
    margin of a few hundred sites is far wider than the range over which linkage carries
    anything. The margin is discarded on both sides so no posterior is read from a position that
    could see the artificial chunk edge.
    """
    n = len(seg)
    if n <= size:
        yield seg, 0, n
        return
    start = 0
    while start < n:
        lo = max(start - margin, 0)
        hi = min(start + size + margin, n)
        yield seg[lo:hi], start - lo, min(start + size, n) - lo
        start += size


def _tally(interior, totals, changed, new_gt):
    for rec, p in interior:
        if not p:
            continue
        best = max(p.items(), key=lambda kv: kv[1])
        stratum = ("GQI<10" if rec.gqi < 10 else
                   "GQI 10-40" if rec.gqi < 40 else "GQI>=40")
        totals[stratum] += 1
        called = tuple(sorted(int(t) for t in rec.gt.replace("|", "/").split("/")
                              if t.isdigit()))
        if len(called) == 2 and called != best[0]:
            changed[stratum] += 1
            new_gt[(rec.chrom, rec.pos, rec.snarl)] = (best[0], best[1])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--dataset", help="score the rewritten VCF with score_vcf.py")
    ap.add_argument("--wt", nargs="+", type=float, default=[1.0],
                    help="transition weights to sweep; 0 must reproduce the input")
    ap.add_argument("--scale", type=float, default=10000.0)
    ap.add_argument("--beta", type=float, default=0.0)
    ap.add_argument("--block", type=float, default=10000.0)
    ap.add_argument("--rho-min", type=float, default=1e-3)
    ap.add_argument("--escape", type=float, default=1e-2)
    ap.add_argument("--freq-prior", type=float, default=1.0,
                    help="1 keeps the panel allele-frequency prior the state space implies, "
                         "0 removes it; (0, w_t=0) must reproduce the input exactly")
    ap.add_argument("--max-gap", type=int, default=200000)
    ap.add_argument("--chunk", type=int, default=2000, help="sites per exact-inference chunk")
    ap.add_argument("--margin", type=int, default=250, help="discarded sites either side")
    ap.add_argument("--threads", type=int, default=5)
    # This script needs numpy, which lives in the truvari venv here; score_vcf.py needs the
    # system interpreter and the repo's own deps. So they are deliberately separate.
    ap.add_argument("--scorer-python", default="python3")
    args = ap.parse_args()

    stats = collections.Counter()
    panel, n_hap = load_panel(Path(args.panel))
    recs, header = load_calls(Path(args.calls), panel, stats)
    print(f"=== {args.label}: {len(recs):,} records, {n_hap} panel haplotypes")
    for k, v in stats.most_common():
        print(f"    {k}: {v:,}")

    for w_t in args.wt:
        changed = collections.Counter()
        totals = collections.Counter()
        new_gt = {}
        for seg in segments(recs, args.max_gap):
            if len(seg) < 2:
                continue
            for block_recs, lo, hi in chunks(seg, args.chunk, args.margin):
                if len(block_recs) < 2:
                    continue
                post = forward_backward(block_recs, n_hap, w_t, args.scale, args.rho_min,
                                        args.beta, args.block, args.escape, args.freq_prior)
                interior = list(zip(block_recs, post))[lo:hi]
                _tally(interior, totals, changed, new_gt)

        print(f"\n  w_t = {w_t}")
        print(f"    {'stratum':>10s} {'records':>9s} {'changed':>9s} {'%':>7s}")
        for k in ("GQI<10", "GQI 10-40", "GQI>=40"):
            if totals[k]:
                print(f"    {k:>10s} {totals[k]:>9,} {changed[k]:>9,} "
                      f"{changed[k] / totals[k] * 100:>6.2f}%")
        harm = changed["GQI>=40"] / max(totals["GQI>=40"], 1) * 100
        print(f"    harm metric (GQI>=40 changed): {harm:.3f}%  "
              f"{'PASS' if harm <= 0.1 else 'over the 0.1% budget'}")

        out = WORK / "sv-atlas" / f"linkage-{args.label}-fp{args.freq_prior}-wt{w_t}.vcf"
        with open(out, "w") as fh:
            for h in header:
                fh.write(h)
            for rec in recs:
                key = (rec.chrom, rec.pos, rec.snarl)
                if key in new_gt:
                    (i, j), prob = new_gt[key]
                    f = list(rec.fields)
                    keys = f[8].split(":")
                    vals = f[9].split(":")
                    d = dict(zip(keys, vals))
                    d["GT"] = f"{i}/{j}"
                    q = 60.0 if prob >= 1.0 else -10.0 * math.log10(max(1.0 - prob, 1e-6))
                    if "GQ" in d:
                        d["GQ"] = str(int(min(256, max(0, round(q)))))
                    f[9] = ":".join(d[k] for k in keys)
                    fh.write("\t".join(f) + "\n")
                else:
                    fh.write("\t".join(rec.fields) + "\n")
        subprocess.run(["bgzip", "-f", str(out)], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(out) + ".gz"], check=True)
        print(f"    wrote {out}.gz")

        if args.dataset:
            tag = f"linkage-{args.label}-fp{args.freq_prior}-wt{w_t}"
            subprocess.run([args.scorer_python, str(HERE / "score_vcf.py"),
                            "--vcf", str(out) + ".gz", "--label", tag,
                            "--dataset", args.dataset, "--threads", str(args.threads)],
                           check=True, capture_output=True)
            s = json.loads((WORK / "sv-atlas" / f"score-{tag}.json").read_text())
            sv = s.get("sv") or {}
            gt = next((float(r["metric_f1"]) for r in (s.get("smallvar") or [])
                       if (r["comparison"], r["region_label"], r["filter"],
                           r["variant_type"]) == ("GT", "ALL", "ALL", "ALL")), 0.0)
            print(f"    scored: small-variant GT F1 {gt:.4f}, SV F1 {sv.get('f1', 0):.4f}")


if __name__ == "__main__":
    main()
