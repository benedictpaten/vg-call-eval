#!/usr/bin/env python3
"""Put the 4-haplotype and 32-haplotype chr20 runs side by side.

Both runs use the same reads sample, the same truth slices, the same confident regions
and the same reference sequence (prep_hap32_chr20.sh refuses to proceed unless the two
graphs' CHM13 paths are byte-identical). What differs is the graph and, unavoidably, the
alignments: reads mapped to one graph cannot be scored against the other because the node
ID spaces differ. So a row here is "richer graph, and reads remapped to it" -- which is
what adopting such a graph actually involves, but it is not a single-variable experiment
and the tables say so.

The interesting column is `-z`. Those arms enumerate alleles from the GBWT haplotypes, so
going from 4 to 34 haplotypes changes what alleles are *available to call* rather than how
they are scored. The design doc's tier-2 finding was that enumeration matters more than the
genotyper, especially for SVs; this is the direct test of that.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

ARM_ORDER = ["poisson", "poisson-z", "readlik", "readlik-nomismap", "readlik-z"]
SMALL_TYPES = [("ALL", "ALL"), ("Snv", "SNV"), ("Insertion", "Insertion (<50 bp)"),
               ("Deletion", "Deletion (<50 bp)")]
SV_TYPES = [("SvInsertion", "SV insertion"), ("SvDeletion", "SV deletion"),
            ("JointStructuralVariant", "SV (joint)")]


def load_arms(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {e["arm"]: e for e in json.loads(path.read_text())}


def sv_precision(results: Path, arm: str) -> float | None:
    """Recompute SV precision from aardvark's per-variant decisions.

    Its summary leaves query_total/query_tp/query_fp at zero for the Sv* categories, so
    the published precision and F1 come out as 0/0. Recall is fine and is used as
    published. Without this, a run that called far more SVs would look like a pure
    recall win when it had in fact traded precision away.
    """
    p = results / f"aardvark-sv-{arm}" / "query.vcf.gz"
    if not p.exists():
        return None
    counts: Counter = Counter()
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            alt = f[4].split(",")[0]
            if alt.startswith("<") or alt == "*":
                continue
            if abs(len(alt) - len(f[3])) < 50:
                continue
            fmt, val = f[8].split(":"), f[9].split(":")
            counts[dict(zip(fmt, val)).get("BD", "?")] += 1
    tp, fp = counts.get("TP", 0), counts.get("FP", 0)
    return tp / (tp + fp) if (tp + fp) else None


def harmonic(r, p):
    return 2 * r * p / (r + p) if (r and p) else None


def load_sv(results: Path) -> dict[str, list[dict]]:
    out = {}
    for arm in ARM_ORDER:
        p = results / f"aardvark-sv-{arm}" / "summary.tsv"
        if p.exists():
            with open(p) as fh:
                out[arm] = list(csv.DictReader(fh, delimiter="\t"))
    return out


def pick(rows, comparison, vtype):
    for r in rows or []:
        if r.get("comparison", "").upper() == comparison and r.get("variant_type") == vtype:
            return r
    return None


def fnum(r, key):
    if not r or r.get(key) in ("", None):
        return None
    try:
        return float(r[key])
    except ValueError:
        return None


def fmt(x, nd=4):
    return "—" if x is None else f"{x:.{nd}f}"


def delta(a, b):
    """b - a, rendered with a sign, or an em dash if either side is missing."""
    if a is None or b is None:
        return "—"
    d = b - a
    return f"{d:+.4f}"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--old", default=str(REPO / "work/tier2-chr20/results"))
    p.add_argument("--new", default=str(REPO / "work/tier2-chr20-hap32/results"))
    p.add_argument("--out", default=str(REPO / "docs/tier2-chr20-hap32.md"))
    args = p.parse_args()

    old_res, new_res = Path(args.old), Path(args.new)
    old = load_arms(old_res / "arms.json")
    new = load_arms(new_res / "arms.json")
    old_sv_rows, new_sv_rows = load_sv(old_res), load_sv(new_res)

    if not new:
        raise SystemExit(f"no arms.json under {new_res}; run run_hap32_chr20.sh first")

    L: list[str] = []
    L.append("# chr20: 4-haplotype vs 32-haplotype graph")
    L.append("")
    L.append("Same sample, same reads, same truth, same confident regions, same reference "
             "sequence. What changes is the graph — and, unavoidably, the alignments.")
    L.append("")
    L.append("| | 4-haplotype | 32-haplotype |")
    L.append("|---|---|---|")
    L.append("| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz` | `…HG002.hap32.gbz` |")
    L.append("| haplotypes | 4 (CHM13, GRCh38, 2 recombinants) | **34** (CHM13, GRCh38, "
             "32 recombinants) |")
    L.append("| HG002 present? | no | **no** — samples are `CHM13`, `GRCh38`, `recombination` |")
    L.append("| alignments | `…HG002.gaf.gz` | `…HG002.hap32.gaf.gz` (remapped) |")
    L.append("")
    L.append("**This is not a single-variable experiment.** Reads mapped to one graph cannot be "
             "scored against the other, because the node ID spaces differ — so the 32-haplotype "
             "arm necessarily uses its own alignments. Graph and alignment move together. That is "
             "what adopting a richer graph actually involves, but it means a difference below "
             "cannot be attributed to the graph alone.")
    L.append("")
    L.append("The rows to watch are the **`-z` arms**, which enumerate alleles from the GBWT "
             "haplotypes. Going from 4 to 34 changes which alleles are *available to call* rather "
             "than how they are scored, and the tier-2 finding was that enumeration matters more "
             "than the genotyper — most of all for SVs. This is the direct test.")
    L.append("")

    # Computed, never transcribed. An earlier version of this section hardcoded the
    # deltas and they went stale the moment the --mismap-max default changed.
    def gtf1(d, arm, vtype="ALL"):
        return fnum(pick(d.get(arm, {}).get("metrics", {}).get("summary"), "GT", vtype),
                    "metric_f1")

    def smf1(arm, vtype="ALL"):
        o = fnum(pick(old_sm.get(arm, {}).get("metrics", {}).get("summary"),
                      "BASEPAIR", vtype), "metric_f1")
        n = fnum(pick(new_sm.get(arm, {}).get("metrics", {}).get("summary"),
                      "BASEPAIR", vtype), "metric_f1")
        return o, n

    old_sm = load_arms(old_res / "arms-size-matched.json")
    new_sm = load_arms(new_res / "arms-size-matched.json")

    L.append("## What this says")
    L.append("")
    L.append("**The read-likelihood caller is better on the richer graph; the Poisson caller is "
             "much worse on it.** That split is the result. More haplotypes offer more true alleles "
             "and more wrong ones, and what decides the outcome is whether the genotyper can tell "
             "them apart read by read.")
    L.append("")
    L.append("| arm | 4-hap GT F1 | 32-hap GT F1 | Δ |")
    L.append("|---|---|---|---|")
    for a in ("poisson-z", "readlik-z"):
        o, n = gtf1(old, a), gtf1(new, a)
        L.append(f"| `{a}` | {fmt(o)} | {fmt(n)} | **{delta(o, n)}** |")
    L.append("")
    gap_old = (gtf1(old, "readlik-z") or 0) - (gtf1(old, "poisson-z") or 0)
    gap_new = (gtf1(new, "readlik-z") or 0) - (gtf1(new, "poisson-z") or 0)
    L.append(f"The read-likelihood caller's margin over the Poisson caller goes from "
             f"**{gap_old:+.4f}** on the 4-haplotype graph to **{gap_new:+.4f}** on the "
             f"32-haplotype one"
             + (f" — {gap_new/gap_old:.1f}x wider." if gap_old else "."))
    L.append("")
    L.append("**This depended on a default that was wrong for graphs like this.** With "
             "`--mismap-max` at its old 0.1, `readlik-z` on the 32-haplotype graph carried 1,597 "
             "false-positive SNVs against the 4-haplotype graph's 375, and looked like a "
             "precision-for-recall trade. The cap was overriding the mapper: at those sites 23.3% "
             "of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the "
             "current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the "
             "derivation; the point for this page is that a caller-level default, not the graph, "
             "was the difference between the two readings.")
    L.append("")
    L.append("**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so "
             "the cap cannot reach it — and on the richer graph it still carries "
             f"{int(pick(new.get('readlik-nomismap', {}).get('metrics', {}).get('summary'), 'GT', 'Snv')['query_fp']):,} "
             "spurious SNVs. The term is what does the work.")
    L.append("")
    o_sm, n_sm = smf1("sm50-readlik-z")
    if o_sm and n_sm:
        L.append(f"Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — "
                 f"`readlik-z` goes {fmt(o_sm)} to {fmt(n_sm)}.")
        L.append("")
    L.append("**One caveat this data cannot settle.** Some of the remaining false positives may not "
             "be error: a graph carrying 32 haplotypes will call real variation a draft benchmark "
             "does not cover, and that scores as a false positive. Separating them needs a more "
             "complete truth set, not a different metric.")
    L.append("")

    L.append("## Cost")
    L.append("")
    L.append("| arm | 4-hap wall | 32-hap wall | 4-hap RSS | 32-hap RSS | 4-hap variants | "
             "32-hap variants |")
    L.append("|---|---|---|---|---|---|---|")
    for a in ARM_ORDER:
        o, n = old.get(a), new.get(a)
        if not n:
            continue
        L.append(f"| `{a}` | {o['seconds']:.0f} s | {n['seconds']:.0f} s | "
                 f"{o['peak_rss_gb']:.1f} GB | {n['peak_rss_gb']:.1f} GB | "
                 f"{o['variants']:,} | {n['variants']:,} |" if o else
                 f"| `{a}` | — | {n['seconds']:.0f} s | — | {n['peak_rss_gb']:.1f} GB | — | "
                 f"{n['variants']:,} |")
    L.append("")

    for comparison in ("GT", "BASEPAIR"):
        L.append(f"## Small variants — {comparison} F1")
        L.append("")
        L.append("| arm | class | 4-hap | 32-hap | Δ |")
        L.append("|---|---|---|---|---|")
        for a in ARM_ORDER:
            if a not in new:
                continue
            for vtype, label in SMALL_TYPES:
                ov = fnum(pick(old.get(a, {}).get("metrics", {}).get("summary"),
                               comparison, vtype), "metric_f1")
                nv = fnum(pick(new[a]["metrics"]["summary"], comparison, vtype), "metric_f1")
                L.append(f"| `{a}` | {label} | {fmt(ov)} | {fmt(nv)} | {delta(ov, nv)} |")
        L.append("")

    L.append("## Structural variants (GIAB `stvar`)")
    L.append("")
    L.append("Recall is aardvark's published value. **Precision is recomputed** from its per-variant "
             "`BD` decisions, because its summary leaves the query columns at zero for the `Sv*` "
             "categories — without that, a run calling far more SVs would read as a pure recall win "
             "when it had traded precision away. F1 is derived from the two.")
    L.append("")
    L.append("| arm | 4-hap recall | 32-hap recall | Δ | 4-hap prec | 32-hap prec | Δ | "
             "4-hap F1 | 32-hap F1 | **Δ F1** |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for a in ARM_ORDER:
        if a not in new_sv_rows:
            continue
        orr = fnum(pick(old_sv_rows.get(a), "GT", "JointStructuralVariant"), "metric_recall")
        nrr = fnum(pick(new_sv_rows[a], "GT", "JointStructuralVariant"), "metric_recall")
        opp, npp = sv_precision(old_res, a), sv_precision(new_res, a)
        of, nf = harmonic(orr, opp), harmonic(nrr, npp)
        L.append(f"| `{a}` | {fmt(orr)} | {fmt(nrr)} | {delta(orr, nrr)} | {fmt(opp)} | {fmt(npp)} | "
                 f"{delta(opp, npp)} | {fmt(of)} | {fmt(nf)} | **{delta(of, nf)}** |")
    L.append("")
    L.append("Per class, recall only:")
    L.append("")
    L.append("| arm | class | 4-hap | 32-hap | Δ |")
    L.append("|---|---|---|---|---|")
    for a in ARM_ORDER:
        if a not in new_sv_rows:
            continue
        for vtype, label in SV_TYPES[:2]:
            ov = fnum(pick(old_sv_rows.get(a), "GT", vtype), "metric_recall")
            nv = fnum(pick(new_sv_rows[a], "GT", vtype), "metric_recall")
            L.append(f"| `{a}` | {label} | {fmt(ov)} | {fmt(nv)} | {delta(ov, nv)} |")
    L.append("")

    # Size-matched: the only apples-to-apples read of the BASEPAIR insertion numbers.
    if new_sm:
        L.append("## Small variants restricted to <50 bp — BASEPAIR")
        L.append("")
        L.append("The `smvar` truth set holds no record >=50 bp, so a large insertion called inside "
                 "its confident region scores FP on every base however right it is. A richer graph "
                 "calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any "
                 "record with a called allele >=50 bp from *both* sides is the only like-for-like "
                 "read of these numbers.")
        L.append("")
        L.append("| arm | class | 4-hap recall | 32-hap recall | 4-hap prec | 32-hap prec | "
                 "4-hap F1 | 32-hap F1 | **Δ F1** |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for a in ("sm50-poisson-z", "sm50-readlik-z"):
            if a not in new_sm:
                continue
            for vtype, label in (("Insertion", "Insertion"), ("Deletion", "Deletion"), ("ALL", "ALL")):
                o_rows = old_sm.get(a, {}).get("metrics", {}).get("summary")
                n_rows = new_sm[a]["metrics"]["summary"]
                orr = fnum(pick(o_rows, "BASEPAIR", vtype), "metric_recall")
                nrr = fnum(pick(n_rows, "BASEPAIR", vtype), "metric_recall")
                opp = fnum(pick(o_rows, "BASEPAIR", vtype), "metric_precision")
                npp = fnum(pick(n_rows, "BASEPAIR", vtype), "metric_precision")
                of = fnum(pick(o_rows, "BASEPAIR", vtype), "metric_f1")
                nf = fnum(pick(n_rows, "BASEPAIR", vtype), "metric_f1")
                L.append(f"| `{a}` | {label} | {fmt(orr)} | {fmt(nrr)} | {fmt(opp)} | {fmt(npp)} | "
                         f"{fmt(of)} | {fmt(nf)} | **{delta(of, nf)}** |")
        L.append("")

    Path(args.out).write_text("\n".join(L) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
