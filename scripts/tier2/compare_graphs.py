#!/usr/bin/env python3
"""Put a contig's 4-haplotype and 34-haplotype runs side by side.

Both runs use the same reads sample, the same truth slices, the same confident regions and
the same reference sequence -- prep_contig.sh extracts the FASTA from *each* graph and
refuses to proceed unless the two are byte-identical, so a difference here can never be a
coordinate mismatch. What differs is the graph and, unavoidably, the alignments: reads
mapped to one graph cannot be scored against the other because the node ID spaces differ.
A row is therefore "richer graph, and reads remapped to it" -- which is what adopting such
a graph actually involves, but it is not a single-variable experiment and the page says so.

The interesting column is `-z`. Those arms enumerate alleles from the GBWT haplotypes, so
going from 4 to 34 haplotypes changes what alleles are *available to call* rather than how
they are scored.

**Structural variants come from truvari, not aardvark.** aardvark is scored against the
small-variant benchmark, which contains no record at all above 50 bp, so its Sv* categories
have almost nothing to score against and its summary leaves the query columns at zero
besides. Plan §9.22 established truvari as the SV metric; the aardvark SV block is still
emitted when the artefacts exist, clearly marked as secondary.
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

ARM_ORDER = ["poisson", "poisson-z", "readlik", "readlik-nomismap",
             "readlik-z-nolink", "readlik-z"]
SMALL_TYPES = [("ALL", "ALL"), ("Snv", "SNV"), ("Insertion", "Insertion (<50 bp)"),
               ("Deletion", "Deletion (<50 bp)")]
SV_TYPES = [("SvInsertion", "SV insertion"), ("SvDeletion", "SV deletion"),
            ("JointStructuralVariant", "SV (joint)")]


def load_arms(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    return {e["arm"]: e for e in json.loads(path.read_text())}


def load_truvari(results: Path) -> dict[str, dict]:
    """Per-arm truvari summaries, the primary SV metric."""
    out = {}
    for arm in ARM_ORDER:
        p = results / f"truvari-{arm}" / "summary.json"
        if p.exists():
            out[arm] = json.loads(p.read_text())
    return out


def sv_precision(results: Path, arm: str) -> float | None:
    """Recompute SV precision from aardvark's per-variant decisions.

    Its summary leaves query_total/query_tp/query_fp at zero for the Sv* categories, so the
    published precision and F1 come out as 0/0. Recall is fine and is used as published.
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
    return f"{b - a:+.4f}"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--contig", default="chr20")
    p.add_argument("--old")
    p.add_argument("--new")
    p.add_argument("--out")
    args = p.parse_args()

    c = args.contig
    old_res = Path(args.old or REPO / f"work/tier2-{c}/results")
    new_res = Path(args.new or REPO / f"work/tier2-{c}-hap32/results")
    out_path = Path(args.out or REPO / f"docs/tier2-{c}-hap32.md")

    old = load_arms(old_res / "arms.json")
    new = load_arms(new_res / "arms.json")
    old_sv_rows, new_sv_rows = load_sv(old_res), load_sv(new_res)
    old_tv, new_tv = load_truvari(old_res), load_truvari(new_res)
    old_sm = load_arms(old_res / "arms-size-matched.json")
    new_sm = load_arms(new_res / "arms-size-matched.json")

    if not new:
        raise SystemExit(f"no arms.json under {new_res}; run the arms first")

    L: list[str] = []
    L.append(f"# {c}: 4-haplotype vs 34-haplotype graph")
    L.append("")
    L.append("Same sample, same reads, same truth, same confident regions, same reference "
             "sequence. What changes is the graph — and, unavoidably, the alignments.")
    L.append("")
    L.append("| | 4-haplotype | 34-haplotype |")
    L.append("|---|---|---|")
    L.append("| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz` | `…HG002.hap32.gbz` |")
    L.append("| haplotypes | 4 (CHM13, GRCh38, 2 recombinants) | **34** (CHM13, GRCh38, "
             "**32 recombinants** — the file is named `hap32` after the recombinant count, "
             "not the total) |")
    L.append("| HG002 present? | no | **no** — samples are `CHM13`, `GRCh38`, `recombination` |")
    L.append("| alignments | `…HG002.gaf.gz` | `…HG002.hap32.gaf.gz` (remapped) |")
    L.append("")
    L.append("**This is not a single-variable experiment.** Reads mapped to one graph cannot be "
             "scored against the other, because the node ID spaces differ — so the 34-haplotype "
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

    L.append("## What this says")
    L.append("")
    L.append("**The read-likelihood caller is better on the richer graph; the Poisson caller is "
             "much worse on it.** That split is the result. More haplotypes offer more true alleles "
             "and more wrong ones, and what decides the outcome is whether the genotyper can tell "
             "them apart read by read.")
    L.append("")
    L.append("| arm | 4-hap GT F1 | 34-hap GT F1 | Δ |")
    L.append("|---|---|---|---|")
    for a in ("poisson-z", "readlik-z"):
        o, n = gtf1(old, a), gtf1(new, a)
        L.append(f"| `{a}` | {fmt(o)} | {fmt(n)} | **{delta(o, n)}** |")
    L.append("")
    gap_old = (gtf1(old, "readlik-z") or 0) - (gtf1(old, "poisson-z") or 0)
    gap_new = (gtf1(new, "readlik-z") or 0) - (gtf1(new, "poisson-z") or 0)
    L.append(f"The read-likelihood caller's margin over the Poisson caller goes from "
             f"**{gap_old:+.4f}** on the 4-haplotype graph to **{gap_new:+.4f}** on the "
             f"34-haplotype one"
             + (f" — {gap_new/gap_old:.1f}x wider." if gap_old else "."))
    L.append("")
    L.append("**Two directions, and they are not the same direction.** GT F1 rises on the richer "
             "graph for the read-likelihood caller; BASEPAIR and SV F1 fall for both callers. The "
             "SV fall is **entirely precision** — recall is flat on chr6 and slightly better on "
             "chr20 — and most of it is not the caller getting worse. Two thirds to all of it is "
             "records that are not structural variants plus the cost of scoring unfiltered; at "
             "matched sensitivity the residual is 0.021 on chr6 and zero on chr20. "
             "[tier2-sv-errors.md](tier2-sv-errors.md) has the decomposition.")
    L.append("")
    L.append("Exposure to multi-allelic sites was the earlier explanation and it does not survive "
             "measurement: precision falls within the biallelic stratum, which is 78-82% of "
             "records, by nearly the whole amount. Multi-allelic records do grow (17.6% to 22.1% "
             "of SV-sized records) and are harder, but they are a minor term rather than the "
             "mechanism.")
    L.append("")
    L.append("**This depended on a default that was wrong for graphs like this.** With "
             "`--mismap-max` at its old 0.1, `readlik-z` on the 34-haplotype graph looked like a "
             "precision-for-recall trade" +
             (" — 1,597 false-positive SNVs against the 4-haplotype graph's 375" if c == "chr20"
              else " (measured on chr20: 1,597 false-positive SNVs against 375)") +
             ". The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, "
             "meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that "
             "excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is "
             "that a caller-level default, not the graph, was the difference between the two "
             "readings.")
    L.append("")
    nm = new.get("readlik-nomismap", {}).get("metrics", {}).get("summary")
    nm_fp = pick(nm, "GT", "Snv")
    if nm_fp and nm_fp.get("query_fp"):
        L.append("**`readlik-nomismap` is the control.** It disables the mismapping term entirely, "
                 "so the cap cannot reach it — and on the richer graph it still carries "
                 f"{int(nm_fp['query_fp']):,} spurious SNVs. The term is what does the work.")
        L.append("")
    o_sm, n_sm = smf1("sm50-readlik-z")
    if o_sm and n_sm:
        L.append(f"Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — "
                 f"`readlik-z` goes {fmt(o_sm)} to {fmt(n_sm)}.")
        L.append("")
    L.append("**One caveat this data cannot settle.** Some of the remaining false positives may not "
             "be error: a graph carrying 32 haplotypes will call real variation a draft benchmark "
             "does not cover, and that scores as a false positive. Separating them needs a more "
             "complete truth set, not a different metric. It has since been *bounded* rather than "
             "settled: false calls made by both callers on both graphs with no truth candidate "
             "anywhere nearby number 44 on chr6 and 40 on chr20, which is a lower bound on the "
             "benchmark's share of them.")
    L.append("")

    L.append("## Cost")
    L.append("")
    L.append("| arm | 4-hap wall | 34-hap wall | 4-hap RSS | 34-hap RSS | 4-hap variants | "
             "34-hap variants |")
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
        L.append("| arm | class | 4-hap | 34-hap | Δ |")
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

    if new_tv:
        L.append("## Structural variants — truvari (GIAB `stvar`)")
        L.append("")
        L.append("The SV metric. Reciprocal-overlap matching against the structural benchmark, "
                 "`--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored "
                 "against the *small-variant* truth set and therefore have essentially no truth "
                 "to match above 50 bp (plan §9.22).")
        L.append("")
        L.append("**These errors are broken down per record in "
                 "[tier2-sv-errors.md](tier2-sv-errors.md)**, and three findings there bear "
                 "directly on this table. The 34-haplotype false-positive rise is not the same "
                 "errors plus more — only about two thirds of the 4-haplotype false calls survive "
                 "the graph change, and the new ones are disproportionately calls with no truth "
                 "candidate at all. A quarter of all false positives are placement or bookkeeping "
                 "artefacts of the metric. And harmonising representation with `truvari refine` "
                 "lifts every arm by roughly 0.05 F1.")
        L.append("")
        L.append("| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | "
                 "34-hap F1 | **Δ F1** |")
        L.append("|---|---|---|---|---|---|---|---|")
        for a in ARM_ORDER:
            if a not in new_tv:
                continue
            o, n = old_tv.get(a), new_tv[a]
            g = lambda d, k: (None if d is None else d.get(k))  # noqa: E731
            L.append(f"| `{a}` | {fmt(g(o, 'recall'))} | {fmt(g(n, 'recall'))} | "
                     f"{fmt(g(o, 'precision'))} | {fmt(g(n, 'precision'))} | "
                     f"{fmt(g(o, 'f1'))} | {fmt(g(n, 'f1'))} | "
                     f"**{delta(g(o, 'f1'), g(n, 'f1'))}** |")
        L.append("")

    if new_sv_rows:
        L.append("## Structural variants — aardvark (secondary)")
        L.append("")
        L.append("Kept for continuity with earlier runs. Recall is aardvark's published value; "
                 "**precision is recomputed** from its per-variant `BD` decisions, because its "
                 "summary leaves the query columns at zero for the `Sv*` categories. Prefer the "
                 "truvari table above: these categories are scored against a truth set with no "
                 "record over 50 bp.")
        L.append("")
        L.append("| arm | 4-hap recall | 34-hap recall | Δ | 4-hap prec | 34-hap prec | Δ | "
                 "4-hap F1 | 34-hap F1 | **Δ F1** |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for a in ARM_ORDER:
            if a not in new_sv_rows:
                continue
            orr = fnum(pick(old_sv_rows.get(a), "GT", "JointStructuralVariant"), "metric_recall")
            nrr = fnum(pick(new_sv_rows[a], "GT", "JointStructuralVariant"), "metric_recall")
            opp, npp = sv_precision(old_res, a), sv_precision(new_res, a)
            of, nf = harmonic(orr, opp), harmonic(nrr, npp)
            L.append(f"| `{a}` | {fmt(orr)} | {fmt(nrr)} | {delta(orr, nrr)} | {fmt(opp)} | "
                     f"{fmt(npp)} | {delta(opp, npp)} | {fmt(of)} | {fmt(nf)} | "
                     f"**{delta(of, nf)}** |")
        L.append("")

    if new_sm:
        L.append("## Small variants restricted to <50 bp — BASEPAIR")
        L.append("")
        L.append("The `smvar` truth set holds no record >=50 bp, so a large insertion called inside "
                 "its confident region scores FP on every base however right it is. A richer graph "
                 "calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any "
                 "record with a called allele >=50 bp from *both* sides is the only like-for-like "
                 "read of these numbers.")
        L.append("")
        L.append("| arm | class | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | "
                 "4-hap F1 | 34-hap F1 | **Δ F1** |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for a in ("sm50-poisson-z", "sm50-readlik-z"):
            if a not in new_sm:
                continue
            for vtype, label in (("Insertion", "Insertion"), ("Deletion", "Deletion"),
                                 ("ALL", "ALL")):
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

    L.append("## Quality fields")
    L.append("")
    L.append("Every arm above is scored at **every** GQ, so nothing on this page depends on the "
             "quality field. It matters for how the calls rank, which is a separate page: see "
             "[tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits "
             "`AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the "
             "called genotype explains. That rescales a quality and does not change a genotype, so "
             "**the unfiltered numbers on this page are unaffected by it**.")
    L.append("")

    out_path.write_text("\n".join(L) + "\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
