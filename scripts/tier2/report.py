#!/usr/bin/env python3
"""Render the full tier-2 chr20 results: SNVs, small indels, and structural variants.

Reads the aardvark output directories written by run_arms.py (small-variant benchmark)
and compare_sv.py (structural-variant benchmark). Kept separate from both so the tables
can be regenerated without re-running any calling or comparison.

One number here is computed rather than read. Aardvark's summary.tsv reports
`truth_total` and `truth_tp` for the SvInsertion / SvDeletion / JointStructuralVariant
categories, but leaves `query_total`, `query_tp` and `query_fp` at zero -- so its own
precision and F1 columns come out as 0/0 and are unusable for those rows. The per-variant
decisions *are* present in aardvark's annotated query VCF, so SV precision is recomputed
here by counting BD=TP against BD=FP over query variants of >=50 bp. Recall is taken from
the summary as published.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

ARM_ORDER = ["poisson", "poisson-z", "readlik", "readlik-nomismap", "readlik-z"]

SMALL_TYPES = [("Snv", "SNV"), ("Insertion", "Insertion (<50 bp)"),
               ("Deletion", "Deletion (<50 bp)"), ("Indel", "Indel"),
               ("JointIndel", "Indel (joint)"), ("ALL", "ALL")]
SV_TYPES = [("SvInsertion", "SV insertion (>=50 bp)"),
            ("SvDeletion", "SV deletion (>=50 bp)"),
            ("JointStructuralVariant", "SV (joint)")]

META = {
    "poisson": ("support (Flow)", "yes"),
    "poisson-z": ("haplotype (`-z`)", "yes"),
    "readlik": ("support (Flow)", "yes"),
    "readlik-nomismap": ("support (Flow)", "yes"),
    "readlik-z": ("haplotype (`-z`)", "**no**"),
}


def load_merged(res: Path, pattern: str) -> dict[str, dict]:
    by_name: dict[str, dict] = {}
    for f in sorted(res.glob(pattern), key=lambda f: f.stat().st_mtime):
        for entry in json.loads(f.read_text()):
            by_name[entry["arm"]] = entry
    return by_name


def pick(rows: list[dict], comparison: str, vtype: str) -> dict | None:
    for r in rows:
        if r.get("comparison", "").upper() == comparison and r.get("variant_type") == vtype:
            return r
    return None


def sv_query_decisions(query_vcf: Path) -> Counter:
    """Count aardvark's BD decisions over query variants of >=50 bp."""
    counts: Counter = Counter()
    if not query_vcf.exists():
        return counts
    with gzip.open(query_vcf, "rt") as fh:
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
    return counts


def f(x, nd: int = 4) -> str:
    return "—" if x is None else f"{x:.{nd}f}"


def num(r, key: str):
    if not r or r.get(key) in ("", None):
        return None
    try:
        return float(r[key])
    except ValueError:
        return None


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default=str(REPO / "work/tier2-chr20/results"))
    p.add_argument("--out", default=str(REPO / "docs/tier2-chr20-results.md"))
    args = p.parse_args()
    res = Path(args.results)

    # arms*.json would also match arms-sv.json, whose entries would then overwrite the
    # small-variant ones (same arm names, newer mtime). Load the small-variant batches
    # explicitly instead.
    small = load_merged(res, "arms.json")
    small.update(load_merged(res, "arms.readlik-z.json"))
    sv = load_merged(res, "arms-sv.json")

    L: list[str] = []
    L.append("# Tier 2 results: HG002 chr20 on HPRC v2.1 MC CHM13")
    L.append("")
    L.append("Real reads, real benchmark, run on a 32 GB laptop.")
    L.append("")
    L.append("| | |")
    L.append("|---|---|")
    L.append("| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz`, 100,179,277 nodes, **4 haplotypes** "
             "(CHM13, GRCh38, 2 recombinants). HG002 itself is **absent** — no circularity |")
    L.append("| chromosome | chr20 component, 2,382,533 nodes, IDs 114,818,865–121,250,404 |")
    L.append("| reads | 596,017,764 alignments genome-wide (~28.6×), 13,279,246 on chr20; "
             "151 bp paired Illumina |")
    L.append("| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |")
    L.append("| regions | small variants 58.9 Mb (88.9% of chr20); SVs 59.4 Mb (89.6%) |")
    L.append("| engine | `aardvark compare`; SV runs use `--min-variant-gap 1000` + record-basepair |")
    L.append("")
    L.append("**Read the caveats before the numbers.** The benchmark is a *draft*: its own README "
             "reports known errors in highly homozygous regions, homopolymers and tandem repeats, and "
             "excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm "
             "comparison is what this table is for.")
    L.append("")

    L.append("## Cost")
    L.append("")
    L.append("| arm | enumeration | pack? | variants | wall | peak RSS |")
    L.append("|---|---|---|---|---|---|")
    for a in ARM_ORDER:
        if a not in small:
            continue
        e = small[a]
        enum, pack = META[a]
        L.append(f"| `{a}` | {enum} | {pack} | {e['variants']:,} | {e['seconds']:.0f} s | "
                 f"{e['peak_rss_gb']:.1f} GB |")
    L.append("")

    L.append("## Small variants (GIAB `smvar` benchmark)")
    L.append("")
    L.append("`GT` is the genotype-aware comparison — the one that matters for a genotyper. "
             "`BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the "
             "wrong sequence. Bold marks the best GT F1 in each class.")
    L.append("")
    for vtype, label in SMALL_TYPES:
        rows = [(a, pick(small[a]["metrics"]["summary"], "GT", vtype),
                 pick(small[a]["metrics"]["summary"], "BASEPAIR", vtype))
                for a in ARM_ORDER if a in small]
        if not any(gt for _, gt, _ in rows):
            continue
        L.append(f"### {label}")
        L.append("")
        L.append("| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | "
                 "BP precision | BP F1 |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        best = max((num(gt, "metric_f1") or 0) for _, gt, _ in rows)
        for a, gt, bp in rows:
            if not gt:
                continue
            f1 = num(gt, "metric_f1")
            mark = "**" if f1 and abs(f1 - best) < 1e-9 else ""
            L.append(
                f"| `{a}` | {f(num(gt,'metric_recall'))} | {f(num(gt,'metric_precision'))} | "
                f"{mark}{f(f1)}{mark} | {int(gt['truth_tp']):,} | {int(gt['truth_fn']):,} | "
                f"{int(gt['query_fp']):,} | {f(num(bp,'metric_recall'))} | "
                f"{f(num(bp,'metric_precision'))} | {f(num(bp,'metric_f1'))} |")
        L.append("")

    L.append("## Structural variants (GIAB `stvar` benchmark)")
    L.append("")
    L.append("Of 176,623 chr20 truth records only **2,052 are >=50 bp** — the rest is the local "
             "sequence context an SV-aware haplotype comparison needs to place the SV. The rows below "
             "are the SV-specific categories, not the whole benchmark.")
    L.append("")
    L.append("**Precision here is recomputed, not read from aardvark.** Its summary leaves "
             "`query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own "
             "precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated "
             "query VCF, so precision is counted from those over query variants of >=50 bp; recall is "
             "the published summary value; F1 is derived from the two.")
    L.append("")
    for vtype, label in SV_TYPES:
        rows = [(a, pick(sv[a]["metrics"]["summary"], "GT", vtype)) for a in ARM_ORDER if a in sv]
        if not any(r for _, r in rows):
            continue
        L.append(f"### {label}")
        L.append("")
        L.append("| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\\* | F1\\* |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for a, r in rows:
            if not r:
                continue
            dec = sv_query_decisions(res / f"aardvark-sv-{a}" / "query.vcf.gz")
            tp, fp = dec.get("TP", 0), dec.get("FP", 0)
            prec = tp / (tp + fp) if tp + fp else None
            rec = num(r, "metric_recall")
            f1 = (2 * rec * prec / (rec + prec)) if (rec and prec) else None
            L.append(f"| `{a}` | {f(rec)} | {int(r['truth_tp']):,} | {int(r['truth_fn']):,} | "
                     f"{tp+fp:,} | {tp:,} | {fp:,} | {f(prec)} | {f(f1)} |")
        L.append("")
    L.append("\\* recomputed as described above. The per-variant counts are shared across the three "
             "SV rows because they are counted over all >=50 bp query variants, not split by "
             "insertion/deletion; only recall is category-specific.")
    L.append("")

    L.append("## Raw aardvark summary rows")
    L.append("")
    for title, src in [("small variants", small), ("structural variants", sv)]:
        for a in ARM_ORDER:
            if a not in src:
                continue
            rows = src[a]["metrics"]["summary"]
            L.append(f"<details><summary><code>{a}</code> — {title}</summary>")
            L.append("")
            if rows:
                cols = list(rows[0].keys())
                L.append("| " + " | ".join(cols) + " |")
                L.append("|" + "---|" * len(cols))
                for r in rows:
                    L.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
            L.append("")
            L.append("</details>")
            L.append("")

    Path(args.out).write_text("\n".join(L) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
