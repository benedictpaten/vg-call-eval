#!/usr/bin/env python3
"""Render the tier-2 chr20 arm results as a markdown table.

Reads the aardvark output directories written by run_arms.py. Kept separate so
the table can be regenerated without re-running any calling.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Which aardvark metric rows to show, and in what order. Aardvark emits several
# comparison types; GT is the genotype-aware one that matters for a genotyper,
# BASEPAIR gives partial credit (see plan §2.2).
WANTED_TYPES = ["GT", "BASEPAIR"]

# aardvark's actual summary.tsv schema, confirmed against output rather than assumed:
#   compare_label comparison region_label filter variant_type
#   truth_total truth_tp truth_fn query_total query_tp query_fp
#   metric_recall metric_precision metric_f1 truth_fn_gt query_fp_gt
# comparison is GT | BASEPAIR; variant_type is ALL | Snv | Insertion | Deletion |
# Indel | JointIndel.
VARIANT_TYPES = ["ALL", "Snv", "Indel", "Insertion", "Deletion"]


def pick(rows: list[dict], comparison: str, variant_type: str) -> dict | None:
    for r in rows:
        if (r.get("comparison", "").upper() == comparison.upper()
                and r.get("variant_type", "") == variant_type):
            return r
    return None


def fnum(row: dict | None, *names: str) -> str:
    if not row:
        return "—"
    for n in names:
        if n in row and row[n] not in ("", None):
            try:
                return f"{float(row[n]):.4f}"
            except ValueError:
                return str(row[n])
    return "—"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default=str(HERE / "results"))
    p.add_argument("--out", default=str(HERE / "RESULTS.md"))
    args = p.parse_args()

    res = Path(args.results)
    # Arms were run in two batches (the pack-free one first, before the pack existed),
    # and each batch overwrites arms.json with only its own arms. Merge every
    # arms*.json, keeping the newest entry per arm name, so no arm is silently missing
    # from the table.
    by_name: dict[str, dict] = {}
    for f in sorted(res.glob("arms*.json"), key=lambda f: f.stat().st_mtime):
        for entry in json.loads(f.read_text()):
            by_name[entry["arm"]] = entry
    order = ["poisson", "poisson-z", "readlik", "readlik-nomismap", "readlik-z"]
    payload = ([by_name[n] for n in order if n in by_name]
               + [v for k, v in by_name.items() if k not in order])

    lines: list[str] = []
    lines.append("# Tier 2: HG002 chr20 on HPRC v2.1 MC CHM13")
    lines.append("")
    lines.append("Truth: GIAB HG002 draft benchmark, defrabb V0.019-20241113, CHM13v2.0 small")
    lines.append("variants, restricted to the benchmark BED (58.9 Mb, 88.9% of chr20).")
    lines.append("")
    lines.append("**The benchmark is a draft** — its README notes known errors in homozygous")
    lines.append("regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute")
    lines.append("numbers should be read with that in mind; the arm-to-arm comparison is the point.")
    lines.append("")
    lines.append("**Recall is bounded by the graph, not only by the caller.** This graph carries 4")
    lines.append("haplotypes (CHM13, GRCh38, 2 recombinants), so `-z` arms can only propose alleles")
    lines.append("present in those walks. Compare down a column for the caller effect, across a row")
    lines.append("for what the sampled graph costs.")
    lines.append("")

    lines.append("## Cost")
    lines.append("")
    lines.append("| arm | variants | wall | peak RSS |")
    lines.append("|---|---|---|---|")
    for a in payload:
        lines.append(f"| `{a['arm']}` | {a['variants']:,} | {a['seconds']:.0f} s | "
                     f"{a['peak_rss_gb']:.1f} GB |")
    lines.append("")

    for vtype in VARIANT_TYPES:
        rows_present = any(
            pick(a.get("metrics", {}).get("summary", []), ct, vtype)
            for a in payload for ct in WANTED_TYPES
        )
        if not rows_present:
            continue
        lines.append(f"## {vtype}")
        lines.append("")
        lines.append("| arm | comparison | recall | precision | F1 | TP | FN | FP |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for a in payload:
            summary = a.get("metrics", {}).get("summary", [])
            for ct in WANTED_TYPES:
                r = pick(summary, ct, vtype)
                if not r:
                    continue
                lines.append(
                    f"| `{a['arm']}` | {ct} | {fnum(r, 'metric_recall')} | "
                    f"{fnum(r, 'metric_precision')} | {fnum(r, 'metric_f1')} | "
                    f"{int(r.get('truth_tp', 0)):,} | {int(r.get('truth_fn', 0)):,} | "
                    f"{int(r.get('query_fp', 0)):,} |"
                )
        lines.append("")

    # Raw summary rows, so nothing is hidden behind the curated view above.
    lines.append("## Raw aardvark summary rows")
    lines.append("")
    for a in payload:
        summary = a.get("metrics", {}).get("summary", [])
        lines.append(f"<details><summary><code>{a['arm']}</code></summary>")
        lines.append("")
        if summary:
            cols = list(summary[0].keys())
            lines.append("| " + " | ".join(cols) + " |")
            lines.append("|" + "---|" * len(cols))
            for r in summary:
                lines.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
        else:
            lines.append("_no summary rows_")
        lines.append("")
        lines.append("</details>")
        lines.append("")

    Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
