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

# The mismapping floor (--mismap-min) was raised from 1e-8 to 0.01 after measurement
# showed the old value let a single MAPQ 60 read veto an allele by -13.8 nats. Arms run
# before that change are kept and labelled rather than discarded: the comparison between
# them *is* the result. poisson and poisson-z do not use the read-likelihood model at
# all, so the change cannot affect them and they are not re-run.
FLOOR_UNAFFECTED = {"poisson", "poisson-z"}
CALIB_ARMS = [("fl0.05", "readlik-z, floor 0.05"),
              ("mm0.2", "readlik-z, cap 0.2"),
              ("mm0.4", "readlik-z, cap 0.4")]

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
    old = load_merged(res, "arms.floor-1e-8.json")
    old.update(load_merged(res, "arms.readlik-z.json"))
    new = load_merged(res, "arms.floor-0.01.json")
    # Current-default table: poisson arms are floor-independent, read-likelihood arms
    # come from the re-run at 0.01.
    small = {k: v for k, v in old.items() if k in FLOOR_UNAFFECTED}
    small.update(new)
    # SV metrics come from the aardvark output directories rather than arms-sv.json:
    # compare_sv.py writes only the arms it was asked for, so the JSON is whatever the
    # last invocation happened to cover, while the directories accumulate.
    import csv as _csv
    sv = {}
    for a_ in ARM_ORDER:
        sp = res / f"aardvark-sv-{a_}" / "summary.tsv"
        if sp.exists():
            sv[a_] = {"arm": a_, "metrics": {"summary": list(_csv.DictReader(open(sp), delimiter="\t"))}}

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
    L.append("**All read-likelihood arms below use `--mismap-min 0.01`**, the current default. That floor "
             "caps how much one read can veto an allele; it was raised from 1e-8 after measurement, and "
             "the before/after comparison is in the calibration section at the end. `poisson` and "
             "`poisson-z` do not use the read-likelihood model, so the change cannot affect them.")
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

    sized = load_merged(res, "arms-size-matched.json")
    if sized:
        L.append("## Reading the insertion BASEPAIR numbers")
        L.append("")
        L.append("The insertion `BASEPAIR` precision above understates the read-likelihood caller, "
                 "and the reason is a property of the benchmark rather than of either caller.")
        L.append("")
        L.append("**The `smvar` truth set contains no record >=50 bp** — that size class lives in the "
                 "separate `stvar` benchmark. But the two confident regions overlap almost completely "
                 "(58.9 Mb vs 59.4 Mb). So a >=50 bp insertion called inside the small-variant "
                 "confident region has every one of its bases scored FP, however right the call is. "
                 "It cannot be scored correct.")
        L.append("")
        L.append("That is exactly where the gap lives. 246 `readlik-z` calls carry a >=200 bp "
                 "insertion allele; they contribute **27,951 FP bases and zero TP bases**, which is "
                 "the whole of the precision difference. The Poisson caller scores better there "
                 "because it does not emit them — at the two largest sites it emits nothing at all.")
        L.append("")
        L.append("Restricting **both** callers to the range the benchmark can adjudicate (dropping any "
                 "record with a called allele >=50 bp from REF, applied identically to each) gives the "
                 "size-matched comparison:")
        L.append("")
        L.append("| arm | class | BP recall | BP precision | **BP F1** |")
        L.append("|---|---|---|---|---|")
        for label in ("sm50-poisson-z", "sm50-readlik-z"):
            if label not in sized:
                continue
            rows = sized[label]["metrics"]["summary"]
            for vtype, vlabel in (("Insertion", "Insertion"), ("Deletion", "Deletion"), ("ALL", "ALL")):
                bp = pick(rows, "BASEPAIR", vtype)
                if not bp:
                    continue
                L.append(f"| `{label}` | {vlabel} | {f(num(bp,'metric_recall'))} | "
                         f"{f(num(bp,'metric_precision'))} | **{f(num(bp,'metric_f1'))}** |")
        L.append("")
        L.append("The insertion BASEPAIR precision gap collapses from **0.139 to 0.008**, and "
                 "insertion BASEPAIR F1 flips from a 0.047 loss into a 0.047 win. There is no "
                 "insertion-sequence defect in the likelihood model; what the unrestricted number "
                 "measures is that one caller emits large insertions and the other does not.")
        L.append("")
        L.append("Whether those large calls are *correct* is a separate question, and the `stvar` "
                 "comparison below is what answers it: they are a net win there (SV insertion recall "
                 "0.4976 vs 0.4263), but of the 246, only **35 are confirmed true**, **73 are "
                 "confirmed false**, and **138 fall outside the SV confident region** and cannot be "
                 "judged at all. See *Known bad output* for the worst of the unjudged ones.")
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

    L.append("## Calibration: the mismapping floor")
    L.append("")
    L.append("MAPQ measures confidence that a read is in the right *place*, not that its path through a "
             "given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term "
             "cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats "
             "from one read** at the old floor of 1e-8. Raising the floor caps that veto.")
    L.append("")
    L.append("The *upper* clamp (`--mismap-max`) is inert here: it binds only where `e_r` is already "
             "large, i.e. the 6.3% of reads at MAPQ ≤ 9, while 90% are MAPQ 60.")
    L.append("")
    L.append("| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |")
    L.append("|---|---|---|---|---|---|")
    def calib_row(tag: str, label: str, from_json: dict | None = None) -> str | None:
        """Build one calibration row.

        Sweep arms were produced by ad-hoc scripts and have no arms*.json entry, so their
        aardvark output directory is the only record. The old-default rows are the
        exception and must come from the preserved JSON instead: re-running at the new
        default **overwrote** `aardvark-readlik-z/`, so reading that directory would
        silently report the new numbers under the old label -- which is exactly the
        before/after comparison this table exists to make.
        """
        if from_json is not None:
            rows_ = from_json.get(tag, {}).get("metrics", {}).get("summary", [])
            if not rows_:
                return None
        else:
            summary_path = res / f"aardvark-{tag}" / "summary.tsv"
            if not summary_path.exists():
                return None
            import csv as _csv
            rows_ = list(_csv.DictReader(open(summary_path), delimiter="\t"))

        def g(comparison: str, vtype: str) -> str:
            return f(num(pick(rows_, comparison, vtype), "metric_f1"))

        return (f"| {label} | {g('GT','ALL')} | {g('GT','Snv')} | {g('GT','Insertion')} | "
                f"{g('GT','Deletion')} | {g('BASEPAIR','ALL')} |")

    for tag, label, src in [("readlik-z", "floor 1e-8 (old default)", old),
                            ("readlik-z", "**floor 0.01 (current default)**", new),
                            ("fl0.05", "floor 0.05", None),
                            ("mm0.2", "cap 0.2, floor 1e-8", None),
                            ("mm0.4", "cap 0.4, floor 1e-8", None)]:
        row = calib_row(tag, label, src)
        if row:
            L.append(row)
    L.append("")
    L.append("Raising the floor to 0.01 changed **1,493 genotypes (1.41%)**, of which **94% were "
             "heterozygous → homozygous** (1/0→1/1: 614, 0/1→1/1: 606, 1/2→1/1: 184), and dropped 1,251 "
             "spurious non-reference calls. The failure it corrects is spurious heterozygosity: a few "
             "locally misaligned reads, each able to veto the homozygous hypothesis almost without "
             "bound, conjuring a second allele that is not there.")
    L.append("")
    L.append("Calibrated on one chromosome of one sample. 0.05 is better on indel `GT` but costs SNVs "
             "and BASEPAIR, so the optimum lies between and is not worth over-fitting here.")
    L.append("")

    L.append("## Known bad output")
    L.append("")
    L.append("Neither benchmark scores these, so they appear in no metric on this page. They are "
             "recorded because they are plainly wrong and would mislead anyone reading the VCF.")
    L.append("")
    L.append("`readlik-z` emits a small number of enormous homozygous insertions in and around the "
             "chr20 pericentromere, at depths that are physically impossible:")
    L.append("")
    L.append("| position | called insertion | GT | DP | GQ |")
    L.append("|---|---|---|---|---|")
    for pos, ln, gt, dp, gq in [(25849044, 61958, "1/1", 7873, 256),
                                (32179077, 57716, "1/1", 5337, 256),
                                (1629728, 33050, "1/1", 291, 256),
                                (25873453, 28685, "1/2", 5498, 256),
                                (25792993, 23450, "1/1", 932, 256)]:
        L.append(f"| chr20:{pos:,} | {ln:,} bp | {gt} | {dp:,} | {gq} |")
    L.append("")
    L.append("Chromosome-median DP is **29**, and the Poisson caller's expected depth (`XD`) never "
             "exceeds **167** anywhere on chr20. Median DP rises monotonically with called insertion "
             "length — 28 for 1 bp, 28 for 2–15 bp, 35 for 50–199 bp, **330 for >=1 kb** — so these "
             "are collapsed-repeat pile-ups, not haplotypes.")
    L.append("")
    L.append("The read-likelihood model cannot reject them, and the reason is structural rather than "
             "a tuning failure: it computes P(reads | genotype) **conditioned on the reads it is "
             "given**, and never asks whether that many reads should be there. The Poisson caller gets "
             "this for free, because an observed-vs-expected depth term is the whole of its model. A "
             "depth-plausibility guard is the obvious remedy, and the expected depth is already "
             "reachable — the read-likelihood caller subclasses `SupportBasedSnarlCaller` and holds a "
             "`TraversalSupportFinder` for allele enumeration.")
    L.append("")
    L.append("Filtering on depth is **not** that remedy, and the measurement says so plainly: "
             "dropping every call above DP 200 removes 195 records including all of the giants above, "
             "and moves insertion BASEPAIR precision by 0.0001 (0.6226 → 0.6227). Dropping above DP 58 "
             "removes 1,202 records and does help (+0.087), but costs SV insertion recall "
             "0.4976 → 0.4167 — it is a blunt proxy for length that discards real SVs. The giants are "
             "bad output that no metric charges for; they should be fixed because they are wrong, not "
             "because they cost a score.")
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
