#!/usr/bin/env python3
"""Render one contig's full tier-2 results: SNVs, small indels, and structural variants.

Reads the aardvark output directories written by run_arms.py (small-variant benchmark) and
the truvari directories written by truvari_sv.py (structural benchmark). Kept separate from
both so the tables can be regenerated without re-running any calling or comparison.

`--contig` selects the dataset; per-contig facts in the header are measured from the work
directory rather than typed in. A handful of narrative blocks quote per-site measurements
made on chr20 and are labelled as such when the page is built for another contig.

**SVs come from truvari.** aardvark's Sv* categories are scored against the *small-variant*
truth set, which holds no record over 50 bp, so there is next to nothing for them to match;
its summary also leaves query_total/query_tp/query_fp at zero for those rows, making its own
precision and F1 come out 0/0. The aardvark SV block is still emitted for continuity with
earlier runs, with precision recomputed from the per-variant BD decisions in its annotated
query VCF, but it is secondary and the page says so.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

ARM_ORDER = ["poisson", "poisson-z", "readlik-support", "readlik-nomismap",
             "readlik-nolink", "readlik"]

# The mismapping floor (--mismap-min) went 1e-8 -> 0.01 -> 0.05 -> 0.02, and the cap
# (--mismap-max) 0.1 -> 0.5 -> 0.7; the current defaults are 0.02 and 0.7. Arms run at the older
# values are kept and labelled rather than discarded, because the comparison between them
# *is* the result. poisson and poisson-z do not use the read-likelihood model at all, so
# neither clamp can reach them and they are not re-run.
FLOOR_UNAFFECTED = {"poisson", "poisson-z"}
CALIB_ARMS = [("fl0.05", "readlik, floor 0.05"),
              ("mm0.2", "readlik, cap 0.2"),
              ("mm0.4", "readlik, cap 0.4")]

SMALL_TYPES = [("Snv", "SNV"), ("Insertion", "Insertion (<50 bp)"),
               ("Deletion", "Deletion (<50 bp)"), ("Indel", "Indel"),
               ("JointIndel", "Indel (joint)"), ("ALL", "ALL")]
SV_TYPES = [("SvInsertion", "SV insertion (>=50 bp)"),
            ("SvDeletion", "SV deletion (>=50 bp)"),
            ("JointStructuralVariant", "SV (joint)")]

# Enumeration source, and whether the arm needs a pack file. The two callers reach panel
# enumeration by different routes and the labels say which: `-z` where the flag is what
# selects it, "default" where the caller now does it unasked. Only readlik-support carries a
# flag to get support enumeration, because for that caller support is no longer the default.
META = {
    "poisson": ("support (Flow)", "yes"),
    "poisson-z": ("panel (`-z`)", "yes"),
    "readlik-support": ("support (`--enumerate-support`)", "yes"),
    "readlik-nomismap": ("panel (default)", "**no**"),
    "readlik-nolink": ("panel (default)", "**no**"),
    "readlik": ("panel (default)", "**no**"),
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
    p.add_argument("--contig", default="chr20")
    p.add_argument("--results")
    p.add_argument("--out")
    args = p.parse_args()
    c = args.contig
    # The 34-haplotype graph is the default subject: it is the configuration the caller is tuned
    # for and the one that performs better, so it is what the headline page should describe. The
    # 4-haplotype graph gets its own page from this same script via --results/--out, and the
    # side-by-side comparison lives in compare_graphs.py.
    res = Path(args.results or REPO / f"work/tier2-{c}-hap32/results")
    out_path = Path(args.out or REPO / f"docs/tier2-{c}-results.md")
    # Derived from the results path rather than passed as a label, so a page cannot be titled with
    # the wrong panel size -- the one error that would make these two pages indistinguishable.
    is_hap32 = "hap32" in res.parent.name
    panel_label = "34-haplotype" if is_hap32 else "4-haplotype"
    # Several narrative blocks below quote per-site measurements made on chr20 -- the
    # pericentromeric pile-ups, the 246 large-insertion calls, the floor-change genotype
    # counts. Those are chr20 facts, not general ones, so they are gated rather than
    # reprinted under another contig's heading with chr20's numbers.
    is_chr20 = (c == "chr20")

    # arms*.json would also match arms-sv.json, whose entries would then overwrite the
    # small-variant ones (same arm names, newer mtime). Load the small-variant batches
    # explicitly instead.
    old = load_merged(res, "arms.floor-1e-8.json")
    old.update(load_merged(res, "arms.readlik.json"))
    # Current results: all five arms re-run together at the present defaults, so the
    # wall-clock column compares runs made on the same machine in the same session.
    # arms.floor-0.01.json is the earlier re-run at the new floor and has the same
    # calls, but timings from before the read path was optimised.
    small = load_merged(res, "arms.json")
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
    L.append(f"# Tier 2 results: HG002 {c} on HPRC v2.1 MC CHM13, {panel_label} graph")
    L.append("")
    L.append("Real reads, real benchmark, run on a 32 GB laptop.")
    L.append("")
    if is_hap32:
        L.append("This is the **34-haplotype** graph: CHM13, GRCh38 and 32 recombinants from "
                 "haplotype sampling. It is the primary subject because it is what the caller is "
                 "tuned for -- both the linkage transition and the panel frequency prior are "
                 "panel-size effects and have little to work with on a thin panel -- and because "
                 "it is the better-performing configuration. The 4-haplotype graph has its own "
                 f"page at [tier2-{c}-4hap-results.md](tier2-{c}-4hap-results.md), and the two "
                 f"are put side by side in "
                 f"[tier2-{c}-graph-comparison.md](tier2-{c}-graph-comparison.md).")
    else:
        L.append("This is the **4-haplotype** graph: CHM13, GRCh38 and 2 recombinants. It is kept "
                 "as a thin-panel reference rather than the headline configuration -- the caller "
                 f"is tuned on the 34-haplotype graph, whose page is "
                 f"[tier2-{c}-results.md](tier2-{c}-results.md). The two are compared directly in "
                 f"[tier2-{c}-graph-comparison.md](tier2-{c}-graph-comparison.md).")
    L.append("")
    # Per-contig facts are measured from the prepared work directory rather than typed in.
    # The chr20 page previously carried them as literals, which is how the node count and
    # region sizes would have silently followed the file around to another contig.
    work = res.parent
    def _wc(path: Path) -> int | None:
        try:
            return sum(1 for _ in path.open())
        except OSError:
            return None
    def _bed_span(path: Path) -> float | None:
        try:
            return sum(int(f[2]) - int(f[1])
                       for f in (l.split() for l in path.open()) if len(f) >= 3) / 1e6
        except OSError:
            return None
    nodes = _wc(work / f"{c}_all_nodes.txt")
    sm_mb, sv_mb = _bed_span(work / f"truth.{c}.smvar.bed"), _bed_span(work / f"truth.{c}.stvar.bed")

    L.append("| | |")
    L.append("|---|---|")
    # Was hardcoded to the 4-haplotype graph, which put "4 haplotypes, 2 recombinants" at the top
    # of the 34-haplotype page the moment this script grew a second subject. Keyed off the results
    # path instead, like the title.
    if is_hap32:
        L.append("| graph | `hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz`, 101,366,693 nodes, "
                 "**34 haplotypes** (CHM13, GRCh38, 32 recombinants from haplotype sampling; the "
                 "file is named for the recombinant count, not the total). HG002 itself is "
                 "**absent** — no circularity |")
    else:
        L.append("| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz`, 100,179,277 nodes, "
                 "**4 haplotypes** (CHM13, GRCh38, 2 recombinants). HG002 itself is **absent** — "
                 "no circularity |")
    L.append(f"| chromosome | {c} component"
             + (f", {nodes:,} nodes" if nodes else "") + " |")
    L.append("| reads | 596,017,764 alignments genome-wide (~28.6×); 151 bp paired Illumina |")
    L.append("| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |")
    L.append("| regions | "
             + (f"small variants {sm_mb:.1f} Mb" if sm_mb else "small variants —")
             + (f"; SVs {sv_mb:.1f} Mb" if sv_mb else "") + " |")
    L.append("| engine | `aardvark compare` for small variants; `truvari bench --sizemin 50` "
             "for SVs |")
    L.append("")
    L.append("**All read-likelihood arms below run at the current clamp defaults, "
             "`--mismap-min 0.02` and `--mismap-max 0.5`.** The floor caps how much one read can veto "
             "an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by "
             "measurement — the floor from 1e-8, the cap down from an original 0.1 that was actively "
             "wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. "
             "`poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.")
    L.append("")
    L.append("**Read the caveats before the numbers.** The benchmark is a *draft*: its own README "
             "reports known errors in highly homozygous regions, homopolymers and tandem repeats, and "
             "excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm "
             "comparison is what this table is for.")
    L.append("")

    L.append("## Cost")
    L.append("")
    L.append("Every arm on this page was re-run together on one build, so the wall-clock column "
             "compares runs made on the same machine in the same session rather than a mixture of "
             "vintages.")
    L.append("")
    L.append("Two changes since the accuracy results were first produced left the calls untouched. "
             "The read path was optimised (vg `44fd008`)" +
             (" — on chr20 `readlik` went **506 s to under 100 s**, so the read-likelihood caller "
              "is now near parity with the Poisson caller at matched enumeration rather than 5.9x, "
              "and `readlik-support` is *faster* than `poisson`" if is_chr20 else "") +
             ". Then `AD`, `BL`, `GQI` and the explained-share scaling of `GQ` were added — which "
             "rescales a quality and does not change a genotype. Both are confirmed by the variant "
             "counts below, which are unchanged to the record.")
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
        overlap = (f" ({sm_mb:.1f} Mb vs {sv_mb:.1f} Mb)" if sm_mb and sv_mb else "")
        L.append("**The `smvar` truth set contains no record >=50 bp** — that size class lives in the "
                 f"separate `stvar` benchmark. But the two confident regions overlap almost completely"
                 f"{overlap}. So a >=50 bp insertion called inside the small-variant confident region "
                 "has every one of its bases scored FP, however right the call is. It cannot be "
                 "scored correct.")
        L.append("")
        if is_chr20:
            L.append("That is exactly where the gap lives. 246 `readlik` calls carry a >=200 bp "
                     "insertion allele; they contribute **27,951 FP bases and zero TP bases**, which "
                     "is the whole of the precision difference. The Poisson caller scores better "
                     "there because it does not emit them — at the two largest sites it emits "
                     "nothing at all.")
        else:
            L.append("The same mechanism was traced site by site on chr20, where 246 `readlik` "
                     "calls carrying a >=200 bp insertion allele contribute 27,951 FP bases and zero "
                     "TP bases — the whole of the precision difference there. The size-matched "
                     "control below is the general test.")
        L.append("")
        L.append("Restricting **both** callers to the range the benchmark can adjudicate (dropping any "
                 "record with a called allele >=50 bp from REF, applied identically to each) gives the "
                 "size-matched comparison:")
        L.append("")
        L.append("| arm | class | BP recall | BP precision | **BP F1** |")
        L.append("|---|---|---|---|---|")
        for label in ("sm50-poisson-z", "sm50-readlik"):
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
        # Computed from this contig's own numbers. Hard-coding chr20's 0.139 -> 0.008 is how
        # the claim would have followed the page to another chromosome and been wrong there.
        def _sm(arm, key, vtype="Insertion"):
            return num(pick(sized.get(arm, {}).get("metrics", {}).get("summary", []),
                            "BASEPAIR", vtype), key)
        gap_un = None
        up = pick(small.get("poisson-z", {}).get("metrics", {}).get("summary", []),
                  "BASEPAIR", "Insertion")
        ur = pick(small.get("readlik", {}).get("metrics", {}).get("summary", []),
                  "BASEPAIR", "Insertion")
        if up and ur:
            gap_un = (num(up, "metric_precision") or 0) - (num(ur, "metric_precision") or 0)
        gp, gr = _sm("sm50-poisson-z", "metric_precision"), _sm("sm50-readlik", "metric_precision")
        gap_sm = (gp - gr) if (gp is not None and gr is not None) else None
        fp_, fr_ = _sm("sm50-poisson-z", "metric_f1"), _sm("sm50-readlik", "metric_f1")
        if gap_un is not None and gap_sm is not None:
            L.append(f"The insertion BASEPAIR precision gap collapses from **{gap_un:.3f} to "
                     f"{gap_sm:.3f}**"
                     + (f", and insertion BASEPAIR F1 goes from {fp_:.4f} for `poisson-z` against "
                        f"{fr_:.4f} for `readlik`" if (fp_ and fr_) else "") + ".")
        L.append("There is no insertion-sequence defect in the likelihood model; what the "
                 "unrestricted number measures is that one caller emits large insertions and the "
                 "other does not.")
        L.append("")
        L.append("Whether those large calls are *correct* is a separate question, and the truvari "
                 "comparison below is what answers it." +
                 (" On chr20, of the 246, only **35 are confirmed true**, **73 are confirmed "
                  "false**, and **138 fall outside the SV confident region** and cannot be judged "
                  "at all. See *Known bad output* for the worst of the unjudged ones."
                  if is_chr20 else ""))
        L.append("")

    tv = {}
    for a_ in ARM_ORDER:
        tp_ = res / f"truvari-{a_}" / "summary.json"
        if tp_.exists():
            tv[a_] = json.loads(tp_.read_text())
    if tv:
        L.append("## Structural variants — truvari (GIAB `stvar` benchmark)")
        L.append("")
        L.append("The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's "
                 "`Sv*` categories as the primary measure: those are scored against the "
                 "*small-variant* truth set, which contains no record over 50 bp at all, so they have "
                 "almost nothing to match (plan §9.22). The aardvark block below is kept for "
                 "continuity with earlier runs.")
        L.append("")
        L.append("**What these errors are made of, per record, is in "
                 "[tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a "
                 "quarter of all false positives are the metric rather than the caller, and that "
                 "harmonising representation with `truvari refine` moves every arm up by roughly "
                 "0.05 F1. Read the ranking between arms here; treat the absolute level as "
                 "benchmark-relative.")
        L.append("")
        L.append("| arm | recall | precision | **F1** | TP-base | FP | FN |")
        L.append("|---|---|---|---|---|---|---|")
        best_tv = max((v.get("f1") or 0) for v in tv.values())
        for a_ in ARM_ORDER:
            if a_ not in tv:
                continue
            v = tv[a_]
            mk = "**" if v.get("f1") and abs(v["f1"] - best_tv) < 1e-9 else ""
            L.append(f"| `{a_}` | {f(v.get('recall'))} | {f(v.get('precision'))} | "
                     f"{mk}{f(v.get('f1'))}{mk} | {int(v.get('TP-base', 0)):,} | "
                     f"{int(v.get('FP', 0)):,} | {int(v.get('FN', 0)):,} |")
        L.append("")

    if sv:
      L.append("## Structural variants — aardvark (secondary)")
    L.append("")
    L.append("Kept for continuity. These categories are scored against the small-variant truth set "
             "and should not be read as the SV result; prefer the truvari table above.")
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

    L.append("## Calibration: the two mismapping clamps")
    L.append("")
    L.append("MAPQ measures confidence that a read is in the right *place*, not that its path through a "
             "given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term "
             "cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats "
             "from one read** at the old floor of 1e-8. The floor caps that veto; the current default "
             "is **0.02**.")
    L.append("")
    if is_hap32:
        L.append("The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the "
                 "4-haplotype graph. There it reaches only reads whose `e_r` is already large — "
                 "6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a "
                 "knob with nothing to act on. On this graph 23.3% of reads sit at MAPQ 1, meaning "
                 "p(wrong) = 0.79, and the old cap of 0.1 was telling the model 0.1: overriding the "
                 "mapper at exactly the sites that matter. Raising it removed 94% of the excess "
                 "false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a "
                 "sparse graph is not thereby harmless — see ")
    else:
        L.append("The *upper* clamp (`--mismap-max`) looks inert on this graph, because it binds "
                 "only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% "
                 "at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the "
                 "old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% "
                 "of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising "
                 "it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. "
                 "A clamp that is inert on a sparse graph is not thereby harmless.")
    L.append("")
    L.append(f"The two graphs are put side by side in "
             f"[tier2-{c}-graph-comparison.md](tier2-{c}-graph-comparison.md); the grids are in "
             "plan §9.20.")
    L.append("")
    L.append("| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |")
    L.append("|---|---|---|---|---|---|")
    def calib_row(tag: str, label: str, from_json: dict | None = None) -> str | None:
        """Build one calibration row.

        Sweep arms were produced by ad-hoc scripts and have no arms*.json entry, so their
        aardvark output directory is the only record. The old-default rows are the
        exception and must come from the preserved JSON instead: re-running at the new
        default **overwrote** `aardvark-readlik/`, so reading that directory would
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

    # `small` is whatever the arms were last run at, so it must be labelled as the current
    # default rather than pinned to a value. It was labelled "floor 0.01" for two default
    # changes after that stopped being true.
    for tag, label, src in [("readlik", "floor 1e-8, cap 0.1 (original defaults)", old),
                            ("readlik", "**floor 0.02, cap 0.7 (current defaults)**", small),
                            ("fl0.05", "floor 0.05, cap 0.1", None),
                            ("mm0.2", "cap 0.2, floor 1e-8", None),
                            ("mm0.4", "cap 0.4, floor 1e-8", None)]:
        row = calib_row(tag, label, src)
        if row:
            L.append(row)
    L.append("")
    if old:
        L.append("Sweep rows other than the current one are historical: they were produced at the "
                 "defaults in force at the time and are kept because the comparison between them is "
                 "the result. The full grids are in plan §9.20-§9.21.")
    else:
        # The old-default arms were only ever preserved for the 4-haplotype runs, so on any other
        # dataset this table has one row. Saying where the history is beats printing a single row
        # under a caption that calls it a comparison -- and rows from two graphs must not share a
        # table, which is the failure the one-build-per-matrix rule exists to prevent.
        L.append("Only the current row is available here: the preserved old-default arms "
                 "(`arms.floor-1e-8.json`, `arms.readlik.json`) exist for the 4-haplotype runs "
                 f"alone, so the before-and-after is on "
                 f"[tier2-{c}-4hap-results.md](tier2-{c}-4hap-results.md). Mixing rows from two "
                 "graphs into one table is exactly what the one-build-per-matrix rule forbids. The "
                 "full grids are in plan §9.20-§9.21.")
    if is_chr20:
        L.append("")
        L.append("Raising the floor off 1e-8 changed **1,493 genotypes (1.41%)** on chr20, of which "
                 "**94% were heterozygous → homozygous** (1/0→1/1: 614, 0/1→1/1: 606, 1/2→1/1: 184), "
                 "and dropped 1,251 spurious non-reference calls. The failure it corrects is spurious "
                 "heterozygosity: a few locally misaligned reads, each able to veto the homozygous "
                 "hypothesis almost without bound, conjuring a second allele that is not there.")
    L.append("")
    L.append("The floor was later re-swept at the corrected cap, on both graphs and both benchmarks, "
             "and settled at **0.02**. 0.05 wins on small-variant `GT` but costs about 0.01 of SV F1 "
             "— which the first sweep never saw, because it was scored on one benchmark only. Plan "
             "§9.21 records that as a process rule: a sweep that sets a default has to be scored on "
             "every benchmark the project runs.")
    L.append("")

    L.append("## Known bad output" + ("" if is_chr20 else f" (measured on chr20)"))
    L.append("")
    L.append("Neither benchmark scores these, so they appear in no metric on this page. They are "
             "recorded because they are plainly wrong and would mislead anyone reading the VCF.")
    L.append("")
    L.append("`readlik` emits a small number of enormous homozygous insertions in and around the "
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
    L.append("The same blindness has a second consequence, found later and now corrected. Because the "
             "model only weighs reads it can see, it had no way to know that a heterozygous deletion "
             "produces *no* reads over the deleted interval, and its flat `1/ploidy` mixture asserted "
             "that both haplotypes contributed equally everywhere. That cost it 94% of heterozygous "
             "deletions above 1 kb and mis-genotyped two thirds of heterozygous insertions above 1 kb. "
             "Weighting each haplotype by the reads it is *expected* to contribute at the site is now "
             "the default and fixes both, without moving small variants at all — see "
             "[tier2-sv-errors.md](tier2-sv-errors.md). It did not remove the need for a depth term: "
             "it corrects the *relative* weight between a genotype's haplotypes, while the pile-ups "
             "above are a statement about *absolute* depth. That term is now also the default, at "
             "`--depth-term 0.1`, and the read arms in the tables on this page carry it — see "
             "[tier2-depth-term.md](tier2-depth-term.md). It does not resolve the pile-ups either: "
             "it detects them emphatically and still cannot outvote the read evidence at them, "
             "which is what the `DR` field and `--depth-quality` are for "
             "([tier2-quality-signals.md](tier2-quality-signals.md)).")
    L.append("")
    L.append("Filtering on depth is **not** that remedy, and that has now been tested properly "
             "rather than by two spot checks. Sweeping a two-sided cut on DP over a rolling local "
             "median, across both chromosomes and both graphs, against the one test a hard filter has "
             "to pass — beat lowering the GQ threshold to the same recall:")
    L.append("")
    L.append("- a **minimum** fails in all eight dataset-by-benchmark cells. Few reads already means "
             "a small likelihood gap, so low depth depresses GQ on its own and a separate cut adds "
             "nothing;")
    L.append("- a **maximum** passes in exactly one configuration — 5x the local median, structural "
             "calls, 34-haplotype graph, worth about +0.025 precision — and is dominated everywhere "
             "else. The two original spot checks (DP 200 moving insertion BASEPAIR precision by "
             "0.0001; DP 58 helping by +0.087 but costing SV insertion recall 0.4976 to 0.4167) were "
             "both right and both too narrow to conclude from.")
    L.append("")
    L.append("What shipped instead attacks the same blindness from the other side: **GQ is now scaled "
             "by the fraction of reads the called genotype explains**, so a pile-up the call does not "
             "account for can no longer carry a saturated quality. The giants remain output that no "
             "metric charges for — they should be fixed because they are wrong, not because they cost "
             "a score — but they no longer look confident. See "
             "[tier2-quality-signals.md](tier2-quality-signals.md).")
    L.append("")

    L.append("## Quality fields")
    L.append("")
    L.append("Every arm above is scored at **every** GQ, so nothing on this page depends on the "
             "quality field. `vg call` emits `AD` (per-allele read support, ties split "
             "fractionally), `BL` (mean absolute fit), `GQI` (the raw likelihood-ratio quality) and "
             "`GQ` (that ratio scaled by the fraction of reads the called genotype explains). The "
             "scaling rescales a quality and does not change a genotype, so **the numbers on this "
             "page are unaffected by it**; what it changes is how the calls rank. See "
             "[tier2-quality-signals.md](tier2-quality-signals.md).")
    L.append("")
    L.append("## The genotype mixture")
    L.append("")
    L.append("The read-likelihood arms on this page use the **length-weighted mixture**, which "
             "became the default after it was found that the flat `1/ploidy` weight breaks "
             "heterozygotes whose alleles differ in length. Unlike the `GQ` scaling above, this "
             "*does* change genotypes, so these numbers are not comparable with runs made before it. "
             "`--flat-mixture` restores the previous model exactly. Derivation and measurements: "
             "[tier2-sv-errors.md](tier2-sv-errors.md).")
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

    out_path.write_text("\n".join(L) + "\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
