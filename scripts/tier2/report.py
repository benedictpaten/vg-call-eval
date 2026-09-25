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
precision and F1 come out 0/0. The aardvark SV block is emitted, for continuity with earlier
runs, only when its aardvark-sv-* output directories exist, with precision recomputed from the
per-variant BD decisions in its annotated query VCF; it is secondary and the page says so.
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
# (--mismap-max) 0.1 -> 0.5 -> 0.7 -> 0.95; the current defaults are read off the binary by
# caller_clamp_defaults(), not kept here. Arms run at the older
# values are kept and labelled rather than discarded, because the comparison between them
# *is* the result. poisson and poisson-z do not use the read-likelihood model at all, so
# neither clamp can reach them and they are not re-run.
FLOOR_UNAFFECTED = {"poisson", "poisson-z"}

# A version string that does not name the commit whose code it is. A binary built from a working
# tree before that tree is committed reports the parent. Each entry was checked: the sources were
# last written before the build and the tree was clean at the commit.
BUILD_NOTES = {
    "vg version v1.4.0-18924-g8acbb43a2":
        "the code of `0cab3fbd4`, built before that commit, so the string names its parent",
}
CALIB_ARMS = [("fl0.05", "readlik, floor 0.05"),
              ("mm0.2", "readlik, cap 0.2"),
              ("mm0.4", "readlik, cap 0.4")]

# aardvark's plain `Indel` row is query-only -- truth_total 0, so no recall and no F1 -- and printing
# it put a table of dashes beside the real one. The indel row is JointIndel.
SMALL_TYPES = [("Snv", "SNV"), ("Insertion", "Insertion (<50 bp)"),
               ("Deletion", "Deletion (<50 bp)"), ("JointIndel", "Indel"), ("ALL", "ALL")]
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


def caller_clamp_defaults(vg: str | None = None) -> tuple[str, str]:
    """(--mismap-min, --mismap-max) as `vg call` currently defaults them.

    Parsed from the help text's trailing `[value]` rather than kept as a literal here. The
    literal version of this went stale silently: the cap moved 0.5 -> 0.7 and every generated
    page went on claiming the arms ran at 0.5, which is precisely the kind of drift between
    what a run did and what the page says it did that this harness exists to catch. Falls back
    to the last known values if the binary is absent, so a page can still be built without it.
    """
    import re
    import subprocess

    fallback = ("0.02", "0.95")
    exe = vg or str(Path.home() / "CLionProjects/vg/bin/vg")
    try:
        help_text = subprocess.run([exe, "call", "--help"], capture_output=True,
                                   text=True, timeout=60).stderr
    except (OSError, subprocess.SubprocessError):
        return fallback
    found = {}
    for flag in ("mismap-min", "mismap-max"):
        # The default sits in [brackets] at the end of a multi-line option description, so
        # scan from the flag to the next option rather than to the end of its first line.
        m = re.search(rf"--{flag}\b.*?\[([0-9.]+)\]", help_text, re.S)
        if m:
            found[flag] = m.group(1)
    if len(found) != 2:
        return fallback
    return found["mismap-min"], found["mismap-max"]


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


# --- figures read off the VCFs ------------------------------------------------------------------
# The large-insertion and pile-up paragraphs used to quote chr20 figures from one old run (246
# calls, 27,951 FP bases, a five-row table of giants). They are computed here instead, so a refresh
# moves them with everything else.
BIG_INS, GIANT_INS = 200, 10_000
INS_BINS = [(1, 1, "1 bp"), (2, 15, "2-15 bp"), (16, 49, "16-49 bp"), (50, 199, "50-199 bp"),
            (200, 999, "200-999 bp"), (1000, None, ">=1 kb")]


def vcf_records(path: Path):
    """(chrom, pos, ref, alts, FORMAT dict) per record of a bgzipped single-sample VCF."""
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 10:
                yield f[0], int(f[1]), f[3], f[4].split(","), dict(zip(f[8].split(":"), f[9].split(":")))


def called_insertions(ref: str, alts: list[str], gt: str) -> list[int]:
    """Net insertion length, len(ALT) - len(REF), of each ALT allele the genotype carries."""
    out = []
    for a in {int(x) for x in gt.replace("|", "/").split("/") if x.isdigit() and x != "0"}:
        alt = alts[a - 1] if a <= len(alts) else "*"
        if not alt.startswith("<") and alt != "*" and len(alt) > len(ref):
            out.append(len(alt) - len(ref))
    return out


def called_insertion(ref: str, alts: list[str], gt: str) -> int:
    """Longest insertion, in bp, among the ALT alleles the genotype carries; 0 if none."""
    return max(called_insertions(ref, alts, gt), default=0)


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def insertion_census(vcf: Path) -> dict | None:
    """Contig median DP over every record, and (chrom, pos, longest insertion, GT, DP, GQ, REF length,
    insertion lengths) for every record carrying an insertion."""
    if not vcf.exists():
        return None
    import statistics
    dps, ins = [], []
    for chrom, pos, ref, alts, fmt in vcf_records(vcf):
        dp = _int(fmt.get("DP"))
        if dp is not None:
            dps.append(dp)
        lens = called_insertions(ref, alts, fmt.get("GT", ""))
        if lens:
            ins.append((chrom, pos, max(lens), fmt.get("GT", ""), dp, _int(fmt.get("GQ")), len(ref), lens))
    return {"median_dp": statistics.median(dps) if dps else None, "ins": ins}


def count_big_insertions(vcf: Path, min_len: int = BIG_INS) -> int | None:
    if not vcf.exists():
        return None
    return sum(1 for _, _, ref, alts, fmt in vcf_records(vcf)
               if called_insertion(ref, alts, fmt.get("GT", "")) >= min_len)


def truvari_big_insertions(res: Path, arm: str = "readlik", min_len: int = BIG_INS):
    """(split alleles, confirmed true, confirmed false) for insertion alleles >= min_len, all counted
    after truvari_sv.py split the multiallelic records, so the three share a unit."""
    norm, tv = res / "truvari-norm" / f"{arm}.norm.vcf.gz", res / f"truvari-{arm}"
    if not (norm.exists() and (tv / "tp-comp.vcf.gz").exists() and (tv / "fp.vcf.gz").exists()):
        return None
    big = lambda path: sum(1 for _, _, ref, alts, fmt in vcf_records(path)
                           if called_insertion(ref, alts, fmt.get("GT", "")) >= min_len)
    return big(norm), big(tv / "tp-comp.vcf.gz"), big(tv / "fp.vcf.gz")


def max_expected_depth(vcf: Path) -> float | None:
    """The Poisson caller's largest expected depth (FORMAT/XD) anywhere in its VCF."""
    if not vcf.exists():
        return None
    xs = [x for x in (_int(fmt.get("XD")) for *_, fmt in vcf_records(vcf)) if x is not None]
    return max(xs) if xs else None


# --- long reads -------------------------------------------------------------------------------
# The ONT ablations run on their own graph (E821: 16 sampled haplotypes plus CHM13 and GRCh38) with
# the same truth and regions, in work/E821-<contig>/results-<MMDD>-{default,preset}. The newest date
# with both configurations present is the one reported.
ONT_ARMS = ["readlik", "readlik-nomismap", "readlik-nolink"]
ONT_CONFIGS = [("default", "short-read defaults"), ("preset", "`--preset ont`")]
# What `--preset ont` sets, read from the preset block in call_main.cpp rather than from its help
# line, which omits --read-min-mapq 5. Update it with that block.
ONT_PRESET = ("--gap-open 1 --gap-extend 1 --mismap-min 0.05 --read-min-mapq 5 "
              "--insertion-nats 0.9 --read-phasing --regenotype")
# Share of chr20 alignments whose e_r sits on the --mismap-max ceiling, keyed by ceiling; at 0.95 it
# binds on MAPQ 0 alone. Counted over the full tier-2 chr20 GAFs: 658,419 of 13,278,140 short-read
# alignments and 542 of 85,373 ONT reads (work/chr20.ont.gaf.gz). A ceiling with no entry here has
# not been measured, and the page says so rather than quoting another's shares.
CEILING_TAIL = {"0.95": ("MAPQ 0", 4.96, 0.63)}
# ONT read lengths on chr20, measured from the alignments; Illumina is 151 bp paired.
ONT_READ_LEN = {"mean": 33_449, "median": 19_222, "p10": 1_845, "p90": 86_190}
SR_READ_LEN = 151


def ont_results(c: str):
    """(date tag, {config: dir}) for the newest ONT ablation with both configurations, or None."""
    tags: dict[str, dict[str, Path]] = {}
    for d in (REPO / f"work/E821-{c}").glob("results-*-*"):
        tag, cfg = d.name[len("results-"):].rsplit("-", 1)
        if (d / "arms.json").exists():
            tags.setdefault(tag, {})[cfg] = d
    full = [t for t, v in tags.items() if all(k in v for k, _ in ONT_CONFIGS)]
    return (max(full), tags[max(full)]) if full else None


def source_depth(c: str, tech: str) -> float | None:
    """Full read depth as the coverage sweep measured it, from subsample_gaf.py's log."""
    for lg in sorted((REPO / f"work/cov/{tech}-{c}").glob("subsample.*.log")):
        for line in lg.read_text().splitlines():
            if line.startswith("source:"):
                return float(line.rsplit("->", 1)[1].strip().rstrip("x"))
    return None


def gt_metric(entry: dict | None, vtype: str, key: str):
    if not entry:
        return None
    return num(pick(entry["metrics"]["summary"], "GT", vtype), key)


def ont_section(c: str, small: dict, hi: str) -> list[str]:
    found = ont_results(c)
    if not found:
        return []
    tag, dirs = found
    data = {cfg: load_merged(d, "arms.json") for cfg, d in dirs.items()}
    L: list[str] = []
    L.append("## Long reads — ONT, 16-haplotype E821 graph")
    L.append("")
    depth = source_depth(c, "ont")
    L.append("Same sample, truth and confident regions as above; different reads and a different "
             "graph. The ONT reads"
             + (f" ({depth:.1f}x on {c}, " if depth else " (")
             + f"mean length {ONT_READ_LEN['mean']:,} bp on chr20) are aligned to the E821 graph: 16 "
             "haplotypes from haplotype sampling plus CHM13 and GRCh38, 18 in the panel against 34 "
             "above. So a figure here is not comparable with one in the short-read tables. The pair "
             "of columns is comparable, and it answers what the long-read preset buys on long-read "
             f"data. The preset sets `{ONT_PRESET}`.")
    L.append("")
    builds = sorted({e.get("vg_version", "") for v in data.values() for e in v.values()
                     if e.get("vg_version")})
    sr_builds = sorted({small[a].get("vg_version", "") for a in ARM_ORDER
                        if a in small and small[a].get("vg_version")})
    if len(builds) == 1:
        same = " — the same build as the short-read arms" if builds == sr_builds else ""
        note = BUILD_NOTES.get(builds[0])
        L.append(f"Build: `{builds[0]}`" + (f" ({note})" if note else "") + f"{same}. Directories: "
                 + ", ".join(f"`{dirs[cfg].relative_to(REPO)}`" for cfg, _ in ONT_CONFIGS) + ".")
    else:
        L.append("**These ONT rows come from more than one build and are not comparable**: "
                 + ", ".join(f"`{b}`" for b in builds) + ". Re-run them.")
    L.append("")

    rows = [("readlik", "ALL", "metric_f1", "`readlik` ALL F1"),
            ("readlik", "Snv", "metric_f1", "`readlik` SNV F1"),
            ("readlik", "JointIndel", "metric_f1", "`readlik` indel F1"),
            ("readlik", "JointIndel", "metric_recall", "`readlik` indel recall"),
            ("readlik", "JointIndel", "metric_precision", "`readlik` indel precision"),
            ("readlik-nomismap", "ALL", "metric_f1", "`readlik-nomismap` ALL F1"),
            ("readlik-nolink", "ALL", "metric_f1", "`readlik-nolink` ALL F1")]
    L.append("| | " + " | ".join(lab for _, lab in ONT_CONFIGS) + " | Δ |")
    L.append("|---|---|---|---|")
    for arm, vtype, key, label in rows:
        a, b = (gt_metric(data[cfg].get(arm), vtype, key) for cfg, _ in ONT_CONFIGS)
        # the difference of the printed values, so the column reads off the table
        d = f"{round(b, 4) - round(a, 4):+.4f}" if a is not None and b is not None else "—"
        L.append(f"| {label} | {f(a)} | {f(b)} | {d} |")
    sv = {cfg: dirs[cfg] / "truvari-readlik" / "summary.json" for cfg, _ in ONT_CONFIGS}
    if all(pth.exists() for pth in sv.values()):
        a, b = (json.loads(sv[cfg].read_text()).get("f1") for cfg, _ in ONT_CONFIGS)
        L.append(f"| `readlik` SV F1 (truvari) | {f(a)} | {f(b)} | {round(b, 4) - round(a, 4):+.4f} |")
    L.append("")
    L.append("Small variants are aardvark's GT comparison, with indels from its joint indel row.")
    L.append("")

    L.append("| arm | configuration | variants | wall | CPU | peak RSS |")
    L.append("|---|---|---|---|---|---|")
    for arm in ONT_ARMS:
        for cfg, lab in ONT_CONFIGS:
            e = data[cfg].get(arm)
            if not e:
                continue
            cpu = e.get("cpu_seconds")
            L.append(f"| `{arm}` | {lab} | {e['variants']:,} | {e['seconds']:.0f} s | "
                     + (f"{cpu:,.0f} s |" if cpu else "— |") + f" {e['peak_rss_gb']:.1f} GB |")
    L.append("")
    L.append("Run serially on the same machine as the short-read arms, `--threads 5`.")
    L.append("")

    def delta(src, arm_b):
        a, b = gt_metric(src.get("readlik"), "ALL", "metric_f1"), gt_metric(src.get(arm_b), "ALL", "metric_f1")
        return None if a is None or b is None else a - b
    dm = {cfg: delta(data[cfg], "readlik-nomismap") for cfg, _ in ONT_CONFIGS}
    dl = {cfg: delta(data[cfg], "readlik-nolink") for cfg, _ in ONT_CONFIGS}
    dm_sr, dl_sr = delta(small, "readlik-nomismap"), delta(small, "readlik-nolink")
    sgn = lambda x: "—" if x is None else ("0.0000" if abs(x) < 5e-5 else f"{x:+.4f}")
    L.append(f"*The MAPQ mismapping term.* `readlik` minus `readlik-nomismap` is {sgn(dm['default'])} "
             f"ALL F1 on ONT at short-read defaults and {sgn(dm['preset'])} under the preset, against "
             f"{sgn(dm_sr)} on short reads. `--no-mismap-term` puts every read on the floor, so the "
             "term changes only a read whose `e_r = clamp(10^(−MAPQ/10), --mismap-min, "
             "--mismap-max)` is above it (MAPQ below 17 at the 0.02 floor), and it reduces a read on "
             "the ceiling to near-silence.")
    if hi in CEILING_TAIL:
        q, sr_pct, ont_pct = CEILING_TAIL[hi]
        L[-1] += (f" At the default ceiling of {hi} that is {q} alone: {sr_pct}% of chr20 short-read "
                 f"alignments against {ont_pct}% of ONT, {sr_pct / ont_pct:.0f}x fewer. Long reads "
                 "anchor uniquely, so there is almost no ambiguous-placement class for the term to act "
                 "on. Under the preset, `--read-min-mapq 5` drops those reads before scoring, so no "
                 "read reaches the ceiling and the term acts only between MAPQ 5 and 13, where "
                 "10^(−MAPQ/10) lies above the preset's 0.05 floor.")
    else:
        L[-1] += (f" The share of reads on the ceiling has not been measured at the current ceiling "
                  f"of {hi}.")
    L.append("")
    L.append(f"*The linkage layer.* `readlik` minus `readlik-nolink` is {sgn(dl['default'])} on ONT at "
             f"short-read defaults and {sgn(dl['preset'])} under the preset, against {sgn(dl_sr)} on "
             "short reads. `--linkage-weight 0` turns the whole layer off, and read phasing and "
             "phase-driven re-genotyping run inside it, so under the preset that arm loses both as "
             "well as the transition model. Switch error across depth, and ONT's SV figures against "
             "short reads with confidence intervals, are in [coverage.md](coverage.md).")
    L.append("")
    return L


def mixture_at_long_reads() -> list[str]:
    """The length-weighted mixture's weights at Illumina and ONT read lengths, computed rather than
    quoted: w_h is proportional to L_h + R - 1, so a deletion of D bp against its reference allele
    has a weight ratio of (D + R - 1) / (R - 1)."""
    L = ["### At long read length", ""]
    L.append("`w_h ∝ L_h + R − 1` counts the read START positions that yield an overlap, so R is the "
             "full read length; the *scoring* window is the read's overlap with the site, but the "
             "weight answers a sampling question, not a fit question. At ONT read length the "
             f"correction all but vanishes. With the measured chr20 mean of {ONT_READ_LEN['mean']:,} bp "
             f"against Illumina's {SR_READ_LEN}:")
    L.append("")
    L.append(f"| | Illumina R={SR_READ_LEN} | ONT R={ONT_READ_LEN['mean']:,} |")
    L.append("|---|---|---|")
    ratio = lambda dlen, r: (dlen + r - 1) / (r - 1)
    fx = lambda x: f"{x:.2f}x" if x < 10 else (f"{x:.1f}x" if x < 100 else f"{x:.0f}x")
    for label, dlen in (("50 bp deletion", 50), ("300 bp Alu", 300), ("5 kb deletion", 5000),
                        ("30 kb deletion", 30000)):
        cells = [f"{fx(ratio(dlen, r))}, w = {ratio(dlen, r) / (1 + ratio(dlen, r)):.3f}"
                 for r in (SR_READ_LEN, ONT_READ_LEN["mean"])]
        L.append(f"| {label} | " + " | ".join(cells) + " |")
    L.append("")
    L.append("That is the model being right rather than distorted: a haplotype carrying a 5 kb "
             f"deletion really does yield only {100 * (1 - 1 / ratio(5000, ONT_READ_LEN['mean'])):.0f}% "
             f"fewer {ONT_READ_LEN['mean'] / 1000:.0f} kb reads over the site, where it yields "
             f"{ratio(5000, SR_READ_LEN):.0f}x fewer {SR_READ_LEN} bp reads. The flat mixture that `--flat-mixture` restores, and that "
             "loses large heterozygous deletions on short reads, is very nearly what the length "
             "weighting already computes for ONT.")
    L.append("")
    L.append(f"The caveat is the mean, not the term. ONT read lengths are heavily skewed (mean "
             f"{ONT_READ_LEN['mean']:,}, median {ONT_READ_LEN['median']:,}, p10 {ONT_READ_LEN['p10']:,}, "
             f"p90 {ONT_READ_LEN['p90']:,}), so one mean R summarises a distribution that spans two "
             "orders of magnitude. A per-read R is the obvious refinement and has not been measured.")
    L.append("")
    return L


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
        vcf_ = res / f"{a_}.vcf.gz"
        # Only when it scores the VCF on disk. These directories are not re-made by a refresh, so
        # an old one would put a previous build's numbers under this page's build line.
        if sp.exists() and vcf_.exists() and sp.stat().st_mtime >= vcf_.stat().st_mtime:
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
    # Read the clamps off the binary rather than restating them. The arms pass no clamp flags,
    # so "the current defaults" is true by construction -- but the *values* drift, and this
    # sentence sat at 0.5 for as long as the real default was 0.7, silently mis-describing how
    # every arm on all four pages had been run.
    lo, hi = caller_clamp_defaults()
    L.append(f"**All read-likelihood arms below run at the current clamp defaults, "
             f"`--mismap-min {lo}` and `--mismap-max {hi}`.** The floor caps how much one read can "
             "veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by "
             "measurement — the floor from 1e-8, the cap raised from an original 0.1 that was actively "
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
    # The build that produced these numbers, read from the results rather than asserted. "One
    # build, one pass" is the rule this harness exists to enforce, and it used to rest on procedure
    # alone -- nothing in arms.json said which binary ran, so a page could claim it and be wrong.
    # Over the arms the page shows: arms.json keeps entries for retired arms (chr6 still holds a
    # readlik-nomismap-support from an old build), and those are not rows of this table.
    builds = sorted({small[a].get("vg_version", "") for a in ARM_ORDER
                     if a in small and small[a].get("vg_version")})
    L.append("Every number on this page — accuracy and cost alike — comes from one `vg` build in "
             "one pass, which is what the refresh harness exists to guarantee: a table whose rows "
             "come from different builds is not a comparison, it is a mixture of vintages.")
    L.append("")
    if len(builds) == 1:
        note = BUILD_NOTES.get(builds[0])
        L.append(f"Build: `{builds[0]}`" + (f" ({note})." if note else "."))
    elif len(builds) > 1:
        L.append("**These rows come from more than one build and are not comparable**: "
                 + ", ".join(f"`{b}`" for b in builds) + ". Re-run the matrix.")
    else:
        L.append("The build is not recorded in these results, so the guarantee above rests on "
                 "procedure rather than on a recorded fact. A refresh from here on records it.")
    L.append("")
    L.append("The wall column is what the caller costs unaided, and the repeatability note below "
             "applies to it harder than to the memory column. It includes snarl decomposition, "
             "which is single-threaded — 46 s of a 197 s chr20 run — and which `vg call -r` skips "
             "for byte-identical output given `vg snarls -T -P <ref path>`. The whole-genome "
             "harness caches one snarl file per contig for exactly that reason; this matrix does "
             "not, so these figures include it.")
    L.append("")
    # A cost delta where the previous refresh was a different build, because "is it getting
    # slower" is a question this table gets asked and could not answer without diffing two git
    # revisions of a markdown file by eye. Only shown where every arm has a comparison point, so
    # the column is either meaningful for the whole table or absent from it.
    prev = {a: small[a].get("previous") for a in ARM_ORDER if a in small}
    show_delta = bool(prev) and all(p and p.get("seconds") for p in prev.values())
    have_cpu = any(small[a].get("cpu_seconds") for a in ARM_ORDER if a in small)
    head = "| arm | enumeration | pack? | variants | wall |"
    if have_cpu:
        head += " CPU |"
    if show_delta:
        head += " Δ wall |"
        if have_cpu:
            head += " Δ CPU |"
    head += " peak RSS |"
    L.append(head)
    L.append("|---" * head.count("|") + "|" if False else "|" + "---|" * (head.count("|") - 1))
    for a in ARM_ORDER:
        if a not in small:
            continue
        e = small[a]
        enum, pack = META[a]
        row = f"| `{a}` | {enum} | {pack} | {e['variants']:,} | {e['seconds']:.0f} s |"
        if have_cpu:
            cpu = e.get("cpu_seconds") or 0.0
            mult = cpu / e["seconds"] if e["seconds"] else 0.0
            row += f" {cpu:,.0f} s ({mult:.1f}x) |" if cpu else " — |"
        if show_delta:
            row += f" {e['seconds'] - prev[a]['seconds']:+.0f} s |"
            if have_cpu:
                was_cpu = prev[a].get("cpu_seconds")
                cpu = e.get("cpu_seconds")
                row += (f" {cpu - was_cpu:+,.0f} s |" if (was_cpu and cpu) else " — |")
        row += f" {e['peak_rss_gb']:.1f} GB |"
        L.append(row)
    L.append("")
    if have_cpu:
        L.append("`CPU` is user+sys, with the multiple of wall clock beside it. It is the column "
                 "that separates work from waiting: this caller has phases that run on one thread "
                 "and phases that block on a subprocess, so a wall-clock change can come from "
                 "either doing less or waiting less, and only CPU distinguishes them. A multiple "
                 "well under `--threads` means the run spent its time parked rather than computing.")
        L.append("")
    if show_delta:
        was_builds = sorted({p.get("vg_version") for p in prev.values() if p.get("vg_version")})
        L.append("`Δ wall` is against " + (", ".join(f"`{b}`" for b in was_builds) if was_builds
                                           else "the previous refresh, whose build was not recorded")
                 + ". Read it with the repeatability note below: run-to-run variance on this "
                   "measurement is larger than most of these deltas, and the Poisson arms are the "
                   "control -- their code is untouched by any read-likelihood change, so a Δ on "
                   "those rows is the machine and not the caller.")
        L.append("")
    # The "repeatability note below" that the wall and Δ paragraphs cite. It was hand-added to the
    # pages rather than emitted here, so a regeneration dropped it and left both citations dangling.
    L.append("**Peak RSS in this table is repeatable to about ±0.35 GB, so read it accordingly.** "
             "Three back-to-back runs of one binary on chr6-4hap, identical parameters and a warm "
             "cache, gave 7.3, 6.6 and 7.0 GB -- a 0.7 GB spread on a 7 GB measurement. Differences "
             "smaller than that are not evidence of anything, and a single measurement of each of two "
             "arms cannot resolve one. Thread count matters too: the same run at `--threads 6` "
             "instead of 5 measured 8.7 GB, because the read and GBWT caches are per thread. Wall "
             "clock is worse still -- a run immediately after a full rebuild took 956 s against "
             "260 s warm, purely from page cache.")
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
        L.append("The insertion `BASEPAIR` precision above understates both callers, and the reason "
                 "is a property of the benchmark rather than of either caller.")
        L.append("")
        overlap = (f" ({sm_mb:.1f} Mb vs {sv_mb:.1f} Mb)" if sm_mb and sv_mb else "")
        L.append("**The `smvar` truth set contains no record >=50 bp** — that size class lives in the "
                 f"separate `stvar` benchmark. But the two confident regions overlap almost completely"
                 f"{overlap}. So a >=50 bp insertion called inside the small-variant confident region "
                 "has every one of its bases scored FP, however right the call is. It cannot be "
                 "scored correct.")
        L.append("")
        n_rl = count_big_insertions(res / "readlik.vcf.gz")
        n_pz = count_big_insertions(res / "poisson-z.vcf.gz")
        if n_rl is not None:
            L.append(f"{n_rl:,} `readlik` calls carry an insertion allele of {BIG_INS} bp or more"
                     + (f", against {n_pz:,} from `poisson-z`" if n_pz is not None else "")
                     + ". Every base of those alleles inside the small-variant confident region is "
                     "scored FP, and the size-matched control below measures what that does to each "
                     "caller's precision.")
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
        # Computed from this contig's own numbers, direction included. The old sentence read "the gap
        # collapses" and argued that the unrestricted figure understates readlik's lead; on the
        # 2026-09-23 build the restriction raises poisson-z more than readlik, so the unrestricted
        # figure OVERSTATES the lead. Only data can say which way round it is.
        def _bp(src, arm, key):
            return num(pick(src.get(arm, {}).get("metrics", {}).get("summary", []),
                            "BASEPAIR", "Insertion"), key)
        p0z, p0r = _bp(small, "poisson-z", "metric_precision"), _bp(small, "readlik", "metric_precision")
        p1z, p1r = _bp(sized, "sm50-poisson-z", "metric_precision"), _bp(sized, "sm50-readlik", "metric_precision")
        f0z, f0r = _bp(small, "poisson-z", "metric_f1"), _bp(small, "readlik", "metric_f1")
        f1z, f1r = _bp(sized, "sm50-poisson-z", "metric_f1"), _bp(sized, "sm50-readlik", "metric_f1")
        if None not in (p0z, p0r, p1z, p1r):
            lead0, lead1 = p0r - p0z, p1r - p1z
            if abs(lead1) < abs(lead0):
                how = "so the unrestricted comparison overstates the difference between them"
            else:
                how = "so the unrestricted comparison understates the difference between them"
            L.append(f"Restricting raises insertion BASEPAIR precision from {p0r:.4f} to {p1r:.4f} for "
                     f"`readlik` and from {p0z:.4f} to {p1z:.4f} for `poisson-z`. `readlik` minus "
                     f"`poisson-z` goes from {lead0:+.3f} to {lead1:+.3f}, {how}."
                     + (f" Insertion BASEPAIR F1 is {f1z:.4f} for `poisson-z` against {f1r:.4f} for "
                        f"`readlik` restricted, and {f0z:.4f} against {f0r:.4f} unrestricted."
                        if None not in (f0z, f0r, f1z, f1r) else ""))
        L.append("")
        judged = truvari_big_insertions(res)
        tail = ""
        if judged:
            tot, tp_, fp_n = judged
            tail = (f" Split into single alleles there are {tot:,} insertions of {BIG_INS} bp or more; "
                    f"truvari confirms **{tp_:,}** and rejects **{fp_n:,}**, and the other "
                    f"**{tot - tp_ - fp_n:,}** fall outside the SV confident region or above "
                    "truvari's 50 kb size cap, so it does not judge them. *Known bad output* lists "
                    "the largest.")
        L.append("Whether those large calls are *correct* is a separate question, and the truvari "
                 "comparison below is what answers it." + tail)
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
                 "almost nothing to match (plan §9.22)."
                 + (" The aardvark block below is kept for continuity with earlier runs."
                    if sv else ""))
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

    # Gated as a whole, prose and footnote included, on the aardvark-sv-* directories existing:
    # chr6 has none, and the prose must not describe a table that is not on the page.
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

    L.extend(ont_section(c, small, hi))

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
                 "knob with nothing to act on. On this graph the old cap of 0.1 was overriding the "
                 "mapper at exactly the sites that matter: 23.3% of reads at those sites sit at "
                 "MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising it to 0.5 "
                 "removed 94% of the excess false-positive SNVs, and the default is now "
                 f"**{hi}**. A clamp that is inert on a sparse graph is not thereby harmless.")
    else:
        L.append("The *upper* clamp (`--mismap-max`) looks inert on this graph, because it binds "
                 "only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% "
                 "at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the "
                 "old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% "
                 "of reads at those sites sit at MAPQ 1, meaning p(wrong) = 0.79, and were being "
                 "told 0.1. Raising it to 0.5 removed 94% of the excess false-positive SNVs, and the "
                 f"default is now **{hi}**. A clamp that is inert on a sparse graph is not thereby "
                 "harmless.")
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
                            ("readlik", f"**floor {lo}, cap {hi} (current defaults)**", small),
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
    elif is_chr20 and not is_hap32:
        # This is the page the history belongs on, so there is nowhere to point: say what is missing.
        L.append("Only the current row is available: the preserved old-default arms "
                 "(`arms.floor-1e-8.json`, `arms.readlik.json`) are not in this results directory. "
                 "The full grids are in plan §9.20-§9.21.")
    else:
        # The old-default arms were only ever preserved for the chr20 4-haplotype run, so on any
        # other dataset this table has one row. Saying where the history is beats printing a single
        # row under a caption that calls it a comparison -- and rows from two datasets must not
        # share a table, which is the failure the one-build-per-matrix rule exists to prevent.
        # The link is fixed to chr20: pointing at tier2-{c}-4hap-results.md made the chr6 4-hap
        # page cite itself.
        L.append("Only the current row is available here: the preserved old-default arms "
                 "(`arms.floor-1e-8.json`, `arms.readlik.json`) exist for the chr20 4-haplotype run "
                 "alone, so the before-and-after is on "
                 "[tier2-chr20-4hap-results.md](tier2-chr20-4hap-results.md). Mixing rows from two "
                 "datasets into one table is exactly what the one-build-per-matrix rule forbids. The "
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

    census = insertion_census(res / "readlik.vcf.gz")
    giants = sorted((r for r in (census or {}).get("ins", []) if r[2] >= GIANT_INS),
                    key=lambda r: -r[2])
    L.append("## Known bad output")
    L.append("")
    L.append("Neither benchmark scores these, so they appear in no metric on this page. They are "
             "recorded because they are plainly wrong and would mislead anyone reading the VCF.")
    L.append("")
    if census and giants:
        import statistics
        med, gdp = census["median_dp"], statistics.median(r[4] for r in giants if r[4] is not None)
        n_alleles = sum(sum(1 for x in r[7] if x >= GIANT_INS) for r in giants)
        # A giant that starts inside another giant's REF span is not an independent event (checked by
        # sequence, some are the same insertion written twice); count them rather than present the
        # list as that many events.
        nested = sum(1 for r in giants if any(o is not r and o[0] == r[0] and o[1] < r[1] < o[1] + o[6]
                                              for o in giants))
        fmt_n = lambda x: f"{x:,.0f}" if float(x).is_integer() else f"{x:,.1f}"
        L.append(f"`readlik` calls **{len(giants):,} records carrying an insertion of "
                 f"{GIANT_INS // 1000} kb or more** on {c}"
                 + (f" ({n_alleles:,} such alleles)" if n_alleles != len(giants) else "")
                 + f", with median DP **{fmt_n(gdp)}** against a median of **{fmt_n(med)}** over all "
                 "of the contig's records."
                 + (f" {nested:,} of them start inside another giant's reference span, so the records "
                    "overstate the number of independent events." if nested else "")
                 + " Length is the net insertion, ALT minus REF. The five largest:")
        L.append("")
        L.append("| position | net insertion | REF length | GT | DP | GQ |")
        L.append("|---|---|---|---|---|---|")
        for chrom, pos, n, gt, dp, gq, rl, _ in giants[:5]:
            L.append(f"| {chrom}:{pos:,} | {n:,} bp | {rl:,} bp | `{gt.replace('|', chr(92) + '|')}` | "
                     f"{'—' if dp is None else f'{dp:,}'} | {'—' if gq is None else gq} |")
        L.append("")
        xd = max_expected_depth(res / "poisson.vcf.gz")
        bins = []
        for lo_, hi_, lab in INS_BINS:
            d = [r[4] for r in census["ins"] if r[4] is not None and r[2] >= lo_
                 and (hi_ is None or r[2] <= hi_)]
            if d:
                bins.append(f"{statistics.median(d):,.0f} for {lab}")
        # Quoted only when it makes the point: on chr6 the Poisson caller's own expected depth
        # reaches 26,019 somewhere, above the giants' median, and then the maximum shows nothing.
        L.append((f"The Poisson caller's expected depth (`XD`) never exceeds **{xd:,}** anywhere on "
                  f"{c}. " if xd and xd < gdp else "")
                 + "Median DP by called insertion length: " + ", ".join(bins) + "."
                 + (" So the giants are collapsed-repeat pile-ups, not haplotypes."
                    if gdp > 3 * med else ""))
        L.append("")
    elif census:
        L.append(f"`readlik` calls no insertion allele of {GIANT_INS // 1000} kb or more on {c}.")
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
    sat = sum(1 for r in giants if r[5] is not None and r[5] >= 256) if census else 0
    L.append("What shipped instead attacks the same blindness from the other side: **GQ is scaled "
             "by the fraction of reads the called genotype explains**, which lowers the quality of a "
             "pile-up the call does not account for."
             + (f" It does not reach all of them: {sat:,} of the {len(giants):,} giants above still "
                "carry GQ 256." if sat else "")
             + " The giants remain output that no metric charges for, and they should be fixed "
             "because they are wrong, not because they cost a score. See "
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
    L.extend(mixture_at_long_reads())

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
