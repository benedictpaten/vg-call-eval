#!/usr/bin/env python3
"""Do vg's quality signals separate true structural calls from false ones, and does gating on
them help?

Two measurements over truvari's own output, autosomes only:

1. **Signal separation.** Median GQ, GQN and DR, and the heterozygous share, for the matched
   (TP-comp) and unmatched (FP) call populations. A signal that separates them is a candidate
   filter.

2. **Gate sweep.** For each threshold, how many true and false calls survive and what the SV F1
   becomes. Recall is accounted on the *base* side, as truvari scores it: gating a comp call moves
   its matched truth record from TP-base to FN rather than leaving the truth count untouched.
   Counting the comp side for both halves of F1 -- the obvious shortcut -- mixes two denominators
   and reads a percent low, because truvari matches more truth records than it has matching calls
   (13,526 against 13,302 on this arm, under `--pick ac`).

   Base and comp are joined on the **second** field of INFO/MatchId, the comp-side id. The first
   field, the base-side id, looks like the natural key and is not: where one call matches several
   truth records only one of them carries that call's base id, so joining on it drops 14 of
   chr20's 410 matched truth records and quietly understates the ungated F1 by 0.009. The ungated
   row is asserted against truvari's own summary.json to catch exactly this.

These numbers were hand-measured once and quoted in sv-residual-errors.md, which is how they went
stale: the arm was recalled and the prose was not. Run this instead.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

# (label, FORMAT field, minimum). None means "no gate", the anchor row.
GATES = [("none", None, 0.0),
         ("DR >= 0.3", "DR", 0.3), ("DR >= 0.5", "DR", 0.5),
         ("GQ >= 3", "GQ", 3.0), ("GQ >= 10", "GQ", 10.0), ("GQ >= 20", "GQ", 20.0),
         ("GQN >= 0.02", "GQN", 0.02), ("GQN >= 0.05", "GQN", 0.05)]


def records(vcf: Path):
    """One tuple per record: MatchId, GT, GQ, GQN, DR, PctSeqSimilarity, size change in bp.
    None for any field written as '.'.

    PctSeqSimilarity is truvari's, and its presence is the discriminator between the two
    false-positive populations: annotated means truvari found a candidate truth SV in range and
    rejected the call on similarity; absent means there was nothing nearby to compare against."""
    if not vcf.exists():
        return []
    q = subprocess.run(["bcftools", "query", "-f",
                        "%INFO/MatchId\t%INFO/PctSeqSimilarity\t%REF\t%ALT\t"
                        "[%GT\t%GQ\t%GQN\t%DR]\n", str(vcf)],
                       capture_output=True, text=True)
    out = []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 8:
            continue

        def num(x):
            try:
                return float(x)
            except ValueError:
                return None

        # Size change of the longest ALT against REF -- enough to bucket the long tail, and it
        # does not need the graph.
        alts = [a for a in f[3].split(",") if a not in (".", "*")]
        size = max((abs(len(a) - len(f[2])) for a in alts), default=0)
        out.append((f[0], f[4], num(f[5]), num(f[6]), num(f[7]), num(f[1]), size))
    return out


def match_ids(vcf: Path):
    """INFO/MatchId only. The truth-side VCF carries none of the caller's FORMAT fields, so the
    record reader above filters every one of its lines away."""
    if not vcf.exists():
        return []
    q = subprocess.run(["bcftools", "query", "-f", "%INFO/MatchId\n", str(vcf)],
                       capture_output=True, text=True)
    return [l.strip() for l in q.stdout.splitlines()]


def fp_classes(vcf: Path):
    """(compared-and-rejected?, size change) per record, using no FORMAT field at all.

    The reader above needs vg's GQ/GQN/DR and so returns nothing for a call set that does not
    write them -- which is every other tool's, and is why the comparison below silently read zero
    records the first time."""
    if not vcf.exists():
        return []
    q = subprocess.run(["bcftools", "query", "-f",
                        "%INFO/PctSeqSimilarity\t%REF\t%ALT\n", str(vcf)],
                       capture_output=True, text=True)
    out = []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        alts = [a for a in f[2].split(",") if a not in (".", "*")]
        out.append((f[0] != ".", max((abs(len(a) - len(f[1])) for a in alts), default=0)))
    return out


def median(xs):
    xs = sorted(x for x in xs if x is not None)
    return xs[len(xs) // 2] if xs else float("nan")


def het(gts):
    """Share of genotypes that are heterozygous, over those whose zygosity is actually known.

    A half-called genotype like `.|1` is dropped rather than counted het: its two alleles do
    differ as strings, so the naive test calls it heterozygous when what it really says is that
    one strand was not determined."""
    known = []
    for g in gts:
        a = g.replace("|", "/").split("/")
        if a and all(x != "." for x in a):
            known.append(a)
    h = sum(1 for a in known if len(set(a)) > 1)
    return 100.0 * h / len(known) if known else float("nan")


def comp_key(match_id: str) -> str:
    """The comp-side half of truvari's MatchId pair -- the half that joins base to comp exactly."""
    f = match_id.split(",")
    return f[1] if len(f) > 1 else ""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score", default="work/wgs-single/score")
    ap.add_argument("--pangenie", default="work/pangenie/score",
                    help="second call set to decompose the same way, so the false-positive excess "
                         "can be attributed to one population or the other; skipped if absent")
    ap.add_argument("--contigs", nargs="*", default=AUTOSOMES)
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    score = Path(args.score)
    tp, fp = [], []
    truth_total = truvari_tp_base = truvari_fp = 0
    # One entry per matched truth record, holding the calls that matched it: gating them all is
    # what turns that truth record back into a false negative.
    truth_calls = []
    for c in args.contigs:
        d = score / f"{c}.truvari"
        summary = d / "summary.json"
        if not summary.exists():
            continue
        js = json.loads(summary.read_text())
        truth_total += int(js.get("TP-base", 0)) + int(js.get("FN", 0))
        truvari_tp_base += int(js.get("TP-base", 0))
        truvari_fp += int(js.get("FP", 0))
        by_comp = {}
        for r in records(d / "tp-comp.vcf.gz"):
            tp.append(r)
            by_comp.setdefault(comp_key(r[0]), []).append(r)
        for mid in match_ids(d / "tp-base.vcf.gz"):
            truth_calls.append(by_comp.get(comp_key(mid), []))
        fp.extend(records(d / "fp.vcf.gz"))

    # The join is only trustworthy if it accounts for every truth record truvari matched.
    joined = sum(1 for calls in truth_calls if calls)
    if joined != truvari_tp_base:
        raise SystemExit(f"MatchId join accounted for {joined} of {truvari_tp_base} matched truth "
                         f"records; the gate sweep would understate recall")

    L = ["# Do vg's quality signals separate true structural calls from false ones?", "",
         "Generated by `scripts/wgs/sv_quality_gates.py`. Autosomes only, SVs >=50 bp, truth",
         "T2T-Q100, over truvari's own matched and unmatched record sets.", "",
         f"Matched calls (TP-comp) **{len(tp):,}**, unmatched (FP) **{len(fp):,}**, "
         f"truth SVs **{truth_total:,}**.", "",
         "## Signal separation", "",
         "| | matched | unmatched |", "|---|---|---|"]
    for label, idx in (("median GQ", 2), ("median GQN", 3), ("median DR", 4)):
        a, b = median([r[idx] for r in tp]), median([r[idx] for r in fp])
        fmt = "{:.0f}" if idx == 2 else "{:.3f}" if idx == 3 else "{:.2f}"
        L.append(f"| {label} | {fmt.format(a)} | {fmt.format(b)} |")
    lo_tp = 100.0 * sum(1 for r in tp if r[4] is not None and r[4] < 0.75) / max(len(tp), 1)
    lo_fp = 100.0 * sum(1 for r in fp if r[4] is not None and r[4] < 0.75) / max(len(fp), 1)
    L += [f"| DR < 0.75 | {lo_tp:.1f}% | {lo_fp:.1f}% |",
          f"| heterozygous | {het([r[1] for r in tp]):.1f}% | {het([r[1] for r in fp]):.1f}% |",
          "", "## Gate sweep", "",
          "| gate | TP kept | FP kept | SV F1 |", "|---|---|---|---|"]

    rows = []
    for label, field, thresh in GATES:
        idx = {"GQ": 2, "GQN": 3, "DR": 4}.get(field)

        def keeps(r):
            if idx is None:
                return True
            v = r[idx]
            return v is not None and v >= thresh

        kept_tp = sum(1 for r in tp if keeps(r))
        kept_fp = sum(1 for r in fp if keeps(r))
        # Base-side recall: a truth SV stays a TP only while some call matching it survives.
        base_tp = sum(1 for calls in truth_calls if any(keeps(r) for r in calls))
        fn = truth_total - base_tp
        f1 = 2 * base_tp / (2 * base_tp + kept_fp + fn) if base_tp else float("nan")
        rows.append((f1, label, kept_tp, kept_fp))
    anchor = next(r for r in rows if r[1] == "none")
    for f1, label, kept_tp, kept_fp in [anchor] + sorted(
            (r for r in rows if r[1] != "none"), key=lambda r: -r[0]):
        mark = "**" if label == "none" else ""
        L.append(f"| {label} | {kept_tp:,} | {kept_fp:,} | {mark}{f1:.4f}{mark} |")

    # The two false-positive populations, which are not alike: a call truvari compared against a
    # real truth SV and rejected is a near miss, and a call with nothing nearby is a different
    # failure. Quoted in sv-residual-errors.md and pangenie-comparison.md.
    rejected = [r for r in fp if r[5] is not None]
    phantom = [r for r in fp if r[5] is None]
    L += ["", "## The two false-positive populations", "",
          "A call truvari compared against a real truth SV in range and rejected on similarity is a",
          "near miss. A call with no truth SV nearby to compare against is a different failure, and",
          "the presence of truvari's `PctSeqSimilarity` annotation separates them.", "",
          "| | no truth SV nearby | compared and rejected |", "|---|---|---|",
          f"| n | {len(phantom):,} | {len(rejected):,} |",
          f"| share | {100 * len(phantom) / max(len(fp), 1):.1f}% | "
          f"{100 * len(rejected) / max(len(fp), 1):.1f}% |",
          f"| median GQ | {median([r[2] for r in phantom]):.0f} | "
          f"{median([r[2] for r in rejected]):.0f} |",
          f"| median DR | {median([r[4] for r in phantom]):.2f} | "
          f"{median([r[4] for r in rejected]):.2f} |",
          f"| heterozygous | {het([r[1] for r in phantom]):.1f}% | "
          f"{het([r[1] for r in rejected]):.1f}% |",
          f"| 700 bp and over | "
          f"{100 * sum(1 for r in phantom if r[6] >= 700) / max(len(phantom), 1):.1f}% | "
          f"{100 * sum(1 for r in rejected if r[6] >= 700) / max(len(rejected), 1):.1f}% |",
          ""]

    # Where the false-positive excess over the other tool actually sits. Quoted in
    # pangenie-comparison.md, which had it hand-written and stale.
    pg = Path(args.pangenie)
    if (pg / f"{args.contigs[0]}.truvari" / "fp.vcf.gz").exists():
        pg_fp = []
        for c in args.contigs:
            pg_fp.extend(fp_classes(pg / f"{c}.truvari" / "fp.vcf.gz"))
        pg_phantom = [r for r in pg_fp if not r[0]]
        if not pg_fp:
            raise SystemExit(f"read no false positives from {pg}; the comparison would be vacuous")
        excess = len(fp) - len(pg_fp)
        phantom_excess = len(phantom) - len(pg_phantom)
        u300 = 100 * sum(1 for r in phantom if r[6] < 300) / max(len(phantom), 1)
        L += ["## Where the excess over the other call set sits", "",
              f"| | this arm | {pg.parent.name} |", "|---|---|---|",
              f"| false positives | {len(fp):,} | {len(pg_fp):,} |",
              f"| no truth SV nearby | {len(phantom):,} | {len(pg_phantom):,} |",
              f"| compared and rejected | {len(rejected):,} | "
              f"{len(pg_fp) - len(pg_phantom):,} |", "",
              f"Excess {excess:,} false positives, of which {phantom_excess:,} "
              f"({100 * phantom_excess / excess:.0f}%) are calls with no truth SV in reach rather "
              f"than near misses.",
              f"{u300:.1f}% of those are under 300 bp.", ""]

    # The ungated row must be truvari's own F1, or the accounting above is wrong somewhere else.
    truvari_f1 = (2 * truvari_tp_base
                  / (2 * truvari_tp_base + truvari_fp + (truth_total - truvari_tp_base)))
    if abs(anchor[0] - truvari_f1) > 5e-5:
        raise SystemExit(f"ungated F1 {anchor[0]:.4f} does not reproduce truvari's "
                         f"{truvari_f1:.4f}")

    best = max(r[0] for r in rows if r[1] != "none")
    L += ["", f"Best gate {best:.4f} against {anchor[0]:.4f} ungated, which reproduces truvari's "
              f"own SV F1: {'no gate helps' if best <= anchor[0] else 'a gate helps'}.", ""]
    Path(args.out).write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
