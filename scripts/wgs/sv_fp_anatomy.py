#!/usr/bin/env python3
"""Take vg's structural-variant false positives apart, allele by allele.

The question this exists to answer is whether the caller is still writing *one* long
insertion or deletion where the sample really carries *several* small changes. That was the
pathology symbolic-allele nested calling was built to remove, and a summary F1 cannot say
whether it is gone: a record that bundles four 200 bp changes and a record that carries one
genuine 800 bp deletion are the same row in a size histogram.

So each false positive is decomposed rather than counted. REF and the called ALT are stripped
of their common prefix and suffix, and what is left is aligned; runs of exact match at least
--min-anchor long separate one change from the next. A record whose core resolves to a single
change is one variant, honestly written, right or wrong. A record that resolves to several is
several variants written as one, and the interesting sub-case is the one where *none* of the
individual changes reaches 50 bp -- such a record is scored as structural only because of what
it was bundled with, and both the false positive and the missed small variants underneath it
are artefacts of the writing rather than of the evidence.

The anchor length is the one judgement call. Too short and chance matches inside a tandem
repeat split a single event into fragments; too long and genuinely separate changes merge. 30 bp
is well past the point where an exact match between two unrelated sequences is plausible, and
well short of the 50 bp structural threshold, so a "change" here can never be an artefact of
splitting one structural event in two.

Run with --baseline to decompose a second call set the same way; the difference between the two
populations is what says whether nested calling did what it was meant to.
"""

from __future__ import annotations

import argparse
import difflib
import json
import subprocess
from collections import defaultdict
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]

# The largest core alignment worth attempting. SequenceMatcher is quadratic, and a handful of
# multi-kilobase records would otherwise dominate the runtime. Skipped records are reported
# rather than dropped, because a silent cap reads as "we looked at everything".
MAX_CELLS = 400_000_000

BANDS = [(50, 100), (100, 300), (300, 700), (700, 2000), (2000, 10000), (10000, 10 ** 9)]

FMT = "%CHROM\t%POS\t%FILTER\t%REF\t%ALT\t%INFO/PctSeqSimilarity\t%INFO/SizeDiff\n"


def band(n: int) -> str:
    for lo, hi in BANDS:
        if lo <= n < hi:
            return f"{lo}-{hi}" if hi < 10 ** 9 else f"{lo}+"
    return "<50"


def query(path: Path, fmt: str) -> list[list[str]]:
    if not path.exists():
        return []
    r = subprocess.run(["bcftools", "query", "-f", fmt, str(path)],
                       capture_output=True, text=True)
    return [ln.split("\t") for ln in r.stdout.splitlines() if ln]


def trim(ref: str, alt: str) -> tuple[str, str]:
    """Strip the common prefix and suffix. What remains is the part that actually differs."""
    i, n = 0, min(len(ref), len(alt))
    while i < n and ref[i] == alt[i]:
        i += 1
    j = 0
    while j < n - i and ref[len(ref) - 1 - j] == alt[len(alt) - 1 - j]:
        j += 1
    return ref[i:len(ref) - j], alt[i:len(alt) - j]


def changes(ref: str, alt: str, min_anchor: int) -> tuple[list[int], bool]:
    """Independent changes between REF and ALT, as their signed length deltas.

    Returns (deltas, decomposed). `decomposed` is False when the pair was too large to align,
    in which case the record is counted but not classified.
    """
    cr, ca = trim(ref, alt)
    if not cr or not ca:
        return ([len(ca) - len(cr)] if (cr or ca) else []), True
    if len(cr) * len(ca) > MAX_CELLS:
        return [len(ca) - len(cr)], False

    sm = difflib.SequenceMatcher(None, cr, ca, autojunk=False)
    out: list[int] = []
    open_r = open_a = 0
    pending = False
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            if i2 - i1 >= min_anchor:
                # A long exact match separates one change from the next.
                if pending:
                    out.append(open_a - open_r)
                    open_r = open_a = 0
                    pending = False
            else:
                # Too short to separate anything: it is interior to the change being built.
                if pending:
                    open_r += i2 - i1
                    open_a += j2 - j1
        else:
            open_r += i2 - i1
            open_a += j2 - j1
            pending = True
    if pending:
        out.append(open_a - open_r)
    return out, True


def classify(ref: str, alt: str, min_anchor: int) -> dict:
    deltas, ok = changes(ref, alt, min_anchor)
    total = len(alt) - len(ref)
    largest = max((abs(d) for d in deltas), default=0)
    return {
        "n_changes": len(deltas),
        "largest": largest,
        "total_delta": total,
        "span": max(len(ref), len(alt)),
        "decomposed": ok,
        # Scored structural, but no single change inside it is: the record is structural only
        # by virtue of what it was written together with.
        "no_dominant": ok and len(deltas) >= 2 and largest < 50,
        "bundled": ok and len(deltas) >= 2,
    }


def sv_type(ref: str, alt: str) -> str:
    d = len(alt) - len(ref)
    return "SUB" if d == 0 else ("INS" if d > 0 else "DEL")


def collect(score: Path, min_anchor: int) -> list[dict]:
    out = []
    for c in AUTOSOMES:
        for chrom, pos, filt, ref, alt, seqsim, sizediff in query(
                score / f"{c}.truvari" / "fp.vcf.gz", FMT):
            # truvari's fp.vcf.gz is biallelic-split, so ALT is a single allele.
            rec = classify(ref, alt, min_anchor)
            rec.update(chrom=chrom, pos=int(pos), filter=filt,
                       type=sv_type(ref, alt),
                       size=max(abs(len(alt) - len(ref)), 0) or max(len(ref), len(alt)),
                       seqsim=None if seqsim in (".", "") else float(seqsim))
            rec["band"] = band(rec["size"])
            out.append(rec)
    return out


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f}%" if b else "--"


def report(name: str, recs: list[dict], L: list[str]) -> None:
    n = len(recs)
    skipped = sum(1 for r in recs if not r["decomposed"])
    dec = [r for r in recs if r["decomposed"]]
    bundled = [r for r in dec if r["bundled"]]
    nodom = [r for r in dec if r["no_dominant"]]
    L.append(f"### {name}")
    L.append("")
    L.append(f"- false positives decomposed: **{len(dec):,}** of {n:,} "
             f"({skipped} too large to align, counted but unclassified)")
    L.append(f"- one change: **{len(dec) - len(bundled):,}** ({pct(len(dec) - len(bundled), len(dec))})")
    L.append(f"- several changes: **{len(bundled):,}** ({pct(len(bundled), len(dec))})")
    L.append(f"- several changes, **none of them 50 bp or larger**: "
             f"**{len(nodom):,}** ({pct(len(nodom), len(dec))})")
    L.append("")
    L.append("| type | size | records | several changes | none >=50 bp |")
    L.append("|---|---|---|---|---|")
    cells = defaultdict(lambda: [0, 0, 0])
    for r in dec:
        c = cells[(r["type"], r["band"])]
        c[0] += 1
        c[1] += r["bundled"]
        c[2] += r["no_dominant"]
    for (t, b), (tot, bun, nd) in sorted(cells.items(), key=lambda kv: -kv[1][0]):
        if tot < 20:
            continue
        L.append(f"| {t} | {b} | {tot:,} | {bun:,} ({pct(bun, tot)}) | {nd:,} ({pct(nd, tot)}) |")
    L.append("")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--score", default="work/wgs-nested8/score")
    ap.add_argument("--baseline", default=None,
                    help="a second score dir to decompose the same way, for comparison")
    ap.add_argument("--pangenie", default="work/pangenie/score")
    ap.add_argument("--min-anchor", type=int, default=30)
    ap.add_argument("--out", default="docs/sv-fp-anatomy.md")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    L: list[str] = []
    L.append("# What vg's structural-variant false positives are made of")
    L.append("")
    L.append("Generated by `scripts/wgs/sv_fp_anatomy.py`. Autosomes only, truth T2T-Q100.")
    L.append("")
    L.append(f"Each false positive is stripped of its common prefix and suffix with the reference "
             f"and what is left is aligned; a run of {args.min_anchor} bp or more of exact match "
             f"separates one change from the next. A record resolving to several changes is "
             f"several variants written as one. Where none of those changes reaches 50 bp, the "
             f"record is scored as structural only because of the bundling.")
    L.append("")

    sets = {"vg call, current default": Path(args.score)}
    if args.baseline:
        sets["vg call, --no-nested baseline"] = Path(args.baseline)
    if args.pangenie:
        sets["PanGenie"] = Path(args.pangenie)

    dumped = {}
    for name, path in sets.items():
        recs = collect(path, args.min_anchor)
        dumped[name] = recs
        report(name, recs, L)

    Path(args.out).write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"\nwrote {args.out}")
    if args.json:
        Path(args.json).write_text(json.dumps(
            {k: [{kk: vv for kk, vv in r.items()} for r in v] for k, v in dumped.items()}))
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
