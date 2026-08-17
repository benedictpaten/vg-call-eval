#!/usr/bin/env python3
"""For truth SVs vg emitted nothing for, ask which stage dropped them.

sv_fn_mechanism.py isolates 1,230 autosomal truth SVs that PanGenie called, vg missed, and vg wrote
no record within 100 bp of. "No record" is where the trail goes cold in a normal run, because
`vg call` emits only non-reference calls: a site genotyped confidently hom-ref and a site the caller
never looked at produce identical output, namely nothing.

`-a/--genotype-snarls` emits every snarl including reference calls, which separates them. Re-called
that way, each locus falls into one of four states, and they implicate different code:

  * **no snarl** -- nothing within the match window even with reference calls emitted. The variant is
    not expressible as a snarl the decomposition found, so this is topology, not scoring.
  * **no comparable allele offered** -- a snarl is there, but none of its alleles changes length in
    the truth variant's direction by a comparable amount. The panel carries the allele, since
    PanGenie called it, so the traversal enumeration did not put it in front of the model.
  * **offered and called, truvari rejected it** -- the allele was there, vg called non-reference, and
    the comparison still refused the match. A representation disagreement, not a scoring one.
  * **offered and called hom-ref** -- the allele was there and the model preferred the reference.
    This is the only bucket a change to the likelihood can fix.

The order of those tests is load-bearing. Asking "did vg call anything non-reference nearby?" before
"was the right allele even on offer?" scores a locus where vg called some unrelated 60 bp event as
though the model had weighed the 4 kb deletion and decided against it. Enumeration is checked first
for that reason, and the first draft of this script got it wrong and reported zero enumeration
failures.

For the hom-ref bucket the quality fields are summarised too, because "the model preferred the
reference" has two very different versions: it had reads and read them wrongly, or it had almost no
reads and defaulted.
"""

from __future__ import annotations

import argparse
import bisect
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import median

WINDOW = 500       # truvari's default refdist
# An allele is "comparable" to the truth variant when its signed length change is the same direction
# and within this factor either way. Both halves matter: matching on the *largest* change at the
# snarl would let a bubble offering a 5 kb allele count as covering a 64 bp deletion, and matching
# without direction would let an insertion stand in for a deletion.
SIZE_LO, SIZE_HI = 0.5, 2.0


def query(path: str, fmt: str, region: str | None = None) -> list[str]:
    cmd = ["bcftools", "query", "-f", fmt] + (["-r", region] if region else []) + [path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True, help="TSV: contig pos type size zyg")
    ap.add_argument("--allvcf", required=True, help="VCF from a -a/--genotype-snarls run")
    ap.add_argument("--contig", required=True)
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    loci = []
    for ln in Path(args.loci).read_text().splitlines()[1:]:
        c, pos, t, size, zyg = ln.split("\t")
        if c == args.contig:
            loci.append((int(pos), t, int(size), zyg))

    # Every snarl the -a run emitted on this contig, with what it decided and how sure it was.
    recs = []
    for ln in query(args.allvcf, "%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AD\t%DR\t%GQN]\n"):
        f = ln.split("\t")
        pos, ref, alt, gt, gq, dp, ad, dr, gqn = f[:9]
        alts = [a for a in alt.split(",") if a not in (".", "")]
        recs.append({
            "pos": int(pos),
            # Signed length change of every allele on offer, so a size match can require the same
            # direction and a comparable magnitude rather than just "something big is here".
            "changes": [len(a) - len(ref) for a in alts],
            "nalt": len(alts),
            "gt": gt, "gq": gq, "dp": dp, "ad": ad, "dr": dr, "gqn": gqn,
            "nonref": any(x not in ("0", ".") for x in gt.replace("|", "/").split("/")),
        })
    recs.sort(key=lambda r: r["pos"])
    positions = [r["pos"] for r in recs]

    tally = defaultdict(int)
    by_band = defaultdict(lambda: defaultdict(int))
    homref_q = defaultdict(list)

    for pos, t, size, zyg in loci:
        lo = bisect.bisect_left(positions, pos - WINDOW)
        hi = bisect.bisect_right(positions, pos + WINDOW)
        near = recs[lo:hi]
        band = ("50-100" if size < 100 else "100-300" if size < 300 else
                "300-700" if size < 700 else "700+")
        # Same direction as the truth event, and within a factor of two of its size.
        want = -1 if t == "DEL" else 1
        def comparable(r):
            return any((c * want > 0) and (SIZE_LO * size <= abs(c) <= SIZE_HI * size)
                       for c in r["changes"])

        # Order matters, and getting it wrong hides the question. Asking "did vg call anything
        # non-reference nearby?" before "was the right allele even on offer?" classifies a locus
        # where vg called some unrelated 60 bp event as though the model had considered the 4 kb
        # deletion and made a choice. The enumeration test has to come first.
        comp = [r for r in near if comparable(r)]
        if not near:
            state = "no snarl"
        elif not comp:
            state = "no comparable allele offered"
        elif any(r["nonref"] for r in comp):
            state = "offered and called, truvari rejected it"
        else:
            state = "offered and called hom-ref"
            best = max(comp, key=lambda r: r["nalt"])
            for k in ("gq", "dp", "gqn", "dr"):
                try:
                    homref_q[k].append(float(best[k]))
                except (TypeError, ValueError):
                    pass
            homref_q["nalt"].append(best["nalt"])
            # AD is per-allele; the share landing on any non-reference allele is the direct question.
            try:
                parts = [float(x) for x in best["ad"].split(",") if x not in (".", "")]
                if len(parts) > 1 and sum(parts) > 0:
                    homref_q["alt_share"].append(sum(parts[1:]) / sum(parts))
            except ValueError:
                pass
        tally[state] += 1
        by_band[band][state] += 1

    total = sum(tally.values()) or 1
    STATES = ["no snarl", "no comparable allele offered",
              "offered and called, truvari rejected it", "offered and called hom-ref"]
    L = [f"# What happened to the SVs vg wrote nothing for ({args.contig})", "",
         f"{total} loci, re-called with `-a/--genotype-snarls` so reference calls are emitted too.",
         f"Match window {WINDOW} bp (truvari's refdist). An allele counts as comparable when it",
         f"changes length in the same direction as the truth event and by {SIZE_LO:g}x to {SIZE_HI:g}x",
         "its size.", "",
         "| what the caller did | n | share |", "|---|---|---|"]
    for s in STATES:
        L.append(f"| {s} | {tally[s]} | {100*tally[s]/total:.1f}% |")
    L += ["", "## By truth variant size", "",
          "| size | " + " | ".join(STATES) + " |",
          "|---|" + "---|" * len(STATES)]
    for b in ("50-100", "100-300", "300-700", "700+"):
        L.append(f"| {b} | " + " | ".join(str(by_band[b][s]) for s in STATES) + " |")

    if homref_q["gq"]:
        L += ["", "## The hom-ref calls: was there evidence to read?", "",
              "| quantity | median | min | max |", "|---|---|---|---|"]
        for k, label in (("dp", "DP (reads at the site)"),
                         ("nalt", "alleles offered"),
                         ("alt_share", "share of AD on a non-ref allele"),
                         ("gq", "GQ of the hom-ref call"),
                         ("gqn", "GQN"),
                         ("dr", "DR")):
            v = homref_q.get(k)
            if v:
                L.append(f"| {label} | {median(v):.3g} | {min(v):.3g} | {max(v):.3g} |")
        L += ["",
              "A hom-ref call with ample DP and a real share of AD on the alternate allele is the",
              "model mis-weighing evidence it had. One with DP near zero is a site the reads never",
              "covered, which no scoring change reaches."]
    Path(args.out).write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
