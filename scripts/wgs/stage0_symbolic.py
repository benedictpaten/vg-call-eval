#!/usr/bin/env python3
"""Stage 0 of the nested-calling plan: would symbolic alleles actually collapse these records?

doc/nested-calling-design.md proposes rewriting each traversal of a non-leaf snarl as a symbolic
allele, replacing each excursion through a nested chain with a symbol for that chain. A traversal
whose symbolic form equals the reference traversal's is then the reference allele at that level, and
its differences descend to the nested sites that own them.

That is only worth building if the records doing the damage are non-leaf snarls. If the 4,710 bp
allele differing at three bases turns out to be a *leaf* snarl -- one large complex bubble with no
nested structure -- there is nothing to symbolise and the design is aimed at the wrong population.
This measures that before any C++ is written, from artefacts that already exist:

  * the snarl inventory, as the ID field (`>start>end`) of every record from a `-a -A` run, which
    enumerates nested snarls as well as top-level ones -- 222,797 against 192,207 on chr20;
  * the production run's records, whose `AT` INFO field already gives each allele's traversal as an
    oriented node path.

Projection: walk a traversal left to right; wherever a node opens a snarl other than this record's
own whose closing node appears later in the same traversal, replace that whole span with one symbol
and continue from the closing node. What remains is the symbolic allele.

Child-snarl granularity is used rather than chain granularity. The two agree for this test: a chain
of several child snarls symbolises to several symbols instead of one, but the nodes between them are
shared by every traversal that crosses the chain, so two traversals collapse together under one
scheme exactly when they do under the other.
"""

from __future__ import annotations

import argparse
import re
import subprocess
from collections import defaultdict
from pathlib import Path

STEP = re.compile(r"[<>]\d+")


def query(path: str, fmt: str, region: str | None = None) -> list[str]:
    cmd = ["bcftools", "query", "-f", fmt] + (["-r", region] if region else []) + [path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def steps(path: str) -> list[str]:
    """'>1>2<3' -> ['>1', '>2', '<3']"""
    return STEP.findall(path)


def node_of(step: str) -> int:
    return int(step[1:])


def symbolic(trav: list[str], own: tuple[int, int], opens: dict) -> tuple:
    """Replace each excursion into a child snarl with a single symbol.

    `opens` maps a node id to the set of node ids that close a snarl starting there. A step is the
    start of a child excursion when some snarl opens at that node, is not this record's own snarl,
    and its closing node occurs later in this traversal. The longest such span is taken, which picks
    the outermost child rather than one of its descendants.
    """
    # Where each node id appears, so "does it close later" is a lookup rather than a scan.
    later = defaultdict(list)
    for i, s in enumerate(trav):
        later[node_of(s)].append(i)

    out = []
    i = 0
    n = len(trav)
    while i < n:
        a = node_of(trav[i])
        best_j = -1
        best_end = None
        for b in opens.get(a, ()):
            if (a, b) == own or (b, a) == own:
                continue
            for j in later.get(b, ()):
                if j > i and j > best_j:
                    best_j, best_end = j, b
        if best_j > 0:
            out.append(("S", a, best_end))     # one symbol for the whole child excursion
            i = best_j                          # continue *from* the closing node, which is shared
        else:
            out.append(trav[i])
            i += 1
    return tuple(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contig", default="chr20")
    ap.add_argument("--nested-vcf", required=True, help="VCF from a -a -A run: the snarl inventory")
    ap.add_argument("--calls", required=True, help="the production run's VCF")
    ap.add_argument("--vg-score", default="work/wgs/score")
    ap.add_argument("--pg-score", default="work/pangenie/score")
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    # --- snarl inventory ------------------------------------------------------------------
    opens = defaultdict(set)
    n_snarls = 0
    for sid in query(args.nested_vcf, "%ID\n"):
        s = steps(sid)
        if len(s) == 2:
            opens[node_of(s[0])].add(node_of(s[1]))
            n_snarls += 1

    # --- truth SNVs vg missed, to cross-tabulate ------------------------------------------
    missed = set()
    pg_tp = set()
    for ln in query(f"{args.vg_score}/{args.contig}.aardvark/truth.vcf.gz",
                    "%POS\t%REF\t%ALT\t[%BD]\n"):
        f = ln.split("\t")
        if len(f) >= 4 and f[3] == "FN" and len(f[1]) == 1 and len(f[2].split(",")[0]) == 1:
            missed.add(int(f[0]))
    for ln in query(f"{args.pg_score}/{args.contig}.aardvark/truth.vcf.gz",
                    "%POS\t%REF\t%ALT\t[%BD]\n"):
        f = ln.split("\t")
        if len(f) >= 4 and f[3] == "TP":
            pg_tp.add(int(f[0]))

    # --- project every emitted record -----------------------------------------------------
    tally = defaultdict(int)
    collapse_swallow = defaultdict(int)
    examples = []

    for ln in query(args.calls, "%POS\t%ID\t%REF\t%ALT\t%INFO/AT\n"):
        f = ln.split("\t")
        if len(f) < 5:
            continue
        pos, sid, ref, alt, at = int(f[0]), f[1], f[2], f[3], f[4]
        sids = steps(sid)
        if len(sids) != 2:
            continue
        own = (node_of(sids[0]), node_of(sids[1]))
        travs = [steps(t) for t in at.split(",")]
        if len(travs) < 2:
            continue
        syms = [symbolic(t, own, opens) for t in travs]
        nonleaf = any(any(isinstance(x, tuple) for x in s) for s in syms)
        # Every ALT symbolically identical to REF: the whole record is reference at this level.
        collapses = all(s == syms[0] for s in syms[1:])

        reflen = len(ref)
        same_len = any(len(a) == reflen for a in alt.split(",")) and reflen >= 50
        swallowed_here = [p for p in range(pos + 1, pos + reflen) if p in missed] if reflen >= 50 else []

        tally["records"] += 1
        tally["non-leaf"] += nonleaf
        tally["collapse"] += collapses
        if reflen >= 50:
            tally["large records"] += 1
            tally["large & non-leaf"] += nonleaf
            tally["large & collapse"] += collapses
        if same_len:
            tally["same-length >=50bp"] += 1
            tally["same-length & non-leaf"] += nonleaf
            tally["same-length & collapse"] += collapses
        if swallowed_here:
            tally["records swallowing a missed SNV"] += 1
            tally["  of those, non-leaf"] += nonleaf
            tally["  of those, collapse"] += collapses
            collapse_swallow["swallowed SNVs total"] += len(swallowed_here)
            if collapses:
                collapse_swallow["swallowed SNVs in a collapsing record"] += len(swallowed_here)
            if nonleaf:
                collapse_swallow["swallowed SNVs in a non-leaf record"] += len(swallowed_here)
            collapse_swallow["  of those PanGenie called"] += sum(1 for p in swallowed_here if p in pg_tp)
            if len(examples) < 5 and collapses and reflen >= 200:
                examples.append((pos, reflen, len(travs), len(swallowed_here)))

    L = [f"# Stage 0: would symbolic alleles collapse the records that are losing variants? ({args.contig})",
         "",
         f"Snarl inventory from the `-a -A` run: {n_snarls:,} snarls. Production records projected:",
         f"{tally['records']:,}.",
         "",
         "A record *collapses* when every ALT's symbolic allele equals the reference traversal's, so the",
         "whole record is reference at this level and its differences belong to nested sites. A record is",
         "*non-leaf* when any allele's traversal crosses a child snarl at all -- the precondition for",
         "symbolising anything.",
         "",
         "| | n | share of its row group |", "|---|---|---|"]
    r = tally["records"] or 1
    L.append(f"| all emitted records | {tally['records']:,} | |")
    L.append(f"| ... non-leaf | {tally['non-leaf']:,} | {100*tally['non-leaf']/r:.1f}% |")
    L.append(f"| ... would collapse to reference | {tally['collapse']:,} | {100*tally['collapse']/r:.1f}% |")
    for grp, keys in (("large records (REF >=50 bp)",
                       ("large records", "large & non-leaf", "large & collapse")),
                      ("same-length alleles >=50 bp (the substitution FPs)",
                       ("same-length >=50bp", "same-length & non-leaf", "same-length & collapse")),
                      ("records swallowing a missed SNV",
                       ("records swallowing a missed SNV", "  of those, non-leaf", "  of those, collapse"))):
        tot = tally[keys[0]] or 1
        L.append(f"| **{grp}** | {tally[keys[0]]:,} | |")
        L.append(f"| ... non-leaf | {tally[keys[1]]:,} | {100*tally[keys[1]]/tot:.1f}% |")
        L.append(f"| ... would collapse | {tally[keys[2]]:,} | {100*tally[keys[2]]/tot:.1f}% |")
    L += ["", "## Missed SNVs recovered by collapsing", "", "| | n |", "|---|---|"]
    for k, v in collapse_swallow.items():
        L.append(f"| {k} | {v:,} |")
    if examples:
        L += ["", "## Examples of collapsing records", "",
              "| POS | REF bp | alleles | missed SNVs inside |", "|---|---|---|---|"]
        for p, rl, nt, ns in examples:
            L.append(f"| {p:,} | {rl:,} | {nt} | {ns} |")
    Path(args.out).write_text("\n".join(L) + "\n")
    for k in ("records", "non-leaf", "collapse", "large records", "large & non-leaf",
              "large & collapse", "same-length >=50bp", "same-length & collapse",
              "records swallowing a missed SNV", "  of those, non-leaf", "  of those, collapse"):
        print(f"  {k:34s} {tally[k]:>8,}")
    for k, v in collapse_swallow.items():
        print(f"  {k:34s} {v:>8,}")


if __name__ == "__main__":
    main()
