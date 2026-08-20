#!/usr/bin/env python3
"""For the vg-only FNs that DO have a comparable call nearby, ask why truvari refused the match.

sv_fn_mechanism.py finds that 33.8% of the truth SVs vg misses and PanGenie finds have a vg record
of comparable size within 100 bp. Truvari's default refdist is 500 bp, so proximity was not the
obstacle. Three things can still block the match, and they call for different work:

  (a) The call was already matched to a *different* truth variant. Truvari's matching is one-to-one,
      so where several truth SVs cluster, one call can only satisfy one of them and the rest are
      false negatives however good the call is. This is a truth-representation artefact, not a
      caller error -- and it is visible directly, because such a call appears in tp-comp.vcf.gz.
  (b) The call is unmatched -- it is in fp.vcf.gz -- so truvari compared the two and rejected them
      on sequence or size similarity. That is a real disagreement about what the variant is.
  (c) The call is in neither, which would mean it fell outside the scored region.

The chr8 example from the previous script is case (a) in miniature: two truth deletions at 62,875
and 62,942 both point at one vg record at 62,860.

Distinguishing them matters because (a) is not fixable in the caller and should be excluded from any
statement about how much recall is recoverable.
"""

from __future__ import annotations

import argparse
import bisect
import subprocess
from collections import defaultdict
from pathlib import Path

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
WINDOW = 100
CLUSTER = 500   # truvari's default refdist: another truth SV this close competes for the same call


def query(path: str, fmt: str, region: str | None = None) -> list[str]:
    cmd = ["bcftools", "query", "-f", fmt] + (["-r", region] if region else []) + [path]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg-score", default="work/wgs-current/score")
    ap.add_argument("--pg-score", default="work/pangenie/score")
    ap.add_argument("--vg-vcf", default="work/wgs-current/HG002.vcf.gz")
    ap.add_argument("--out", default="/dev/stdout")
    args = ap.parse_args()

    tally = defaultdict(int)
    clustered = defaultdict(int)

    for c in AUTOSOMES:
        def truthset(d):
            out = {}
            for ln in query(f"{d}/{c}.truvari/fn.vcf.gz", "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n"):
                chrom, pos, t, sl = ln.split("\t")
                out[f"{chrom}:{pos}:{t}:{sl}"] = (int(pos), abs(int(sl)))
            return out

        vg_fn, pg_fn = truthset(args.vg_score), truthset(args.pg_score)
        vg_only = {k: v for k, v in vg_fn.items() if k not in pg_fn}
        if not vg_only:
            continue

        # All truth SVs on this contig -- fn plus tp-base -- so a cluster can be detected.
        truth_pos = sorted(p for p, _ in vg_fn.values())
        for ln in query(f"{args.vg_score}/{c}.truvari/tp-base.vcf.gz", "%POS\n"):
            truth_pos.append(int(ln))
        truth_pos.sort()

        # Positions of vg calls truvari matched (tp-comp) and rejected (fp), on the call side.
        matched, rejected = [], []
        for which, dest in (("tp-comp", matched), ("fp", rejected)):
            for ln in query(f"{args.vg_score}/{c}.truvari/{which}.vcf.gz", "%POS\t%REF\t%ALT\n"):
                pos, ref, alt = ln.split("\t")
                dest.append((int(pos), max(abs(len(a) - len(ref)) for a in alt.split(","))))
        matched.sort(); rejected.sort()
        mpos = [p for p, _ in matched]
        rpos = [p for p, _ in rejected]

        for pos, n in vg_only.values():
            def near(lst, plist):
                lo = bisect.bisect_left(plist, pos - WINDOW)
                hi = bisect.bisect_right(plist, pos + WINDOW)
                return [e for e in lst[lo:hi] if e[1] >= 0.5 * n]
            in_m, in_r = near(matched, mpos), near(rejected, rpos)
            if not (in_m or in_r):
                continue          # handled by sv_fn_mechanism.py; not this question
            # Is another truth SV close enough to compete for the same call?
            lo = bisect.bisect_left(truth_pos, pos - CLUSTER)
            hi = bisect.bisect_right(truth_pos, pos + CLUSTER)
            competing = (hi - lo) > 1
            if in_m:
                tally["call was matched to a different truth variant"] += 1
                clustered["matched-elsewhere"] += competing
            else:
                tally["call was compared and rejected on similarity"] += 1
                clustered["rejected"] += competing

    total = sum(tally.values())
    L = ["# Why truvari refused a comparable nearby call", "",
         f"Restricted to vg-only false negatives with a vg record of comparable size within {WINDOW} bp",
         f"({total:,} of them). Proximity was not the obstacle -- truvari's refdist default is {CLUSTER} bp.",
         "", "| reason | n | share | of which another truth SV within 500 bp |", "|---|---|---|---|"]
    for k, v in sorted(tally.items(), key=lambda kv: -kv[1]):
        tag = "matched-elsewhere" if k.startswith("call was matched") else "rejected"
        L.append(f"| {k} | {v:,} | {100*v/total:.1f}% | {clustered[tag]:,} ({100*clustered[tag]/v:.0f}%) |")
    L += ["",
          "A call matched to a *different* truth variant is not a caller error at all. Truvari matches",
          "one-to-one, so where truth SVs cluster, one call can satisfy only one of them and the",
          "remainder are false negatives no matter how good the call is. That share of the recall gap",
          "is not recoverable by changing the model, and should be excluded before quoting how much",
          "is.",
          "",
          "A call compared and rejected on similarity is a genuine disagreement about the variant's",
          "sequence or size, and is the population where a better model or a better traversal could",
          "change the outcome."]
    Path(args.out).write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
