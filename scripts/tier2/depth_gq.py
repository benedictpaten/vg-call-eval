#!/usr/bin/env python3
"""Can the depth ratio go into GQ itself, the way the explained share did?

`DR` = reads observed at a site over the number the called genotype predicts. It is the
one signal the read-likelihood model is structurally blind to: `P(reads | G)` is
conditioned on the reads it was handed and never asks whether that many reads should be
there, which is why collapsed-repeat pile-ups survive it. The depth *term* acts on the
genotype; this asks the separate question of whether the same observable should also
lower the emitted *quality* at a site whose read count the call cannot account for.

Same frame as `share_gq.py`, deliberately, because that is the precedent this has to
clear: a **single fixed formula with nothing fitted** that improves the ranking in all
eight cells -- two chromosomes, two graphs, two benchmarks -- and makes none of them
worse. Per-dataset weights are not shippable and are not measured here.

Three families, and the second is the one with a mechanism behind it:

    two-sided   GQ' = GQ * exp(-a * |ln DR|)
    excess      GQ' = GQ * exp(-a * max(0, ln DR))
    phred       GQ' = min(GQ, -10 log10(excess share of reads), with a tolerance)

Why `excess` deserves separate treatment rather than being folded into `|ln DR|`: the two
tails are not symmetric in what they mean. `DR` above 1 is reads the call cannot explain,
which is evidence against the call in the same sense the unexplained share is. `DR` below
1 is a genotype claiming more sequence than the reads cover -- also suspicious, but it is
the *same* quantity the depth term already acts on when it is armed, and it is the tail
where a missed heterozygous deletion sits. Discounting the low tail can therefore punish
the calls the rest of the work is trying to recover.

Scored as a ranking -- AUC, then precision and surviving false calls at matched recall,
recall on the base side (see filter_lib.truth_counts). No train/test split: nothing here
is fitted, so there is nothing to hold out.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from share_gq import at_recall, auc  # noqa: E402

REPO = HERE.parent.parent
WORK = REPO / "work"

DATASETS = {
    "chr20-4hap": "tier2-chr20",
    "chr20-34hap": "tier2-chr20-hap32",
    "chr6-4hap": "tier2-chr6",
    "chr6-34hap": "tier2-chr6-hap32",
}

MAX_Q = 60.0


def labels(work: Path, ds: str, tag: str, kind: str) -> dict:
    bd = {}
    if kind == "truvari":
        d = work / "results" / f"truvari-{ds}-{tag}"
        for fn, lab in (("tp-comp.vcf.gz", "TP"), ("fp.vcf.gz", "FP")):
            with gzip.open(d / fn, "rt") as fh:
                for line in fh:
                    if not line.startswith("#"):
                        bd[int(line.split("\t", 2)[1])] = lab
    else:
        with gzip.open(work / "results" / f"aardvark-{ds}-{tag}" / "query.vcf.gz", "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                v = dict(zip(f[8].split(":"), f[9].split(":"))).get("BD")
                if v in ("TP", "FP"):
                    bd[int(f[1])] = v
    return bd


def truth_counts(ds: str, tag: str, kind: str) -> tuple[int, int]:
    """(base-side true positives, truth total) from this run's own score file.

    Base-side, not query-side, and the distinction matters: truvari matches several base
    records to one query record wherever the caller emitted one structural variant that
    the benchmark decomposed, so dividing query true positives by the truth total
    understates recall by nearly a factor of two. filter_lib.truth_counts has the full
    note; this reads the same numbers out of score_vcf.py's json instead of arms.json so
    it works for an arbitrary tag.
    """
    d = json.loads((WORK / "sv-atlas" / f"score-{ds}-{tag}.json").read_text())
    if kind == "truvari":
        sv = d["sv"]
        return int(sv["TP-base"]), int(sv["TP-base"]) + int(sv["FN"])
    for r in d.get("smallvar") or []:
        if (r["comparison"], r["region_label"], r["filter"], r["variant_type"]) == \
                ("GT", "ALL", "ALL", "ALL"):
            return int(r["truth_tp"]), int(r["truth_total"])
    raise SystemExit(f"no small-variant truth counts for {ds} {tag}")


def collect(work: Path, ds: str, tag: str, kind: str) -> list[dict]:
    bd = labels(work, ds, tag, kind)
    vcf = str(work / "results" / f"sweep-{ds}-{tag}.vcf.gz")
    q = subprocess.run(["bcftools", "query", "-f",
                        "%POS\t%REF\t%ALT[\t%DP\t%AD\t%GQ\t%GQI\t%DR\t%GT]\n", vcf],
                       capture_output=True, text=True)
    if q.returncode != 0:
        raise SystemExit(q.stderr.strip()[:400])
    out = []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        try:
            pos = int(f[0])
            ref, alts = f[1], f[2].split(",")
            dp = float(f[3])
            ad = [int(x) for x in f[4].split(",")]
            gq, gqi, dr = float(f[5]), float(f[6]), float(f[7])
            gt = f[8]
        except (ValueError, IndexError):
            continue
        if dp <= 0 or pos not in bd or dr < 0:
            continue
        # Largest length change among the *called* alleles, which is what decides
        # whether DR's geometry is dominated by the read length or by the allele.
        svlen = 0
        for tok in gt.replace("|", "/").split("/"):
            if tok.isdigit() and 0 < int(tok) <= len(alts):
                svlen = max(svlen, abs(len(alts[int(tok) - 1]) - len(ref)))
        out.append({"pos": pos, "label": bd[pos], "dp": dp, "gq": gq, "gqi": gqi,
                    "dr": dr, "svlen": svlen, "share": min(1.0, sum(ad) / dp)})
    return out


def forms(centre: float):
    """Candidate qualities. `centre` is where DR is treated as unremarkable.

    Passing anything but 1.0 is a diagnostic, not a candidate: a per-run empirical centre
    cannot ship, and if it were needed that would be an argument for fixing lambda rather
    than for calibrating around it.
    """
    def ln(dr):
        return math.log(max(dr, 1e-3) / centre)

    out = [("GQ (unchanged)", lambda r: r["gq"])]
    for a in (0.25, 0.5, 1.0):
        out.append((f"two-sided a={a:g}",
                    lambda r, a=a: r["gq"] * math.exp(-a * abs(ln(r["dr"])))))
    for a in (0.25, 0.5, 1.0):
        out.append((f"excess   a={a:g}",
                    lambda r, a=a: r["gq"] * math.exp(-a * max(0.0, ln(r["dr"])))))

    # Phred reading of the excess. If a site holds twice the reads the call predicts,
    # half of them are unaccounted for; treat that as an error probability. A tolerance
    # t is needed for the same reason it was needed for share -- without one, ordinary
    # Poisson scatter caps good calls -- and t is the fraction of excess forgiven.
    def phred(r, t):
        dr = max(r["dr"], 1e-3) / centre
        if dr <= 1.0:
            return r["gq"]
        excess = (1.0 - 1.0 / dr - t) / (1.0 - t)
        if excess <= 0.0:
            return r["gq"]
        return min(r["gq"], -10.0 * math.log10(min(1.0, excess)))

    for t in (0.1, 0.25, 0.5):
        out.append((f"phred    t={t:g}", lambda r, t=t: phred(r, t)))

    # Size-gated. Ungated, the two-sided discount gains on structural variants
    # everywhere and *loses* on small variants on the 34-haplotype graphs -- the same
    # sign reversal that has blocked every previous depth signal. There is a mechanism
    # for gating rather than a threshold to fit: at a SNV, lambda's geometry is
    # dominated by the read length, so DR mostly reports local coverage scatter, and
    # ranking on it adds noise to a GQ that is already a good ranker there. At a large
    # event the geometry is dominated by the allele, so DR reports whether the called
    # sequence is actually present. 50 bp is not tuned here either: it is the boundary
    # the two benchmarks already draw, and the small-variant truth set holds no record
    # above it.
    for lo in (50, 300):
        for a in (0.5, 1.0):
            out.append((f"gated>={lo} a={a:g}",
                        lambda r, a=a, lo=lo: r["gq"] * (
                            math.exp(-a * abs(ln(r["dr"]))) if r["svlen"] >= lo else 1.0)))
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tag", default="depthcount-off",
                   help="run tag under results/, e.g. depthcount-off or "
                        "depthcount-eff0.25")
    p.add_argument("--only", nargs="*")
    p.add_argument("--centre", type=float, default=1.0,
                   help="diagnostic: treat this DR as unremarkable instead of 1.0")
    args = p.parse_args()

    for ds, sub in DATASETS.items():
        if args.only and ds not in args.only:
            continue
        W = WORK / sub
        for kind, targets in (("aardvark", [0.90, 0.93]), ("truvari", [0.35, 0.42])):
            if not (W / "results" / f"{kind}-{ds}-{args.tag}").exists():
                print(f"\n=== {ds} {kind}: not scored for tag {args.tag}")
                continue
            rows = collect(W, ds, args.tag, kind)
            if not rows:
                print(f"\n=== {ds} {kind}: no labelled rows")
                continue
            base_tp, total = truth_counts(ds, args.tag, kind)
            ys = [1 if r["label"] == "TP" else 0 for r in rows]
            bench = "small variants" if kind == "aardvark" else "SVs"
            drs = sorted(r["dr"] for r in rows)
            print(f"\n=== {ds}, {bench}: {sum(ys):,} TP / {len(ys) - sum(ys):,} FP, "
                  f"median DR {drs[len(drs) // 2]:.3f} ===")
            print(f"  {'form':<18}{'AUC':>8}" +
                  "".join(f"{'P@R' + str(t):>10}{'FP':>7}" for t in targets))
            base = None
            for name, fn in forms(args.centre):
                sc = [fn(r) for r in rows]
                a = auc([s for s, y in zip(sc, ys) if y == 1],
                        [s for s, y in zip(sc, ys) if y == 0])
                pr = at_recall(sc, ys, base_tp, total, targets)
                cells = ""
                for t in targets:
                    v, fp = pr.get(t, (float("nan"), -1))
                    cells += f"{v:>10.4f}{fp:>7d}"
                mark = ""
                if base is None:
                    base = a
                else:
                    mark = f"  {a - base:+.4f}"
                print(f"  {name:<18}{a:>8.4f}" + cells + mark)


if __name__ == "__main__":
    main()
