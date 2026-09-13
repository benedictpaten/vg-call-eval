#!/usr/bin/env python3
"""Decompose 1 bp indel calls by direction and reference homopolymer context.

The flat affine gap penalty is directionally biased: in a two-allele vote it demands
66.7% of reads to call a contraction where the ONT data demand 55-61%, and only 33.3%
to call an expansion where they demand 40-46%. The basecaller's own error runs the same
way (43.6% insertion against 25.5% deletion miscount at HP>=13), so the two COMPOUND for
insertions and CANCEL for deletions.

The signature is therefore an insertion excess among false positives with recall
untouched. This prints it. Usage:

    python3 scripts/tier2/indel_direction.py [--arm ont-chr20] <tag> [<tag> ...]

Each tag is an arm scored by work/ont-preset/arm.py, i.e. it has a results directory
work/<dataset>/results/aardvark-<tag>/ holding query.vcf.gz and truth.vcf.gz with
per-record BD. FN comes from the TRUTH file, never by subtracting query TPs from the
truth total -- aardvark's truth_tp and query_tp differ by design.
"""
import sys, gzip, collections
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "work/ont-preset"))
import arm  # noqa: E402

HPMIN = 5


def classify(path, want, ref):
    """(kind, hp) counter over 1 bp indel records whose BD is in `want`."""
    out = collections.Counter()
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            bd = dict(zip(c[8].split(":"), c[9].split(":"))).get("BD")
            if bd not in want:
                continue
            r, a = c[3], c[4]
            if abs(len(a) - len(r)) != 1 or min(len(r), len(a)) != 1:
                continue
            ins = len(a) > len(r)
            base = a[1] if ins else r[1]
            if arm.hp_run(int(c[1]) + 1, base) < HPMIN:
                continue
            out[("INS" if ins else "DEL", bd)] += 1
    return out


def report(tag):
    adir = arm.W / "results" / f"aardvark-{tag}"
    if not adir.exists():
        print(f"{tag}: no results dir {adir}")
        return None
    ref = arm.ref_seq()
    q = classify(adir / "query.vcf.gz", {"TP", "FP"}, ref)
    t = classify(adir / "truth.vcf.gz", {"FN"}, ref)
    rows = {}
    for kind in ("INS", "DEL"):
        tp, fp, fn = q[(kind, "TP")], q[(kind, "FP")], t[(kind, "FN")]
        rows[kind] = (tp, fp, fn,
                      tp / (tp + fp) if tp + fp else 0.0,
                      tp / (tp + fn) if tp + fn else 0.0)
    print(f"\n=== {tag} : 1 bp indels, reference homopolymer run >= {HPMIN} ===")
    print(f"{'':<6}{'TP':>7}{'FP':>7}{'FN':>7}{'precision':>11}{'recall':>9}")
    for kind in ("INS", "DEL"):
        tp, fp, fn, p, r = rows[kind]
        print(f"{kind:<6}{tp:>7}{fp:>7}{fn:>7}{p:>11.4f}{r:>9.4f}")
    ti, td = rows["INS"], rows["DEL"]
    fp_share = ti[1] / (ti[1] + td[1]) if ti[1] + td[1] else 0.0
    tp_share = ti[0] / (ti[0] + td[0]) if ti[0] + td[0] else 0.0
    print(f"  FP insertion share {fp_share:.1%} against a TP share of {tp_share:.1%}"
          f"   (excess {fp_share - tp_share:+.1%})")
    print(f"  precision gap INS vs DEL: {ti[3] - td[3]:+.4f}    recall gap: {ti[4] - td[4]:+.4f}")
    # What matching insertion precision to deletion precision would be worth.
    if td[3] > 0:
        target_fp = ti[0] * (1 - td[3]) / td[3]
        print(f"  matching INS precision to DEL would remove {ti[1] - target_fp:.0f} FPs")
    return rows


if __name__ == "__main__":
    args = sys.argv[1:]
    which = "ont-chr20"
    if args and args[0] == "--arm":
        which, args = args[1], args[2:]
    if not args:
        sys.exit(__doc__)
    arm.select_arm(which)
    for tag in args:
        report(tag)
