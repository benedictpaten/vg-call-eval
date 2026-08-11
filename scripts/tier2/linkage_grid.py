#!/usr/bin/env python3
"""Stage 2 of the linkage HMM: search the transition weight against the distance scale.

`w_t` and `scale` are crossed rather than swept one at a time, because they interact by
construction. A longer scale lowers the switch probability at a given gap, which strengthens the
transition; a larger `w_t` tempers the switch probability downward, which also strengthens it.
Coordinate descent over that surface finds whichever corner it started nearest -- the same
argument that put the depth weight and the mismapping floor on one grid.

The block-switch rate `beta` used to be the second axis, and it is gone. Two findings retired it,
both recorded in `subchain_linkage.py` and in the `scale` comment in `linkage_model.hpp`. It was
not a second axis: smeared over `gap / block_length` it is exactly a shorter `scale`, so
`beta = 0.57` at 10 kb blocks was `--linkage-scale 5423` and the grid was measuring one axis
twice. And the premise failed on its own terms -- panel linkage across a real subchain boundary is
weaker by only 0.008 NMI at the gaps where adjacent calls actually sit (z = 1.1, gap-matched with
a permutation control), because the sampler switches to another assembly in the same panel and
human haplotypes mostly agree.

`--linkage-freq-prior` stays at its default of 0 and is *not* an axis here. Over a panel that
haplotype sampling selected against these same reads, panel allele frequency is already
conditioned on the data being genotyped, so using it as a prior counts the same evidence twice.
Measured offline it is worth about half the total gain, which is precisely why it should not be
switched on without a graph whose panel was chosen independently.

Reported per point: both benchmarks, the heterozygous class breakdown, and the two numbers this
model is judged on -- the share of undecided genotypes it moves, and the share of confident ones
it moves. The second is a diagnostic rather than a veto: it was written as a 0.1% budget, but
that budget was invented here rather than measured, and at `w = 2` it flags points whose F1 is
better on every cell of every dataset. It tracks something real -- it first crosses 0.1% at the
same weight where thin-panel SNV F1 starts to regress -- so it is still printed, and it is not
allowed to disqualify a measured outcome on its own.
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import param_sweep as ps  # noqa: E402

WORK = ps.WORK


def genotype_shift(baseline_vcf: Path, arm_vcf: Path) -> dict:
    """How many genotypes moved, split by how sure the per-site model was.

    Keyed on GQI, the per-site likelihood ratio, which the linkage pass leaves untouched -- so the
    strata are the same before and after and the comparison is not circular.
    """
    def load(path: Path) -> dict:
        out = {}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                s = dict(zip(f[8].split(":"), f[9].split(":")))
                try:
                    gqi = float(s.get("GQI", "nan"))
                except ValueError:
                    gqi = float("nan")
                out[(f[0], f[1], f[2])] = (s.get("GT"), gqi)
        return out

    base, arm = load(baseline_vcf), load(arm_vcf)
    totals, changed = collections.Counter(), collections.Counter()
    for key in set(base) & set(arm):
        gt_b, gqi = base[key]
        gt_a, _ = arm[key]
        if gqi != gqi:
            continue
        stratum = "lo" if gqi < 10 else ("mid" if gqi < 40 else "hi")
        totals[stratum] += 1
        if gt_a != gt_b:
            changed[stratum] += 1
    return {"totals": dict(totals), "changed": dict(changed),
            "lo_pct": 100.0 * changed["lo"] / max(totals["lo"], 1),
            "hi_pct": 100.0 * changed["hi"] / max(totals["hi"], 1)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["chr20-34hap"])
    ap.add_argument("--weights", nargs="+", default=["0", "0.5", "1", "2"])
    ap.add_argument("--scales", nargs="+", default=["10000"])
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    ap.add_argument("--threads", type=int, default=5)
    args = ap.parse_args()

    rows = []
    for ds in args.datasets:
        # The w = 0 control is run, not carried over: it must come from this binary and these
        # flags, or the comparison is against a different experiment.
        baseline_vcf = None
        for scale in args.scales:
            for w in args.weights:
                # Scale is inert at zero weight -- the transition is uniform -- so one control
                # per dataset, not one per scale.
                if float(w) == 0.0 and scale != args.scales[0]:
                    continue
                params = {"linkage-scale": scale}
                if float(w) > 0:
                    params["linkage-weight"] = w
                # Scale belongs in the tag. It was fixed when this was written, so leaving it out
                # was harmless; with it swept, two scales would collide on one cached VCF and the
                # second point would silently report the first one's numbers.
                tag = (f"{ds}-lgrid-w{w}" if float(w) == 0.0
                       else f"{ds}-lgrid-w{w}-s{scale}")
                print(f"=== {ds} w={w} scale={scale}", flush=True)
                vcf = ps.call(args.vg, ds, tag, params, args.threads)
                sc = ps.score(vcf, ds, tag, args.threads)
                if float(w) == 0.0:
                    baseline_vcf = vcf
                shift = genotype_shift(baseline_vcf, vcf) if baseline_vcf else {}
                rows.append((ds, w, scale, sc, shift))

    hdr = (f"{'dataset':12s} {'w':>4s} {'scale':>6s} {'SV F1':>7s} {'SVrec':>7s} {'SVprec':>7s} "
           f"{'smallGT':>8s} {'SNV':>7s} {'hDEL1k':>7s} {'moved@GQI<10':>13s} "
           f"{'moved@GQI>=40':>14s}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for ds, w, scale, s, shift in rows:
        sv = s.get("sv") or {}
        hi = shift.get("hi_pct")
        flag = "" if hi is None or hi <= 0.1 else "  OVER BUDGET"
        print(f"{ds:12s} {w:>4s} {scale:>6s} {sv.get('f1', 0):7.4f} {sv.get('recall', 0):7.4f} "
              f"{sv.get('precision', 0):7.4f} {(ps.smallvar_f1(s) or 0):8.4f} "
              f"{(ps.smallvar_f1(s, 'Snv') or 0):7.4f} {ps.cls(s, 'DEL 1k+ het'):>7s} "
              f"{shift.get('lo_pct', 0):12.2f}% {shift.get('hi_pct', 0):13.3f}%{flag}")

    dest = WORK / "sv-atlas" / "sweep-linkage-grid.json"
    dest.write_text(json.dumps(
        [{"dataset": ds, "weight": w, "scale": scale, "sv": s.get("sv"),
          "smallvar_all_f1": ps.smallvar_f1(s), "smallvar_snv_f1": ps.smallvar_f1(s, "Snv"),
          "sv_by_class": s.get("sv_by_class"), "shift": shift}
         for ds, w, scale, s, shift in rows], indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
