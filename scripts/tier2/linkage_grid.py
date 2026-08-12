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

`--linkage-freq-prior` is the third axis and, measured, the most important one. It was excluded
from the first version of this search on the argument that panel allele frequency over a
read-sampled panel counts the same evidence twice. That is true and irrelevant: nothing connects
the panel to the *truth set*, and reusing one's own reads is what mapping and calling already do.
Double counting can leave `GQ` overconfident while the genotype improves, which is a calibration
question, not a reason to hold a parameter at zero.

It also had a ceiling. Both application sites guarded on `freq_prior < 1.0`, which is a no-op at
exactly 1 and therefore looked free, while silently making every larger value behave as 1. Above 1
the exponent goes negative and the prior is amplified past the state space's own multiplicity.
Uncapped it peaks near 5, is worth more than the transition weight, and inverts past 8 -- so pass
it explicitly rather than relying on any default when searching.

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
    ap.add_argument("--freq-priors", nargs="+", default=["0"])
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    ap.add_argument("--threads", type=int, default=5)
    args = ap.parse_args()

    rows = []
    for ds in args.datasets:
        # The w = 0 control is run, not carried over: it must come from this binary and these
        # flags, or the comparison is against a different experiment.
        baseline_vcf = None
        for scale in args.scales:
            for fp in args.freq_priors:
                for w in args.weights:
                    # Both axes are inert at zero weight -- the layer is switched off entirely --
                    # so one control per dataset, not one per cell.
                    if float(w) == 0.0 and (scale != args.scales[0]
                                            or fp != args.freq_priors[0]):
                        continue
                    params = {"linkage-scale": scale}
                    if float(w) > 0:
                        params["linkage-weight"] = w
                        params["linkage-freq-prior"] = fp
                    # Every swept axis belongs in the tag. Each was fixed when it was written, so
                    # leaving it out was harmless then; swept, two cells would collide on one
                    # cached VCF and the second would silently report the first one's numbers.
                    tag = (f"{ds}-lgrid-w{w}" if float(w) == 0.0
                           else f"{ds}-lgrid-w{w}-s{scale}-f{fp}")
                    print(f"=== {ds} w={w} scale={scale} f={fp}", flush=True)
                    vcf = ps.call(args.vg, ds, tag, params, args.threads)
                    sc = ps.score(vcf, ds, tag, args.threads)
                    if float(w) == 0.0:
                        baseline_vcf = vcf
                    shift = genotype_shift(baseline_vcf, vcf) if baseline_vcf else {}
                    rows.append((ds, w, scale, fp, sc, shift))

    hdr = (f"{'dataset':12s} {'w':>4s} {'scale':>6s} {'f':>5s} {'SV F1':>7s} {'SVrec':>7s} {'SVprec':>7s} "
           f"{'smallGT':>8s} {'SNV':>7s} {'hDEL1k':>7s} {'moved@GQI<10':>13s} "
           f"{'moved@GQI>=40':>14s}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for ds, w, scale, fp, s, shift in rows:
        sv = s.get("sv") or {}
        hi = shift.get("hi_pct")
        flag = "" if hi is None or hi <= 0.1 else "  OVER BUDGET"
        print(f"{ds:12s} {w:>4s} {scale:>6s} {fp:>5s} {sv.get('f1', 0):7.4f} {sv.get('recall', 0):7.4f} "
              f"{sv.get('precision', 0):7.4f} {(ps.smallvar_f1(s) or 0):8.4f} "
              f"{(ps.smallvar_f1(s, 'Snv') or 0):7.4f} {ps.cls(s, 'DEL 1k+ het'):>7s} "
              f"{shift.get('lo_pct', 0):12.2f}% {shift.get('hi_pct', 0):13.3f}%{flag}")

    dest = WORK / "sv-atlas" / "sweep-linkage-grid.json"
    dest.write_text(json.dumps(
        [{"dataset": ds, "weight": w, "scale": scale, "freq_prior": fp, "sv": s.get("sv"),
          "smallvar_all_f1": ps.smallvar_f1(s), "smallvar_snv_f1": ps.smallvar_f1(s, "Snv"),
          "sv_by_class": s.get("sv_by_class"), "shift": shift}
         for ds, w, scale, fp, s, shift in rows], indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
