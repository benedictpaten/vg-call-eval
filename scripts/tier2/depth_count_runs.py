#!/usr/bin/env python3
"""A/B the depth term's read count: one read apiece, or `1 - e_r` of a read?

The Stage 1 depth term compared a raw row count `N` against `lambda_G`. That asserts
something the rest of the model explicitly declines to assert. The read term already
believes each read only to the extent of `1 - e_r`; a MAPQ 0 read enters it at 0.3 of
its weight and then, one line later, was counted as a whole read of depth.

The correction is `N_eff = sum_r (1 - e_r)` -- and, critically, the *same* weighting
applied when the local rate is measured, because the rate and the count have to be in
the same units. Weighting one side only would put a constant factor between `N` and
`lambda` and push every `DR` in the same direction, which is a bias, not a signal.

That symmetry is also why the expected effect is small and local. A site whose mapping
quality matches its neighbourhood's sees the factor cancel exactly. The term moves only
where a site is more or less ambiguously mapped than the sequence around it -- which is
the case worth testing, not a case the previous arm could express at all.

Run both arms from one script so the command lines are identical apart from the flag
under test. `--depth-count-raw` restores the Stage 1 behaviour exactly.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import param_sweep as ps  # noqa: E402

WORK = ps.WORK


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["chr20-4hap"])
    ap.add_argument("--weight", default="0.25")
    ap.add_argument("--arms", nargs="+", default=["raw", "eff"],
                    choices=["raw", "eff", "off"])
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    ap.add_argument("--threads", type=int, default=5)
    args = ap.parse_args()

    # vg's own defaults for everything else, so this measures the flag and nothing
    # else. mismap-max is 0.7 and mismap-min 0.02 in the shipped build.
    arms = {
        "off": {},
        "raw": {"depth-term": args.weight, "depth-count-raw": ""},
        "eff": {"depth-term": args.weight},
    }

    rows = []
    for ds in args.datasets:
        for arm in args.arms:
            params = {k: v for k, v in arms[arm].items()}
            tag = f"{ds}-depthcount-{arm}{'' if arm == 'off' else args.weight}"
            print(f"=== {ds} {arm}", flush=True)
            vcf = ps.call(args.vg, ds, tag, params, args.threads)
            rows.append((ds, arm, ps.score(vcf, ds, tag, args.threads)))

    hdr = (f"{'dataset':12s} {'arm':>6s} {'SV F1':>7s} {'SV TP':>6s} {'SV FP':>6s} "
           f"{'smallGT':>8s} {'SNV':>7s} {'hetDEL1k':>9s} {'hetDEL3-9':>10s} "
           f"{'hetINS1k':>9s} {'het frac':>9s}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for ds, arm, s in rows:
        sv = s.get("sv") or {}
        print(f"{ds:12s} {arm:>6s} {sv.get('f1', 0):7.4f} {sv.get('TP-base', 0):6d} "
              f"{sv.get('FP', 0):6d} "
              f"{(ps.smallvar_f1(s) or 0):8.4f} {(ps.smallvar_f1(s, 'Snv') or 0):7.4f} "
              f"{ps.cls(s, 'DEL 1k+ het'):>9s} {ps.cls(s, 'DEL 300-999 het'):>10s} "
              f"{ps.cls(s, 'INS 1k+ het'):>9s} "
              f"{(s.get('genotype_mix') or {}).get('het_frac', 0):9.4f}")

    dest = WORK / "sv-atlas" / "sweep-depth-count.json"
    dest.write_text(json.dumps(
        [{"dataset": d, "arm": a, "sv": s.get("sv"),
          "smallvar_all_f1": ps.smallvar_f1(s), "smallvar_snv_f1": ps.smallvar_f1(s, "Snv"),
          "sv_by_class": s.get("sv_by_class"), "genotype_mix": s.get("genotype_mix")}
         for d, a, s in rows], indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
