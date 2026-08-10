#!/usr/bin/env python3
"""Re-search the depth term's weight jointly with the mismapping floor.

Why now. Every number that supported `--depth-term` was measured against a lambda that
counted the snarl's two boundary nodes -- sequence that recruits no read into the
likelihood matrix, because a read lying entirely inside a boundary node cannot
discriminate and is dropped as uninformative. That put the median DR at 0.59 instead of
1, so typical sites were being scored on the Poisson's steep low-count flank rather than
near its mode. Both the size of the per-genotype differences and where on the curve they
are read off have changed, so `w_d = 0.25` is no longer a searched value.

Why these two axes together. `e_r` now enters the objective in three places -- the
per-read background term, the observation `N_eff = sum_r (1 - e_r)`, and the local rate
`c(s)`, which is measured under the same weighting so the correction stays relative.
`--mismap-min` therefore moves the depth term as well as the read term, and by more than
one route. Coordinate descent on a surface like that finds whichever corner it started
nearest, which is exactly the mistake the earlier cap sweep refused to make when it
declined to carry "floor 0.05 is optimal" through a cap change.

What is deliberately *not* here. `--depth-quality` is a ranking discount: it cannot
change a genotype, so putting it in this grid would multiply vg runs for nothing. It is
swept offline from the winning point's VCF by depth_gq.py -- but it does have to be
re-swept, because a non-zero `w_d` changes which genotypes are called and so changes the
DR distribution it acts on. `--depth-ploidy` is structural, not free. `--depth-count-raw`
was settled on all four datasets. `--depth-window` gets its own 1-D scan afterwards.

No single objective, on purpose: every point is scored on both benchmarks plus the
heterozygous class breakdown, and the surface is printed whole.

Two of the reported columns are a *harm* check rather than a gain check, and can veto a
point whatever its F1 does. lambda grows with allele length, so an anomalously large N
mechanically favours whichever genotype presents the most sequence -- a preference for
long alleles, not a rejection of the site, and the way this term can do damage. Stage 0
checked that and it passed at w_d = 0.25, but the anchor constant used to damp lambda's
length dependence and no longer does. `footprint@DR>3` is the median called footprint at
pile-up sites, against the w_d = 0 control on the same dataset.
"""

from __future__ import annotations

import argparse
import gzip
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import param_sweep as ps  # noqa: E402

WORK = ps.WORK


def dr_stats(vcf: Path) -> dict:
    """Centring, tail, and the long-allele-preference check, from the emitted VCF."""
    drs, foot_hi = [], []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            s = dict(zip(f[8].split(":"), f[9].split(":")))
            dr = s.get("DR")
            if dr in (None, ".") or float(dr) < 0:
                continue
            dr = float(dr)
            drs.append(dr)
            if dr > 3.0:
                # Total sequence the call claims: the summed length of the called
                # alleles. If the depth term is buying its gains by preferring longer
                # alleles at pile-ups, this rises with w_d.
                alts = f[4].split(",")
                total = 0
                for tok in s.get("GT", "").replace("|", "/").split("/"):
                    if tok == "0":
                        total += len(f[3])
                    elif tok.isdigit() and 0 < int(tok) <= len(alts):
                        total += len(alts[int(tok) - 1])
                foot_hi.append(total)
    if not drs:
        return {}
    drs.sort()
    return {"median_dr": statistics.median(drs),
            "n_dr_gt3": sum(1 for v in drs if v > 3.0),
            "n": len(drs),
            "footprint_hi": statistics.median(foot_hi) if foot_hi else 0.0}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="+", default=["chr20-4hap"])
    ap.add_argument("--depth-weights", nargs="+",
                    default=["0", "0.1", "0.25", "0.5", "1.0"])
    ap.add_argument("--floors", nargs="+", default=["0.01", "0.02", "0.05"])
    ap.add_argument("--caps", nargs="+", default=["0.7"])
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    ap.add_argument("--threads", type=int, default=5)
    args = ap.parse_args()

    rows = []
    for ds in args.datasets:
        for cap in args.caps:
            for floor in args.floors:
                for w in args.depth_weights:
                    params = {"mismap-max": cap, "mismap-min": floor}
                    # w = 0 is the control and must be *run*, not carried over from an
                    # earlier arm: the floor moves underneath it.
                    if float(w) > 0:
                        params["depth-term"] = w
                    tag = f"{ds}-dgrid-w{w}-f{floor}-c{cap}"
                    print(f"=== {ds} w_d={w} floor={floor} cap={cap}", flush=True)
                    vcf = ps.call(args.vg, ds, tag, params, args.threads)
                    s = ps.score(vcf, ds, tag, args.threads)
                    rows.append((ds, w, floor, cap, s, dr_stats(vcf)))

    hdr = (f"{'dataset':12s} {'w_d':>5s} {'floor':>6s} {'cap':>5s} {'SV F1':>7s} "
           f"{'SV TP':>6s} {'SV FP':>6s} {'smallGT':>8s} {'SNV':>7s} {'hDEL1k':>7s} "
           f"{'hDEL3-9':>8s} {'hINS1k':>7s} {'hetfr':>6s} {'medDR':>6s} {'DR>3':>5s} "
           f"{'foot@hi':>8s}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for ds, w, floor, cap, s, d in rows:
        sv = s.get("sv") or {}
        print(f"{ds:12s} {w:>5s} {floor:>6s} {cap:>5s} {sv.get('f1', 0):7.4f} "
              f"{sv.get('TP-base', 0):6d} {sv.get('FP', 0):6d} "
              f"{(ps.smallvar_f1(s) or 0):8.4f} {(ps.smallvar_f1(s, 'Snv') or 0):7.4f} "
              f"{ps.cls(s, 'DEL 1k+ het'):>7s} {ps.cls(s, 'DEL 300-999 het'):>8s} "
              f"{ps.cls(s, 'INS 1k+ het'):>7s} "
              f"{(s.get('genotype_mix') or {}).get('het_frac', 0):6.4f} "
              f"{d.get('median_dr', 0):6.3f} {d.get('n_dr_gt3', 0):5d} "
              f"{d.get('footprint_hi', 0):8.0f}")

    dest = WORK / "sv-atlas" / "sweep-depth-grid.json"
    dest.write_text(json.dumps(
        [{"dataset": ds, "depth_weight": w, "mismap_min": floor, "mismap_max": cap,
          "sv": s.get("sv"), "smallvar_all_f1": ps.smallvar_f1(s),
          "smallvar_snv_f1": ps.smallvar_f1(s, "Snv"), "sv_by_class": s.get("sv_by_class"),
          "genotype_mix": s.get("genotype_mix"), "dr": d}
         for ds, w, floor, cap, s, d in rows], indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
