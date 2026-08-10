#!/usr/bin/env python3
"""Re-sweep the read-likelihood caller's tuned parameters, scoring both benchmarks.

Why re-sweep. `--mismap-max` and `--mismap-min` were tuned against a model whose
mixture weights were flat. That model has been replaced, and not marginally: the
weighted mixture changes which genotype wins at every site whose alleles differ in
length. `--mismap-min` in particular was *directly* implicated in the failure the new
mixture fixes -- lowering it recovered a third of the lost heterozygous deletions on its
own -- so its old optimum was partly compensating for a defect that no longer exists.
Carrying those settings across unexamined is the same mistake the original cap sweep
called out when it refused to carry "floor 0.05 is optimal" across a cap change.

Method. One parameter at a time by default, and a 2-D grid where the parameters are
known to interact. `--mismap-max` and `--mismap-min` do interact, and the evidence is in
the previous sweep's own output: on the 34-haplotype graph at floor 0.01 the cap ranked
0.5 > 0.9 (0.9520 against 0.9517), while at floor 0.05 cap 0.9 gave 0.9567 -- the best
point in that whole sweep, and one never carried forward. Coordinate descent through a
surface like that finds whichever corner it started nearest. Pass --param2/--values2 for
the grid.

Objective. Deliberately none. Every point is scored on *both* benchmarks plus the
heterozygous SV class breakdown and the genotype mix, and the surface is printed whole.
A setting that buys structural-variant F1 with small-variant genotype F1 is a judgement
call about what the caller is for, not a maximisation -- and picking a single number to
maximise is how the mismapping cap ended up at 0.1 in the first place.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"

DATASETS = {
    "chr6-4hap": ("tier2-chr6", "chr6", "graph.gbz.db", "reads.gaf.db"),
    "chr6-34hap": ("tier2-chr6-hap32", "chr6", "graph.hap32.gbz.db", "reads.hap32.gaf.db"),
    "chr20-4hap": ("tier2-chr20", "chr20", "graph.gbz.db", "reads.gaf.db"),
    "chr20-34hap": ("tier2-chr20-hap32", "chr20", "graph.hap32.gbz.db", "reads.hap32.gaf.db"),
}

# The shipped operating point. Sweeping one parameter holds the others here, so this
# has to track vg's own defaults or "holding the others fixed" quietly means holding
# them somewhere vg no longer is: the cap moved to 0.7 as a *result* of this sweep and
# this line kept the pre-sweep value.
# --read-weight was removed from vg after this sweep showed it cannot change a
# genotype; passing it to a current build is an argument error, not a no-op.
DEFAULTS = {"mismap-max": "0.7", "mismap-min": "0.02"}


def gbz_base_binary() -> str:
    found = shutil.which("gbz-base")
    if found:
        return found
    hits = sorted(Path("/private/tmp/claude-501").rglob("gbz-base"))
    for h in hits:
        if h.is_file():
            return str(h)
    sys.exit("gbz-base not found")


def call(vg: str, ds: str, tag: str, params: dict, threads: int) -> Path:
    sub, contig, gbzdb, gafdb = DATASETS[ds]
    w = WORK / sub
    out = w / "results" / f"sweep-{tag}.vcf.gz"
    # Size, not just existence. A failed `vg call` still leaves a 28-byte empty
    # bgzip and a valid index behind, and an existence check happily caches that
    # forever -- the point re-reads as "already measured" on every later run.
    if out.exists() and out.stat().st_size > 4096 and out.with_suffix(".gz.tbi").exists():
        print(f"  {ds} {tag}: cached", flush=True)
        return out
    cmd = [vg, "call", str(w / f"{contig}_0_{contig}.gbz"), "-p", f"CHM13#0#{contig}",
           "-t", str(threads), "--read-likelihood", "-z",
           "--gaf-base", str(WORK / gafdb), "--gbz-base", str(WORK / gbzdb),
           "--gaf-base-binary", gbz_base_binary()]
    for k, v in params.items():
        # An empty value means a no-argument flag, not a flag with an empty argument.
        cmd += [f"--{k}"] + ([] if v == "" else [str(v)])
    log = w / "results" / f"sweep-{tag}.log"
    with open(out.with_suffix(""), "wb") as fh, open(log, "wb") as errfh:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=errfh)
        if proc.returncode != 0 or not proc.stdout:
            out.with_suffix("").unlink(missing_ok=True)
            sys.exit(f"vg call failed for {tag} (rc={proc.returncode}); see {log}")
        subprocess.run(["bgzip", "-c"], input=proc.stdout, stdout=fh, check=True)
    shutil.move(str(out.with_suffix("")), str(out))
    subprocess.run(["tabix", "-f", "-p", "vcf", str(out)], check=True)
    return out


def score(vcf: Path, ds: str, tag: str, threads: int) -> dict:
    dest = WORK / "sv-atlas" / f"score-{tag}.json"
    if dest.exists():
        return json.loads(dest.read_text())
    subprocess.run([sys.executable, str(HERE / "score_vcf.py"), "--vcf", str(vcf),
                    "--label", tag, "--dataset", ds, "--threads", str(threads)],
                   check=True, capture_output=True)
    return json.loads(dest.read_text())


def smallvar_f1(s: dict, vtype: str = "ALL") -> float | None:
    for r in s.get("smallvar") or []:
        if (r.get("comparison") == "GT" and r.get("region_label") == "ALL"
                and r.get("filter") == "ALL" and r.get("variant_type") == vtype):
            try:
                return float(r["metric_f1"])
            except (KeyError, ValueError):
                return None
    return None


def cls(s: dict, key: str) -> str:
    v = (s.get("sv_by_class") or {}).get(key)
    if not v:
        return "-"
    tp, n = v
    return f"{tp/n:.3f}" if n else "-"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--param", required=True, choices=sorted(DEFAULTS))
    ap.add_argument("--values", nargs="+", required=True)
    ap.add_argument("--param2", choices=sorted(DEFAULTS),
                    help="second axis, for a 2-D grid over interacting parameters")
    ap.add_argument("--values2", nargs="+", default=None)
    ap.add_argument("--datasets", nargs="+", default=["chr6-4hap"])
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    ap.add_argument("--threads", type=int, default=5)
    args = ap.parse_args()

    axis2 = args.values2 if args.param2 else [None]
    rows = []
    for ds in args.datasets:
        for val in args.values:
            for val2 in axis2:
                params = dict(DEFAULTS)
                params[args.param] = val
                label = f"{args.param}{val}"
                if args.param2:
                    # vg rejects --mismap-min above --mismap-max, so a grid over both
                    # necessarily has invalid corners. Skip them rather than logging a
                    # failed run that later reads as a missing measurement.
                    params[args.param2] = val2
                    label += f"-{args.param2}{val2}"
                    if float(params["mismap-min"]) > float(params["mismap-max"]):
                        print(f"=== {ds} {label}: skipped (min > max)", flush=True)
                        continue
                tag = f"{ds}-{label}"
                print(f"=== {ds} {label}", flush=True)
                vcf = call(args.vg, ds, tag, params, args.threads)
                s = score(vcf, ds, tag, args.threads)
                rows.append((ds, label, s))

    hdr = (f"{'dataset':12s} {'setting':>26s} {'SV F1':>7s} {'SV TP':>6s} {'SV FP':>6s} "
           f"{'smallGT':>8s} {'SNV':>7s} {'hetDEL1k':>9s} {'hetDEL3-9':>10s} "
           f"{'hetINS1k':>9s} {'het frac':>9s}")
    print("\n" + hdr)
    print("-" * len(hdr))
    for ds, val, s in rows:
        sv = s.get("sv") or {}
        print(f"{ds:12s} {val:>26s} {sv.get('f1', 0):7.4f} {sv.get('TP-base', 0):6d} "
              f"{sv.get('FP', 0):6d} "
              f"{(smallvar_f1(s) or 0):8.4f} {(smallvar_f1(s, 'Snv') or 0):7.4f} "
              f"{cls(s, 'DEL 1k+ het'):>9s} {cls(s, 'DEL 300-999 het'):>10s} "
              f"{cls(s, 'INS 1k+ het'):>9s} "
              f"{(s.get('genotype_mix') or {}).get('het_frac', 0):9.4f}")

    suffix = args.param + (f"-x-{args.param2}" if args.param2 else "")
    dest = WORK / "sv-atlas" / f"sweep-{suffix}.json"
    dest.write_text(json.dumps(
        [{"dataset": d, "value": v, "sv": s.get("sv"),
          "smallvar_all_f1": smallvar_f1(s), "smallvar_snv_f1": smallvar_f1(s, "Snv"),
          "sv_by_class": s.get("sv_by_class"), "genotype_mix": s.get("genotype_mix")}
         for d, v, s in rows], indent=2))
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
