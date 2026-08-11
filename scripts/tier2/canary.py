#!/usr/bin/env python3
"""Verify that a cached arm is still valid for this build, instead of assuming it is not.

`refresh_all.sh` re-runs all five arms from one build because a table whose rows come from
different builds is the failure this harness exists to prevent -- it has happened here once.
But blind re-running only *assumes* the rows would agree. The Poisson arms cannot be touched
by a read-likelihood change, and they are 43% of the matrix's run time (36 minutes of 86).

So verify instead: re-run one cheap Poisson arm, compare the VCF body -- header excluded,
since it carries no genotype information -- against the cached copy, and skip the rest only
if they are identical. That is strictly stronger than re-running, because a mismatch tells
you something touched shared code, which re-running would have quietly absorbed.

Exit 0 if identical (the cached Poisson rows are usable), 1 if not (re-run everything).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from param_sweep import DATASETS, WORK  # noqa: E402


def body_hash(path: Path) -> str:
    h = hashlib.sha256()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                h.update(line.encode())
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    ap.add_argument("--dataset", default="chr20-4hap", help="cheapest dataset by default")
    ap.add_argument("--arm", default="poisson-z", help="cheapest unaffected arm by default")
    ap.add_argument("--threads", type=int, default=5)
    args = ap.parse_args()

    sub, contig, _, _ = DATASETS[args.dataset]
    w = WORK / sub
    cached = w / "results" / f"{args.arm}.vcf.gz"
    if not cached.exists():
        print(f"canary: no cached {args.arm} for {args.dataset}; nothing to verify against")
        sys.exit(1)

    extra = {"poisson": ["-k", str(w / f"{contig}.pack")],
             "poisson-z": ["-z", "-k", str(w / f"{contig}.pack")]}[args.arm]
    cmd = [args.vg, "call", str(w / f"{contig}_0_{contig}.gbz"), "-p", f"CHM13#0#{contig}",
           "-t", str(args.threads)] + extra

    with tempfile.TemporaryDirectory() as td:
        raw = Path(td) / "canary.vcf.gz"
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        if proc.returncode != 0 or not proc.stdout:
            print(f"canary: vg call failed (rc={proc.returncode})")
            sys.exit(1)
        with open(raw, "wb") as fh:
            subprocess.run(["bgzip", "-c"], input=proc.stdout, stdout=fh, check=True)
        fresh, old = body_hash(raw), body_hash(cached)

    if fresh == old:
        print(f"canary PASS: {args.arm} on {args.dataset} is byte-identical to the cached "
              f"arm ({fresh[:16]}), so the Poisson rows carry over")
        sys.exit(0)
    print(f"canary FAIL: {args.arm} on {args.dataset} changed ({old[:16]} -> {fresh[:16]}). "
          f"Something touched shared code; re-run every arm.")
    sys.exit(1)


if __name__ == "__main__":
    main()
