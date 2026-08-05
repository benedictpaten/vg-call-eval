#!/usr/bin/env python3
"""Compare the already-called arm VCFs against the GIAB *structural variant* benchmark.

Separate from run_arms.py because no re-calling is needed: the callers emit one VCF
containing everything, and the small-variant and structural-variant benchmarks are two
different truth sets over (nearly) the same confident region -- 58.9 Mb for smvar,
59.4 Mb for stvar, both ~89% of chr20.

Aardvark's SV preset is used (`--min-variant-gap 1000`, record-basepair metrics on),
which is what the harness plan §2.1 specifies for this size class: wider clustering, so
an SV and the small variants around it are resolved together rather than as independent
calls.

Note on what "the SV benchmark" contains. Of 176,623 chr20 truth records only **2,052
are >=50 bp**; the rest are the local sequence context an SV-aware haplotype comparison
needs in order to place the SV correctly. So these numbers are dominated by the same
small variants scored in the small-variant run -- what changes is the clustering, and
therefore how a large event and its neighbours are credited. Read the SV-specific rows,
not the ALL row, for anything about SV calling.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))

from vgcalleval.engines import aardvark  # noqa: E402

ARMS = ["poisson", "poisson-z", "readlik", "readlik-nomismap", "readlik-z"]


def sv_counts(vcf: Path) -> dict:
    """Count >=50 bp events in a VCF, so the SV content is stated rather than assumed."""
    out = subprocess.run(["bcftools", "query", "-f", "%REF\t%ALT\n", str(vcf)],
                         capture_output=True, text=True)
    ins = dele = 0
    for line in out.stdout.splitlines():
        try:
            ref, alts = line.split("\t")
        except ValueError:
            continue
        for alt in alts.split(","):
            if alt.startswith("<") or alt == "*":
                continue
            d = len(alt) - len(ref)
            if d >= 50:
                ins += 1
            elif d <= -50:
                dele += 1
    return {"sv_insertions": ins, "sv_deletions": dele, "sv_total": ins + dele}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default=str(REPO / "work/tier2-chr20/results"))
    p.add_argument("--truth-vcf", default=str(REPO / "work/tier2-chr20/truth.chr20.stvar.vcf.gz"))
    p.add_argument("--truth-bed", default=str(REPO / "work/tier2-chr20/truth.chr20.stvar.bed"))
    p.add_argument("--reference", default=str(REPO / "work/tier2-chr20/chr20.fa"))
    p.add_argument("--aardvark", default="aardvark")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--arms", nargs="*", default=ARMS)
    args = p.parse_args()

    res = Path(args.results)
    payload = []

    for name in args.arms:
        gz = res / f"{name}.vcf.gz"
        if not gz.exists():
            print(f"skip {name}: {gz} missing", flush=True)
            continue
        counts = sv_counts(gz)
        print(f"\n=== {name}: {counts['sv_total']:,} called SVs "
              f"({counts['sv_insertions']:,} ins / {counts['sv_deletions']:,} del)", flush=True)

        adir = res / f"aardvark-sv-{name}"
        started = time.time()
        try:
            aardvark.compare(
                aardvark=args.aardvark,
                reference=Path(args.reference),
                truth_vcf=Path(args.truth_vcf),
                query_vcf=gz,
                regions_bed=Path(args.truth_bed),
                out_dir=adir,
                truth_sample="HG002",
                query_sample="SAMPLE",
                label=f"sv-{name}",
                options=aardvark.AardvarkOptions.for_structural_variants(threads=args.threads),
            )
            summary = aardvark.read_summary(adir)
        except Exception as exc:  # noqa: BLE001
            print(f"  FAILED: {exc}", flush=True)
            summary = []
        elapsed = time.time() - started
        print(f"  aardvark {elapsed:.0f}s, {len(summary)} summary rows", flush=True)
        payload.append({"arm": name, "seconds": round(elapsed, 1),
                        **counts, "metrics": {"summary": summary}})

    (res / "arms-sv.json").write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {res / 'arms-sv.json'}")


if __name__ == "__main__":
    main()
