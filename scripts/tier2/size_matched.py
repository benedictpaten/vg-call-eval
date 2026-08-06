#!/usr/bin/env python3
"""Score the arms against the small-variant benchmark on the size range it can score.

Why this exists. The GIAB draft small-variant benchmark for chr20 contains **no truth
record >=50 bp** -- that size class lives in the separate structural-variant benchmark.
But the two benchmarks' *confident regions* overlap almost completely (58.9 Mb vs
59.4 Mb, both ~89% of chr20). So a caller that emits a 300 bp insertion inside the
small-variant confident region has every one of those 300 bases scored as a false
positive, however right the call is. It is unscoreable-as-correct by construction.

That is not a hypothetical. The read-likelihood caller emits large insertions the
Poisson caller does not, which made its insertion BASEPAIR precision look 0.139 worse
while its recall was 0.094 *better*. Restricting both callers to <50 bp -- the range the
benchmark can actually adjudicate -- collapses that gap to 0.008 and flips insertion
BASEPAIR F1 from a 0.047 loss into a 0.047 win. See docs/tier2-chr20-results.md.

The restriction drops a whole VCF record when any *called* allele differs from REF by
>=50 bp, applied identically to every arm. Dropping the record rather than the allele
loses the record's other alleles too, which is why it is applied to both sides: the
comparison is symmetric even though neither side is complete.

This is a diagnostic, not the headline metric. It answers "is the read-likelihood model
worse at insertion sequence?" (no) and deliberately says nothing about whether the large
calls it removes are correct -- for that, read the structural-variant comparison in
compare_sv.py, which scores exactly the class this script discards.
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

ARMS = ["poisson-z", "readlik-z"]
MAX_LEN = 50


def restrict(src: Path, dst: Path, max_len: int) -> tuple[int, int]:
    """Copy src to dst, dropping records with a called allele >= max_len from REF."""
    reader = subprocess.Popen(["bcftools", "view", str(src)], stdout=subprocess.PIPE, text=True)
    writer = subprocess.Popen(["bcftools", "view", "-Oz", "-o", str(dst), "-"],
                              stdin=subprocess.PIPE, text=True)
    kept = dropped = 0
    for line in reader.stdout:
        if line.startswith("#"):
            writer.stdin.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        ref, alts, gt = fields[3], fields[4].split(","), fields[9].split(":")[0]
        called = {int(t) for t in gt.replace("|", "/").split("/") if t.isdigit() and int(t) > 0}
        big = any(i <= len(alts) and not alts[i - 1].startswith("<") and alts[i - 1] != "*"
                  and abs(len(alts[i - 1]) - len(ref)) >= max_len for i in called)
        if big:
            dropped += 1
        else:
            kept += 1
            writer.stdin.write(line)
    writer.stdin.close()
    writer.wait()
    reader.wait()
    subprocess.run(["bcftools", "index", "-t", "-f", str(dst)], check=True)
    return kept, dropped


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--results", default=str(REPO / "work/tier2-chr20/results"))
    p.add_argument("--truth-vcf", default=str(REPO / "work/tier2-chr20/truth.chr20.smvar.vcf.gz"))
    p.add_argument("--truth-bed", default=str(REPO / "work/tier2-chr20/truth.chr20.smvar.bed"))
    p.add_argument("--reference", default=str(REPO / "work/tier2-chr20/chr20.fa"))
    p.add_argument("--aardvark", default="aardvark")
    p.add_argument("--threads", type=int, default=6)
    p.add_argument("--max-len", type=int, default=MAX_LEN)
    p.add_argument("--arms", nargs="*", default=ARMS)
    args = p.parse_args()

    res = Path(args.results)
    payload = []

    for name in args.arms:
        src = res / f"{name}.vcf.gz"
        if not src.exists():
            print(f"skip {name}: {src} missing", flush=True)
            continue
        label = f"sm{args.max_len}-{name}"
        dst = res / f"{label}.vcf.gz"
        kept, dropped = restrict(src, dst, args.max_len)
        print(f"\n=== {label}: kept {kept:,}, dropped {dropped:,} "
              f"records with a >={args.max_len} bp called allele", flush=True)

        adir = res / f"aardvark-{label}"
        started = time.time()
        aardvark.compare(
            aardvark=args.aardvark,
            reference=Path(args.reference),
            truth_vcf=Path(args.truth_vcf),
            query_vcf=dst,
            regions_bed=Path(args.truth_bed),
            out_dir=adir,
            truth_sample="HG002",
            query_sample="SAMPLE",
            label=label,
            options=aardvark.AardvarkOptions(enable_record_basepair=True, threads=args.threads),
        )
        summary = aardvark.read_summary(adir)
        elapsed = time.time() - started
        print(f"  aardvark {elapsed:.0f}s, {len(summary)} summary rows", flush=True)
        payload.append({"arm": label, "source_arm": name, "max_len": args.max_len,
                        "kept": kept, "dropped": dropped, "seconds": round(elapsed, 1),
                        "metrics": {"summary": summary}})

    # Merge, for the same reason compare_sv.py does: a run restricted to --arms must not
    # silently delete the arms it was not asked about.
    out_path = res / "arms-size-matched.json"
    merged = {}
    if out_path.exists():
        for entry in json.loads(out_path.read_text()):
            merged[entry["arm"]] = entry
    for entry in payload:
        merged[entry["arm"]] = entry
    out_path.write_text(json.dumps(list(merged.values()), indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
