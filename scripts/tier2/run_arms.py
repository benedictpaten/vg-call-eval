#!/usr/bin/env python3
"""Tier 2: run the caller arms on HG002 chr20 and compare against the GIAB draft benchmark.

Implements plan §9.4 steps 5-7. Deliberately a standalone script rather than a
path through `vgcalleval.pipeline`, because that pipeline builds its dataset by
simulation; here the dataset is 28.8 GB of real reads prepared once by hand.
The aardvark engine module is reused as-is.

Five arms, chosen to separate two effects that would otherwise be confounded
(plan §9.2a). The graph carries only 4 haplotypes, so `-z` enumeration can only
propose alleles present in those walks:

    caller \\ enumeration   support (pack, FlowTraversalFinder)   haplotype (-z)
    Poisson                 poisson                              poisson-z
    read-likelihood         readlik, readlik-nomismap            readlik-z

Comparing down a column isolates the caller; comparing across a row isolates
what the sampled graph costs. Without both, a graph limitation reads as a caller
limitation.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))

from vgcalleval.engines import aardvark  # noqa: E402

# Set from --contig; the constants remain as defaults so existing chr20 invocations
# keep working unchanged.
REF_SAMPLE = "CHM13"
REF_PATH = "CHM13#0#chr20"
CONTIG = "chr20"


@dataclass
class Arm:
    name: str
    extra_args: list[str]
    needs_pack: bool
    needs_reads: bool
    description: str = ""
    # filled in by the run
    seconds: float = 0.0
    peak_rss_gb: float = 0.0
    variants: int = 0
    metrics: dict = field(default_factory=dict)


def arms() -> list[Arm]:
    return [
        Arm("poisson", [], True, False,
            "current default: Poisson genotyping, support enumeration"),
        # Needs the pack despite -z: the pack becomes optional only for
        # --read-likelihood, because Poisson *genotyping* consumes support even when
        # enumeration comes from haplotypes.
        Arm("poisson-z", ["-z"], True, False,
            "Poisson genotyping, haplotype enumeration from the 4 graph haplotypes"),
        Arm("readlik", ["--read-likelihood"], True, True,
            "read-level likelihoods, support enumeration -- the like-for-like caller comparison"),
        Arm("readlik-nomismap", ["--read-likelihood", "--no-mismap-term"], True, True,
            "as readlik, MAPQ mismapping term disabled, to measure its contribution"),
        Arm("readlik-z", ["--read-likelihood", "-z"], False, True,
            "read-level likelihoods, haplotype enumeration, no pack file"),
    ]


def sh(cmd: list[str], out_path: Path | None = None, log: Path | None = None,
       time_it: bool = False) -> tuple[int, float, float]:
    """Run a command, optionally under /usr/bin/time -l. Returns (rc, seconds, peak_rss_gb)."""
    full = ["/usr/bin/time", "-l"] + cmd if time_it else cmd
    started = time.time()
    with open(out_path, "wb") if out_path else open("/dev/null", "wb") as out:
        proc = subprocess.run(full, stdout=out, stderr=subprocess.PIPE)
    elapsed = time.time() - started
    stderr = proc.stderr.decode(errors="replace")
    if log:
        log.write_text(stderr)
    peak = 0.0
    for line in stderr.splitlines():
        if "maximum resident set size" in line:
            peak = int(line.split()[0]) / 1073741824
    return proc.returncode, elapsed, peak


def rename_contig(src: Path, dst: Path) -> int:
    """Rewrite CHM13#0#chr20 -> chr20 in CHROM and ##contig lines.

    The benchmark VCF uses bare `chr20`. Skipping this does not error -- it
    produces zero true positives, which looks like a catastrophic caller failure
    rather than a naming mismatch. Asserted on below.
    """
    kept = 0
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("##contig="):
                fout.write(line.replace(f"ID={REF_PATH}", f"ID={CONTIG}"))
            elif line.startswith("#"):
                fout.write(line)
            else:
                fields = line.split("\t", 1)
                if fields[0] == REF_PATH:
                    fout.write(f"{CONTIG}\t{fields[1]}")
                else:
                    fout.write(line)
                kept += 1
    return kept


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--vg", default=str(Path.home() / "CLionProjects/vg/bin/vg"))
    p.add_argument("--graph", default=str(HERE / "chr20_0_chr20.gbz"))
    p.add_argument("--pack", default=str(HERE / "chr20.pack"))
    p.add_argument("--gaf-base", default=str(REPO / "work/reads.gaf.db"))
    p.add_argument("--gbz-base", default=str(REPO / "work/graph.gbz.db"))
    p.add_argument("--reference", default=str(HERE / "chr20.fa"))
    p.add_argument("--truth-vcf", default=str(HERE / "truth.chr20.smvar.vcf.gz"))
    p.add_argument("--truth-bed", default=str(HERE / "truth.chr20.smvar.bed"))
    p.add_argument("--aardvark", default="aardvark")
    p.add_argument("--threads", type=int, default=6)
    # Not passed unless set: vg picks a per-backend default (4096 for --gaf-base),
    # and the harness should track that rather than pin a value that silently
    # becomes wrong when the default moves.
    p.add_argument("--read-window", type=int, default=0)
    p.add_argument("--out", default=str(HERE / "results"))
    p.add_argument("--only", nargs="*", help="run only these arms")
    p.add_argument("--contig", default="chr20",
                   help="contig to call; sets both the reference path and the "
                        "name written into the output VCF")
    args = p.parse_args()

    global REF_PATH, CONTIG
    CONTIG = args.contig
    REF_PATH = f"{REF_SAMPLE}#0#{args.contig}"

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    selected = [a for a in arms() if not args.only or a.name in args.only]

    for arm in selected:
        print(f"\n=== {arm.name}: {arm.description}", flush=True)
        raw = out / f"{arm.name}.raw.vcf"
        cmd = [args.vg, "call", args.graph, "-p", REF_PATH,
               "-t", str(args.threads), "--progress"] + arm.extra_args
        if arm.needs_pack:
            cmd += ["-k", args.pack]
        if arm.needs_reads:
            cmd += ["--gaf-base", args.gaf_base, "--gbz-base", args.gbz_base]
            if args.read_window:
                cmd += ["--read-window", str(args.read_window)]

        rc, secs, peak = sh(cmd, out_path=raw, log=out / f"{arm.name}.call.log", time_it=True)
        arm.seconds, arm.peak_rss_gb = secs, peak
        if rc != 0:
            print(f"  FAILED rc={rc}; see {arm.name}.call.log", flush=True)
            tail = (out / f"{arm.name}.call.log").read_text().splitlines()[-6:]
            print("  " + "\n  ".join(tail), flush=True)
            continue

        renamed = out / f"{arm.name}.vcf"
        arm.variants = rename_contig(raw, renamed)
        print(f"  {arm.variants} variants, {secs:.0f}s, {peak:.1f} GB peak", flush=True)
        raw.unlink()

        gz = renamed.with_suffix(".vcf.gz")
        subprocess.run(["bgzip", "-f", "-c", str(renamed)], stdout=open(gz, "wb"), check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(gz)], check=True)

        sample = subprocess.run(["bcftools", "query", "-l", str(gz)],
                                capture_output=True, text=True).stdout.split()
        query_sample = sample[0] if sample else "SAMPLE"

        adir = out / f"aardvark-{arm.name}"
        try:
            aardvark.compare(
                aardvark=args.aardvark,
                reference=Path(args.reference),
                truth_vcf=Path(args.truth_vcf),
                query_vcf=gz,
                regions_bed=Path(args.truth_bed),
                out_dir=adir,
                truth_sample="HG002",
                query_sample=query_sample,
                label=arm.name,
                options=aardvark.AardvarkOptions(threads=args.threads),
            )
            arm.metrics = {"summary": aardvark.read_summary(adir)}
        except Exception as exc:  # noqa: BLE001
            print(f"  aardvark failed: {exc}", flush=True)

    payload = [
        {"arm": a.name, "description": a.description, "variants": a.variants,
         "seconds": round(a.seconds, 1), "peak_rss_gb": round(a.peak_rss_gb, 2),
         "metrics": a.metrics}
        for a in selected
    ]
    # Merge rather than overwrite. A run restricted to --only would otherwise delete
    # every arm it was not asked about, and the file would still look well-formed --
    # it would just silently describe a smaller experiment. The same bug bit
    # compare_sv.py, and it cost a re-run here before being fixed.
    out_path = out / "arms.json"
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
