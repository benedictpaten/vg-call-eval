"""Build a tier-0 dataset, run the caller matrix over it, and compare to truth.

Every step writes into a run directory and is skipped if its output already
exists and is newer than its inputs, so re-running is cheap and a failed run can
be resumed. Nothing here imports vg; it drives the binary given in the config, so
two builds can be compared in one matrix.
"""

from __future__ import annotations

import json
import os
import re
import resource
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path

from .simulate import (
    SimParams,
    generate_truth,
    read_count_for_depth,
    write_confident_bed,
    write_reference,
    write_truth_vcf,
)


@dataclass
class StepResult:
    command: str
    seconds: float
    peak_rss_mb: float
    returncode: int


def _peak_rss_mb() -> float:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    # Linux reports KB, macOS reports bytes.
    return usage / (1024 * 1024) if usage > 10**7 else usage / 1024


def run(cmd: list[str] | str, *, stdout: Path | None = None, cwd: Path | None = None,
        env: dict | None = None) -> StepResult:
    """Run a command, capturing wall-clock and peak RSS. Raises on failure."""
    shell = isinstance(cmd, str)
    printable = cmd if shell else " ".join(str(c) for c in cmd)

    before = _peak_rss_mb()
    start = time.time()
    out_handle = open(stdout, "wb") if stdout else None
    try:
        proc = subprocess.run(
            cmd, shell=shell, stdout=out_handle or subprocess.PIPE,
            stderr=subprocess.PIPE, cwd=cwd, env=env,
        )
    finally:
        if out_handle:
            out_handle.close()
    elapsed = time.time() - start
    peak = max(0.0, _peak_rss_mb() - before)

    if proc.returncode != 0:
        stderr = proc.stderr.decode(errors="replace")[-4000:]
        raise RuntimeError(f"command failed ({proc.returncode}):\n  {printable}\n{stderr}")

    return StepResult(printable, elapsed, peak, proc.returncode)


def _stale(output: Path, inputs: list[Path]) -> bool:
    """make-style staleness: rebuild if missing or older than any input."""
    if not output.exists():
        return True
    out_mtime = output.stat().st_mtime
    return any(i.exists() and i.stat().st_mtime > out_mtime for i in inputs)


# ---------------------------------------------------------------------------
# Dataset construction
# ---------------------------------------------------------------------------


@dataclass
class Dataset:
    directory: Path
    reference: Path
    truth_vcf: Path
    confident_bed: Path
    gbz: Path
    mapped_gam: Path
    simulated_gam: Path
    pack: Path
    params: SimParams

    def write_manifest(self) -> None:
        (self.directory / "manifest.json").write_text(
            json.dumps({"params": asdict(self.params)}, indent=2, default=str)
        )


def build_tier0_dataset(directory: Path, params: SimParams, vg: str, threads: int = 4) -> Dataset:
    directory.mkdir(parents=True, exist_ok=True)

    reference = directory / "ref.fa"
    truth_plain = directory / "truth.vcf"
    truth_vcf = directory / "truth.vcf.gz"
    bed = directory / "confident.bed"

    if _stale(truth_vcf, []):
        truth = generate_truth(params)
        write_reference(truth, reference, params.contig)
        write_truth_vcf(truth, truth_plain, params)
        write_confident_bed(bed, params)
        run(["bgzip", "-f", str(truth_plain)])
        run(["tabix", "-f", "-p", "vcf", str(truth_vcf)])
        run(["samtools", "faidx", str(reference)])
        counts = truth.counts()
        (directory / "truth_counts.json").write_text(json.dumps(counts, indent=2))

    # Graph + indexes. autoindex gives a GBZ carrying the sample's two phased
    # haplotypes as paths, which is what we simulate from and also what -g/-z
    # enumeration will later use as alleles.
    gbz = directory / "idx.giraffe.gbz"
    if _stale(gbz, [truth_vcf, reference]):
        run(
            [vg, "autoindex", "-r", str(reference), "-v", str(truth_vcf),
             "-w", "giraffe", "-p", str(directory / "idx"), "-t", str(threads)],
        )

    simulated = directory / "simulated.gam"
    mapped = directory / "mapped.gam"
    pack = directory / "reads.pack"

    if _stale(simulated, [gbz]):
        hap_paths = _haplotype_paths(vg, gbz, params.sample)
        if len(hap_paths) < 2:
            raise RuntimeError(
                f"expected 2 haplotype paths for sample {params.sample}, found {hap_paths}"
            )
        total_reads = read_count_for_depth(params)
        per_hap = max(1, total_reads // 2)
        parts = []
        for i, path_name in enumerate(hap_paths[:2]):
            part = directory / f"sim_hap{i}.gam"
            run(
                [vg, "sim", "-x", str(gbz), "-P", path_name, "-n", str(per_hap),
                 "-l", str(params.read_length), "-a", "-s", str(params.seed + 100 + i),
                 "-e", str(params.error_rate)],
                stdout=part,
            )
            parts.append(part)
        with open(simulated, "wb") as out:
            for part in parts:
                with open(part, "rb") as fh:
                    shutil.copyfileobj(fh, out)

    # Re-map. Calling on vg sim's own alignments would assume away every mapping
    # error, including the ones the mismapping term exists to handle.
    if _stale(mapped, [simulated, gbz]):
        reads_fq = directory / "reads.fq"
        run([vg, "view", "-X", str(simulated)], stdout=reads_fq)
        run([vg, "giraffe", "-Z", str(gbz), "-f", str(reads_fq), "-t", str(threads)],
            stdout=mapped)

    if _stale(pack, [mapped, gbz]):
        run([vg, "pack", "-x", str(gbz), "-g", str(mapped), "-o", str(pack), "-t", str(threads)])

    dataset = Dataset(
        directory=directory, reference=reference, truth_vcf=truth_vcf, confident_bed=bed,
        gbz=gbz, mapped_gam=mapped, simulated_gam=simulated, pack=pack, params=params,
    )
    dataset.write_manifest()
    return dataset


def _haplotype_paths(vg: str, gbz: Path, sample: str) -> list[str]:
    """Find the sample's haplotype paths in the GBZ, rather than assuming names."""
    proc = subprocess.run([vg, "paths", "-L", "-x", str(gbz)], capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.decode(errors="replace"))
    names = [n for n in proc.stdout.decode().splitlines() if n.strip()]
    # PanSN-ish: SAMPLE#HAP#CONTIG#FRAGMENT
    matching = [n for n in names if n.split("#")[0] == sample]
    matching.sort()
    return matching


# ---------------------------------------------------------------------------
# Caller matrix
# ---------------------------------------------------------------------------


@dataclass
class Arm:
    name: str
    vg: str
    extra_args: list[str]
    needs_reads: bool = False
    needs_pack: bool = True
    description: str = ""


def default_arms(vg: str, vg_depthfix: str | None = None) -> list[Arm]:
    arms = [
        Arm("poisson", vg, [], description="current default, as shipped"),
        Arm("readlik", vg, ["--read-likelihood"], needs_reads=True,
            description="read-level likelihood model"),
        Arm("readlik-nomismap", vg, ["--read-likelihood", "--no-mismap-term"], needs_reads=True,
            description="mismapping term disabled, to measure its contribution"),
        Arm("readlik-gbwt-nopack", vg, ["--read-likelihood", "-z"], needs_reads=True,
            needs_pack=False, description="haplotype enumeration, no pack file"),
    ]
    if vg_depthfix:
        arms.insert(
            1,
            Arm("poisson-depthfix", vg_depthfix, [],
                description="Poisson with the depth_err one-liner patched"),
        )
    return arms


def run_arm(arm: Arm, dataset: Dataset, out_dir: Path, threads: int = 4) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    vcf = out_dir / f"{arm.name}.vcf"

    cmd = [arm.vg, "call", str(dataset.gbz), "-s", dataset.params.sample, "-t", str(threads)]
    if arm.needs_pack:
        cmd += ["-k", str(dataset.pack)]
    if arm.needs_reads:
        cmd += ["--gam", str(dataset.mapped_gam)]
    cmd += arm.extra_args

    inputs = [dataset.gbz, dataset.pack, dataset.mapped_gam]
    if _stale(vcf, inputs):
        result = run(cmd, stdout=vcf)
        perf = {"seconds": result.seconds, "peak_rss_mb": result.peak_rss_mb}
        (out_dir / f"{arm.name}.perf.json").write_text(json.dumps(perf, indent=2))
    else:
        perf = json.loads((out_dir / f"{arm.name}.perf.json").read_text())

    gz = out_dir / f"{arm.name}.norm.vcf.gz"
    if _stale(gz, [vcf]):
        # Normalise before comparison: split multi-allelics and left-align, or
        # representation differences show up as caller differences.
        run(f"bcftools norm -f {dataset.reference} -m -any -Oz -o {gz} {vcf} 2>/dev/null")
        run(["tabix", "-f", "-p", "vcf", str(gz)])

    return {
        "arm": arm.name,
        "description": arm.description,
        "vg_version": vg_version(arm.vg),
        "vcf": str(gz),
        "command": " ".join(cmd),
        **perf,
    }


def vg_version(vg: str) -> str:
    proc = subprocess.run([vg, "version"], capture_output=True)
    first = proc.stdout.decode(errors="replace").splitlines()
    return first[0].strip() if first else "unknown"
