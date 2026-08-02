"""Aardvark adapter: the primary comparison engine.

Aardvark builds truth and query haplotype sequences within each sub-region and
searches for the phased zygosity assignment minimising edit distance between
them, so representation differences (left-alignment, anchor bases, multi-allelic
grouping) resolve correctly and the comparison is variant-type agnostic.

Two of its properties matter enough to call out:

* `--truth-sample` / `--query-sample` mean no bcftools renaming step.
* `--max-branch-factor` bounds the optimiser on dense variant regions, which is
  the guard that makes it plausible on pangenome-derived callsets.

There is no macOS or ARM release binary; on this platform aardvark is built from
source (`cargo build --release`). See docs/install.md.
"""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

from ..pipeline import run


@dataclass
class AardvarkOptions:
    min_variant_gap: int = 50
    max_branch_factor: int = 50
    enable_record_basepair: bool = False
    threads: int = 4

    @classmethod
    def for_structural_variants(cls, threads: int = 4) -> "AardvarkOptions":
        """Recommended SV settings: wider clustering, record-basepair scoring."""
        return cls(min_variant_gap=1000, enable_record_basepair=True, threads=threads)


def compare(
    *,
    aardvark: str,
    reference: Path,
    truth_vcf: Path,
    query_vcf: Path,
    regions_bed: Path,
    out_dir: Path,
    truth_sample: str,
    query_sample: str,
    label: str,
    options: AardvarkOptions | None = None,
) -> Path:
    """Run `aardvark compare`. Returns the output directory."""
    options = options or AardvarkOptions()
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        aardvark, "compare",
        "-r", str(reference),
        "-t", str(truth_vcf),
        "-q", str(query_vcf),
        "-b", str(regions_bed),
        "-o", str(out_dir),
        "--truth-sample", truth_sample,
        "--query-sample", query_sample,
        "--compare-label", label,
        "--min-variant-gap", str(options.min_variant_gap),
        "--max-branch-factor", str(options.max_branch_factor),
        "--threads", str(options.threads),
    ]
    if options.enable_record_basepair:
        cmd.append("--enable-record-basepair-metrics")

    run(cmd)
    return out_dir


def read_summary(out_dir: Path) -> list[dict]:
    """Parse aardvark's summary.tsv into a list of metric rows."""
    summary = out_dir / "summary.tsv"
    if not summary.exists():
        raise FileNotFoundError(f"no summary.tsv in {out_dir}")
    with open(summary) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_query_decisions(out_dir: Path) -> list[dict]:
    """Per-variant TP/FP decisions from the labelled query VCF.

    Aardvark annotates each query variant with a benchmark decision. GQ is NOT
    assumed to survive into this file -- the documentation does not promise it --
    so the caller joins GQ back from the original VCF on (chrom, pos, ref, alt).
    """
    candidates = list(out_dir.glob("*query*.vcf.gz")) + list(out_dir.glob("query.vcf.gz"))
    if not candidates:
        raise FileNotFoundError(f"no query VCF in {out_dir}: {[p.name for p in out_dir.iterdir()]}")
    path = candidates[0]

    records = []
    with gzip.open(path, "rt") as fh:
        sample_cols: list[str] = []
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_cols = line.rstrip("\n").split("\t")[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            fmt = fields[8].split(":")
            values = fields[9].split(":")
            annotations = dict(zip(fmt, values))
            records.append(
                {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "ref": fields[3],
                    "alt": fields[4],
                    "decision": _decision_of(annotations),
                    "annotations": annotations,
                }
            )
    return records


def _decision_of(annotations: dict) -> str:
    """Aardvark's benchmark decision field. Field name varies by version, so try
    the documented ones and fall back to UNKNOWN rather than guessing wrong."""
    for key in ("BD", "BENCHMARK_DECISION", "AD_BD"):
        if key in annotations:
            return annotations[key]
    return "UNKNOWN"
