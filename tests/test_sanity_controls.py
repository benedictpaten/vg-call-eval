"""The controls that make concordance numbers believable.

A comparison harness that is subtly wrong is worse than no harness, because it
produces confident numbers. These two tests are the gate: the engine must score
identical inputs perfectly, and must detect deliberate damage. Neither is
optional, and a run that has not passed them should not be quoted.
"""

from __future__ import annotations

import csv
import gzip
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from vgcalleval.engines import aardvark
from vgcalleval.simulate import (
    SimParams,
    generate_truth,
    write_confident_bed,
    write_reference,
    write_truth_vcf,
)

pytestmark = pytest.mark.skipif(
    not (shutil.which("aardvark") and shutil.which("bgzip")),
    reason="needs aardvark and bgzip on PATH",
)


def _gt_all(out_dir: Path) -> dict:
    for row in aardvark.read_summary(out_dir):
        if (row["comparison"], row["region_label"], row["filter"], row["variant_type"]) == (
            "GT", "ALL", "ALL", "ALL",
        ):
            return row
    raise AssertionError("no GT/ALL row in summary")


@pytest.fixture(scope="module")
def truth_set(tmp_path_factory):
    d = tmp_path_factory.mktemp("truth")
    params = SimParams(ref_length=20_000, seed=5)
    truth = generate_truth(params)
    write_reference(truth, d / "ref.fa", params.contig)
    write_truth_vcf(truth, d / "truth.vcf", params)
    write_confident_bed(d / "confident.bed", params)
    subprocess.run(["bgzip", "-f", str(d / "truth.vcf")], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(d / "truth.vcf.gz")], check=True)
    subprocess.run(["samtools", "faidx", str(d / "ref.fa")], check=True)
    return d, params


def test_identical_inputs_score_perfectly(truth_set, tmp_path):
    """Positive control. If this fails, the engine is misconfigured."""
    d, params = truth_set
    out = aardvark.compare(
        aardvark=shutil.which("aardvark"),
        reference=d / "ref.fa", truth_vcf=d / "truth.vcf.gz", query_vcf=d / "truth.vcf.gz",
        regions_bed=d / "confident.bed", out_dir=tmp_path / "pos",
        truth_sample=params.sample, query_sample=params.sample, label="pos",
    )
    row = _gt_all(out)
    assert float(row["metric_recall"]) == 1.0
    assert float(row["metric_precision"]) == 1.0
    assert int(row["truth_fn"]) == 0


def test_corrupted_calls_are_detected(truth_set, tmp_path):
    """Negative control. If this passes at 1.0, the harness cannot see errors and
    every number it has ever produced is meaningless."""
    d, params = truth_set
    rng = random.Random(0)
    corrupt = tmp_path / "corrupt.vcf"
    dropped = flipped = 0
    with gzip.open(d / "truth.vcf.gz", "rt") as fh, open(corrupt, "w") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            roll = rng.random()
            if roll < 0.25:
                dropped += 1
                continue
            if roll < 0.55:
                f = line.rstrip("\n").split("\t")
                f[9] = "1|1" if f[9] != "1|1" else "0|1"
                line = "\t".join(f) + "\n"
                flipped += 1
            out.write(line)
    assert dropped > 0 and flipped > 0, "corruption did nothing; test is vacuous"

    subprocess.run(["bgzip", "-f", str(corrupt)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(corrupt) + ".gz"], check=True)

    out = aardvark.compare(
        aardvark=shutil.which("aardvark"),
        reference=d / "ref.fa", truth_vcf=d / "truth.vcf.gz",
        query_vcf=Path(str(corrupt) + ".gz"),
        regions_bed=d / "confident.bed", out_dir=tmp_path / "neg",
        truth_sample=params.sample, query_sample=params.sample, label="neg",
    )
    row = _gt_all(out)
    assert float(row["metric_recall"]) < 0.95, "dropped calls were not counted as FN"
    assert float(row["metric_precision"]) < 0.95, "flipped genotypes were not counted as FP"
    assert int(row["truth_fn"]) >= dropped
