"""Tier-0 simulation: a diploid sample with exactly known phased genotypes.

See docs/simulation.md for the design and, more importantly, for what tier 0
cannot tell you. The short version: reads are simulated from the graph, so
mapping is unrealistically easy and absolute precision/recall will flatter every
caller. Tier 0 exists to compare callers to each other and to calibrate GQ.
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path

BASES = "ACGT"


@dataclass
class SimParams:
    """Everything that shapes the simulated dataset. Seeded for reproducibility."""

    ref_length: int = 200_000
    seed: int = 1

    # Variant density, as a probability per eligible position.
    snp_rate: float = 0.002
    indel_rate: float = 0.0004
    sv_rate: float = 0.00002

    indel_size_range: tuple[int, int] = (1, 20)
    sv_size_range: tuple[int, int] = (50, 2000)

    # Proportion of variants that are heterozygous rather than homozygous alt.
    het_fraction: float = 0.6

    # Minimum gap between variants. Variants close enough to interact make the
    # truth ambiguous: the pair can be written several ways and a haplotype-aware
    # comparison engine may legitimately resolve it differently from how we did.
    # This is enforced, not hoped for.
    min_variant_gap: int = 20

    contig: str = "sim"
    sample: str = "SIMSAMPLE"

    # Read simulation.
    depth: float = 30.0
    read_length: int = 150
    # vg sim's substitution error rate; indel error left at vg's default.
    error_rate: float = 0.01


@dataclass
class Variant:
    pos: int  # 1-based VCF position
    ref: str
    alt: str
    hap0: int  # 0 or 1
    hap1: int

    @property
    def size_change(self) -> int:
        return len(self.alt) - len(self.ref)

    @property
    def kind(self) -> str:
        d = abs(self.size_change)
        if d == 0 and len(self.ref) == 1:
            return "SNP"
        if d >= 50:
            return "SV"
        return "INDEL"


@dataclass
class SimulatedTruth:
    reference: str
    variants: list[Variant] = field(default_factory=list)

    def counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for v in self.variants:
            out[v.kind] = out.get(v.kind, 0) + 1
        return out


def _random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def _different_base(rng: random.Random, base: str) -> str:
    choices = [b for b in BASES if b != base.upper()]
    return rng.choice(choices)


def generate_truth(params: SimParams) -> SimulatedTruth:
    """Build a reference and a set of non-overlapping phased variants."""
    rng = random.Random(params.seed)
    reference = _random_sequence(rng, params.ref_length)

    variants: list[Variant] = []
    # Leave room at both ends so no variant runs off the contig and every variant
    # has flanking sequence for the mapper to anchor on.
    pos = 100
    limit = params.ref_length - 100

    while pos < limit:
        roll = rng.random()
        variant = None

        if roll < params.sv_rate:
            size = rng.randint(*params.sv_size_range)
            if pos + size + 1 >= limit:
                pos += 1
                continue
            if rng.random() < 0.5:
                # Deletion: REF spans the deleted bases plus one anchor base.
                ref_allele = reference[pos - 1 : pos + size]
                alt_allele = reference[pos - 1]
            else:
                # Insertion, anchored on the preceding base.
                ref_allele = reference[pos - 1]
                alt_allele = ref_allele + _random_sequence(rng, size)
            variant = (ref_allele, alt_allele)

        elif roll < params.sv_rate + params.indel_rate:
            size = rng.randint(*params.indel_size_range)
            if pos + size + 1 >= limit:
                pos += 1
                continue
            if rng.random() < 0.5:
                ref_allele = reference[pos - 1 : pos + size]
                alt_allele = reference[pos - 1]
            else:
                ref_allele = reference[pos - 1]
                alt_allele = ref_allele + _random_sequence(rng, size)
            variant = (ref_allele, alt_allele)

        elif roll < params.sv_rate + params.indel_rate + params.snp_rate:
            ref_allele = reference[pos - 1]
            alt_allele = _different_base(rng, ref_allele)
            variant = (ref_allele, alt_allele)

        if variant is None:
            pos += 1
            continue

        ref_allele, alt_allele = variant

        # Assign a phased genotype. Homozygous reference is not emitted: a truth
        # set records what the sample has, and 0|0 rows would just be noise.
        if rng.random() < params.het_fraction:
            hap0, hap1 = (1, 0) if rng.random() < 0.5 else (0, 1)
        else:
            hap0, hap1 = 1, 1

        variants.append(
            Variant(pos=pos, ref=ref_allele, alt=alt_allele, hap0=hap0, hap1=hap1)
        )

        # Skip past this variant plus the mandated gap.
        pos += max(len(ref_allele), len(alt_allele)) + params.min_variant_gap

    truth = SimulatedTruth(reference=reference, variants=variants)
    _validate(truth, params)
    return truth


def _validate(truth: SimulatedTruth, params: SimParams) -> None:
    """Assert the invariants that produce silently wrong truth if violated."""
    previous_end = 0
    for v in truth.variants:
        # REF must actually match the contig, or every downstream comparison is
        # against a truth set that does not describe this reference.
        actual = truth.reference[v.pos - 1 : v.pos - 1 + len(v.ref)]
        if actual != v.ref:
            raise AssertionError(
                f"REF mismatch at {v.pos}: VCF says {v.ref!r}, reference has {actual!r}"
            )
        # Variants must not overlap or abut.
        if v.pos <= previous_end + params.min_variant_gap:
            raise AssertionError(
                f"variant at {v.pos} is within {params.min_variant_gap}bp of the previous one"
            )
        previous_end = v.pos + len(v.ref) - 1
        if not (v.hap0 or v.hap1):
            raise AssertionError(f"variant at {v.pos} is homozygous reference")


def write_reference(truth: SimulatedTruth, path: Path, contig: str, width: int = 60) -> None:
    with open(path, "w") as fh:
        fh.write(f">{contig}\n")
        for i in range(0, len(truth.reference), width):
            fh.write(truth.reference[i : i + width] + "\n")


def write_truth_vcf(truth: SimulatedTruth, path: Path, params: SimParams) -> None:
    """Write the phased truth VCF. Written uncompressed; the caller bgzips it."""
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={params.contig},length={params.ref_length}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", params.sample]
        ),
    ]
    for v in truth.variants:
        lines.append(
            "\t".join(
                [
                    params.contig,
                    str(v.pos),
                    ".",
                    v.ref,
                    v.alt,
                    ".",
                    "PASS",
                    ".",
                    "GT",
                    f"{v.hap0}|{v.hap1}",
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n")


def write_confident_bed(path: Path, params: SimParams) -> None:
    """For tier 0 the whole contig is confident, by construction."""
    path.write_text(f"{params.contig}\t0\t{params.ref_length}\n")


def read_count_for_depth(params: SimParams) -> int:
    """Convert target depth to a read count, so configs specify depth not counts."""
    return max(1, int(params.depth * params.ref_length / params.read_length))
