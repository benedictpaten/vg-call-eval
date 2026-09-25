"""Recall, precision and F1 from benchmark counts, each rate in its own side's units.

aardvark and truvari both count true positives once per side, because one truth record can be
matched by several query records and one query record can match several truth records. The two
counts differ -- chr1's JointIndel row has 65,232 on the truth side and 72,989 on the query side --
and each rate takes the count from its own side:

    recall    = truth_tp / (truth_tp + truth_fn)
    precision = query_tp / (query_tp + query_fp)
    F1        = 2 * precision * recall / (precision + recall)

That is aardvark's `metric_f1` (truth_tp, query_tp) and truvari's `f1` (TP-base, TP-comp), and
tests/test_bench_metrics.py holds every per-contig score on disk to both, bit for bit. One TP used
for both halves pairs a count from one side with a denominator from the other. The mismatch does
not cancel between call sets either: the ratio of the two TP counts is a property of how a caller
writes its records, so it differs between callers and between read technologies.

An aggregate sums the four counts over contigs and then takes the rates. Averaging per-contig F1s
would weight chr21 like chr1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Counts:
    truth_tp: int = 0
    truth_fn: int = 0
    query_tp: int = 0
    query_fp: int = 0

    def __add__(self, other: "Counts") -> "Counts":
        return Counts(self.truth_tp + other.truth_tp, self.truth_fn + other.truth_fn,
                      self.query_tp + other.query_tp, self.query_fp + other.query_fp)

    @property
    def recall(self) -> float:
        n = self.truth_tp + self.truth_fn
        return self.truth_tp / n if n else math.nan

    @property
    def precision(self) -> float:
        n = self.query_tp + self.query_fp
        return self.query_tp / n if n else math.nan

    @property
    def f1(self) -> float:
        # NaN, not 0, when both rates are 0: aardvark writes NaN there, and this must match it.
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if p + r > 0 else math.nan


def f1(recall: float, precision: float) -> float:
    return 2 * precision * recall / (precision + recall) if precision + recall > 0 else math.nan


def aardvark_counts(row: dict | None) -> Counts:
    """One row of aardvark's summary.tsv, as read_summary returns it."""
    if not row:
        return Counts()
    return Counts(*(int(row.get(k, 0) or 0)
                    for k in ("truth_tp", "truth_fn", "query_tp", "query_fp")))


def truvari_counts(summary: dict | None) -> Counts:
    """truvari bench's summary.json."""
    if not summary:
        return Counts()
    return Counts(*(int(summary.get(k, 0) or 0) for k in ("TP-base", "FN", "TP-comp", "FP")))


def pick(rows, comparison: str, vtype: str) -> dict | None:
    """The aardvark row for one comparison and variant type, over all regions and filters."""
    for r in rows or []:
        if (r.get("comparison", "").upper() == comparison and r.get("variant_type") == vtype
                and r.get("region_label", "ALL") == "ALL" and r.get("filter", "ALL") == "ALL"):
            return r
    return None


def small(results, contigs, vtype: str, comparison: str = "GT") -> Counts:
    """Summed aardvark counts over bench_wgs.py's per-contig results."""
    total = Counts()
    for r in results:
        if r["contig"] in contigs:
            total += aardvark_counts(pick(r.get("aardvark"), comparison, vtype))
    return total


def sv(results, contigs) -> Counts:
    """Summed truvari counts over bench_wgs.py's per-contig results."""
    total = Counts()
    for r in results:
        if r["contig"] in contigs:
            total += truvari_counts(r.get("truvari"))
    return total
