"""Per-contig rates from bench_metrics must be aardvark's and truvari's own, bit for bit.

The genome-wide figures are sums of per-contig counts, so the only way to know the sum means what
the tools mean is to hold each contig's rates to theirs. Data-driven, like test_switch_bed.py:
every bench_wgs.py score under work/ is checked, and a newly scored arm is picked up without
editing this file.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import bench_metrics as bm  # noqa: E402


def scored(path: Path) -> bool:
    """False for a run whose every contig failed to score -- it has nothing to check."""
    return any(r.get("aardvark") or r.get("truvari") for r in json.loads(path.read_text()))


CORPUS = [p for p in sorted((ROOT / "work").glob("*/score/per-contig.json")) if scored(p)]


def same(ours: float, theirs: str | float | None) -> bool:
    """Exact equality, with aardvark's empty or NaN cell matching an undefined rate."""
    theirs = math.nan if theirs in ("", None) else float(theirs)
    return ours == theirs or (math.isnan(ours) and math.isnan(theirs))


def test_each_rate_takes_its_own_sides_tp():
    # chr1 JointIndel from work/wgs-mq0, where the two TP counts are 65,232 and 72,989.
    c = bm.Counts(truth_tp=65232, truth_fn=4608, query_tp=72989, query_fp=5710)
    assert c.recall == 65232 / 69840
    assert c.precision == 72989 / 78699
    single_tp = 2 * 65232 / (2 * 65232 + 5710 + 4608)
    assert round(c.f1, 4) == 0.9307 and round(single_tp, 4) == 0.9267


def test_aggregate_sums_counts_not_rates():
    results = [{"contig": "a", "aardvark": [{"comparison": "GT", "region_label": "ALL",
                                             "filter": "ALL", "variant_type": "ALL",
                                             "truth_tp": "9", "truth_fn": "1",
                                             "query_tp": "8", "query_fp": "2"}],
                "truvari": {"TP-base": 3, "FN": 1, "TP-comp": 2, "FP": 2}},
               {"contig": "b", "aardvark": [{"comparison": "GT", "region_label": "ALL",
                                             "filter": "ALL", "variant_type": "ALL",
                                             "truth_tp": "1", "truth_fn": "9",
                                             "query_tp": "1", "query_fp": "0"}],
                "truvari": {"TP-base": 1, "FN": 3, "TP-comp": 1, "FP": 0}}]
    assert bm.small(results, {"a", "b"}, "ALL") == bm.Counts(10, 10, 9, 2)
    assert bm.sv(results, {"a", "b"}) == bm.Counts(4, 4, 3, 2)
    assert bm.small(results, {"a"}, "ALL") == bm.Counts(9, 1, 8, 2)


def test_both_rates_zero_is_nan_as_aardvark_writes_it():
    assert math.isnan(bm.Counts(0, 55, 0, 1456).f1)


@pytest.mark.skipif(not CORPUS, reason="no bench_wgs.py scores in work/")
@pytest.mark.parametrize("path", CORPUS, ids=lambda p: str(p.relative_to(ROOT)))
def test_reproduces_the_tools(path):
    for r in json.loads(path.read_text()):
        for row in r.get("aardvark") or []:
            c = bm.aardvark_counts(row)
            where = f"{r['contig']} {row['comparison']} {row['variant_type']}"
            assert same(c.f1, row["metric_f1"]), where
            assert same(c.recall, row["metric_recall"]), where
            assert same(c.precision, row["metric_precision"]), where
        t = r.get("truvari")
        if t and t.get("f1") is not None:
            c = bm.truvari_counts(t)
            assert (c.f1, c.recall, c.precision) == (t["f1"], t["recall"], t["precision"]), \
                r["contig"]


@pytest.mark.skipif(not CORPUS, reason="no bench_wgs.py scores in work/")
def test_the_single_tp_formula_would_have_failed():
    """The check above can fire: under one TP for both halves most scored rows disagree."""
    wrong = total = 0
    for path in CORPUS:
        for r in json.loads(path.read_text()):
            for row in r.get("aardvark") or []:
                c = bm.aardvark_counts(row)
                if row["metric_f1"] in ("", "NaN") or not c.truth_tp:
                    continue
                total += 1
                single = 2 * c.truth_tp / (2 * c.truth_tp + c.query_fp + c.truth_fn)
                wrong += single != float(row["metric_f1"])
    assert total and wrong > total // 2
