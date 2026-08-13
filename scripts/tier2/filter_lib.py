#!/usr/bin/env python3
"""Shared loading for the filter experiments.

Pulled out of coverage_model.py so the depth sweep and the GQ adjustment see exactly the
same calls, labels and local-depth baseline. Two things here are easy to get wrong and
both were:

  - the local depth median must be taken over *all* calls, not only labelled ones. With
    truvari labels only a few hundred calls carry a label, and a rolling median over those
    is a median over megabases rather than a neighbourhood.
  - the truth totals must come from the original run, not from the labelled calls, or a
    filter's recall is measured against a denominator the filter itself shrank.
"""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path


def labels(work: Path, kind: str) -> dict:
    bd = {}
    if kind == "truvari":
        for fn, lab in (("tp-comp.vcf.gz", "TP"), ("fp.vcf.gz", "FP")):
            with gzip.open(work / "results/truvari-readlik" / fn, "rt") as fh:
                for line in fh:
                    if not line.startswith("#"):
                        bd[int(line.split("\t", 2)[1])] = lab
    else:
        with gzip.open(work / "results/aardvark-readlik/query.vcf.gz", "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                v = dict(zip(f[8].split(":"), f[9].split(":"))).get("BD")
                if v in ("TP", "FP"):
                    bd[int(f[1])] = v
    return bd


def truth_counts(work: Path, kind: str) -> tuple[int, int]:
    """(base-side true positives, truth total) from the *unfiltered* run.

    Both numbers are needed, and the base-side count is the one that is easy to get wrong.
    A filter can only be scored on the query side -- it keeps or drops query records -- but
    recall is a base-side quantity, and the two counts are not equal: truvari matches
    several base records to one query record when the caller emitted a structural variant
    as a single record and the benchmark decomposed it. On chr6 4-hap that is 792 base
    matches against 425 query records. Dividing the query count by the truth total
    understates recall by a factor of nearly two, which makes every filter look far more
    damaging than it is. The fix is to scale: a filter that retains a fraction f of query
    true positives is credited with f of the base-side matches.
    """
    if kind == "truvari":
        d = json.loads((work / "results/truvari-readlik/summary.json").read_text())
        return int(d["TP-base"]), int(d["TP-base"]) + int(d["FN"])
    for arm in json.loads((work / "results/arms.json").read_text()):
        if arm["arm"] != "readlik":
            continue
        for m in (arm["metrics"].get("summary") or []):
            if (m["region_label"], m["variant_type"], m["filter"]) == ("ALL", "ALL", "ALL") \
                    and m["comparison"] == "GT":
                return int(m["truth_tp"]), int(m["truth_total"])
    raise SystemExit(f"no truth counts for {work} {kind}")


def rolling_median(vals: list[float], window: int = 201) -> list[float]:
    half = window // 2
    out = []
    for i in range(len(vals)):
        w = sorted(vals[max(0, i - half):min(len(vals), i + half + 1)])
        out.append(w[len(w) // 2] if w else 0.0)
    return out


def collect(work: Path, kind: str) -> list[dict]:
    """One record per labelled call: label, DP, DP/local median, explained share, GQ."""
    bd = labels(work, kind)
    vcf = str(work / "results/readlik.vcf.gz")
    # GQI only exists in files from a build that emits it. Ask for it, and fall back
    # rather than failing, so this works against older result sets too.
    has_gqi = subprocess.run(["bcftools", "view", "-h", vcf],
                             capture_output=True, text=True).stdout.find("ID=GQI") >= 0
    fmt = "%POS[\t%DP\t%AD\t%GQ" + ("\t%GQI" if has_gqi else "") + "]\n"
    q = subprocess.run(["bcftools", "query", "-f", fmt, vcf],
                       capture_output=True, text=True)
    if q.returncode != 0:
        raise SystemExit(q.stderr.strip()[:400])

    pos, dps, shares, gqs, gqis = [], [], [], [], []
    for line in q.stdout.splitlines():
        f = line.split("\t")
        try:
            dp = float(f[1])
            ad = [int(x) for x in f[2].split(",")]
            gq = float(f[3])
            gqi = float(f[4]) if has_gqi else gq
        except (ValueError, IndexError):
            continue
        if dp <= 0:
            continue
        pos.append(int(f[0]))
        dps.append(dp)
        shares.append(min(1.0, sum(ad) / dp))
        gqs.append(gq)
        gqis.append(gqi)

    med = rolling_median(dps)
    return [{"pos": p, "label": bd[p], "dp": d, "ratio": d / m if m > 0 else 1.0,
             "share": s, "gq": g, "gqi": gi}
            for p, d, s, g, gi, m in zip(pos, dps, shares, gqs, gqis, med) if p in bd]


# Absolute, so these scripts work from any directory. A relative path here silently
# produced four empty result tables rather than an error, because the caller was catching
# the missing-file exception in order to skip datasets that genuinely lack one benchmark.
REPO = Path(__file__).resolve().parent.parent.parent

DATASETS = [("chr6 4-hap", REPO / "work/tier2-chr6"),
            ("chr6 34-hap", REPO / "work/tier2-chr6-hap32"),
            ("chr20 4-hap", REPO / "work/tier2-chr20"),
            ("chr20 34-hap", REPO / "work/tier2-chr20-hap32")]


def prf(tp: int, fp: int, tp_all: int, base_tp: int, total: int):
    """Precision from the kept query calls; recall from the base-side count scaled by the
    fraction of query true positives kept. See truth_counts for why the scaling is there."""
    p = tp / (tp + fp) if tp + fp else 0.0
    r = (base_tp * (tp / tp_all) / total) if tp_all and total else 0.0
    return p, r, (2 * p * r / (p + r) if p + r else 0.0)
