#!/usr/bin/env python3
"""Tabulate the mismapping floor/cap sweeps for both graphs.

Reads whatever `aardvark-sweep-*` and `aardvark-cap*` directories exist under each
results dir, so it can be run while a sweep is still going and will simply show fewer
rows. Reports GT and BASEPAIR together: §9.15 found the two knobs pull against each
other, with the floor that is best for indel GT costing SNVs and basepair accuracy, so a
single F1 column would hide the trade rather than settle it.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

FLOOR_RE = re.compile(r"^aardvark-sweep-([\d.]+)$")
CAP_RE = re.compile(r"^aardvark-cap([\d.]+)-fl([\d.]+)$")


def metrics(d: Path) -> dict:
    p = d / "summary.tsv"
    if not p.exists():
        return {}
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r.get("comparison", "").upper(), r.get("variant_type"))
            if r.get("region_label") == "ALL" and r.get("filter") == "ALL":
                out[key] = r
    return out


def val(rows, comparison, vtype, key):
    r = rows.get((comparison, vtype))
    if not r or r.get(key) in ("", None):
        return None
    try:
        return float(r[key])
    except ValueError:
        return None


def fmt(x, nd=4):
    return "  --  " if x is None else f"{x:.{nd}f}"


def emit(title: str, rows: list[tuple[str, dict]], default_tag: str | None) -> None:
    if not rows:
        return
    print(f"\n=== {title} ===")
    print(f"{'setting':<16}{'GT rec':>9}{'GT prec':>9}{'GT F1':>9}"
          f"{'SNV F1':>9}{'INS F1':>9}{'DEL F1':>9}{'BP F1':>9}")
    best = None
    for tag, m in rows:
        gt_f1 = val(m, "GT", "ALL", "metric_f1")
        if gt_f1 is not None and (best is None or gt_f1 > best[1]):
            best = (tag, gt_f1)
    for tag, m in rows:
        mark = ""
        if best and tag == best[0]:
            mark = "  <- best GT F1"
        if default_tag and tag == default_tag:
            mark += "  (current default)"
        print(f"{tag:<16}"
              f"{fmt(val(m,'GT','ALL','metric_recall')):>9}"
              f"{fmt(val(m,'GT','ALL','metric_precision')):>9}"
              f"{fmt(val(m,'GT','ALL','metric_f1')):>9}"
              f"{fmt(val(m,'GT','Snv','metric_f1')):>9}"
              f"{fmt(val(m,'GT','Insertion','metric_f1')):>9}"
              f"{fmt(val(m,'GT','Deletion','metric_f1')):>9}"
              f"{fmt(val(m,'BASEPAIR','ALL','metric_f1')):>9}{mark}")


def collect(res: Path):
    floors, caps = [], []
    for d in sorted(res.glob("aardvark-*")):
        if not d.is_dir():
            continue
        m = FLOOR_RE.match(d.name)
        if m:
            floors.append((float(m.group(1)), f"floor {m.group(1)}", metrics(d)))
            continue
        m = CAP_RE.match(d.name)
        if m:
            caps.append(((float(m.group(1)), float(m.group(2))),
                         f"cap {m.group(1)} fl {m.group(2)}", metrics(d)))
    floors.sort(key=lambda x: x[0])
    caps.sort(key=lambda x: x[0])
    return ([(t, m) for _, t, m in floors if m],
            [(t, m) for _, t, m in caps if m])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--repo", default=str(Path(__file__).resolve().parent.parent.parent))
    args = p.parse_args()
    repo = Path(args.repo)

    for label, sub in [("4-haplotype graph", "work/tier2-chr20/results"),
                       ("32-haplotype graph", "work/tier2-chr20-hap32/results")]:
        res = repo / sub
        floors, caps = collect(res)
        emit(f"{label}: --mismap-min sweep (cap at default 0.1)", floors, "floor 0.01")
        emit(f"{label}: --mismap-max sweep", caps, None)


if __name__ == "__main__":
    main()
