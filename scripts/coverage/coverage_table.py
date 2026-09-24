#!/usr/bin/env python3
"""Tabulate the coverage sweep: one row per (technology, contig, depth), four-count F1 throughout.

Reads what sweep.py wrote: work/cov/<tech>-<contig>/score.<tag>.json (SV from truvari, phasing from
whatshap) and the aardvark summary.tsv behind each label. Median DP comes from the called VCF itself,
and the measured source depth from subsample_gaf.py's log. Emits TSV to stdout, or with --md the
per-technology tables that docs/coverage.md prints.
"""
import csv
import json
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path.home() / "PycharmProjects/vg-call-eval"


def f1(tpb, fn, tpc, fp):
    r = tpb / (tpb + fn) if tpb + fn else 0.0
    p = tpc / (tpc + fp) if tpc + fp else 0.0
    return 2 * p * r / (p + r) if p + r else 0.0


def aardvark(contig, label):
    p = REPO / f"work/tier2-{contig}-hap32/results/aardvark-{label}/summary.tsv"
    out = {}
    if not p.exists():
        return out
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["comparison"] == "GT" and r["region_label"] == "ALL" and r["filter"] == "ALL":
            c = [int(r[k]) for k in ("truth_tp", "truth_fn", "query_tp", "query_fp")]
            # aardvark's plain 'Indel' row is query-only (truth_total 0); the indel row is JointIndel
            rec = c[0] / (c[0] + c[1]) if c[0] + c[1] else float("nan")
            prec = c[2] / (c[2] + c[3]) if c[2] + c[3] else float("nan")
            out[r["variant_type"]] = (f1(*c), rec, prec)
    return out


def median_dp(vcf):
    q = subprocess.run(["bcftools", "query", "-f", "[%DP]\n", str(vcf)], capture_output=True, text=True).stdout
    v = [int(x) for x in q.split() if x.isdigit()]
    return statistics.median(v) if v else float("nan")


def source_depth(d):
    for lg in sorted(d.glob("subsample.*.log")):
        for line in lg.read_text().splitlines():
            if line.startswith("source:"):
                return float(line.rsplit("->", 1)[1].strip().rstrip("x"))
    return float("nan")


def rows():
    for d in sorted((REPO / "work/cov").glob("*-*")):
        tech, contig = d.name.split("-", 1)
        src = source_depth(d)
        tags = sorted((p.name[len("score."):-len(".json")] for p in d.glob("score.*.json")),
                      key=lambda t: float("inf") if t == "full" else float(t.rstrip("x")))
        for tag in tags:
            s = json.loads((d / f"score.{tag}.json").read_text())
            a = aardvark(contig, s["label"])
            sv = s.get("sv") or {}
            ph = s.get("phase") or {}
            depth = src if tag == "full" else float(tag.rstrip("x"))
            yield {
                "tech": tech, "contig": contig, "depth": f"{depth:.1f}" + ("*" if tag == "full" else ""),
                "median_DP": median_dp(d / f"{tag}.vcf.gz"),
                "ALL_F1": a.get("ALL", (float("nan"),) * 3)[0],
                "SNV_F1": a.get("Snv", (float("nan"),) * 3)[0],
                "SNV_R": a.get("Snv", (float("nan"),) * 3)[1],
                "SNV_P": a.get("Snv", (float("nan"),) * 3)[2],
                "INDEL_F1": a.get("JointIndel", (float("nan"),) * 3)[0],
                "INDEL_R": a.get("JointIndel", (float("nan"),) * 3)[1],
                "INDEL_P": a.get("JointIndel", (float("nan"),) * 3)[2],
                "SV_F1": sv.get("f1", float("nan")),
                "SV_R": sv.get("recall", float("nan")),
                "SV_P": sv.get("precision", float("nan")),
                "switch_pct": ph.get("switch_pct", "").rstrip("%"),
                "pairs": ph.get("pairs", ""),
            }


DOC_TABLES = [("sr", "Short reads (Illumina, 32-haplotype graph)"),
              ("ont", "ONT (16-haplotype E821 graph, `--preset ont`)")]
DOC_COLS = [("median DP", "median_DP"), ("ALL F1", "ALL_F1"), ("SNV R", "SNV_R"), ("SNV P", "SNV_P"),
            ("SNV F1", "SNV_F1"), ("indel R", "INDEL_R"), ("indel P", "INDEL_P"), ("indel F1", "INDEL_F1"),
            ("SV R", "SV_R"), ("SV P", "SV_P"), ("SV F1", "SV_F1"), ("switch %", "switch_pct"),
            ("pairs", "pairs")]


def doc_md(rs):
    """The per-technology tables in docs/coverage.md, exactly as the page prints them."""
    out = []
    for tech, title in DOC_TABLES:
        out += [f"#### {title}", "",
                "| contig | depth | " + " | ".join(h for h, _ in DOC_COLS) + " |",
                "|" + "---|" * (len(DOC_COLS) + 2)]
        for r in (r for r in rs if r["tech"] == tech):
            full = r["depth"].endswith("*")
            d = float(r["depth"].rstrip("*"))
            depth = f"**{d:.1f}x (full)**" if full else f"{d:g}x"
            cells = []
            for _, k in DOC_COLS:
                v = r[k]
                if k == "median_DP":
                    cells.append(f"{v:.0f}")
                elif k == "switch_pct":
                    cells.append(f"{float(v):.2f}")
                elif k == "pairs":
                    cells.append(f"{int(v):,}")
                else:
                    cells.append(f"{v:.4f}")
            out.append(f"| {r['contig']} | {depth} | " + " | ".join(cells) + " |")
        out.append("")
    return "\n".join(out).rstrip("\n")


def main():
    rs = list(rows())
    cols = list(rs[0]) if rs else []
    fmt = lambda v: f"{v:.4f}" if isinstance(v, float) and v == v and v < 1.01 else (
        f"{v:.0f}" if isinstance(v, float) and v == v else str(v))
    if "--md" in sys.argv:
        print(doc_md(rs))
    else:
        print("\t".join(cols))
        for r in rs:
            print("\t".join(fmt(r[c]) for c in cols))


if __name__ == "__main__":
    main()
