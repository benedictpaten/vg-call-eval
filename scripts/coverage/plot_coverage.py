#!/usr/bin/python3
"""Coverage curves for docs/coverage.md, as two figures of small multiples (columns = contig, short reads vs
ONT in every panel):

  coverage-curves-<theme>.png            SNV F1, indel F1, SV F1, switch error
  coverage-recall-precision-<theme>.png  recall and precision for SNVs, indels and SVs

Within a variant class and contig, the recall and precision panels share one y-range, so the gap between
the two is read off the axis rather than hidden by auto-scaling. Light and dark variants, served through
<picture> so dark mode is selected, not inverted. Palette: dataviz reference slots 1-2, validated (CVD dE
24.7 light / 26.8 dark).

Reads the TSV from coverage_table.py. Needs matplotlib (the system python3 has it).
Usage: coverage_table.py > summary.tsv; plot_coverage.py summary.tsv docs/figures
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

THEMES = {
    "light": dict(surface="#fcfcfb", text="#0b0b0b", text2="#52514e", grid="#e6e5e1",
                  series={"sr": "#2a78d6", "ont": "#eb6834"}),
    "dark": dict(surface="#1a1a19", text="#ffffff", text2="#c3c2b7", grid="#302f2d",
                 series={"sr": "#3987e5", "ont": "#d95926"}),
}
NAMES = {"sr": "Short reads", "ont": "ONT"}
CONTIGS = ["chr20", "chr6"]
TECHS = ("sr", "ont")
F1_ROWS = [("SNV_F1", "SNV F1"), ("INDEL_F1", "Indel F1"), ("SV_F1", "SV F1 (>= 50 bp)"),
           ("switch_pct", "Switch error (%)")]
RP_CLASSES = [("SNV", "SNV"), ("INDEL", "Indel"), ("SV", "SV (>= 50 bp)")]


def load(tsv):
    data = {}
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        pt = {"depth": float(r["depth"].rstrip("*")), "full": r["depth"].endswith("*")}
        for k, v in r.items():
            if k not in ("tech", "contig", "depth"):
                try:
                    pt[k] = float(v)
                except ValueError:
                    pass
        data.setdefault((r["tech"], r["contig"]), []).append(pt)
    for v in data.values():
        v.sort(key=lambda p: p["depth"])
    return data


def style_axes(ax, t):
    ax.set_facecolor(t["surface"])
    ax.grid(True, axis="y", color=t["grid"], linewidth=1, linestyle="-")
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)


def panel(ax, t, data, contig, metric, ylim=None, min_span=None):
    """Two series (short reads, ONT): 2px lines, ringed markers, a diamond at full source depth, and
    direct labels placed where the two lines are furthest apart."""
    style_axes(ax, t)
    series = {}
    for tech in TECHS:
        pts = [p for p in data.get((tech, contig), []) if metric in p]
        if not pts:
            continue
        series[tech] = pts
        c = t["series"][tech]
        ax.plot([p["depth"] for p in pts], [p[metric] for p in pts], color=c, linewidth=2,
                solid_capstyle="round", solid_joinstyle="round", zorder=2)
        sub = [p for p in pts if not p["full"]]
        ax.scatter([p["depth"] for p in sub], [p[metric] for p in sub], s=36, color=c,
                   edgecolor=t["surface"], linewidth=2, zorder=3)
        full = [p for p in pts if p["full"]]
        ax.scatter([p["depth"] for p in full], [p[metric] for p in full], s=64, marker="D", color=c,
                   edgecolor=t["surface"], linewidth=2, zorder=4)
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif min_span is not None:
        lo, hi = ax.get_ylim()
        if hi - lo < min_span:
            mid = (hi + lo) / 2
            ax.set_ylim(mid - min_span / 2, mid + min_span / 2)
    if len(series) == 2:
        a, b = series["sr"], series["ont"]
        shared = sorted(set(p["depth"] for p in a) & set(p["depth"] for p in b))
        at = lambda s, d: next(p[metric] for p in s if p["depth"] == d)
        d = max(shared, key=lambda x: abs(at(a, x) - at(b, x)))
        ya, yb = at(a, d), at(b, d)
        xs_all = sorted(set(p["depth"] for v in series.values() for p in v))
        # anchor inward at the plot edges so a label never runs into the tick labels
        ha, dx = ("left", 4) if d == xs_all[0] else (("right", -4) if d == xs_all[-1] else ("center", 0))
        for tech, y, above in (("sr", ya, ya >= yb), ("ont", yb, ya < yb)):
            ax.annotate(NAMES[tech], (d, y), xytext=(dx, 9 if above else -9), textcoords="offset points",
                        ha=ha, va="bottom" if above else "top", fontsize=9.5, color=t["text2"])
    return series


def legend(fig, t):
    handles = [plt.Line2D([], [], color=t["series"][k], linewidth=2, marker="o", markersize=6,
                          markeredgecolor=t["surface"], label=NAMES[k]) for k in TECHS]
    handles.append(plt.Line2D([], [], color=t["text2"], linewidth=0, marker="D", markersize=6,
                              label="full source depth (short 30x, ONT 43-45x)"))
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, 0.995), labelcolor=t["text2"])


def setup(theme, nrows, height):
    t = THEMES[theme]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": t["grid"],
                         "axes.labelcolor": t["text2"], "xtick.color": t["text2"], "ytick.color": t["text2"],
                         "text.color": t["text"]})
    fig, axes = plt.subplots(nrows, len(CONTIGS), figsize=(10, height), sharex=True)
    fig.patch.set_facecolor(t["surface"])
    return t, fig, axes


def finish(fig, t, axes, out):
    for j, contig in enumerate(CONTIGS):
        axes[0][j].set_title(contig, fontsize=12, color=t["text"], pad=10, loc="left", fontweight="bold")
        axes[-1][j].set_xlabel("Read depth (x)", fontsize=10.5)
        axes[-1][j].set_xticks([5, 10, 15, 20, 25, 30, 40, 45])
    legend(fig, t)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    fig.savefig(out, dpi=150, facecolor=t["surface"])
    plt.close(fig)


def draw_f1(theme, data, out):
    t, fig, axes = setup(theme, len(F1_ROWS), 12.5)
    for j, contig in enumerate(CONTIGS):
        for i, (metric, title) in enumerate(F1_ROWS):
            ax = axes[i][j]
            # SV differences under ~0.01 are not significant here (paired bootstrap); a 0.05 floor on the
            # visible span keeps that noise from reading as a trend.
            panel(ax, t, data, contig, metric, min_span=0.05 if metric == "SV_F1" else None)
            if metric == "INDEL_F1" and ("ont", contig) in data:
                pk = max(data[("ont", contig)], key=lambda p: p[metric])
                ax.annotate(f"ONT peak {pk['depth']:g}x", (pk["depth"], pk[metric]), xytext=(8, 16),
                            textcoords="offset points", fontsize=9, color=t["text2"],
                            arrowprops=dict(arrowstyle="-", color=t["text2"], linewidth=0.8))
            if j == 0:
                ax.set_ylabel(title, fontsize=10.5)
    finish(fig, t, axes, out)


def draw_rp(theme, data, out):
    rows = [(cls, m, f"{label} {'recall' if m == 'R' else 'precision'}")
            for cls, label in RP_CLASSES for m in ("R", "P")]
    t, fig, axes = setup(theme, len(rows), 17)
    for j, contig in enumerate(CONTIGS):
        for cls, _ in RP_CLASSES:
            # one y-range for recall and precision of this class on this contig
            vals = [p[f"{cls}_{m}"] for tech in TECHS for p in data.get((tech, contig), [])
                    for m in ("R", "P") if f"{cls}_{m}" in p]
            lo, hi = min(vals), max(vals)
            pad = max((hi - lo) * 0.08, 0.004)
            ylim = (lo - pad, hi + pad)
            for i, (rc, m, title) in enumerate(rows):
                if rc != cls:
                    continue
                ax = axes[i][j]
                panel(ax, t, data, contig, f"{cls}_{m}", ylim=ylim)
                if j == 0:
                    ax.set_ylabel(title, fontsize=10.5)
    finish(fig, t, axes, out)


if __name__ == "__main__":
    data = load(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    for theme in THEMES:
        for name, fn in (("coverage-curves", draw_f1), ("coverage-recall-precision", draw_rp)):
            fn(theme, data, outdir / f"{name}-{theme}.png")
            print(f"wrote {outdir / f'{name}-{theme}.png'}")
