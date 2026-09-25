#!/usr/bin/python3
"""ONT indel accuracy against depth before and after --hp-prior, for docs/ont-hp-prior.md.

Small multiples: rows indel F1, precision, recall; columns chr20 (fitted) and chr6 (held out). The
two series are the same caller with and without the raised frequency exponent at run-length sites,
so each panel's gap is the change. Light and dark variants, served through <picture>.

Reads work/run-cov/hp-prior-series.tsv (contig, depth, arm, indel_f1, indel_r, indel_p, ...).
Usage: plot_hp_prior.py <series.tsv> <outdir>
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

THEMES = {
    "light": dict(surface="#fcfcfb", text="#0b0b0b", text2="#52514e", grid="#e6e5e1",
                  series={"before": "#2a78d6", "hp-prior": "#eb6834"}),
    "dark": dict(surface="#1a1a19", text="#ffffff", text2="#c3c2b7", grid="#302f2d",
                 series={"before": "#3987e5", "hp-prior": "#d95926"}),
}
NAMES = {"before": "--preset ont before", "hp-prior": "with --hp-prior 20"}
ARMS = ("before", "hp-prior")
CONTIGS = [("chr20", "chr20 (fitted)"), ("chr6", "chr6 (held out)")]
ROWS = [("indel_f1", "Indel F1"), ("indel_p", "Indel precision"), ("indel_r", "Indel recall")]


def load(tsv):
    data = {}
    for r in csv.DictReader(open(tsv), delimiter="\t"):
        pt = {"depth": float(r["depth"].rstrip("*")), "full": r["depth"].endswith("*")}
        pt.update({k: float(r[k]) for k, _ in ROWS})
        data.setdefault((r["arm"], r["contig"]), []).append(pt)
    for v in data.values():
        v.sort(key=lambda p: p["depth"])
    return data


def draw(theme, data, out):
    t = THEMES[theme]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.edgecolor": t["grid"],
                         "axes.labelcolor": t["text2"], "xtick.color": t["text2"], "ytick.color": t["text2"],
                         "text.color": t["text"]})
    fig, axes = plt.subplots(len(ROWS), len(CONTIGS), figsize=(10, 9.5), sharex=True)
    fig.patch.set_facecolor(t["surface"])
    for j, (contig, title) in enumerate(CONTIGS):
        for i, (metric, label) in enumerate(ROWS):
            ax = axes[i][j]
            ax.set_facecolor(t["surface"])
            ax.grid(True, axis="y", color=t["grid"], linewidth=1)
            ax.set_axisbelow(True)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
            ax.tick_params(length=0)
            for arm in ARMS:
                pts = data[(arm, contig)]
                c = t["series"][arm]
                ax.plot([p["depth"] for p in pts], [p[metric] for p in pts], color=c, linewidth=2,
                        solid_capstyle="round", zorder=2)
                sub = [p for p in pts if not p["full"]]
                ax.scatter([p["depth"] for p in sub], [p[metric] for p in sub], s=36, color=c,
                           edgecolor=t["surface"], linewidth=2, zorder=3)
                full = [p for p in pts if p["full"]]
                ax.scatter([p["depth"] for p in full], [p[metric] for p in full], s=64, marker="D", color=c,
                           edgecolor=t["surface"], linewidth=2, zorder=4)
            # direct labels on the top row only, at the full-depth end where the lines are furthest apart
            a, b = data[("before", contig)][-1], data[("hp-prior", contig)][-1]
            labels = () if i else (("before", a, a[metric] >= b[metric]),
                                   ("hp-prior", b, b[metric] > a[metric]))
            for arm, p, above in labels:
                ax.annotate(NAMES[arm], (p["depth"], p[metric]), xytext=(-4, 9 if above else -9),
                            textcoords="offset points", ha="right", va="bottom" if above else "top",
                            fontsize=9.5, color=t["text2"])
            if i == 0:
                ax.set_title(title, fontsize=12, color=t["text"], pad=10, loc="left", fontweight="bold")
            if j == 0:
                ax.set_ylabel(label, fontsize=10.5)
            if i == len(ROWS) - 1:
                ax.set_xlabel("ONT read depth (x)", fontsize=10.5)
                ax.set_xticks([5, 10, 15, 20, 25, 30, 40, 45])
    handles = [plt.Line2D([], [], color=t["series"][k], linewidth=2, marker="o", markersize=6,
                          markeredgecolor=t["surface"], label=NAMES[k]) for k in ARMS]
    handles.append(plt.Line2D([], [], color=t["text2"], linewidth=0, marker="D", markersize=6,
                              label="full source depth (chr20 43x, chr6 45x)"))
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, 0.995), labelcolor=t["text2"])
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=150, facecolor=t["surface"])
    plt.close(fig)


if __name__ == "__main__":
    data = load(sys.argv[1])
    outdir = Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    for theme in THEMES:
        draw(theme, data, outdir / f"ont-hp-prior-{theme}.png")
        print(f"wrote {outdir / f'ont-hp-prior-{theme}.png'}")
