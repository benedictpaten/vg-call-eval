#!/usr/bin/env python3
"""Should the explained-share discount also apply to GQ on records the linkage HMM changed?

The caller computes `GQ = GQI * share`, where `share` is the fraction of reads whose
best-fitting allele is one the call contains. Reads outside the call enter every
genotype's likelihood and cancel from the best-versus-second-best gap, so GQ cannot see
them; the discount is what stops a call being Q60 while a tenth of its pile-up argues for
something it does not contain.

When the linkage pass changes a genotype it **replaces** GQ with the phred complement of
the HMM posterior and applies no discount. So a single VCF carries two definitions, and
`GQ <= GQI` -- an invariant everywhere else -- fails on 4.63% of chr20-34hap records.

The argument for discounting those too is that the blindness is identical: the posterior
is built from the same emission, from the same matrix, so reads outside the call are
exactly as invisible to it. The argument against is that a posterior is already a
probability and multiplying it by a heuristic makes it neither.

That is an empirical question, and this answers it without rebuilding vg, because the
proposed change is a pure post-hoc transform of an emitted field:

    changed records:  GQ' = GQ * share        (share from AD/DP, as the caller computes it)
    other records:    GQ' = GQ                (already discounted)

Records the HMM changed are identified by diffing `readlik-z` against `readlik-z-nolink`,
which differ only in `--linkage-weight`.

Scored the way every other quality signal in this project has been: AUC over true against
false calls, and false calls surviving at matched recall. GQ does not change any genotype,
so accuracy metrics are blind to this -- ranking is the only thing that can move.
"""

from __future__ import annotations

import argparse
import gzip
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from filter_lib import DATASETS, labels, truth_counts  # noqa: E402


def auc(pos: list[float], neg: list[float]) -> float:
    """P(a true call outranks a false one), ties counted as half."""
    if not pos or not neg:
        return float("nan")
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum, i = 0.0, 0
    while i < len(merged):
        j = i
        while j < len(merged) and merged[j][0] == merged[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2.0
        rank_sum += sum(avg_rank for k in range(i, j) if merged[k][1] == 1)
        i = j
    n_pos, n_neg = len(pos), len(neg)
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def genotypes(path: Path) -> dict:
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out[int(f[1])] = f[9].split(":")[0]
    return out


def query(vcf: Path) -> dict:
    """pos -> (dp, share, gq, gqi), share from AD/DP as the caller computes it."""
    fmt = "%POS[\t%DP\t%AD\t%GQ\t%GQI]\n"
    q = subprocess.run(["bcftools", "query", "-f", fmt, str(vcf)],
                       capture_output=True, text=True)
    if q.returncode != 0:
        raise SystemExit(q.stderr.strip()[:400])
    out = {}
    for line in q.stdout.splitlines():
        f = line.split("\t")
        try:
            dp = float(f[1])
            ad = [int(x) for x in f[2].split(",")]
            gq, gqi = float(f[3]), float(f[4])
        except (ValueError, IndexError):
            continue
        if dp <= 0:
            continue
        out[int(f[0])] = (dp, min(1.0, sum(ad) / dp), gq, gqi)
    return out


def at_recall(scores: list[float], ys: list[int], base_tp: int, total: int,
              target: float) -> tuple[float, int]:
    """Precision and surviving false calls at the first threshold reaching base-side `target`.

    Recall is scaled onto the base side: the query-side true-positive fraction times the
    unfiltered base-side count over the truth total. Copied from share_gq.py rather than
    rewritten, because the base-side conversion is the part that is easy to get wrong and
    every other quality signal here was scored with exactly this.
    """
    order = sorted(range(len(ys)), key=lambda i: -scores[i])
    tp_all = sum(ys)
    tp = 0
    for rank, i in enumerate(order, 1):
        tp += ys[i]
        rec = base_tp * (tp / tp_all) / total if tp_all and total else 0.0
        if rec >= target:
            return tp / rank, rank - tp
    return float("nan"), -1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", default="smvar", choices=["smvar", "stvar"])
    ap.add_argument("--recall", type=float, default=0.95)
    args = ap.parse_args()

    print(f"benchmark: {args.kind}, matched recall {args.recall}\n")
    hdr = (f"{'dataset':14s} {'chg':>6s} {'AUC now':>9s} {'AUC xshare':>11s} {'AUC cap':>9s} "
           f"{'FP now':>7s} {'FP xshare':>10s} {'FP cap':>7s}")
    print(hdr)
    print("-" * len(hdr))

    for name, work in DATASETS:
        vcf = work / "results/readlik-z.vcf.gz"
        nolink = work / "results/readlik-z-nolink.vcf.gz"
        if not vcf.exists() or not nolink.exists():
            print(f"{name:14s} (missing readlik-z or readlik-z-nolink)")
            continue
        try:
            bd = labels(work, args.kind)
            base_tp, total = truth_counts(work, args.kind)
        except Exception as exc:  # noqa: BLE001
            print(f"{name:14s} ({exc})")
            continue

        with_link, without = genotypes(vcf), genotypes(nolink)
        changed = {p for p in with_link if p in without and with_link[p] != without[p]}
        rows = query(vcf)

        now_s, new_s, cap_s, ys = [], [], [], []
        for pos, (dp, share, gq, gqi) in rows.items():
            if pos not in bd:
                continue
            # labels() yields the strings "TP"/"FP"; both are truthy, and treating
            # them as booleans made every record a positive and every AUC nan.
            ys.append(1 if bd[pos] == "TP" else 0)
            now_s.append(gq)
            # Only records the HMM touched are transformed; the rest are already discounted.
            new_s.append(gq * share if pos in changed else gq)
            # Third variant: also cap at GQI, so the linkage pass may lower confidence but
            # never raise it above what the per-site evidence alone supported. This is the
            # only one of the three that restores `GQ <= GQI` -- multiplying by the share
            # does not, because the posterior-based quality is not derived from GQI and
            # bears no arithmetic relation to it.
            cap_s.append(min(gq * share, gqi) if pos in changed else gq)

        n_changed = sum(1 for p in rows if p in bd and p in changed)
        def score(vals):
            return auc([v for v, y in zip(vals, ys) if y],
                       [v for v, y in zip(vals, ys) if not y])
        a_now, a_new, a_cap = score(now_s), score(new_s), score(cap_s)
        _, fp_now = at_recall(now_s, ys, base_tp, total, args.recall)
        _, fp_new = at_recall(new_s, ys, base_tp, total, args.recall)
        _, fp_cap = at_recall(cap_s, ys, base_tp, total, args.recall)
        print(f"{name:14s} {n_changed:6d} {a_now:9.5f} {a_new:11.5f} {a_cap:9.5f} "
              f"{fp_now:7d} {fp_new:10d} {fp_cap:7d}")

    print("\nGQ changes no genotype, so accuracy metrics cannot see this. A higher AUC or "
          "fewer surviving false calls is the whole of the difference.")


if __name__ == "__main__":
    main()
