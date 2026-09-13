#!/usr/bin/env python3
"""Tabulate arm.py result JSONs side by side, optionally as deltas against a baseline.

    python3 scripts/tier2/arm_table.py [--base TAG] [--label TAG=text ...] TAG [TAG ...]

Reads work/ont-preset/<tag>.json, which arm.py writes. Prints SNV / indel / ALL F1 with
the counts that F1 hides -- a genotype-quality regression scores identically, so the
counters are the gate, not the F1.

ALL TAGS MUST BE THE SAME CONTIG AND READ SET. Nothing in a json records which arm produced
it, so this cannot check for you, and a cross-contig delta looks exactly like a real one.
That mistake has been made in this project twice: an ONT whole-genome TP count compared
against a short-read autosome count, and an indel precision quoted chrX-inclusive under an
autosome table.
"""
import sys, json
from pathlib import Path

OUT = Path(__file__).resolve().parents[2] / "work/ont-preset"


def load(tag):
    p = OUT / f"{tag}.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def main():
    args = sys.argv[1:]
    base, labels = None, {}
    rest = []
    i = 0
    while i < len(args):
        if args[i] == "--base":
            base = args[i + 1]; i += 2
        elif args[i] == "--label":
            k, _, v = args[i + 1].partition("="); labels[k] = v; i += 2
        else:
            rest.append(args[i]); i += 1
    if not rest:
        sys.exit(__doc__)

    b = load(base) if base else None
    w = max([len(labels.get(t, t)) for t in rest] + [14])
    print(f"{'arm':<{w}}{'SNV F1':>10}{'SNVfp':>7}{'Indel F1':>10}{'Indfp':>7}{'Indfn':>7}"
          f"{'ALL F1':>10}{'ALLfp':>7}{'HP>=5fp':>9}")
    for t in rest:
        j = load(t)
        if j is None:
            print(f"{labels.get(t, t):<{w}}  -- no json --"); continue
        h5 = j["hp"]["rec"].get(">=5", [0, 0])[0]
        print(f"{labels.get(t, t):<{w}}{j['snv_f1']:>10.5f}{j['snv_fp']:>7d}"
              f"{j['indel_f1']:>10.5f}{j['indel_fp']:>7d}{j['indel_fn']:>7d}"
              f"{j['all_f1']:>10.5f}{j['all_fp']:>7d}{h5:>9d}")
    if b:
        print(f"\ndeltas against {base}:")
        for t in rest:
            j = load(t)
            if j is None or t == base:
                continue
            print(f"  {labels.get(t, t):<{w}} SNV {j['snv_f1']-b['snv_f1']:+.5f}   "
                  f"Indel {j['indel_f1']-b['indel_f1']:+.5f}   "
                  f"ALL {j['all_f1']-b['all_f1']:+.5f}   "
                  f"ALLfp {j['all_fp']-b['all_fp']:+d}")


if __name__ == "__main__":
    main()
