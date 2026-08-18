#!/usr/bin/env python3
"""Stage 4 of the mosaic output: how good is the phasing, measured against a phased truth?

The phase here does not come from reads spanning consecutive sites, the way a read-based phaser
works. It comes from the haplotype panel: the linkage layer's Viterbi path says which pair of panel
haplotypes best explains the calls, and the order of that pair is the phase. So a phase block is a
whole chain rather than a read-length island, and switch error is the metric that has to carry the
weight -- a phaser can trivially win on switch error by emitting short blocks, and this does the
opposite. Block length is therefore reported beside it, always.

**The switch error reported here already excludes genotype errors, and that is not a choice made
here.** This script was written to report two numbers -- all het sites, and only correctly
genotyped ones -- on the reasoning that the phasing is constrained to our own calls, so a
mis-genotyped site forces the path through it and may cost a switch that is really a calling
error. Measured, the two are identical to the digit, because `whatshap compare` assesses only
variants that are heterozygous *and* identically genotyped in both files. On chr20 that drops 815
het-in-both sites whose genotypes disagree, and removing them by hand changes neither the assessed
pair count nor the switch count.

So there is one number, and it flatters the phasing relative to what a user experiences: sites we
call wrongly do not count against it, they simply disappear from the denominator. The
correct-genotype filter is kept because it makes that explicit rather than implicit, and because a
future scorer that does not intersect the same way would need it.

Hamming distance is reported but should not be read as a quality: over a single chromosome-length
block every switch flips the relative phase of everything downstream, so any non-zero switch rate
drives it toward 50%. It measures block length here, not phasing.

Truth is the T2T-Q100 HG002 benchmark, which is assembly-derived and fully phased. HG002 is
excluded from the HPRC graph, so this measures *imputation* phasing against a panel that does not
contain the sample -- the fair test, and the one whose result will not generalise to a sample the
panel represents poorly.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WHATSHAP = "/tmp/wsenv/bin/whatshap"


def sh(cmd, **kw):
    proc = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if proc.returncode != 0:
        print(" ".join(str(c) for c in cmd), file=sys.stderr)
        print(proc.stderr[-3000:], file=sys.stderr)
        raise SystemExit(f"command failed: {cmd[0]}")
    return proc.stdout


def prepare_calls(calls: Path, out: Path, sample: str) -> Path:
    """Rename the sample to match the truth and index, since whatshap pairs by sample name."""
    names = out.with_suffix(".sample.txt")
    names.write_text(f"{sample}\n")
    sh(["bcftools", "reheader", "-s", str(names), "-o", str(out), str(calls)])
    sh(["bcftools", "index", "-f", "-t", str(out)])
    return out


def gt_key(gt: str):
    """Unordered allele set, or None where the genotype is missing or half-called."""
    if not gt or "." in gt:
        return None
    toks = gt.replace("|", "/").split("/")
    if len(toks) != 2:
        return None
    return tuple(sorted(toks))


def correct_only(truth: Path, calls: Path, out: Path) -> int:
    """Write the calls restricted to sites whose genotype matches the truth's, unordered.

    Keeps phase error separable from genotype error. Matching is on the *unordered* pair on
    purpose: a call that is right but phased the other way round is a phase observation, which is
    what we came to measure, not a genotype error to be filtered out.
    """
    tsv = sh(["bcftools", "query", "-f", "%CHROM\\t%POS\\t[%GT]\\n", str(truth)])
    want = {}
    for line in tsv.splitlines():
        c, pos, gt = line.split("\t")
        k = gt_key(gt)
        if k:
            want[(c, pos)] = k

    keep = []
    tsv = sh(["bcftools", "query", "-f", "%CHROM\\t%POS\\t[%GT]\\n", str(calls)])
    for line in tsv.splitlines():
        c, pos, gt = line.split("\t")
        k = gt_key(gt)
        if k and want.get((c, pos)) == k:
            keep.append(f"{c}\t{pos}")
    regions = out.with_suffix(".sites.txt")
    regions.write_text("\n".join(keep) + "\n")
    sh(["bcftools", "view", "-T", str(regions), "-o", str(out), "-O", "z", str(calls)])
    sh(["bcftools", "index", "-f", "-t", str(out)])
    return len(keep)


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _num(row, key):
    try:
        return float(row.get(key, "") or 0)
    except ValueError:
        return 0.0


def compare(truth: Path, calls: Path, out_prefix: Path, sample: str) -> dict:
    """Switch statistics, summed across chromosomes.

    whatshap writes **one row per chromosome**, and this used to return `rows[0]`. On a
    single-contig run that was right; on a genome it silently reported chr1 and called it the
    genome -- 191,676 pairs where the autosomes carry 2,442,552, and a rate that happened to look
    plausible. Counts are summed and the rate recomputed, never averaged over chromosomes, which
    would weight chr21 like chr1.
    """
    tsv = out_prefix.with_suffix(".compare.tsv")
    sh([WHATSHAP, "compare", "--sample", sample, "--names", "truth,calls",
        "--tsv-pairwise", str(tsv), str(truth), str(calls)])
    rows = read_tsv(tsv)
    if not rows:
        return {}
    pairs = sum(_num(r, "all_assessed_pairs") for r in rows)
    switches = sum(_num(r, "all_switches") for r in rows)
    out = dict(rows[0])
    out["all_assessed_pairs"] = f"{pairs:.0f}"
    out["all_switches"] = f"{switches:.0f}"
    out["all_switchflip_rate"] = f"{switches / pairs:.6f}" if pairs else "nan"
    out["blockwise_hamming"] = f"{sum(_num(r, 'blockwise_hamming') for r in rows):.0f}"
    out["intersection_blocks"] = f"{sum(_num(r, 'intersection_blocks') for r in rows):.0f}"
    out["chromosomes"] = str(len(rows))
    return out


def stats(calls: Path, out_prefix: Path, sample: str) -> dict:
    """Block statistics, summed across chromosomes -- same per-chromosome trap as compare()."""
    tsv = out_prefix.with_suffix(".stats.tsv")
    sh([WHATSHAP, "stats", "--sample", sample, "--tsv", str(tsv), str(calls)])
    rows = read_tsv(tsv)
    if not rows:
        return {}
    # whatshap emits an ALL row alongside the per-chromosome ones in some versions; drop it so it
    # is not counted twice.
    per = [r for r in rows if (r.get("chromosome") or "").lower() not in ("all", "")]
    out = dict(per[0]) if per else dict(rows[0])
    out["phased"] = f"{sum(_num(r, 'phased') for r in per):.0f}"
    out["blocks"] = f"{sum(_num(r, 'blocks') for r in per):.0f}"
    # N50 does not sum. The longest block is the meaningful summary here, since the claim being
    # made is chain-length blocks.
    out["block_n50"] = f"{max((_num(r, 'block_n50') for r in per), default=0):.0f}"
    out["bp_per_block_max"] = f"{max((_num(r, 'bp_per_block_max') for r in per), default=0):.0f}"
    out["chromosomes"] = str(len(per))
    return out


def fmt(row: dict, key: str, default="-"):
    v = row.get(key)
    return default if v in (None, "", "nan") else v


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls", required=True, help="phased VCF from vg call --phased")
    p.add_argument("--truth", required=True, help="phased truth VCF")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--sample", default="HG002")
    args = p.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    truth = Path(args.truth)

    calls_gz = Path(str(out) + ".calls.vcf.gz")
    src = Path(args.calls)
    if src.suffix != ".gz":
        packed = Path(str(out) + ".src.vcf.gz")
        sh(["bcftools", "view", "-o", str(packed), "-O", "z", str(src)])
        sh(["bcftools", "index", "-f", "-t", str(packed)])
        src = packed
    prepare_calls(src, calls_gz, args.sample)

    print("== all het sites ==", flush=True)
    all_cmp = compare(truth, calls_gz, Path(str(out) + ".all"), args.sample)
    all_stats = stats(calls_gz, Path(str(out) + ".all"), args.sample)

    print("== correctly genotyped only ==", flush=True)
    good = Path(str(out) + ".correct.vcf.gz")
    n_kept = correct_only(truth, calls_gz, good)
    good_cmp = compare(truth, good, Path(str(out) + ".correct"), args.sample)

    rows = [("whatshap intersection", all_cmp), ("+ correct-GT filter", good_cmp)]
    print()
    # Both rates as actual percentages.
    #
    # `all_switchflip_rate` is recomputed above as switches/pairs, which is a *fraction*, and it used
    # to be printed under a "switch %" heading -- so chr20's 0.027667 was read and published as
    # 0.0277% when it is 2.77%. Every phasing figure quoted for this caller was a hundred times
    # better than the measurement.
    #
    # The hamming column is here for the same reason. Switch rate and block length together read as
    # near-perfect chromosome-scale phase; hamming says what the orientation is actually worth, and
    # at 2.77% per adjacent pair over a 248 Mb block it is 49.3% -- a coin flip. A long block is a
    # statement about how the sites are grouped into one PS, not about long-range phase accuracy.
    print(f"{'subset':26s} {'pairs':>10s} {'switches':>9s} {'switch %':>9s} "
          f"{'hamming':>9s} {'hamming %':>10s} {'blocks':>7s}")
    for name, r in rows:
        pairs = _num(r, "all_assessed_pairs")
        rate = _num(r, "all_switchflip_rate")
        ham = _num(r, "blockwise_hamming")
        print(f"{name:26s} {fmt(r, 'all_assessed_pairs'):>10s} "
              f"{fmt(r, 'all_switches'):>9s} {100.0 * rate:>8.4f}% "
              f"{fmt(r, 'blockwise_hamming'):>9s} "
              f"{(100.0 * ham / pairs if pairs else float('nan')):>9.2f}% "
              f"{fmt(r, 'intersection_blocks'):>7s}")
    print()
    print(f"phased variants: {fmt(all_stats, 'phased')}   "
          f"blocks: {fmt(all_stats, 'blocks')}   "
          f"block N50: {fmt(all_stats, 'block_n50')}   "
          f"longest: {fmt(all_stats, 'bp_per_block_max')}")
    print(f"sites kept for the correct-genotype subset: {n_kept}")
    print()
    print("Switch error is only meaningful beside block length: short blocks make it small for")
    print("free, and these blocks are chain-length by construction.")
    print()
    print("The two rows are expected to agree: whatshap already assesses only variants that are")
    print("het and identically genotyped in both files, so the switch error is genotype-error")
    print("free before the second filter is applied. A difference between them would mean that")
    print("assumption had stopped holding.")
    print("Hamming is not a quality measure over one chromosome-length block -- each switch flips")
    print("everything downstream, so it tends to 50% at any non-zero switch rate.")


if __name__ == "__main__":
    main()
