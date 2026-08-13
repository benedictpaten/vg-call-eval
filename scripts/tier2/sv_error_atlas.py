#!/usr/bin/env python3
"""Per-record forensics for the truvari SV error sets.

Why this exists. The SV tables report recall, precision and F1 per arm, and two gaps
show up in them -- the read model trails the depth model, and both callers lose
precision on the 34-haplotype graph. Neither gap is interpretable from a summary
statistic, because "FP" and "FN" each cover several unrelated failures:

  * an FP may be a call at a locus with no truth at all, or a *correct* event written
    500 bp away in a tandem repeat, or a correct event whose base record was already
    consumed by another query record under `--pick ac`;
  * an FN may be a truth SV the graph never offered, one the genotyper saw and called
    reference, one emitted as several sub-50 bp records and therefore invisible at
    `--sizemin 50`, or one emitted at full size and simply not matched.

This builds one table per side with enough columns to tell those apart, and every
downstream question in the SV investigation is a query against it.

Two things make it cheap. Truvari already annotates *every* record, including FP and
FN, with its best near-miss (`PctSeqSimilarity`, `PctSizeSimilarity`, `StartDistance`,
`TruScore`, `MatchId`) -- so the "was this nearly right?" question needs no new
alignment. And the GIAB structural truth carries tandem-repeat and low-complexity
annotations (`TRF*`, `LCR`, `REMAP`, `RM_clsfam`), so repeat context is free.

`MatchId` is shared between the base and comp sides of a candidate pair, which makes it
the correct join key between arms and between the FN and FP files. Keying on
CHROM/POS/REF/ALT instead loses records to collisions at repeat loci -- about 12 of 849
on chr6 4-hap when I tried it.

One thing the FN taxonomy would like and cannot have: `vg call` emits **no 0/0
records**, so a missed truth SV looks identical whether the site was never offered or
was genotyped and lost to reference.

`site_probe.sh` was written to supply that by re-running with `-a/--genotype-snarls`,
which does emit reference calls. **It does not work, and the columns derived from it
have been removed.** Adding `-a` changes the snarl decomposition, not just which
genotypes are printed: on chr6 4-hap, `poisson-z` calls 47 of 48 large heterozygous
deletions in its normal run, and its *own* `-a` probe contains a matching allele at
only 26 of them. The loss is concentrated in exactly the large, nested sites the probe
was built to interrogate, so overall agreement looks reassuring -- 286,557 non-reference
records either way, 93.9% position agreement -- while the answer to the question being
asked is wrong by half.

Taken at face value the probe said that only 13% of missed SVs ever had an allele of
the right size offered, which would have made this a graph-content story rather than a
genotyper one. It is kept here as a negative result: `-a` is not a site-existence
oracle, and answering "was this allele enumerated?" needs `-T/--traversals`, which
reports candidate traversals without genotyping and has not been tried yet.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import gzip
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
WORK = REPO / "work"
OUT = WORK / "sv-atlas"

# label -> (work subdirectory, contig, graph richness, chromosome)
DATASETS = {
    "chr6-4hap": ("tier2-chr6", "chr6", "4hap"),
    "chr6-34hap": ("tier2-chr6-hap32", "chr6", "34hap"),
    "chr20-4hap": ("tier2-chr20", "chr20", "4hap"),
    "chr20-34hap": ("tier2-chr20-hap32", "chr20", "34hap"),
}
ARMS = ["poisson", "poisson-z", "readlik-support", "readlik-nomismap", "readlik"]

SIZEMIN = 50
REFDIST = 500          # truvari's --refdist for these runs; see params.json
PCT_OK = 0.7           # truvari's --pctseq / --pctsize

# Truvari annotations copied through to both tables. Present on TP, FP and FN alike:
# on a non-TP they describe the best candidate that was considered and rejected.
MATCH_KEYS = ["PctSeqSimilarity", "PctSizeSimilarity", "PctRecOverlap",
              "SizeDiff", "StartDistance", "EndDistance", "TruScore", "GTMatch"]
# GIAB structural-truth context, present on the base side only.
CONTEXT_KEYS = ["TRF", "TRFperiod", "TRFcopies", "TRFsim", "LCR", "REMAP",
                "RM_clsfam", "SVTYPE", "SVLEN"]


def parse_info(field: str) -> dict:
    out = {}
    for kv in field.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = "1"
    return out


def read_vcf(path: Path):
    """Yield (chrom, pos, ref, alt, qual, filt, info_dict, fmt_keys, sample_fields)."""
    if not path.exists():
        return
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt = f[8].split(":") if len(f) > 9 else []
            smp = f[9].split(":") if len(f) > 9 else []
            yield f[0], int(f[1]), f[3], f[4], f[5], f[6], parse_info(f[7]), fmt, smp


def size_bin(dlen: int) -> str:
    a = abs(dlen)
    if a < 100:
        return "50-99"
    if a < 300:
        return "100-299"
    if a < 1000:
        return "300-999"
    return "1k+"


def svtype_of(ref: str, alt: str, info: dict) -> tuple[str, int]:
    """Prefer the record's own SVTYPE/SVLEN; fall back to the allele lengths.

    The GIAB truth carries both. vg emits neither, and its records are always fully
    resolved sequence, so the fallback is exact rather than approximate.
    """
    dlen = len(alt) - len(ref)
    t = info.get("SVTYPE")
    if t in ("INS", "DEL", "INV", "DUP"):
        return t, int(info.get("SVLEN", dlen) or dlen)
    return ("INS" if dlen > 0 else "DEL" if dlen < 0 else "OTHER"), dlen


def fmt_get(fmt: list[str], smp: list[str], key: str):
    if key in fmt:
        i = fmt.index(key)
        if i < len(smp):
            v = smp[i]
            return None if v == "." else v
    return None


def is_nonref(gt: str | None) -> bool:
    if not gt:
        return False
    return any(a not in ("0", ".") for a in gt.replace("|", "/").split("/"))


class PosIndex:
    """Sorted position index over a VCF, for window queries.

    Built once per (dataset, arm) rather than tabix-ing per FN: about 15,000 lookups
    across the matrix, and a subprocess each would dominate the runtime.
    """

    def __init__(self, path: Path, want_fmt=("GT", "DP", "AD", "GQ", "BL")):
        self.recs = []
        for chrom, pos, ref, alt, _q, _f, info, fmt, smp in read_vcf(path):
            gt = fmt_get(fmt, smp, "GT")
            self.recs.append({
                "pos": pos, "ref": ref, "alt": alt,
                "dlen": len(alt) - len(ref),
                "gt": gt, "nonref": is_nonref(gt),
                "n_at": len(info.get("AT", "").split(",")) if "AT" in info else 0,
                **{k: fmt_get(fmt, smp, k) for k in want_fmt if k != "GT"},
            })
        self.recs.sort(key=lambda r: r["pos"])
        self.keys = [r["pos"] for r in self.recs]

    def window(self, lo: int, hi: int) -> list[dict]:
        i = bisect.bisect_left(self.keys, lo)
        j = bisect.bisect_right(self.keys, hi)
        return self.recs[i:j]

    def __len__(self):
        return len(self.recs)


def classify_fn(fn_pos: int, fn_dlen: int, calls: PosIndex,
                info: dict,
                tp_pos: set[int], fp_pos: set[int]) -> tuple[str, str, int]:
    """Return (call_class, match_class, net_dlen) for a missed truth SV.

    call_class -- what the caller did at the locus:
        no-call       nothing non-reference emitted anywhere in the window. Either the
                      site was never offered or reference won; see the module docstring
                      for why `-a/--genotype-snarls` cannot separate those two.
        small-only    non-reference calls present but none resembles the event: the
                      net length change nearby is under half the truth SV's. Nearby
                      SNVs are not evidence that the caller found a 2 kb deletion,
                      and without this check they would be counted as if they were.
        fragmented    likewise all below --sizemin, but the net length change nearby
                      recovers at least half the event with the same sign. The caller
                      found it and decomposed it, so truvari never sees it at all.
        split-merge   a >=50 bp call is right there, but every one of them is already
                      a true positive matched to some *other* truth record. This is
                      the many-to-one case: one emitted event covering two benchmark
                      records, the second of which can only ever be an FN.
        called-large  a >=50 bp call is right there, is itself an FP, and the pair was
                      still not matched. The genuine representation failure.

    match_class -- why truvari rejected the best candidate it did find:
        placement     sequence and size both pass, but beyond --refdist. The event is
                      correct and written elsewhere in the repeat.
        consumed      sequence, size and distance all pass, yet still not a TP --
                      the base record was matched to a different query record.
        dissimilar    a candidate was considered and failed sequence or size.
        none          no candidate at all within the chunk.
    """
    span = max(abs(fn_dlen), 1)
    lo, hi = fn_pos - REFDIST, fn_pos + span + REFDIST
    near = calls.window(lo, hi)
    nonref = [r for r in near if r["nonref"]]
    large = [r for r in nonref if abs(r["dlen"]) >= SIZEMIN]
    net = sum(r["dlen"] for r in nonref)

    if large:
        # Is any of the large neighbours unmatched? If they are all TPs credited to
        # other truth records, this FN is an accounting artefact, not a miss.
        if any(r["pos"] in fp_pos for r in large):
            call_class = "called-large"
        elif any(r["pos"] in tp_pos for r in large):
            call_class = "split-merge"
        else:
            call_class = "called-large"
    elif nonref:
        recovered = net * (1 if fn_dlen > 0 else -1) >= 0.5 * abs(fn_dlen)
        call_class = "fragmented" if recovered else "small-only"
    else:
        call_class = "no-call"

    seq = info.get("PctSeqSimilarity")
    if seq in (None, "."):
        match_class = "none"
    else:
        ok = (float(seq) >= PCT_OK
              and float(info.get("PctSizeSimilarity") or 0) >= PCT_OK)
        if not ok:
            match_class = "dissimilar"
        else:
            d = abs(float(info.get("StartDistance") or 0))
            match_class = "placement" if d > REFDIST else "consumed"
    return call_class, match_class, net


def classify_fp(info: dict) -> str:
    """Same match taxonomy as the FN side, seen from the query record."""
    seq = info.get("PctSeqSimilarity")
    if seq in (None, "."):
        return "none"
    if not (float(seq) >= PCT_OK and float(info.get("PctSizeSimilarity") or 0) >= PCT_OK):
        return "dissimilar"
    return "placement" if abs(float(info.get("StartDistance") or 0)) > REFDIST else "consumed"


def ad_share(ad: str | None, dp: str | None) -> tuple[float | None, float | None]:
    """sum(AD) and sum(AD)/DP -- the explained share, recomputed here.

    These VCFs predate the GQI/share change in the caller, so their GQ is the raw
    likelihood ratio and `share` has to be derived. That is fine for this purpose and
    is stated on the output page: nothing here depends on the emitted GQ being the
    discounted one.
    """
    if ad is None or dp in (None, "0"):
        return None, None
    try:
        s = sum(int(x) for x in ad.split(",") if x != ".")
        d = int(dp)
    except ValueError:
        return None, None
    return (s, min(1.0, s / d)) if d else (s, None)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=ARMS)
    ap.add_argument("--datasets", nargs="*", default=list(DATASETS))
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    truth_rows, call_rows = [], []

    for ds in args.datasets:
        sub, contig, richness = DATASETS[ds]
        res = WORK / sub / "results"
        for arm in args.arms:
            tdir = res / f"truvari-{arm}"
            if not (tdir / "summary.json").exists():
                print(f"skip {ds} {arm}: no truvari output", file=sys.stderr)
                continue

            norm = res / "truvari-norm" / f"{arm}.norm.vcf.gz"
            calls = PosIndex(norm)

            # ---- query side first: tp-comp + fp ----
            # The FN taxonomy needs to know whether a large call beside a missed truth
            # record is itself unmatched, so the query side has to be indexed before
            # the truth side is classified.
            tp_pos = {p for _c, p, *_ in read_vcf(tdir / "tp-comp.vcf.gz")}
            fp_pos = {p for _c, p, *_ in read_vcf(tdir / "fp.vcf.gz")}

            # ---- truth side: tp-base + fn ----
            for side, path in (("TP", tdir / "tp-base.vcf.gz"), ("FN", tdir / "fn.vcf.gz")):
                for chrom, pos, ref, alt, _q, _f, info, fmt, smp in read_vcf(path):
                    svtype, svlen = svtype_of(ref, alt, info)
                    row = {
                        "dataset": ds, "contig": contig, "richness": richness, "arm": arm,
                        "outcome": side, "chrom": chrom, "pos": pos,
                        "svtype": svtype, "svlen": svlen, "sizebin": size_bin(svlen),
                        "gt_truth": fmt_get(fmt, smp, "GT"),
                        "matchid": info.get("MatchId", ""),
                        "call_class": "", "match_class": "", "net_dlen": "",
                    }
                    for k in MATCH_KEYS + CONTEXT_KEYS:
                        row[k] = info.get(k, "")
                    if side == "FN":
                        cc, mc, net = classify_fn(pos, svlen, calls, info,
                                                  tp_pos, fp_pos)
                        row["call_class"], row["match_class"] = cc, mc
                        row["net_dlen"] = net
                    truth_rows.append(row)

            for side, path in (("TP", tdir / "tp-comp.vcf.gz"), ("FP", tdir / "fp.vcf.gz")):
                for chrom, pos, ref, alt, qual, _f, info, fmt, smp in read_vcf(path):
                    svtype, svlen = svtype_of(ref, alt, info)
                    dp = fmt_get(fmt, smp, "DP")
                    ad = fmt_get(fmt, smp, "AD")
                    adsum, share = ad_share(ad, dp)
                    row = {
                        "dataset": ds, "contig": contig, "richness": richness, "arm": arm,
                        "outcome": side, "chrom": chrom, "pos": pos,
                        "svtype": svtype, "svlen": svlen, "sizebin": size_bin(svlen),
                        "gt": fmt_get(fmt, smp, "GT"),
                        "dp": dp, "ad_sum": adsum,
                        "share": round(share, 4) if share is not None else "",
                        "gq": fmt_get(fmt, smp, "GQ"), "bl": fmt_get(fmt, smp, "BL"),
                        "qual": qual,
                        "n_alleles": len(info.get("AT", "").split(",")) if "AT" in info else "",
                        "matchid": info.get("MatchId", ""),
                        "match_class": classify_fp(info) if side == "FP" else "",
                    }
                    for k in MATCH_KEYS:
                        row[k] = info.get(k, "")
                    call_rows.append(row)

            print(f"{ds:12s} {arm:18s} truth={sum(1 for r in truth_rows if r['arm'] == arm and r['dataset'] == ds):5d}"
                  f" calls={sum(1 for r in call_rows if r['arm'] == arm and r['dataset'] == ds):5d}"
                  f" norm={len(calls):6d}", flush=True)

    def write(path: Path, rows: list[dict]) -> None:
        if not rows:
            return
        cols = list(rows[0].keys())
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {path} ({len(rows)} rows)")

    write(out / "truth.tsv", truth_rows)
    write(out / "calls.tsv", call_rows)

    # A short console summary, so a run that produced nothing useful says so rather
    # than leaving a well-formed empty table to be discovered later.
    fn = collections.Counter((r["dataset"], r["arm"], r["call_class"])
                             for r in truth_rows if r["outcome"] == "FN")
    if fn:
        print("\nFN call_class counts:")
        for (ds, arm, cc), n in sorted(fn.items()):
            print(f"  {ds:12s} {arm:18s} {cc:14s} {n:5d}")


if __name__ == "__main__":
    main()
