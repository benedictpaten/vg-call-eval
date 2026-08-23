# chr6: 4-haplotype vs 34-haplotype graph

> **Stale for the caller as of decide-then-render (2026-08).** Every vg figure below was measured
> before genotypes were settled ahead of record construction. That change moved the whole-genome
> autosomal numbers -- ALL F1 0.9703 -> 0.9729, Indel 0.9195 -> 0.9272, SV >=50 bp 0.5488 -> 0.5596,
> with both precision and recall improving in every class -- so the figures here understate the
> current caller by roughly that much, and any *analysis* built on which calls were wrong may have
> picked a different population. Not re-run: these arms use their own reads, truth sets and graphs, and
> re-measuring them is hours of runs that were not spent. Current numbers:
> [wgs-results.md](wgs-results.md), [pangenie-comparison.md](pangenie-comparison.md).
> What changed and what is still open: `planning/decide-then-render.md`.


Same sample, same reads, same truth, same confident regions, same reference sequence. What changes is the graph — and, unavoidably, the alignments.

| | 4-haplotype | 34-haplotype |
|---|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz` | `…HG002.hap32.gbz` |
| haplotypes | 4 (CHM13, GRCh38, 2 recombinants) | **34** (CHM13, GRCh38, **32 recombinants** — the file is named `hap32` after the recombinant count, not the total) |
| HG002 present? | no | **no** — samples are `CHM13`, `GRCh38`, `recombination` |
| alignments | `…HG002.gaf.gz` | `…HG002.hap32.gaf.gz` (remapped) |

**This is not a single-variable experiment.** Reads mapped to one graph cannot be scored against the other, because the node ID spaces differ — so the 34-haplotype arm necessarily uses its own alignments. Graph and alignment move together. That is what adopting a richer graph actually involves, but it means a difference below cannot be attributed to the graph alone.

The rows to watch are the **`-z` arms**, which enumerate alleles from the GBWT haplotypes. Going from 4 to 34 changes which alleles are *available to call* rather than how they are scored, and the tier-2 finding was that enumeration matters more than the genotyper — most of all for SVs. This is the direct test.

## What this says

**The read-likelihood caller is better on the richer graph; the Poisson caller is much worse on it.** That split is the result. More haplotypes offer more true alleles and more wrong ones, and what decides the outcome is whether the genotyper can tell them apart read by read.

| arm | 4-hap GT F1 | 34-hap GT F1 | Δ |
|---|---|---|---|
| `poisson-z` | 0.9466 | 0.9318 | **-0.0148** |
| `readlik` | 0.9598 | 0.9750 | **+0.0152** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0132** on the 4-haplotype graph to **+0.0432** on the 34-haplotype one — 3.3x wider.

**Two directions, and they are not the same direction.** GT F1 rises on the richer graph for the read-likelihood caller; BASEPAIR and SV F1 fall for both callers. The SV fall is **entirely precision** — recall is flat on chr6 and slightly better on chr20 — and most of it is not the caller getting worse. Two thirds to all of it is records that are not structural variants plus the cost of scoring unfiltered; at matched sensitivity the residual is 0.021 on chr6 and zero on chr20. [tier2-sv-errors.md](tier2-sv-errors.md) has the decomposition.

Exposure to multi-allelic sites was the earlier explanation and it does not survive measurement: precision falls within the biallelic stratum, which is 78-82% of records, by nearly the whole amount. Multi-allelic records do grow (17.6% to 22.1% of SV-sized records) and are harder, but they are a minor term rather than the mechanism.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik` on the 34-haplotype graph looked like a precision-for-recall trade (measured on chr20: 1,597 false-positive SNVs against 375). The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 3,132 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik` goes 0.9597 to 0.9714.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric. It has since been *bounded* rather than settled: false calls made by both callers on both graphs with no truth candidate anywhere nearby number 44 on chr6 and 40 on chr20, which is a lower bound on the benchmark's share of them.

## Cost

| arm | 4-hap wall | 34-hap wall | 4-hap RSS | 34-hap RSS | 4-hap variants | 34-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 354 s | 603 s | 5.3 GB | 4.8 GB | 288,849 | 294,626 |
| `poisson-z` | 169 s | 214 s | 5.6 GB | 5.2 GB | 289,002 | 294,835 |
| `readlik-support` | 325 s | 338 s | 6.7 GB | 9.0 GB | 293,247 | 298,098 |
| `readlik-nomismap` | 260 s | 339 s | 7.2 GB | 8.8 GB | 294,788 | 305,480 |
| `readlik-nolink` | 255 s | 272 s | 5.8 GB | 7.7 GB | 293,250 | 297,938 |
| `readlik` | 281 s | 376 s | 6.6 GB | 5.3 GB | 292,762 | 297,484 |

## Small variants — GT F1

| arm | class | 4-hap | 34-hap | Δ |
|---|---|---|---|---|
| `poisson` | ALL | 0.9461 | 0.9297 | -0.0164 |
| `poisson` | SNV | 0.9788 | 0.9685 | -0.0103 |
| `poisson` | Insertion (<50 bp) | 0.8101 | 0.8118 | +0.0017 |
| `poisson` | Deletion (<50 bp) | 0.8356 | 0.7824 | -0.0532 |
| `poisson-z` | ALL | 0.9466 | 0.9318 | -0.0148 |
| `poisson-z` | SNV | 0.9791 | 0.9703 | -0.0088 |
| `poisson-z` | Insertion (<50 bp) | 0.8111 | 0.8157 | +0.0046 |
| `poisson-z` | Deletion (<50 bp) | 0.8367 | 0.7852 | -0.0514 |
| `readlik-support` | ALL | 0.9582 | 0.9673 | +0.0091 |
| `readlik-support` | SNV | 0.9816 | 0.9859 | +0.0043 |
| `readlik-support` | Insertion (<50 bp) | 0.8491 | 0.8945 | +0.0454 |
| `readlik-support` | Deletion (<50 bp) | 0.8911 | 0.9093 | +0.0183 |
| `readlik-nomismap` | ALL | 0.9595 | 0.9704 | +0.0109 |
| `readlik-nomismap` | SNV | 0.9820 | 0.9823 | +0.0003 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8542 | 0.9184 | +0.0642 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8952 | 0.9396 | +0.0444 |
| `readlik-nolink` | ALL | 0.9581 | 0.9674 | +0.0093 |
| `readlik-nolink` | SNV | 0.9816 | 0.9860 | +0.0044 |
| `readlik-nolink` | Insertion (<50 bp) | 0.8486 | 0.8949 | +0.0463 |
| `readlik-nolink` | Deletion (<50 bp) | 0.8908 | 0.9098 | +0.0190 |
| `readlik` | ALL | 0.9598 | 0.9750 | +0.0152 |
| `readlik` | SNV | 0.9821 | 0.9870 | +0.0049 |
| `readlik` | Insertion (<50 bp) | 0.8555 | 0.9214 | +0.0659 |
| `readlik` | Deletion (<50 bp) | 0.8958 | 0.9422 | +0.0465 |

## Small variants — BASEPAIR F1

| arm | class | 4-hap | 34-hap | Δ |
|---|---|---|---|---|
| `poisson` | ALL | 0.9162 | 0.8219 | -0.0943 |
| `poisson` | SNV | 0.9826 | 0.9731 | -0.0095 |
| `poisson` | Insertion (<50 bp) | 0.8199 | 0.7094 | -0.1105 |
| `poisson` | Deletion (<50 bp) | 0.8068 | 0.6626 | -0.1441 |
| `poisson-z` | ALL | 0.9169 | 0.8076 | -0.1093 |
| `poisson-z` | SNV | 0.9827 | 0.9747 | -0.0080 |
| `poisson-z` | Insertion (<50 bp) | 0.8210 | 0.6845 | -0.1365 |
| `poisson-z` | Deletion (<50 bp) | 0.8080 | 0.6278 | -0.1803 |
| `readlik-support` | ALL | 0.9250 | 0.8769 | -0.0481 |
| `readlik-support` | SNV | 0.9841 | 0.9868 | +0.0027 |
| `readlik-support` | Insertion (<50 bp) | 0.8032 | 0.6821 | -0.1211 |
| `readlik-support` | Deletion (<50 bp) | 0.8594 | 0.8390 | -0.0204 |
| `readlik-nomismap` | ALL | 0.9106 | 0.8550 | -0.0556 |
| `readlik-nomismap` | SNV | 0.9843 | 0.9841 | -0.0003 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7382 | 0.6067 | -0.1315 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8796 | 0.8906 | +0.0110 |
| `readlik-nolink` | ALL | 0.9248 | 0.8801 | -0.0447 |
| `readlik-nolink` | SNV | 0.9840 | 0.9868 | +0.0027 |
| `readlik-nolink` | Insertion (<50 bp) | 0.8030 | 0.6815 | -0.1216 |
| `readlik-nolink` | Deletion (<50 bp) | 0.8590 | 0.8410 | -0.0180 |
| `readlik` | ALL | 0.9388 | 0.9033 | -0.0355 |
| `readlik` | SNV | 0.9844 | 0.9876 | +0.0031 |
| `readlik` | Insertion (<50 bp) | 0.8477 | 0.7219 | -0.1258 |
| `readlik` | Deletion (<50 bp) | 0.8639 | 0.8799 | +0.0160 |

## Structural variants — truvari (GIAB `stvar`)

The SV metric. Reciprocal-overlap matching against the structural benchmark, `--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored against the *small-variant* truth set and therefore have essentially no truth to match above 50 bp (plan §9.22).

**These errors are broken down per record in [tier2-sv-errors.md](tier2-sv-errors.md)**, and three findings there bear directly on this table. The 34-haplotype false-positive rise is not the same errors plus more — only about two thirds of the 4-haplotype false calls survive the graph change, and the new ones are disproportionately calls with no truth candidate at all. A quarter of all false positives are placement or bookkeeping artefacts of the metric. And harmonising representation with `truvari refine` lifts every arm by roughly 0.05 F1.

| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5320 | 0.5512 | 0.4618 | 0.5490 | 0.4944 | **-0.0546** |
| `poisson-z` | 0.5488 | 0.5417 | 0.5468 | 0.4442 | 0.5478 | 0.4881 | **-0.0597** |
| `readlik-support` | 0.5714 | 0.5740 | 0.5677 | 0.5057 | 0.5695 | 0.5377 | **-0.0318** |
| `readlik-nomismap` | 0.5811 | 0.5992 | 0.5369 | 0.4657 | 0.5581 | 0.5241 | **-0.0341** |
| `readlik-nolink` | 0.5779 | 0.5941 | 0.5647 | 0.4981 | 0.5712 | 0.5419 | **-0.0294** |
| `readlik` | 0.5766 | 0.5856 | 0.5765 | 0.5499 | 0.5766 | 0.5672 | **-0.0094** |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7938 | 0.8019 | 0.8825 | 0.8645 | 0.8358 | 0.8320 | **-0.0038** |
| `sm50-poisson-z` | Deletion | 0.8945 | 0.9239 | 0.8348 | 0.7428 | 0.8636 | 0.8235 | **-0.0401** |
| `sm50-poisson-z` | ALL | 0.9219 | 0.9271 | 0.9498 | 0.9091 | 0.9356 | 0.9180 | **-0.0176** |
| `sm50-readlik` | Insertion | 0.8898 | 0.9267 | 0.8907 | 0.9024 | 0.8902 | 0.9144 | **+0.0242** |
| `sm50-readlik` | Deletion | 0.8995 | 0.9480 | 0.9097 | 0.9078 | 0.9046 | 0.9275 | **+0.0229** |
| `sm50-readlik` | ALL | 0.9464 | 0.9702 | 0.9733 | 0.9725 | 0.9597 | 0.9714 | **+0.0117** |

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. It matters for how the calls rank, which is a separate page: see [tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits `AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the called genotype explains. That rescales a quality and does not change a genotype, so **the unfiltered numbers on this page are unaffected by it**.

