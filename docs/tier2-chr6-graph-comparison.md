# chr6: 4-haplotype vs 34-haplotype graph

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
| `readlik-z` | 0.9602 | 0.9689 | **+0.0088** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0136** on the 4-haplotype graph to **+0.0372** on the 34-haplotype one — 2.7x wider.

**Two directions, and they are not the same direction.** GT F1 rises on the richer graph for the read-likelihood caller; BASEPAIR and SV F1 fall for both callers. The SV fall is **entirely precision** — recall is flat on chr6 and slightly better on chr20 — and most of it is not the caller getting worse. Two thirds to all of it is records that are not structural variants plus the cost of scoring unfiltered; at matched sensitivity the residual is 0.021 on chr6 and zero on chr20. [tier2-sv-errors.md](tier2-sv-errors.md) has the decomposition.

Exposure to multi-allelic sites was the earlier explanation and it does not survive measurement: precision falls within the biallelic stratum, which is 78-82% of records, by nearly the whole amount. Multi-allelic records do grow (17.6% to 22.1% of SV-sized records) and are harder, but they are a minor term rather than the mechanism.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik-z` on the 34-haplotype graph looked like a precision-for-recall trade (measured on chr20: 1,597 false-positive SNVs against 375). The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 3,298 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik-z` goes 0.9605 to 0.9672.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric. It has since been *bounded* rather than settled: false calls made by both callers on both graphs with no truth candidate anywhere nearby number 44 on chr6 and 40 on chr20, which is a lower bound on the benchmark's share of them.

## Cost

| arm | 4-hap wall | 34-hap wall | 4-hap RSS | 34-hap RSS | 4-hap variants | 34-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 337 s | 605 s | 4.8 GB | 4.4 GB | 288,849 | 294,626 |
| `poisson-z` | 169 s | 226 s | 4.4 GB | 4.1 GB | 289,002 | 294,835 |
| `readlik` | 337 s | 406 s | 5.0 GB | 4.7 GB | 286,465 | 284,466 |
| `readlik-nomismap` | 346 s | 384 s | 4.7 GB | 5.4 GB | 287,943 | 290,599 |
| `readlik-z-nolink` | 296 s | 267 s | 4.6 GB | 6.5 GB | 286,474 | 284,529 |
| `readlik-z` | 286 s | 338 s | 5.0 GB | 4.7 GB | 286,474 | 284,529 |

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
| `readlik` | ALL | 0.9583 | 0.9588 | +0.0006 |
| `readlik` | SNV | 0.9815 | 0.9771 | -0.0044 |
| `readlik` | Insertion (<50 bp) | 0.8502 | 0.8885 | +0.0384 |
| `readlik` | Deletion (<50 bp) | 0.8921 | 0.9038 | +0.0116 |
| `readlik-nomismap` | ALL | 0.9578 | 0.9546 | -0.0033 |
| `readlik-nomismap` | SNV | 0.9813 | 0.9726 | -0.0086 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8489 | 0.8854 | +0.0365 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8913 | 0.9015 | +0.0102 |
| `readlik-z-nolink` | ALL | 0.9585 | 0.9615 | +0.0030 |
| `readlik-z-nolink` | SNV | 0.9818 | 0.9794 | -0.0024 |
| `readlik-z-nolink` | Insertion (<50 bp) | 0.8499 | 0.8913 | +0.0414 |
| `readlik-z-nolink` | Deletion (<50 bp) | 0.8922 | 0.9066 | +0.0144 |
| `readlik-z` | ALL | 0.9602 | 0.9689 | +0.0088 |
| `readlik-z` | SNV | 0.9822 | 0.9804 | -0.0017 |
| `readlik-z` | Insertion (<50 bp) | 0.8564 | 0.9168 | +0.0604 |
| `readlik-z` | Deletion (<50 bp) | 0.8978 | 0.9394 | +0.0417 |

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
| `readlik` | ALL | 0.9355 | 0.8877 | -0.0478 |
| `readlik` | SNV | 0.9841 | 0.9790 | -0.0051 |
| `readlik` | Insertion (<50 bp) | 0.8408 | 0.7274 | -0.1133 |
| `readlik` | Deletion (<50 bp) | 0.8613 | 0.8377 | -0.0236 |
| `readlik-nomismap` | ALL | 0.9147 | 0.8512 | -0.0635 |
| `readlik-nomismap` | SNV | 0.9839 | 0.9758 | -0.0081 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7572 | 0.6246 | -0.1326 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8728 | 0.8555 | -0.0173 |
| `readlik-z-nolink` | ALL | 0.9359 | 0.8941 | -0.0418 |
| `readlik-z-nolink` | SNV | 0.9842 | 0.9802 | -0.0040 |
| `readlik-z-nolink` | Insertion (<50 bp) | 0.8412 | 0.7286 | -0.1126 |
| `readlik-z-nolink` | Deletion (<50 bp) | 0.8611 | 0.8413 | -0.0198 |
| `readlik-z` | ALL | 0.9500 | 0.9135 | -0.0365 |
| `readlik-z` | SNV | 0.9845 | 0.9808 | -0.0037 |
| `readlik-z` | Insertion (<50 bp) | 0.8693 | 0.7575 | -0.1118 |
| `readlik-z` | Deletion (<50 bp) | 0.8866 | 0.8815 | -0.0051 |

## Structural variants — truvari (GIAB `stvar`)

The SV metric. Reciprocal-overlap matching against the structural benchmark, `--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored against the *small-variant* truth set and therefore have essentially no truth to match above 50 bp (plan §9.22).

**These errors are broken down per record in [tier2-sv-errors.md](tier2-sv-errors.md)**, and three findings there bear directly on this table. The 34-haplotype false-positive rise is not the same errors plus more — only about two thirds of the 4-haplotype false calls survive the graph change, and the new ones are disproportionately calls with no truth candidate at all. A quarter of all false positives are placement or bookkeeping artefacts of the metric. And harmonising representation with `truvari refine` lifts every arm by roughly 0.05 F1.

| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5320 | 0.5512 | 0.4618 | 0.5490 | 0.4944 | **-0.0546** |
| `poisson-z` | 0.5488 | 0.5417 | 0.5468 | 0.4442 | 0.5478 | 0.4881 | **-0.0597** |
| `readlik` | 0.5417 | 0.5255 | 0.5684 | 0.4766 | 0.5547 | 0.4999 | **-0.0548** |
| `readlik-nomismap` | 0.5391 | 0.5307 | 0.5383 | 0.4319 | 0.5387 | 0.4762 | **-0.0625** |
| `readlik-z-nolink` | 0.5507 | 0.5456 | 0.5728 | 0.4716 | 0.5616 | 0.5059 | **-0.0557** |
| `readlik-z` | 0.5514 | 0.5352 | 0.5881 | 0.5187 | 0.5691 | 0.5268 | **-0.0423** |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7938 | 0.8019 | 0.8825 | 0.8645 | 0.8358 | 0.8320 | **-0.0038** |
| `sm50-poisson-z` | Deletion | 0.8945 | 0.9239 | 0.8348 | 0.7428 | 0.8636 | 0.8235 | **-0.0401** |
| `sm50-poisson-z` | ALL | 0.9219 | 0.9271 | 0.9498 | 0.9091 | 0.9356 | 0.9180 | **-0.0176** |
| `sm50-readlik-z` | Insertion | 0.8895 | 0.9168 | 0.8940 | 0.9090 | 0.8917 | 0.9129 | **+0.0211** |
| `sm50-readlik-z` | Deletion | 0.8990 | 0.9340 | 0.9150 | 0.9172 | 0.9070 | 0.9255 | **+0.0186** |
| `sm50-readlik-z` | ALL | 0.9461 | 0.9572 | 0.9754 | 0.9774 | 0.9605 | 0.9672 | **+0.0067** |

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. It matters for how the calls rank, which is a separate page: see [tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits `AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the called genotype explains. That rescales a quality and does not change a genotype, so **the unfiltered numbers on this page are unaffected by it**.

