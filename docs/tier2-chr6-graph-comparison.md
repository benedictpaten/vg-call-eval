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
| `readlik` | 0.9575 | 0.9775 | **+0.0200** |

**And the cost side runs the same way.** Going from four haplotypes to thirty-four costs the read-likelihood arms 1.1x to 1.3x more CPU and `poisson` **2.75x** more, so the caller that gets better on the richer graph is also the one whose compute barely grows. The Cost section below has the per-arm figures and the caveats.

The read-likelihood caller's margin over the Poisson caller goes from **+0.0109** on the 4-haplotype graph to **+0.0457** on the 34-haplotype one — 4.2x wider.

**Two directions, and they are not the same direction.** GT F1 rises on the richer graph for the read-likelihood caller; BASEPAIR falls for both callers, and SV F1 falls for every arm **except** `readlik`, which holds flat at +0.0011. The SV fall is **entirely precision** — recall is flat on chr6 and slightly better on chr20 — and most of it is not the caller getting worse. Two thirds to all of it is records that are not structural variants plus the cost of scoring unfiltered; at matched sensitivity the residual is 0.021 on chr6 and zero on chr20. [tier2-sv-errors.md](tier2-sv-errors.md) has the decomposition.

Exposure to multi-allelic sites was the earlier explanation and it does not survive measurement: precision falls within the biallelic stratum, which is 78-82% of records, by nearly the whole amount. Multi-allelic records do grow (17.6% to 22.1% of SV-sized records) and are harder, but they are a minor term rather than the mechanism.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik` on the 34-haplotype graph looked like a precision-for-recall trade (measured on chr20: 1,597 false-positive SNVs against 375). The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 3,085 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik` goes 0.9593 to 0.9772.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric. It has since been *bounded* rather than settled: false calls made by both callers on both graphs with no truth candidate anywhere nearby number 44 on chr6 and 40 on chr20, which is a lower bound on the benchmark's share of them.

## Cost

| arm | 4-hap wall | 34-hap wall | 4-hap CPU | 34-hap CPU | **CPU x** | 4-hap RSS | 34-hap RSS | 4-hap variants | 34-hap variants |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 329 s | 666 s | 981 s | 2,696 s | **2.75x** | 6.1 GB | 6.0 GB | 288,849 | 294,626 |
| `poisson-z` | 161 s | 198 s | 366 s | 524 s | **1.43x** | 5.9 GB | 6.1 GB | 289,002 | 294,835 |
| `readlik-support` | 308 s | 350 s | 1,004 s | 1,232 s | **1.23x** | 8.0 GB | 8.2 GB | 293,606 | 299,877 |
| `readlik-nomismap` | 267 s | 356 s | 912 s | 1,144 s | **1.25x** | 7.3 GB | 8.5 GB | 296,674 | 303,729 |
| `readlik-nolink` | 255 s | 277 s | 875 s | 981 s | **1.12x** | 6.9 GB | 8.0 GB | 293,633 | 299,880 |
| `readlik` | 282 s | 366 s | 892 s | 1,182 s | **1.33x** | 6.7 GB | 7.6 GB | 295,204 | 296,793 |

**`CPU x` is the column to read, and it says something the accuracy tables do not.** CPU is user+sys, so unlike wall clock it measures work rather than elapsed time -- it does not move with how much of the machine a phase manages to use, or with how warm the page cache was. Going from four haplotypes to thirty-four, the read-likelihood arms cost between 1.1x and 1.3x more compute. `poisson` costs **2.75x** more.

So the split this page opens with has a cost side as well as an accuracy side: the caller that gets *better* on the richer graph is also the one whose compute barely grows, and the caller that gets worse is the one that more than doubles. `poisson-z` sits between them at 1.43x, which locates most of the effect in support enumeration rather than in Poisson genotyping.

Read it with the not-a-single-variable caveat above: the two graphs differ in topology and the reads are remapped, so this is not panel size alone. The arm-to-arm contrast across one fixed pair of graphs is what the column supports.

One caveat on the Poisson rows specifically: that path is **not bit-reproducible**. The same binary run twice on the 4-haplotype dataset differs on 20 records of 289,002, in depth-derived fields -- `QUAL`, `GL`, `XD` -- with `GT`, `AD` and `GQ` identical, so no genotype moves and the F1 figures are stable to the digits shown. It does mean two regenerations of this page will not diff clean on those arms, and that byte-identity is not a usable regression gate for them. The read-likelihood arms are exactly reproducible; that was verified across twenty runs in [performance.md](performance.md).

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
| `readlik-support` | ALL | 0.9582 | 0.9672 | +0.0090 |
| `readlik-support` | SNV | 0.9815 | 0.9858 | +0.0042 |
| `readlik-support` | Insertion (<50 bp) | 0.8489 | 0.8941 | +0.0452 |
| `readlik-support` | Deletion (<50 bp) | 0.8910 | 0.9088 | +0.0178 |
| `readlik-nomismap` | ALL | 0.9569 | 0.9726 | +0.0157 |
| `readlik-nomismap` | SNV | 0.9806 | 0.9830 | +0.0025 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8471 | 0.9256 | +0.0785 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8896 | 0.9470 | +0.0574 |
| `readlik-nolink` | ALL | 0.9580 | 0.9674 | +0.0094 |
| `readlik-nolink` | SNV | 0.9815 | 0.9859 | +0.0044 |
| `readlik-nolink` | Insertion (<50 bp) | 0.8485 | 0.8943 | +0.0459 |
| `readlik-nolink` | Deletion (<50 bp) | 0.8908 | 0.9099 | +0.0191 |
| `readlik` | ALL | 0.9575 | 0.9775 | +0.0200 |
| `readlik` | SNV | 0.9809 | 0.9880 | +0.0071 |
| `readlik` | Insertion (<50 bp) | 0.8482 | 0.9287 | +0.0804 |
| `readlik` | Deletion (<50 bp) | 0.8904 | 0.9498 | +0.0595 |

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
| `readlik-support` | ALL | 0.9251 | 0.8776 | -0.0475 |
| `readlik-support` | SNV | 0.9840 | 0.9865 | +0.0025 |
| `readlik-support` | Insertion (<50 bp) | 0.8029 | 0.6817 | -0.1212 |
| `readlik-support` | Deletion (<50 bp) | 0.8593 | 0.8371 | -0.0222 |
| `readlik-nomismap` | ALL | 0.9070 | 0.8620 | -0.0450 |
| `readlik-nomismap` | SNV | 0.9833 | 0.9842 | +0.0009 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7393 | 0.6144 | -0.1249 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8649 | 0.9033 | +0.0384 |
| `readlik-nolink` | ALL | 0.9248 | 0.8804 | -0.0444 |
| `readlik-nolink` | SNV | 0.9840 | 0.9865 | +0.0025 |
| `readlik-nolink` | Insertion (<50 bp) | 0.8027 | 0.6811 | -0.1217 |
| `readlik-nolink` | Deletion (<50 bp) | 0.8588 | 0.8398 | -0.0190 |
| `readlik` | ALL | 0.9442 | 0.9090 | -0.0352 |
| `readlik` | SNV | 0.9836 | 0.9879 | +0.0043 |
| `readlik` | Insertion (<50 bp) | 0.8476 | 0.7300 | -0.1176 |
| `readlik` | Deletion (<50 bp) | 0.8835 | 0.8854 | +0.0019 |

## Structural variants — truvari (GIAB `stvar`)

The SV metric. Reciprocal-overlap matching against the structural benchmark, `--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored against the *small-variant* truth set and therefore have essentially no truth to match above 50 bp (plan §9.22).

**These errors are broken down per record in [tier2-sv-errors.md](tier2-sv-errors.md)**, and three findings there bear directly on this table. The 34-haplotype false-positive rise is not the same errors plus more — only about two thirds of the 4-haplotype false calls survive the graph change, and the new ones are disproportionately calls with no truth candidate at all. A quarter of all false positives are placement or bookkeeping artefacts of the metric. And harmonising representation with `truvari refine` lifts every arm by roughly 0.05 F1.

| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5320 | 0.5512 | 0.4618 | 0.5490 | 0.4944 | **-0.0546** |
| `poisson-z` | 0.5488 | 0.5417 | 0.5468 | 0.4442 | 0.5478 | 0.4881 | **-0.0597** |
| `readlik-support` | 0.5721 | 0.5779 | 0.5712 | 0.5092 | 0.5717 | 0.5414 | **-0.0303** |
| `readlik-nomismap` | 0.5947 | 0.6167 | 0.5421 | 0.4831 | 0.5672 | 0.5418 | **-0.0254** |
| `readlik-nolink` | 0.5779 | 0.5973 | 0.5647 | 0.4928 | 0.5712 | 0.5400 | **-0.0312** |
| `readlik` | 0.5966 | 0.6063 | 0.5691 | 0.5627 | 0.5825 | 0.5837 | **+0.0011** |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7938 | 0.8019 | 0.8825 | 0.8645 | 0.8358 | 0.8320 | **-0.0038** |
| `sm50-poisson-z` | Deletion | 0.8945 | 0.9239 | 0.8348 | 0.7428 | 0.8636 | 0.8235 | **-0.0401** |
| `sm50-poisson-z` | ALL | 0.9219 | 0.9271 | 0.9498 | 0.9091 | 0.9356 | 0.9180 | **-0.0176** |
| `sm50-readlik` | Insertion | 0.8924 | 0.9387 | 0.8845 | 0.9104 | 0.8884 | 0.9244 | **+0.0360** |
| `sm50-readlik` | Deletion | 0.9035 | 0.9517 | 0.9069 | 0.9184 | 0.9052 | 0.9347 | **+0.0295** |
| `sm50-readlik` | ALL | 0.9484 | 0.9747 | 0.9704 | 0.9798 | 0.9593 | 0.9772 | **+0.0180** |

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. It matters for how the calls rank, which is a separate page: see [tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits `AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the called genotype explains. That rescales a quality and does not change a genotype, so **the unfiltered numbers on this page are unaffected by it**.

