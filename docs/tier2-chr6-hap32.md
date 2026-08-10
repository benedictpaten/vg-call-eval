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
| `readlik-z` | 0.9585 | 0.9615 | **+0.0030** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0119** on the 4-haplotype graph to **+0.0297** on the 34-haplotype one — 2.5x wider.

**Two directions, and they are not the same direction.** GT F1 rises on the richer graph for the read-likelihood caller; BASEPAIR and SV F1 fall for both callers. The fall is precision, not recall, and plan §9.24 traces it to exposure: 32 extra haplotypes offer multi-allelic sites the 4-haplotype graph cannot produce at all, and those sites are harder. Multi-allelic records go from about 2.3% of the call set to about 3.4% on both chromosomes tested.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik-z` on the 34-haplotype graph looked like a precision-for-recall trade (measured on chr20: 1,597 false-positive SNVs against 375). The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 3,297 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik-z` goes 0.9575 to 0.9565.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric. It has since been *bounded* rather than settled: false calls made by both callers on both graphs with no truth candidate anywhere nearby number 44 on chr6 and 40 on chr20, which is a lower bound on the benchmark's share of them.

## Cost

| arm | 4-hap wall | 34-hap wall | 4-hap RSS | 34-hap RSS | 4-hap variants | 34-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 346 s | 675 s | 6.0 GB | 7.1 GB | 288,849 | 294,626 |
| `poisson-z` | 160 s | 224 s | 6.4 GB | 6.3 GB | 289,002 | 294,835 |
| `readlik` | 289 s | 587 s | 7.6 GB | 7.8 GB | 286,451 | 284,450 |
| `readlik-nomismap` | 290 s | 461 s | 7.3 GB | 8.4 GB | 287,939 | 290,597 |
| `readlik-z` | 241 s | 642 s | 7.5 GB | 7.4 GB | 286,462 | 284,525 |

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
| `readlik` | Insertion (<50 bp) | 0.8501 | 0.8883 | +0.0382 |
| `readlik` | Deletion (<50 bp) | 0.8922 | 0.9040 | +0.0117 |
| `readlik-nomismap` | ALL | 0.9579 | 0.9546 | -0.0033 |
| `readlik-nomismap` | SNV | 0.9813 | 0.9726 | -0.0086 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8488 | 0.8852 | +0.0363 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8915 | 0.9017 | +0.0102 |
| `readlik-z` | ALL | 0.9585 | 0.9615 | +0.0030 |
| `readlik-z` | SNV | 0.9818 | 0.9794 | -0.0024 |
| `readlik-z` | Insertion (<50 bp) | 0.8498 | 0.8910 | +0.0413 |
| `readlik-z` | Deletion (<50 bp) | 0.8922 | 0.9069 | +0.0147 |

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
| `readlik` | ALL | 0.9359 | 0.8600 | -0.0759 |
| `readlik` | SNV | 0.9840 | 0.9790 | -0.0051 |
| `readlik` | Insertion (<50 bp) | 0.8212 | 0.6293 | -0.1919 |
| `readlik` | Deletion (<50 bp) | 0.8835 | 0.8692 | -0.0142 |
| `readlik-nomismap` | ALL | 0.9150 | 0.8351 | -0.0799 |
| `readlik-nomismap` | SNV | 0.9839 | 0.9758 | -0.0081 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7564 | 0.5847 | -0.1717 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8750 | 0.8580 | -0.0170 |
| `readlik-z` | ALL | 0.9362 | 0.8656 | -0.0706 |
| `readlik-z` | SNV | 0.9842 | 0.9802 | -0.0041 |
| `readlik-z` | Insertion (<50 bp) | 0.8216 | 0.6297 | -0.1919 |
| `readlik-z` | Deletion (<50 bp) | 0.8831 | 0.8735 | -0.0096 |

## Structural variants — truvari (GIAB `stvar`)

The SV metric. Reciprocal-overlap matching against the structural benchmark, `--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored against the *small-variant* truth set and therefore have essentially no truth to match above 50 bp (plan §9.22).

**These errors are broken down per record in [tier2-sv-errors.md](tier2-sv-errors.md)**, and three findings there bear directly on this table. The 34-haplotype false-positive rise is not the same errors plus more — only about two thirds of the 4-haplotype false calls survive the graph change, and the new ones are disproportionately calls with no truth candidate at all. A quarter of all false positives are placement or bookkeeping artefacts of the metric. And harmonising representation with `truvari refine` lifts every arm by roughly 0.05 F1.

| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5320 | 0.5512 | 0.4618 | 0.5490 | 0.4944 | **-0.0546** |
| `poisson-z` | 0.5488 | 0.5417 | 0.5468 | 0.4442 | 0.5478 | 0.4881 | **-0.0597** |
| `readlik` | 0.5301 | 0.5113 | 0.5654 | 0.4689 | 0.5472 | 0.4892 | **-0.0580** |
| `readlik-nomismap` | 0.5339 | 0.5236 | 0.5374 | 0.4283 | 0.5357 | 0.4712 | **-0.0645** |
| `readlik-z` | 0.5385 | 0.5314 | 0.5709 | 0.4655 | 0.5542 | 0.4963 | **-0.0579** |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7938 | 0.8019 | 0.8825 | 0.8645 | 0.8358 | 0.8320 | **-0.0038** |
| `sm50-poisson-z` | Deletion | 0.8945 | 0.9239 | 0.8348 | 0.7428 | 0.8636 | 0.8235 | **-0.0401** |
| `sm50-poisson-z` | ALL | 0.9219 | 0.9271 | 0.9498 | 0.9091 | 0.9356 | 0.9180 | **-0.0176** |
| `sm50-readlik-z` | Insertion | 0.8833 | 0.9102 | 0.8880 | 0.8763 | 0.8856 | 0.8930 | **+0.0073** |
| `sm50-readlik-z` | Deletion | 0.8956 | 0.9277 | 0.9070 | 0.8784 | 0.9012 | 0.9024 | **+0.0011** |
| `sm50-readlik-z` | ALL | 0.9440 | 0.9542 | 0.9714 | 0.9589 | 0.9575 | 0.9565 | **-0.0010** |

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. It matters for how the calls rank, which is a separate page: see [tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits `AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the called genotype explains. That rescales a quality and does not change a genotype, so **the unfiltered numbers on this page are unaffected by it**.

