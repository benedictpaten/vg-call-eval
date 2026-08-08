# chr20: 4-haplotype vs 34-haplotype graph

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
| `poisson-z` | 0.9359 | 0.9124 | **-0.0235** |
| `readlik-z` | 0.9490 | 0.9547 | **+0.0057** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0131** on the 4-haplotype graph to **+0.0423** on the 34-haplotype one — 3.2x wider.

**Two directions, and they are not the same direction.** GT F1 rises on the richer graph for the read-likelihood caller; BASEPAIR and SV F1 fall for both callers. The fall is precision, not recall, and plan §9.24 traces it to exposure: 32 extra haplotypes offer multi-allelic sites the 4-haplotype graph cannot produce at all, and those sites are harder. Multi-allelic records go from about 2.3% of the call set to about 3.4% on both chromosomes tested.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik-z` on the 34-haplotype graph looked like a precision-for-recall trade — 1,597 false-positive SNVs against the 4-haplotype graph's 375. The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 2,481 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik-z` goes 0.9441 to 0.9428.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric.

## Cost

| arm | 4-hap wall | 34-hap wall | 4-hap RSS | 34-hap RSS | 4-hap variants | 34-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 168 s | 294 s | 3.0 GB | 3.5 GB | 106,587 | 124,445 |
| `poisson-z` | 75 s | 112 s | 3.0 GB | 3.1 GB | 106,686 | 124,769 |
| `readlik` | 120 s | 148 s | 3.6 GB | 4.0 GB | 104,462 | 106,100 |
| `readlik-nomismap` | 115 s | 149 s | 3.5 GB | 4.0 GB | 106,295 | 121,427 |
| `readlik-z` | 99 s | 126 s | 3.7 GB | 3.8 GB | 104,470 | 106,172 |

## Small variants — GT F1

| arm | class | 4-hap | 34-hap | Δ |
|---|---|---|---|---|
| `poisson` | ALL | 0.9355 | 0.9107 | -0.0248 |
| `poisson` | SNV | 0.9733 | 0.9558 | -0.0175 |
| `poisson` | Insertion (<50 bp) | 0.7827 | 0.7835 | +0.0008 |
| `poisson` | Deletion (<50 bp) | 0.8148 | 0.7505 | -0.0643 |
| `poisson-z` | ALL | 0.9359 | 0.9124 | -0.0235 |
| `poisson-z` | SNV | 0.9735 | 0.9576 | -0.0159 |
| `poisson-z` | Insertion (<50 bp) | 0.7850 | 0.7866 | +0.0015 |
| `poisson-z` | Deletion (<50 bp) | 0.8148 | 0.7521 | -0.0627 |
| `readlik` | ALL | 0.9489 | 0.9515 | +0.0026 |
| `readlik` | SNV | 0.9756 | 0.9739 | -0.0017 |
| `readlik` | Insertion (<50 bp) | 0.8291 | 0.8729 | +0.0438 |
| `readlik` | Deletion (<50 bp) | 0.8771 | 0.8843 | +0.0073 |
| `readlik-nomismap` | ALL | 0.9488 | 0.9397 | -0.0090 |
| `readlik-nomismap` | SNV | 0.9759 | 0.9615 | -0.0144 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8273 | 0.8647 | +0.0374 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8758 | 0.8809 | +0.0051 |
| `readlik-z` | ALL | 0.9490 | 0.9547 | +0.0057 |
| `readlik-z` | SNV | 0.9757 | 0.9765 | +0.0008 |
| `readlik-z` | Insertion (<50 bp) | 0.8293 | 0.8764 | +0.0472 |
| `readlik-z` | Deletion (<50 bp) | 0.8771 | 0.8883 | +0.0112 |

## Small variants — BASEPAIR F1

| arm | class | 4-hap | 34-hap | Δ |
|---|---|---|---|---|
| `poisson` | ALL | 0.8825 | 0.8068 | -0.0757 |
| `poisson` | SNV | 0.9783 | 0.9641 | -0.0142 |
| `poisson` | Insertion (<50 bp) | 0.7648 | 0.6321 | -0.1327 |
| `poisson` | Deletion (<50 bp) | 0.7630 | 0.6813 | -0.0817 |
| `poisson-z` | ALL | 0.8867 | 0.7861 | -0.1006 |
| `poisson-z` | SNV | 0.9784 | 0.9656 | -0.0127 |
| `poisson-z` | Insertion (<50 bp) | 0.7712 | 0.6306 | -0.1406 |
| `poisson-z` | Deletion (<50 bp) | 0.7782 | 0.6424 | -0.1358 |
| `readlik` | ALL | 0.9038 | 0.8445 | -0.0593 |
| `readlik` | SNV | 0.9799 | 0.9774 | -0.0024 |
| `readlik` | Insertion (<50 bp) | 0.7398 | 0.6079 | -0.1320 |
| `readlik` | Deletion (<50 bp) | 0.8740 | 0.8603 | -0.0138 |
| `readlik-nomismap` | ALL | 0.8912 | 0.8202 | -0.0711 |
| `readlik-nomismap` | SNV | 0.9801 | 0.9681 | -0.0120 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7109 | 0.5764 | -0.1345 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8570 | 0.8395 | -0.0175 |
| `readlik-z` | ALL | 0.9041 | 0.8460 | -0.0581 |
| `readlik-z` | SNV | 0.9800 | 0.9786 | -0.0014 |
| `readlik-z` | Insertion (<50 bp) | 0.7406 | 0.6029 | -0.1377 |
| `readlik-z` | Deletion (<50 bp) | 0.8737 | 0.8626 | -0.0111 |

## Structural variants — truvari (GIAB `stvar`)

The SV metric. Reciprocal-overlap matching against the structural benchmark, `--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored against the *small-variant* truth set and therefore have essentially no truth to match above 50 bp (plan §9.22).

| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|
| `poisson` | 0.4889 | 0.4810 | 0.5021 | 0.4289 | 0.4954 | 0.4535 | **-0.0419** |
| `poisson-z` | 0.4902 | 0.4824 | 0.4959 | 0.4029 | 0.4930 | 0.4391 | **-0.0540** |
| `readlik` | 0.4601 | 0.4484 | 0.5052 | 0.4395 | 0.4816 | 0.4439 | **-0.0377** |
| `readlik-nomismap` | 0.4575 | 0.4614 | 0.4678 | 0.3950 | 0.4626 | 0.4256 | **-0.0369** |
| `readlik-z` | 0.4654 | 0.4719 | 0.5007 | 0.4187 | 0.4824 | 0.4437 | **-0.0387** |

## Structural variants — aardvark (secondary)

Kept for continuity with earlier runs. Recall is aardvark's published value; **precision is recomputed** from its per-variant `BD` decisions, because its summary leaves the query columns at zero for the `Sv*` categories. Prefer the truvari table above: these categories are scored against a truth set with no record over 50 bp.

| arm | 4-hap recall | 34-hap recall | Δ | 4-hap prec | 34-hap prec | Δ | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 0.4643 | -0.0024 | 0.4712 | 0.3915 | -0.0797 | 0.4689 | 0.4248 | **-0.0441** |
| `poisson-z` | 0.5024 | 0.4887 | -0.0137 | 0.5187 | 0.4017 | -0.1170 | 0.5104 | 0.4409 | **-0.0695** |
| `readlik` | 0.4839 | 0.5113 | +0.0274 | 0.4680 | 0.4138 | -0.0541 | 0.4758 | 0.4574 | **-0.0184** |
| `readlik-nomismap` | 0.4810 | 0.5119 | +0.0310 | 0.4586 | 0.3975 | -0.0611 | 0.4695 | 0.4475 | **-0.0220** |
| `readlik-z` | 0.5250 | 0.5690 | +0.0440 | 0.5045 | 0.4321 | -0.0724 | 0.5145 | 0.4912 | **-0.0234** |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7637 | 0.7802 | 0.8700 | 0.8282 | 0.8134 | 0.8035 | **-0.0100** |
| `sm50-poisson-z` | Deletion | 0.8628 | 0.8949 | 0.8094 | 0.6969 | 0.8353 | 0.7836 | **-0.0517** |
| `sm50-poisson-z` | ALL | 0.8993 | 0.9109 | 0.9385 | 0.8746 | 0.9184 | 0.8924 | **-0.0260** |
| `sm50-readlik-z` | Insertion | 0.8632 | 0.8968 | 0.8666 | 0.8425 | 0.8649 | 0.8688 | **+0.0039** |
| `sm50-readlik-z` | Deletion | 0.8639 | 0.9016 | 0.8970 | 0.8528 | 0.8801 | 0.8765 | **-0.0036** |
| `sm50-readlik-z` | ALL | 0.9254 | 0.9421 | 0.9635 | 0.9436 | 0.9441 | 0.9428 | **-0.0012** |

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. It matters for how the calls rank, which is a separate page: see [tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits `AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the called genotype explains. That rescales a quality and does not change a genotype, so **the unfiltered numbers on this page are unaffected by it**.

