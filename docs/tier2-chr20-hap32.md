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
| `readlik-z` | 0.9489 | 0.9547 | **+0.0058** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0130** on the 4-haplotype graph to **+0.0423** on the 34-haplotype one — 3.3x wider.

**Two directions, and they are not the same direction.** GT F1 rises on the richer graph for the read-likelihood caller; BASEPAIR and SV F1 fall for both callers. The fall is precision, not recall, and plan §9.24 traces it to exposure: 32 extra haplotypes offer multi-allelic sites the 4-haplotype graph cannot produce at all, and those sites are harder. Multi-allelic records go from about 2.3% of the call set to about 3.4% on both chromosomes tested.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik-z` on the 34-haplotype graph looked like a precision-for-recall trade — 1,597 false-positive SNVs against the 4-haplotype graph's 375. The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 2,483 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik-z` goes 0.9436 to 0.9423.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric. It has since been *bounded* rather than settled: false calls made by both callers on both graphs with no truth candidate anywhere nearby number 44 on chr6 and 40 on chr20, which is a lower bound on the benchmark's share of them.

## Cost

| arm | 4-hap wall | 34-hap wall | 4-hap RSS | 34-hap RSS | 4-hap variants | 34-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 208 s | 318 s | 2.9 GB | 3.4 GB | 106,587 | 124,445 |
| `poisson-z` | 71 s | 111 s | 2.9 GB | 3.3 GB | 106,686 | 124,769 |
| `readlik` | 127 s | 157 s | 3.4 GB | 4.4 GB | 104,489 | 106,123 |
| `readlik-nomismap` | 129 s | 157 s | 3.4 GB | 4.4 GB | 106,347 | 121,467 |
| `readlik-z` | 108 s | 132 s | 3.8 GB | 3.7 GB | 104,495 | 106,175 |

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
| `readlik` | ALL | 0.9488 | 0.9514 | +0.0026 |
| `readlik` | SNV | 0.9756 | 0.9739 | -0.0017 |
| `readlik` | Insertion (<50 bp) | 0.8297 | 0.8746 | +0.0449 |
| `readlik` | Deletion (<50 bp) | 0.8755 | 0.8814 | +0.0059 |
| `readlik-nomismap` | ALL | 0.9488 | 0.9398 | -0.0089 |
| `readlik-nomismap` | SNV | 0.9759 | 0.9615 | -0.0144 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8281 | 0.8680 | +0.0400 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8749 | 0.8781 | +0.0032 |
| `readlik-z` | ALL | 0.9489 | 0.9547 | +0.0058 |
| `readlik-z` | SNV | 0.9757 | 0.9766 | +0.0009 |
| `readlik-z` | Insertion (<50 bp) | 0.8298 | 0.8783 | +0.0485 |
| `readlik-z` | Deletion (<50 bp) | 0.8753 | 0.8855 | +0.0101 |

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
| `readlik` | ALL | 0.9100 | 0.8524 | -0.0576 |
| `readlik` | SNV | 0.9799 | 0.9774 | -0.0025 |
| `readlik` | Insertion (<50 bp) | 0.7669 | 0.6316 | -0.1353 |
| `readlik` | Deletion (<50 bp) | 0.8560 | 0.8419 | -0.0140 |
| `readlik-nomismap` | ALL | 0.8980 | 0.8291 | -0.0689 |
| `readlik-nomismap` | SNV | 0.9801 | 0.9681 | -0.0120 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7404 | 0.6067 | -0.1337 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8427 | 0.8199 | -0.0227 |
| `readlik-z` | ALL | 0.9103 | 0.8544 | -0.0558 |
| `readlik-z` | SNV | 0.9800 | 0.9786 | -0.0014 |
| `readlik-z` | Insertion (<50 bp) | 0.7675 | 0.6267 | -0.1408 |
| `readlik-z` | Deletion (<50 bp) | 0.8558 | 0.8446 | -0.0111 |

## Structural variants — truvari (GIAB `stvar`)

The SV metric. Reciprocal-overlap matching against the structural benchmark, `--sizemin 50`. This replaced aardvark's `Sv*` categories, which are scored against the *small-variant* truth set and therefore have essentially no truth to match above 50 bp (plan §9.22).

**These errors are broken down per record in [tier2-sv-errors.md](tier2-sv-errors.md)**, and three findings there bear directly on this table. The 34-haplotype false-positive rise is not the same errors plus more — only about two thirds of the 4-haplotype false calls survive the graph change, and the new ones are disproportionately calls with no truth candidate at all. A quarter of all false positives are placement or bookkeeping artefacts of the metric. And harmonising representation with `truvari refine` lifts every arm by roughly 0.05 F1.

| arm | 4-hap recall | 34-hap recall | 4-hap prec | 34-hap prec | 4-hap F1 | 34-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|
| `poisson` | 0.4889 | 0.4810 | 0.5021 | 0.4289 | 0.4954 | 0.4535 | **-0.0419** |
| `poisson-z` | 0.4902 | 0.4824 | 0.4959 | 0.4029 | 0.4930 | 0.4391 | **-0.0540** |
| `readlik` | 0.4784 | 0.4706 | 0.5101 | 0.4440 | 0.4938 | 0.4569 | **-0.0368** |
| `readlik-nomismap` | 0.4771 | 0.4850 | 0.4728 | 0.4015 | 0.4750 | 0.4393 | **-0.0357** |
| `readlik-z` | 0.4797 | 0.4889 | 0.5021 | 0.4298 | 0.4907 | 0.4574 | **-0.0333** |

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
| `sm50-readlik-z` | Insertion | 0.8589 | 0.8939 | 0.8691 | 0.8488 | 0.8640 | 0.8708 | **+0.0068** |
| `sm50-readlik-z` | Deletion | 0.8656 | 0.9042 | 0.8914 | 0.8431 | 0.8783 | 0.8726 | **-0.0058** |
| `sm50-readlik-z` | ALL | 0.9248 | 0.9420 | 0.9630 | 0.9426 | 0.9436 | 0.9423 | **-0.0012** |

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. It matters for how the calls rank, which is a separate page: see [tier2-quality-signals.md](tier2-quality-signals.md). In short, `vg call` now emits `AD`, `BL` and `GQI` alongside `GQ`, and `GQ` is scaled by the fraction of reads the called genotype explains. That rescales a quality and does not change a genotype, so **the unfiltered numbers on this page are unaffected by it**.

