# chr20: 4-haplotype vs 32-haplotype graph

Same sample, same reads, same truth, same confident regions, same reference sequence. What changes is the graph — and, unavoidably, the alignments.

| | 4-haplotype | 32-haplotype |
|---|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz` | `…HG002.hap32.gbz` |
| haplotypes | 4 (CHM13, GRCh38, 2 recombinants) | **34** (CHM13, GRCh38, 32 recombinants) |
| HG002 present? | no | **no** — samples are `CHM13`, `GRCh38`, `recombination` |
| alignments | `…HG002.gaf.gz` | `…HG002.hap32.gaf.gz` (remapped) |

**This is not a single-variable experiment.** Reads mapped to one graph cannot be scored against the other, because the node ID spaces differ — so the 32-haplotype arm necessarily uses its own alignments. Graph and alignment move together. That is what adopting a richer graph actually involves, but it means a difference below cannot be attributed to the graph alone.

The rows to watch are the **`-z` arms**, which enumerate alleles from the GBWT haplotypes. Going from 4 to 34 changes which alleles are *available to call* rather than how they are scored, and the tier-2 finding was that enumeration matters more than the genotyper — most of all for SVs. This is the direct test.

## What this says

**The read-likelihood caller is better on the richer graph; the Poisson caller is much worse on it.** That split is the result. More haplotypes offer more true alleles and more wrong ones, and what decides the outcome is whether the genotyper can tell them apart read by read.

| arm | 4-hap GT F1 | 32-hap GT F1 | Δ |
|---|---|---|---|
| `poisson-z` | 0.9359 | 0.9124 | **-0.0235** |
| `readlik-z` | 0.9490 | 0.9547 | **+0.0057** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0131** on the 4-haplotype graph to **+0.0423** on the 32-haplotype one — 3.2x wider.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik-z` on the 32-haplotype graph carried 1,597 false-positive SNVs against the 4-haplotype graph's 375, and looked like a precision-for-recall trade. The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 2,481 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik-z` goes 0.9441 to 0.9428.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric.

## Cost

| arm | 4-hap wall | 32-hap wall | 4-hap RSS | 32-hap RSS | 4-hap variants | 32-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 156 s | 293 s | 2.9 GB | 3.3 GB | 106,587 | 124,445 |
| `poisson-z` | 72 s | 106 s | 2.9 GB | 3.2 GB | 106,686 | 124,769 |
| `readlik` | 154 s | 144 s | 3.5 GB | 4.5 GB | 104,462 | 106,100 |
| `readlik-nomismap` | 149 s | 142 s | 3.5 GB | 4.1 GB | 106,295 | 121,427 |
| `readlik-z` | 118 s | 124 s | 3.7 GB | 3.9 GB | 104,470 | 106,172 |

## Small variants — GT F1

| arm | class | 4-hap | 32-hap | Δ |
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

| arm | class | 4-hap | 32-hap | Δ |
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

## Structural variants (GIAB `stvar`)

Recall is aardvark's published value. **Precision is recomputed** from its per-variant `BD` decisions, because its summary leaves the query columns at zero for the `Sv*` categories — without that, a run calling far more SVs would read as a pure recall win when it had traded precision away. F1 is derived from the two.

| arm | 4-hap recall | 32-hap recall | Δ | 4-hap prec | 32-hap prec | Δ | 4-hap F1 | 32-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 0.4643 | -0.0024 | 0.4712 | 0.3915 | -0.0797 | 0.4689 | 0.4248 | **-0.0441** |
| `poisson-z` | 0.5024 | 0.4887 | -0.0137 | 0.5187 | 0.4017 | -0.1170 | 0.5104 | 0.4409 | **-0.0695** |
| `readlik` | 0.4839 | 0.5113 | +0.0274 | 0.4680 | 0.4138 | -0.0541 | 0.4758 | 0.4574 | **-0.0184** |
| `readlik-nomismap` | 0.4810 | 0.5119 | +0.0310 | 0.4586 | 0.3975 | -0.0611 | 0.4695 | 0.4475 | **-0.0220** |
| `readlik-z` | 0.5250 | 0.5690 | +0.0440 | 0.5045 | 0.4321 | -0.0724 | 0.5145 | 0.4912 | **-0.0234** |

Per class, recall only:

| arm | class | 4-hap | 32-hap | Δ |
|---|---|---|---|---|
| `poisson` | SV insertion | 0.3877 | 0.3925 | +0.0048 |
| `poisson` | SV deletion | 0.5434 | 0.5340 | -0.0094 |
| `poisson-z` | SV insertion | 0.4263 | 0.4130 | -0.0133 |
| `poisson-z` | SV deletion | 0.5763 | 0.5622 | -0.0141 |
| `readlik` | SV insertion | 0.4577 | 0.4940 | +0.0362 |
| `readlik` | SV deletion | 0.5094 | 0.5282 | +0.0188 |
| `readlik-nomismap` | SV insertion | 0.4589 | 0.4964 | +0.0374 |
| `readlik-nomismap` | SV deletion | 0.5023 | 0.5270 | +0.0246 |
| `readlik-z` | SV insertion | 0.4988 | 0.5556 | +0.0568 |
| `readlik-z` | SV deletion | 0.5505 | 0.5822 | +0.0317 |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 32-hap recall | 4-hap prec | 32-hap prec | 4-hap F1 | 32-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7637 | 0.7802 | 0.8700 | 0.8282 | 0.8134 | 0.8035 | **-0.0100** |
| `sm50-poisson-z` | Deletion | 0.8628 | 0.8949 | 0.8094 | 0.6969 | 0.8353 | 0.7836 | **-0.0517** |
| `sm50-poisson-z` | ALL | 0.8993 | 0.9109 | 0.9385 | 0.8746 | 0.9184 | 0.8924 | **-0.0260** |
| `sm50-readlik-z` | Insertion | 0.8632 | 0.8968 | 0.8666 | 0.8425 | 0.8649 | 0.8688 | **+0.0039** |
| `sm50-readlik-z` | Deletion | 0.8639 | 0.9016 | 0.8970 | 0.8528 | 0.8801 | 0.8765 | **-0.0036** |
| `sm50-readlik-z` | ALL | 0.9254 | 0.9421 | 0.9635 | 0.9436 | 0.9441 | 0.9428 | **-0.0012** |

