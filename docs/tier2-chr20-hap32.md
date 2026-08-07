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
| `readlik-z` | 0.9479 | 0.9520 | **+0.0041** |

The read-likelihood caller's margin over the Poisson caller goes from **+0.0120** on the 4-haplotype graph to **+0.0396** on the 32-haplotype one — 3.3x wider.

**This depended on a default that was wrong for graphs like this.** With `--mismap-max` at its old 0.1, `readlik-z` on the 32-haplotype graph carried 1,597 false-positive SNVs against the 4-haplotype graph's 375, and looked like a precision-for-recall trade. The cap was overriding the mapper: at those sites 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. At the current default of 0.5 that excess is 94% gone. Harness plan §9.20 has the derivation; the point for this page is that a caller-level default, not the graph, was the difference between the two readings.

**`readlik-nomismap` is the control.** It disables the mismapping term entirely, so the cap cannot reach it — and on the richer graph it still carries 2,593 spurious SNVs. The term is what does the work.

Size-matched to <50 bp — the only like-for-like read of the BASEPAIR numbers — `readlik-z` goes 0.9416 to 0.9392.

**One caveat this data cannot settle.** Some of the remaining false positives may not be error: a graph carrying 32 haplotypes will call real variation a draft benchmark does not cover, and that scores as a false positive. Separating them needs a more complete truth set, not a different metric.

## Cost

| arm | 4-hap wall | 32-hap wall | 4-hap RSS | 32-hap RSS | 4-hap variants | 32-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 156 s | 293 s | 2.9 GB | 3.3 GB | 106,587 | 124,445 |
| `poisson-z` | 72 s | 106 s | 2.9 GB | 3.2 GB | 106,686 | 124,769 |
| `readlik` | 118 s | 144 s | 3.9 GB | 4.4 GB | 104,725 | 106,619 |
| `readlik-nomismap` | 118 s | 143 s | 3.8 GB | 4.6 GB | 106,682 | 123,021 |
| `readlik-z` | 101 s | 122 s | 3.8 GB | 3.8 GB | 104,733 | 106,690 |

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
| `readlik` | ALL | 0.9478 | 0.9489 | +0.0011 |
| `readlik` | SNV | 0.9760 | 0.9743 | -0.0018 |
| `readlik` | Insertion (<50 bp) | 0.8235 | 0.8626 | +0.0391 |
| `readlik` | Deletion (<50 bp) | 0.8706 | 0.8707 | +0.0001 |
| `readlik-nomismap` | ALL | 0.9474 | 0.9365 | -0.0109 |
| `readlik-nomismap` | SNV | 0.9761 | 0.9610 | -0.0151 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8215 | 0.8553 | +0.0338 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8689 | 0.8673 | -0.0016 |
| `readlik-z` | ALL | 0.9479 | 0.9520 | +0.0041 |
| `readlik-z` | SNV | 0.9761 | 0.9768 | +0.0007 |
| `readlik-z` | Insertion (<50 bp) | 0.8236 | 0.8662 | +0.0426 |
| `readlik-z` | Deletion (<50 bp) | 0.8706 | 0.8743 | +0.0037 |

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
| `readlik` | ALL | 0.9003 | 0.8397 | -0.0607 |
| `readlik` | SNV | 0.9802 | 0.9776 | -0.0026 |
| `readlik` | Insertion (<50 bp) | 0.7338 | 0.6031 | -0.1307 |
| `readlik` | Deletion (<50 bp) | 0.8665 | 0.8463 | -0.0202 |
| `readlik-nomismap` | ALL | 0.8882 | 0.8204 | -0.0678 |
| `readlik-nomismap` | SNV | 0.9802 | 0.9676 | -0.0126 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7057 | 0.5834 | -0.1223 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8507 | 0.8286 | -0.0221 |
| `readlik-z` | ALL | 0.9003 | 0.8416 | -0.0588 |
| `readlik-z` | SNV | 0.9804 | 0.9789 | -0.0015 |
| `readlik-z` | Insertion (<50 bp) | 0.7344 | 0.5987 | -0.1357 |
| `readlik-z` | Deletion (<50 bp) | 0.8651 | 0.8498 | -0.0153 |

## Structural variants (GIAB `stvar`)

Recall is aardvark's published value. **Precision is recomputed** from its per-variant `BD` decisions, because its summary leaves the query columns at zero for the `Sv*` categories — without that, a run calling far more SVs would read as a pure recall win when it had traded precision away. F1 is derived from the two.

| arm | 4-hap recall | 32-hap recall | Δ | 4-hap prec | 32-hap prec | Δ | 4-hap F1 | 32-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 0.4643 | -0.0024 | 0.4712 | 0.3915 | -0.0797 | 0.4689 | 0.4248 | **-0.0441** |
| `poisson-z` | 0.5024 | 0.4887 | -0.0137 | 0.5187 | 0.4017 | -0.1170 | 0.5104 | 0.4409 | **-0.0695** |
| `readlik` | 0.4833 | 0.5101 | +0.0268 | 0.4716 | 0.4131 | -0.0584 | 0.4774 | 0.4565 | **-0.0208** |
| `readlik-nomismap` | 0.4821 | 0.5113 | +0.0292 | 0.4607 | 0.3988 | -0.0619 | 0.4712 | 0.4481 | **-0.0231** |
| `readlik-z` | 0.5250 | 0.5696 | +0.0446 | 0.5081 | 0.4309 | -0.0772 | 0.5164 | 0.4906 | **-0.0257** |

Per class, recall only:

| arm | class | 4-hap | 32-hap | Δ |
|---|---|---|---|---|
| `poisson` | SV insertion | 0.3877 | 0.3925 | +0.0048 |
| `poisson` | SV deletion | 0.5434 | 0.5340 | -0.0094 |
| `poisson-z` | SV insertion | 0.4263 | 0.4130 | -0.0133 |
| `poisson-z` | SV deletion | 0.5763 | 0.5622 | -0.0141 |
| `readlik` | SV insertion | 0.4553 | 0.4903 | +0.0350 |
| `readlik` | SV deletion | 0.5106 | 0.5293 | +0.0188 |
| `readlik-nomismap` | SV insertion | 0.4589 | 0.4940 | +0.0350 |
| `readlik-nomismap` | SV deletion | 0.5047 | 0.5282 | +0.0235 |
| `readlik-z` | SV insertion | 0.4976 | 0.5543 | +0.0568 |
| `readlik-z` | SV deletion | 0.5516 | 0.5845 | +0.0329 |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 32-hap recall | 4-hap prec | 32-hap prec | 4-hap F1 | 32-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7637 | 0.7802 | 0.8700 | 0.8282 | 0.8134 | 0.8035 | **-0.0100** |
| `sm50-poisson-z` | Deletion | 0.8628 | 0.8949 | 0.8094 | 0.6969 | 0.8353 | 0.7836 | **-0.0517** |
| `sm50-poisson-z` | ALL | 0.8993 | 0.9109 | 0.9385 | 0.8746 | 0.9184 | 0.8924 | **-0.0260** |
| `sm50-readlik-z` | Insertion | 0.8574 | 0.8949 | 0.8639 | 0.8340 | 0.8606 | 0.8634 | **+0.0028** |
| `sm50-readlik-z` | Deletion | 0.8600 | 0.8996 | 0.8874 | 0.8379 | 0.8735 | 0.8676 | **-0.0059** |
| `sm50-readlik-z` | ALL | 0.9233 | 0.9413 | 0.9605 | 0.9371 | 0.9416 | 0.9392 | **-0.0024** |

