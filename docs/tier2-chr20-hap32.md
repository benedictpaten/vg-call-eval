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

**The richer graph trades precision for recall, and every arm ends up worse on F1.** That is the honest headline and it holds on all three measures — small-variant GT, size-matched BASEPAIR, and SV. More haplotypes means more candidate alleles; more of the true ones get offered, and so do more wrong ones.

**But the two callers are not affected equally, and that is the useful result.**

| measure | `poisson-z` Δ | `readlik-z` Δ | ratio |
|---|---|---|---|
| small-variant ALL GT F1 | −0.0235 | **−0.0022** | 11x |
| size-matched ALL BASEPAIR F1 | −0.0260 | **−0.0060** | 4x |
| SV F1 | −0.0695 | **−0.0212** | 3x |

The read-likelihood model absorbs the extra allele ambiguity 3–11x better than the Poisson model on every axis. Scoring each read against each allele degrades gracefully as alleles multiply; aggregating depth does not.

**SV recall is the one outright win, and only for the read-likelihood arms.** `readlik-z` goes 0.5214 → 0.5726 (+0.0512) and `readlik` +0.0292, while `poisson-z` *falls* 0.0137. This is the design's thesis in one line: the extra haplotypes supply better SV alleles, but only a caller that scores alleles individually can use them. It still costs precision, so F1 declines — the operating point moved, it did not simply improve.

**The MAPQ mismapping term earns its keep here.** `readlik` against `readlik-nomismap` on ALL GT F1 is worth +0.0066 on this graph against +0.0007 on the 4-haplotype one, and it suppresses 7,287 calls against 752 — roughly ten times the work. More near-identical haplotypes means more chances for a read to fit an allele it did not come from, which is exactly what that term exists to damp.

**One caveat that this data cannot settle.** Some of the lost precision may not be error: a graph carrying 32 haplotypes will call real variation that a draft benchmark does not cover, and that scores as a false positive. The size-matched insertion row below — recall +0.037, precision −0.034, F1 flat — is equally consistent with "found more truth and more noise in equal measure" and with "found more truth than the benchmark knows about". Separating them needs a more complete truth set, not a different metric.

## Cost

| arm | 4-hap wall | 32-hap wall | 4-hap RSS | 32-hap RSS | 4-hap variants | 32-hap variants |
|---|---|---|---|---|---|---|
| `poisson` | 156 s | 293 s | 2.9 GB | 3.3 GB | 106,587 | 124,445 |
| `poisson-z` | 72 s | 106 s | 2.9 GB | 3.2 GB | 106,686 | 124,769 |
| `readlik` | 115 s | 134 s | 3.8 GB | 4.0 GB | 105,930 | 115,734 |
| `readlik-nomismap` | 115 s | 134 s | 3.5 GB | 4.0 GB | 106,682 | 123,021 |
| `readlik-z` | 97 s | 116 s | 3.5 GB | 3.8 GB | 105,936 | 115,787 |

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
| `readlik` | ALL | 0.9481 | 0.9431 | -0.0050 |
| `readlik` | SNV | 0.9764 | 0.9675 | -0.0089 |
| `readlik` | Insertion (<50 bp) | 0.8230 | 0.8617 | +0.0387 |
| `readlik` | Deletion (<50 bp) | 0.8705 | 0.8702 | -0.0003 |
| `readlik-nomismap` | ALL | 0.9474 | 0.9365 | -0.0109 |
| `readlik-nomismap` | SNV | 0.9761 | 0.9610 | -0.0151 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.8215 | 0.8553 | +0.0338 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8689 | 0.8673 | -0.0016 |
| `readlik-z` | ALL | 0.9482 | 0.9460 | -0.0022 |
| `readlik-z` | SNV | 0.9766 | 0.9700 | -0.0066 |
| `readlik-z` | Insertion (<50 bp) | 0.8231 | 0.8650 | +0.0420 |
| `readlik-z` | Deletion (<50 bp) | 0.8706 | 0.8733 | +0.0027 |

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
| `readlik` | ALL | 0.8969 | 0.8326 | -0.0643 |
| `readlik` | SNV | 0.9805 | 0.9725 | -0.0081 |
| `readlik` | Insertion (<50 bp) | 0.7243 | 0.5980 | -0.1263 |
| `readlik` | Deletion (<50 bp) | 0.8612 | 0.8403 | -0.0209 |
| `readlik-nomismap` | ALL | 0.8882 | 0.8204 | -0.0678 |
| `readlik-nomismap` | SNV | 0.9802 | 0.9676 | -0.0126 |
| `readlik-nomismap` | Insertion (<50 bp) | 0.7057 | 0.5834 | -0.1223 |
| `readlik-nomismap` | Deletion (<50 bp) | 0.8507 | 0.8286 | -0.0221 |
| `readlik-z` | ALL | 0.8973 | 0.8333 | -0.0640 |
| `readlik-z` | SNV | 0.9807 | 0.9737 | -0.0069 |
| `readlik-z` | Insertion (<50 bp) | 0.7247 | 0.5932 | -0.1315 |
| `readlik-z` | Deletion (<50 bp) | 0.8620 | 0.8425 | -0.0195 |

## Structural variants (GIAB `stvar`)

Recall is aardvark's published value. **Precision is recomputed** from its per-variant `BD` decisions, because its summary leaves the query columns at zero for the `Sv*` categories — without that, a run calling far more SVs would read as a pure recall win when it had traded precision away. F1 is derived from the two.

| arm | 4-hap recall | 32-hap recall | Δ | 4-hap prec | 32-hap prec | Δ | 4-hap F1 | 32-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 0.4643 | -0.0024 | 0.4712 | 0.3915 | -0.0797 | 0.4689 | 0.4248 | **-0.0441** |
| `poisson-z` | 0.5024 | 0.4887 | -0.0137 | 0.5187 | 0.4017 | -0.1170 | 0.5104 | 0.4409 | **-0.0695** |
| `readlik` | 0.4810 | 0.5101 | +0.0292 | 0.4671 | 0.4101 | -0.0570 | 0.4739 | 0.4547 | **-0.0193** |
| `readlik-nomismap` | 0.4821 | 0.5113 | +0.0292 | 0.4607 | 0.3988 | -0.0619 | 0.4712 | 0.4481 | **-0.0231** |
| `readlik-z` | 0.5214 | 0.5726 | +0.0512 | 0.5028 | 0.4294 | -0.0734 | 0.5120 | 0.4908 | **-0.0212** |

Per class, recall only:

| arm | class | 4-hap | 32-hap | Δ |
|---|---|---|---|---|
| `poisson` | SV insertion | 0.3877 | 0.3925 | +0.0048 |
| `poisson` | SV deletion | 0.5434 | 0.5340 | -0.0094 |
| `poisson-z` | SV insertion | 0.4263 | 0.4130 | -0.0133 |
| `poisson-z` | SV deletion | 0.5763 | 0.5622 | -0.0141 |
| `readlik` | SV insertion | 0.4553 | 0.4903 | +0.0350 |
| `readlik` | SV deletion | 0.5059 | 0.5293 | +0.0235 |
| `readlik-nomismap` | SV insertion | 0.4589 | 0.4940 | +0.0350 |
| `readlik-nomismap` | SV deletion | 0.5047 | 0.5282 | +0.0235 |
| `readlik-z` | SV insertion | 0.4976 | 0.5580 | +0.0604 |
| `readlik-z` | SV deletion | 0.5446 | 0.5869 | +0.0423 |

## Small variants restricted to <50 bp — BASEPAIR

The `smvar` truth set holds no record >=50 bp, so a large insertion called inside its confident region scores FP on every base however right it is. A richer graph calls more of those, which inflates the apparent BASEPAIR loss above. Dropping any record with a called allele >=50 bp from *both* sides is the only like-for-like read of these numbers.

| arm | class | 4-hap recall | 32-hap recall | 4-hap prec | 32-hap prec | 4-hap F1 | 32-hap F1 | **Δ F1** |
|---|---|---|---|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7637 | 0.7802 | 0.8700 | 0.8282 | 0.8134 | 0.8035 | **-0.0100** |
| `sm50-poisson-z` | Deletion | 0.8628 | 0.8949 | 0.8094 | 0.6969 | 0.8353 | 0.7836 | **-0.0517** |
| `sm50-poisson-z` | ALL | 0.8993 | 0.9109 | 0.9385 | 0.8746 | 0.9184 | 0.8924 | **-0.0260** |
| `sm50-readlik-z` | Insertion | 0.8578 | 0.8953 | 0.8624 | 0.8283 | 0.8601 | 0.8605 | **+0.0004** |
| `sm50-readlik-z` | Deletion | 0.8603 | 0.9009 | 0.8865 | 0.8346 | 0.8732 | 0.8665 | **-0.0067** |
| `sm50-readlik-z` | ALL | 0.9238 | 0.9424 | 0.9596 | 0.9283 | 0.9413 | 0.9353 | **-0.0061** |

