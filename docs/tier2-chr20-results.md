# Tier 2 results: HG002 chr20 on HPRC v2.1 MC CHM13, 34-haplotype graph

Real reads, real benchmark, run on a 32 GB laptop.

This is the **34-haplotype** graph: CHM13, GRCh38 and 32 recombinants from haplotype sampling. It is the primary subject because it is what the caller is tuned for -- both the linkage transition and the panel frequency prior are panel-size effects and have little to work with on a thin panel -- and because it is the better-performing configuration. The 4-haplotype graph has its own page at [tier2-chr20-4hap-results.md](tier2-chr20-4hap-results.md), and the two are put side by side in [tier2-chr20-graph-comparison.md](tier2-chr20-graph-comparison.md).

| | |
|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz`, 101,366,693 nodes, **34 haplotypes** (CHM13, GRCh38, 32 recombinants from haplotype sampling; the file is named for the recombinant count, not the total). HG002 itself is **absent** — no circularity |
| chromosome | chr20 component, 2,781,046 nodes |
| reads | 596,017,764 alignments genome-wide (~28.6×); 151 bp paired Illumina |
| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |
| regions | small variants 58.9 Mb; SVs 59.4 Mb |
| engine | `aardvark compare` for small variants; `truvari bench --sizemin 50` for SVs |

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.95`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap raised from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every number on this page — accuracy and cost alike — comes from one `vg` build in one pass, which is what the refresh harness exists to guarantee: a table whose rows come from different builds is not a comparison, it is a mixture of vintages.

Build: `vg version v1.4.0-18924-g8acbb43a2` (the code of `0cab3fbd4`, built before that commit, so the string names its parent).

The wall column is what the caller costs unaided, and the repeatability note below applies to it harder than to the memory column. It includes snarl decomposition, which is single-threaded — 46 s of a 197 s chr20 run — and which `vg call -r` skips for byte-identical output given `vg snarls -T -P <ref path>`. The whole-genome harness caches one snarl file per contig for exactly that reason; this matrix does not, so these figures include it.

| arm | enumeration | pack? | variants | wall | CPU | Δ wall | Δ CPU | peak RSS |
|---|---|---|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 124,445 | 295 s | 1,226 s (4.2x) | +38 s | — | 3.4 GB |
| `poisson-z` | panel (`-z`) | yes | 124,769 | 104 s | 294 s (2.8x) | -2 s | — | 3.2 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 117,707 | 872 s | 1,492 s (1.7x) | +720 s | — | 8.6 GB |
| `readlik-nomismap` | panel (default) | **no** | 130,969 | 171 s | 719 s (4.2x) | -134 s | — | 7.8 GB |
| `readlik-nolink` | panel (default) | **no** | 117,417 | 154 s | 668 s (4.3x) | +28 s | — | 7.7 GB |
| `readlik` | panel (default) | **no** | 112,207 | 172 s | 716 s (4.2x) | -138 s | — | 7.5 GB |

`CPU` is user+sys, with the multiple of wall clock beside it. It is the column that separates work from waiting: this caller has phases that run on one thread and phases that block on a subprocess, so a wall-clock change can come from either doing less or waiting less, and only CPU distinguishes them. A multiple well under `--threads` means the run spent its time parked rather than computing.

`Δ wall` is against the previous refresh, whose build was not recorded. Read it with the repeatability note below: run-to-run variance on this measurement is larger than most of these deltas, and the Poisson arms are the control -- their code is untouched by any read-likelihood change, so a Δ on those rows is the machine and not the caller.

**Peak RSS in this table is repeatable to about ±0.35 GB, so read it accordingly.** Three back-to-back runs of one binary on chr6-4hap, identical parameters and a warm cache, gave 7.3, 6.6 and 7.0 GB -- a 0.7 GB spread on a 7 GB measurement. Differences smaller than that are not evidence of anything, and a single measurement of each of two arms cannot resolve one. Thread count matters too: the same run at `--threads 6` instead of 5 measured 8.7 GB, because the read and GBWT caches are per thread. Wall clock is worse still -- a run immediately after a full rebuild took 956 s against 260 s warm, purely from page cache.

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9536 | 0.9581 | 0.9558 | 71,533 | 3,484 | 3,050 | 0.9626 | 0.9656 | 0.9641 |
| `poisson-z` | 0.9569 | 0.9583 | 0.9576 | 71,786 | 3,231 | 3,036 | 0.9658 | 0.9655 | 0.9656 |
| `readlik-support` | 0.9732 | 0.9907 | 0.9819 | 73,009 | 2,008 | 666 | 0.9778 | 0.9890 | 0.9834 |
| `readlik-nomismap` | 0.9777 | 0.9685 | 0.9730 | 73,341 | 1,676 | 2,330 | 0.9808 | 0.9723 | 0.9765 |
| `readlik-nolink` | 0.9735 | 0.9908 | 0.9821 | 73,032 | 1,985 | 660 | 0.9779 | 0.9890 | 0.9835 |
| `readlik` | 0.9754 | 0.9950 | **0.9851** | 73,173 | 1,844 | 357 | 0.9792 | 0.9921 | 0.9856 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7516 | 0.8184 | 0.7835 | 7,297 | 2,412 | 1,878 | 0.7817 | 0.5305 | 0.6321 |
| `poisson-z` | 0.7540 | 0.8220 | 0.7866 | 7,321 | 2,388 | 1,836 | 0.7890 | 0.5251 | 0.6306 |
| `readlik-support` | 0.8930 | 0.8636 | 0.8780 | 8,670 | 1,039 | 1,500 | 0.9134 | 0.6749 | 0.7762 |
| `readlik-nomismap` | 0.9194 | 0.8998 | 0.9095 | 8,926 | 783 | 1,072 | 0.9320 | 0.5379 | 0.6821 |
| `readlik-nolink` | 0.8944 | 0.8634 | 0.8787 | 8,684 | 1,025 | 1,503 | 0.9152 | 0.6619 | 0.7682 |
| `readlik` | 0.9201 | 0.9153 | **0.9177** | 8,933 | 776 | 889 | 0.9308 | 0.7636 | 0.8390 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8300 | 0.6849 | 0.7505 | 8,271 | 1,694 | 4,201 | 0.9098 | 0.5445 | 0.6813 |
| `poisson-z` | 0.8341 | 0.6848 | 0.7521 | 8,312 | 1,653 | 4,219 | 0.9096 | 0.4965 | 0.6424 |
| `readlik-support` | 0.9142 | 0.8583 | 0.8854 | 9,110 | 855 | 1,633 | 0.9320 | 0.6388 | 0.7581 |
| `readlik-nomismap` | 0.9417 | 0.9244 | 0.9330 | 9,384 | 581 | 826 | 0.9427 | 0.8094 | 0.8710 |
| `readlik-nolink` | 0.9152 | 0.8593 | 0.8863 | 9,120 | 845 | 1,621 | 0.9316 | 0.7411 | 0.8255 |
| `readlik` | 0.9419 | 0.9369 | **0.9394** | 9,386 | 579 | 679 | 0.9456 | 0.8646 | 0.9033 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7913 | 0.7382 | 0.7639 | 15,568 | 4,106 | 6,568 | 0.8470 | 0.5347 | 0.6555 |
| `poisson-z` | 0.7946 | 0.7392 | 0.7659 | 15,633 | 4,041 | 6,573 | 0.8504 | 0.4980 | 0.6281 |
| `readlik-support` | 0.9037 | 0.8623 | 0.8825 | 17,780 | 1,894 | 3,236 | 0.9228 | 0.6512 | 0.7636 |
| `readlik-nomismap` | 0.9307 | 0.9033 | 0.9168 | 18,310 | 1,364 | 2,208 | 0.9375 | 0.6397 | 0.7604 |
| `readlik-nolink` | 0.9050 | 0.8632 | 0.8836 | 17,804 | 1,870 | 3,216 | 0.9236 | 0.6999 | 0.7963 |
| `readlik` | 0.9311 | 0.9264 | **0.9288** | 18,319 | 1,355 | 1,635 | 0.9383 | 0.8109 | 0.8700 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9198 | 0.9017 | 0.9107 | 87,101 | 7,590 | 9,618 | 0.9140 | 0.7221 | 0.8068 |
| `poisson-z` | 0.9232 | 0.9019 | 0.9124 | 87,419 | 7,272 | 9,609 | 0.9176 | 0.6876 | 0.7861 |
| `readlik-support` | 0.9588 | 0.9590 | 0.9589 | 90,789 | 3,902 | 3,902 | 0.9601 | 0.8189 | 0.8839 |
| `readlik-nomismap` | 0.9679 | 0.9531 | 0.9604 | 91,651 | 3,040 | 4,538 | 0.9695 | 0.8051 | 0.8797 |
| `readlik-nolink` | 0.9593 | 0.9593 | 0.9593 | 90,836 | 3,855 | 3,876 | 0.9606 | 0.8547 | 0.9046 |
| `readlik` | 0.9662 | 0.9788 | **0.9725** | 91,492 | 3,199 | 1,992 | 0.9692 | 0.9294 | 0.9489 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates both callers, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (58.9 Mb vs 59.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

299 `readlik` calls carry an insertion allele of 200 bp or more, against 245 from `poisson-z`. Every base of those alleles inside the small-variant confident region is scored FP, and the size-matched control below measures what that does to each caller's precision.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7802 | 0.8282 | **0.8035** |
| `sm50-poisson-z` | Deletion | 0.8949 | 0.6969 | **0.7836** |
| `sm50-poisson-z` | ALL | 0.9109 | 0.8746 | **0.8924** |
| `sm50-readlik` | Insertion | 0.9247 | 0.8850 | **0.9044** |
| `sm50-readlik` | Deletion | 0.9363 | 0.9003 | **0.9179** |
| `sm50-readlik` | ALL | 0.9652 | 0.9699 | **0.9675** |

Restricting raises insertion BASEPAIR precision from 0.7636 to 0.8850 for `readlik` and from 0.5251 to 0.8282 for `poisson-z`. `readlik` minus `poisson-z` goes from +0.238 to +0.057, so the unrestricted comparison overstates the difference between them. Insertion BASEPAIR F1 is 0.8035 for `poisson-z` against 0.9044 for `readlik` restricted, and 0.6306 against 0.8390 unrestricted.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it. Split into single alleles there are 332 insertions of 200 bp or more; truvari confirms **71** and rejects **40**, and the other **221** fall outside the SV confident region or above truvari's 50 kb size cap, so it does not judge them. *Known bad output* lists the largest.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22).

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.4810 | 0.4289 | 0.4535 | 368 | 478 | 397 |
| `poisson-z` | 0.4824 | 0.4029 | 0.4391 | 369 | 535 | 396 |
| `readlik-support` | 0.5464 | 0.4370 | 0.4856 | 418 | 523 | 347 |
| `readlik-nomismap` | 0.5830 | 0.4061 | 0.4787 | 446 | 642 | 319 |
| `readlik-nolink` | 0.5712 | 0.4329 | 0.4926 | 437 | 558 | 328 |
| `readlik` | 0.5765 | 0.5017 | **0.5365** | 441 | 428 | 324 |

## Long reads — ONT, 16-haplotype E821 graph

Same sample, truth and confident regions as above; different reads and a different graph. The ONT reads (43.1x on chr20, mean length 33,449 bp on chr20) are aligned to the E821 graph: 16 haplotypes from haplotype sampling plus CHM13 and GRCh38, 18 in the panel against 34 above. So a figure here is not comparable with one in the short-read tables. The pair of columns is comparable, and it answers what the long-read preset buys on long-read data. The preset sets `--gap-open 1 --gap-extend 1 --mismap-min 0.05 --read-min-mapq 5 --insertion-nats 0.9 --read-phasing --regenotype`.

Build: `vg version v1.4.0-18924-g8acbb43a2` (the code of `0cab3fbd4`, built before that commit, so the string names its parent) — the same build as the short-read arms. Directories: `work/E821-chr20/results-0923-default`, `work/E821-chr20/results-0923-preset`.

| | short-read defaults | `--preset ont` | Δ |
|---|---|---|---|
| `readlik` ALL F1 | 0.9262 | 0.9582 | +0.0320 |
| `readlik` SNV F1 | 0.9849 | 0.9860 | +0.0011 |
| `readlik` indel F1 | 0.7482 | 0.8632 | +0.1150 |
| `readlik` indel recall | 0.7976 | 0.8735 | +0.0759 |
| `readlik` indel precision | 0.7046 | 0.8531 | +0.1485 |
| `readlik-nomismap` ALL F1 | 0.9262 | 0.9582 | +0.0320 |
| `readlik-nolink` ALL F1 | 0.9119 | 0.9227 | +0.0108 |
| `readlik` SV F1 (truvari) | 0.5548 | 0.5574 | +0.0026 |

Small variants are aardvark's GT comparison, with indels from its joint indel row.

| arm | configuration | variants | wall | CPU | peak RSS |
|---|---|---|---|---|---|
| `readlik` | short-read defaults | 118,416 | 176 s | 638 s | 5.8 GB |
| `readlik` | `--preset ont` | 114,861 | 205 s | 670 s | 5.1 GB |
| `readlik-nomismap` | short-read defaults | 118,856 | 175 s | 640 s | 5.8 GB |
| `readlik-nomismap` | `--preset ont` | 114,898 | 211 s | 688 s | 5.1 GB |
| `readlik-nolink` | short-read defaults | 122,031 | 169 s | 620 s | 5.9 GB |
| `readlik-nolink` | `--preset ont` | 118,610 | 176 s | 630 s | 5.0 GB |

Run serially on the same machine as the short-read arms, `--threads 5`.

*The MAPQ mismapping term.* `readlik` minus `readlik-nomismap` is 0.0000 ALL F1 on ONT at short-read defaults and 0.0000 under the preset, against +0.0120 on short reads. `--no-mismap-term` puts every read on the floor, so the term changes only a read whose `e_r = clamp(10^(−MAPQ/10), --mismap-min, --mismap-max)` is above it (MAPQ below 17 at the 0.02 floor), and it reduces a read on the ceiling to near-silence. At the default ceiling of 0.95 that is MAPQ 0 alone: 4.96% of chr20 short-read alignments against 0.63% of ONT, 8x fewer. Long reads anchor uniquely, so there is almost no ambiguous-placement class for the term to act on. Under the preset, `--read-min-mapq 5` drops those reads before scoring, so no read reaches the ceiling and the term acts only between MAPQ 5 and 13, where 10^(−MAPQ/10) lies above the preset's 0.05 floor.

*The linkage layer.* `readlik` minus `readlik-nolink` is +0.0144 on ONT at short-read defaults and +0.0355 under the preset, against +0.0132 on short reads. `--linkage-weight 0` turns the whole layer off, and read phasing and phase-driven re-genotyping run inside it, so under the preset that arm loses both as well as the transition model. Switch error across depth, and ONT's SV figures against short reads with confidence intervals, are in [coverage.md](coverage.md).

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads at those sites sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising it to 0.5 removed 94% of the excess false-positive SNVs, and the default is now **0.95**. A clamp that is inert on a sparse graph is not thereby harmless.

The two graphs are put side by side in [tier2-chr20-graph-comparison.md](tier2-chr20-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.95 (current defaults)** | 0.9725 | 0.9851 | 0.9177 | 0.9394 | 0.9489 |

Only the current row is available here: the preserved old-default arms (`arms.floor-1e-8.json`, `arms.readlik.json`) exist for the chr20 4-haplotype run alone, so the before-and-after is on [tier2-chr20-4hap-results.md](tier2-chr20-4hap-results.md). Mixing rows from two datasets into one table is exactly what the one-build-per-matrix rule forbids. The full grids are in plan §9.20-§9.21.

Raising the floor off 1e-8 changed **1,493 genotypes (1.41%)** on chr20, of which **94% were heterozygous → homozygous** (1/0→1/1: 614, 0/1→1/1: 606, 1/2→1/1: 184), and dropped 1,251 spurious non-reference calls. The failure it corrects is spurious heterozygosity: a few locally misaligned reads, each able to veto the homozygous hypothesis almost without bound, conjuring a second allele that is not there.

The floor was later re-swept at the corrected cap, on both graphs and both benchmarks, and settled at **0.02**. 0.05 wins on small-variant `GT` but costs about 0.01 of SV F1 — which the first sweep never saw, because it was scored on one benchmark only. Plan §9.21 records that as a process rule: a sweep that sets a default has to be scored on every benchmark the project runs.

## Known bad output

Neither benchmark scores these, so they appear in no metric on this page. They are recorded because they are plainly wrong and would mislead anyone reading the VCF.

`readlik` calls **14 records carrying an insertion of 10 kb or more** on chr20, with median DP **4,034.5** against a median of **29** over all of the contig's records. 4 of them start inside another giant's reference span, so the records overstate the number of independent events. Length is the net insertion, ALT minus REF. The five largest:

| position | net insertion | REF length | GT | DP | GQ |
|---|---|---|---|---|---|
| chr20:32,497,305 | 165,108 bp | 389,082 bp | `1\|0` | 228,151 | 256 |
| chr20:46,000,775 | 72,907 bp | 1 bp | `1\|0` | 6,134 | 256 |
| chr20:25,849,101 | 61,955 bp | 6,057 bp | `0\|1` | 9,902 | 4 |
| chr20:32,152,179 | 57,731 bp | 29,624 bp | `.\|1` | 8,563 | 256 |
| chr20:32,121,004 | 22,462 bp | 4,371 bp | `.\|1` | 3,599 | 139 |

The Poisson caller's expected depth (`XD`) never exceeds **123** anywhere on chr20. Median DP by called insertion length: 28 for 1 bp, 27 for 2-15 bp, 25 for 16-49 bp, 34 for 50-199 bp, 64 for 200-999 bp, 356 for >=1 kb. So the giants are collapsed-repeat pile-ups, not haplotypes.

The read-likelihood model cannot reject them, and the reason is structural rather than a tuning failure: it computes P(reads | genotype) **conditioned on the reads it is given**, and never asks whether that many reads should be there. The Poisson caller gets this for free, because an observed-vs-expected depth term is the whole of its model. A depth-plausibility guard is the obvious remedy, and the expected depth is already reachable — the read-likelihood caller subclasses `SupportBasedSnarlCaller` and holds a `TraversalSupportFinder` for allele enumeration.

The same blindness has a second consequence, found later and now corrected. Because the model only weighs reads it can see, it had no way to know that a heterozygous deletion produces *no* reads over the deleted interval, and its flat `1/ploidy` mixture asserted that both haplotypes contributed equally everywhere. That cost it 94% of heterozygous deletions above 1 kb and mis-genotyped two thirds of heterozygous insertions above 1 kb. Weighting each haplotype by the reads it is *expected* to contribute at the site is now the default and fixes both, without moving small variants at all — see [tier2-sv-errors.md](tier2-sv-errors.md). It did not remove the need for a depth term: it corrects the *relative* weight between a genotype's haplotypes, while the pile-ups above are a statement about *absolute* depth. That term is now also the default, at `--depth-term 0.1`, and the read arms in the tables on this page carry it — see [tier2-depth-term.md](tier2-depth-term.md). It does not resolve the pile-ups either: it detects them emphatically and still cannot outvote the read evidence at them, which is what the `DR` field and `--depth-quality` are for ([tier2-quality-signals.md](tier2-quality-signals.md)).

Filtering on depth is **not** that remedy, and that has now been tested properly rather than by two spot checks. Sweeping a two-sided cut on DP over a rolling local median, across both chromosomes and both graphs, against the one test a hard filter has to pass — beat lowering the GQ threshold to the same recall:

- a **minimum** fails in all eight dataset-by-benchmark cells. Few reads already means a small likelihood gap, so low depth depresses GQ on its own and a separate cut adds nothing;
- a **maximum** passes in exactly one configuration — 5x the local median, structural calls, 34-haplotype graph, worth about +0.025 precision — and is dominated everywhere else. The two original spot checks (DP 200 moving insertion BASEPAIR precision by 0.0001; DP 58 helping by +0.087 but costing SV insertion recall 0.4976 to 0.4167) were both right and both too narrow to conclude from.

What shipped instead attacks the same blindness from the other side: **GQ is scaled by the fraction of reads the called genotype explains**, which lowers the quality of a pile-up the call does not account for. It does not reach all of them: 3 of the 14 giants above still carry GQ 256. The giants remain output that no metric charges for, and they should be fixed because they are wrong, not because they cost a score. See [tier2-quality-signals.md](tier2-quality-signals.md).

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. `vg call` emits `AD` (per-allele read support, ties split fractionally), `BL` (mean absolute fit), `GQI` (the raw likelihood-ratio quality) and `GQ` (that ratio scaled by the fraction of reads the called genotype explains). The scaling rescales a quality and does not change a genotype, so **the numbers on this page are unaffected by it**; what it changes is how the calls rank. See [tier2-quality-signals.md](tier2-quality-signals.md).

## The genotype mixture

The read-likelihood arms on this page use the **length-weighted mixture**, which became the default after it was found that the flat `1/ploidy` weight breaks heterozygotes whose alleles differ in length. Unlike the `GQ` scaling above, this *does* change genotypes, so these numbers are not comparable with runs made before it. `--flat-mixture` restores the previous model exactly. Derivation and measurements: [tier2-sv-errors.md](tier2-sv-errors.md).

### At long read length

`w_h ∝ L_h + R − 1` counts the read START positions that yield an overlap, so R is the full read length; the *scoring* window is the read's overlap with the site, but the weight answers a sampling question, not a fit question. At ONT read length the correction all but vanishes. With the measured chr20 mean of 33,449 bp against Illumina's 151:

| | Illumina R=151 | ONT R=33,449 |
|---|---|---|
| 50 bp deletion | 1.33x, w = 0.571 | 1.00x, w = 0.500 |
| 300 bp Alu | 3.00x, w = 0.750 | 1.01x, w = 0.502 |
| 5 kb deletion | 34.3x, w = 0.972 | 1.15x, w = 0.535 |
| 30 kb deletion | 201x, w = 0.995 | 1.90x, w = 0.655 |

That is the model being right rather than distorted: a haplotype carrying a 5 kb deletion really does yield only 13% fewer 33 kb reads over the site, where it yields 34x fewer 151 bp reads. The flat mixture that `--flat-mixture` restores, and that loses large heterozygous deletions on short reads, is very nearly what the length weighting already computes for ONT.

The caveat is the mean, not the term. ONT read lengths are heavily skewed (mean 33,449, median 19,222, p10 1,845, p90 86,190), so one mean R summarises a distribution that spans two orders of magnitude. A per-read R is the obvious refinement and has not been measured.

## Raw aardvark summary rows

<details><summary><code>poisson</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson | GT | ALL | ALL | ALL | 94691 | 87101 | 7590 | 97839 | 88221 | 9618 | 0.9198445470002429 | 0.9016956428418116 | 0.9106796816751409 | 1684 | 624 |
| poisson | GT | ALL | ALL | Snv | 75017 | 71533 | 3484 | 72747 | 69697 | 3050 | 0.9535571937027607 | 0.9580738724620946 | 0.955810197223596 | 256 | 175 |
| poisson | GT | ALL | ALL | Insertion | 9709 | 7297 | 2412 | 10339 | 8461 | 1878 | 0.7515707075908951 | 0.8183576748234839 | 0.783543585324115 | 825 | 124 |
| poisson | GT | ALL | ALL | Deletion | 9965 | 8271 | 1694 | 13331 | 9130 | 4201 | 0.8300050175614652 | 0.6848698522241392 | 0.7504849741194335 | 603 | 291 |
| poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1422 | 933 | 489 |  | 0.6561181434599156 |  | 0 | 34 |
| poisson | GT | ALL | ALL | JointIndel | 19674 | 15568 | 4106 | 25092 | 18524 | 6568 | 0.7912981600081326 | 0.738243264785589 | 0.7638505601664448 | 1428 | 449 |
| poisson | BASEPAIR | ALL | ALL | ALL | 390888 | 357266 | 33622 | 494750 | 357266 | 137484 | 0.9139855917807659 | 0.7221141990904497 | 0.8067991662507706 |  |  |
| poisson | BASEPAIR | ALL | ALL | Snv | 200440 | 192941 | 7499 | 193848 | 187180 | 6668 | 0.9625873079225703 | 0.9656019149023978 | 0.9640922548348373 |  |  |
| poisson | BASEPAIR | ALL | ALL | Insertion | 95512 | 74664 | 20848 | 139632 | 74080 | 65552 | 0.7817237624591674 | 0.5305374126274779 | 0.6320901817385857 |  |  |
| poisson | BASEPAIR | ALL | ALL | Deletion | 99134 | 90193 | 8941 | 161934 | 88175 | 73759 | 0.9098089454677507 | 0.5445119616633938 | 0.6812827226871123 |  |  |
| poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 19846 | 9605 | 10241 |  | 0.48397661997379826 |  |  |  |
| poisson | BASEPAIR | ALL | ALL | JointIndel | 194646 | 164857 | 29789 | 321412 | 171860 | 149552 | 0.8469580674660666 | 0.5347031224720918 | 0.6555458408689969 |  |  |

</details>

<details><summary><code>poisson-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-z | GT | ALL | ALL | ALL | 94691 | 87419 | 7272 | 97946 | 88337 | 9609 | 0.9232028387069521 | 0.9018949216915443 | 0.9124244958135492 | 1741 | 616 |
| poisson-z | GT | ALL | ALL | Snv | 75017 | 71786 | 3231 | 72747 | 69711 | 3036 | 0.9569297625871469 | 0.9582663202606293 | 0.9575975750522121 | 299 | 174 |
| poisson-z | GT | ALL | ALL | Insertion | 9709 | 7321 | 2388 | 10314 | 8478 | 1836 | 0.7540426408486971 | 0.8219895287958116 | 0.786551400385407 | 835 | 105 |
| poisson-z | GT | ALL | ALL | Deletion | 9965 | 8312 | 1653 | 13384 | 9165 | 4219 | 0.83411941796287 | 0.6847728631201434 | 0.7521038182052947 | 607 | 283 |
| poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1501 | 983 | 518 |  | 0.6548967355096602 |  | 0 | 54 |
| poisson-z | GT | ALL | ALL | JointIndel | 19674 | 15633 | 4041 | 25199 | 18626 | 6573 | 0.7946020128087832 | 0.7391563157268145 | 0.7658769772648958 | 1442 | 442 |
| poisson-z | BASEPAIR | ALL | ALL | ALL | 390838 | 358634 | 32204 | 521548 | 358634 | 162914 | 0.9176026896054119 | 0.687633736492135 | 0.7861453376093013 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Snv | 200440 | 193587 | 6853 | 193852 | 187160 | 6692 | 0.9658102175214528 | 0.9654788188927635 | 0.9656444897740127 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Insertion | 95512 | 75356 | 20156 | 141150 | 74123 | 67027 | 0.788968925370634 | 0.5251363797378675 | 0.6305678602531807 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Deletion | 99134 | 90177 | 8957 | 177622 | 88186 | 89436 | 0.9096475477636331 | 0.49648129173188005 | 0.6423635969182704 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 29162 | 10953 | 18209 |  | 0.375591523215143 |  |  |  |
| poisson-z | BASEPAIR | ALL | ALL | JointIndel | 194646 | 165533 | 29113 | 347934 | 173262 | 174672 | 0.8504310389116653 | 0.4979737536429323 | 0.6281382845858628 |  |  |

</details>

<details><summary><code>readlik-support</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-support | GT | ALL | ALL | ALL | 94691 | 90789 | 3902 | 95180 | 91278 | 3902 | 0.958792282265474 | 0.9590039924353856 | 0.958898125664831 | 764 | 629 |
| readlik-support | GT | ALL | ALL | Snv | 75017 | 73009 | 2008 | 71678 | 71012 | 666 | 0.9732327339136463 | 0.9907084461061971 | 0.9818928380589611 | 95 | 348 |
| readlik-support | GT | ALL | ALL | Insertion | 9709 | 8670 | 1039 | 10995 | 9495 | 1500 | 0.8929858893809867 | 0.8635743519781719 | 0.8780338898609517 | 353 | 107 |
| readlik-support | GT | ALL | ALL | Deletion | 9965 | 9110 | 855 | 11523 | 9890 | 1633 | 0.9141996989463121 | 0.8582834331337326 | 0.8853595749153018 | 316 | 168 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 984 | 881 | 103 |  | 0.8953252032520326 |  | 0 | 6 |
| readlik-support | GT | ALL | ALL | JointIndel | 19674 | 17780 | 1894 | 23502 | 20266 | 3236 | 0.9037308122395039 | 0.8623095906731342 | 0.8825344488106772 | 669 | 281 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 390706 | 375123 | 15583 | 458076 | 375123 | 82953 | 0.9601157903897047 | 0.8189099625389673 | 0.8839089424610794 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 200440 | 195994 | 4446 | 193466 | 191340 | 2126 | 0.9778187986429855 | 0.989010989010989 | 0.983383049402518 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 95512 | 87236 | 8276 | 129434 | 87350 | 42084 | 0.9133512019432113 | 0.6748613192824142 | 0.7762001481211204 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 99134 | 92393 | 6741 | 144314 | 92191 | 52123 | 0.9320011297839288 | 0.638822290283687 | 0.7580522274743333 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11550 | 6247 | 5303 |  | 0.5408658008658008 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 194646 | 179629 | 15017 | 285298 | 185788 | 99510 | 0.9228496860968116 | 0.6512068083197219 | 0.7635888556398048 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 91651 | 3040 | 96738 | 92200 | 4538 | 0.9678955761371196 | 0.9530897889143873 | 0.9604356255228026 | 425 | 201 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 73341 | 1676 | 73902 | 71572 | 2330 | 0.9776583974299159 | 0.9684717598982436 | 0.9730433960677691 | 106 | 67 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 8926 | 783 | 10695 | 9623 | 1072 | 0.9193531774642084 | 0.8997662459093034 | 0.9094542628957619 | 179 | 66 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 9384 | 581 | 10932 | 10106 | 826 | 0.9416959357752133 | 0.9244420051225759 | 0.9329892072876508 | 140 | 66 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1209 | 899 | 310 |  | 0.7435897435897436 |  | 0 | 2 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 18310 | 1364 | 22836 | 20628 | 2208 | 0.9306699196909627 | 0.9033105622700999 | 0.9167861672605967 | 319 | 134 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390690 | 378757 | 11933 | 470436 | 378757 | 91679 | 0.9694566024213571 | 0.8051190810227108 | 0.8796784674948845 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 196597 | 3843 | 197252 | 191783 | 5469 | 0.9808271802035522 | 0.9722740453835702 | 0.9765318846000886 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 89016 | 6496 | 164784 | 88644 | 76140 | 0.9319876036518971 | 0.5379405767550247 | 0.6821475439680932 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 93458 | 5676 | 115168 | 93218 | 21950 | 0.9427441644642605 | 0.8094088635732147 | 0.8710032406860941 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 14802 | 6678 | 8124 |  | 0.4511552492906364 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 182474 | 12172 | 294754 | 188540 | 106214 | 0.9374659638523268 | 0.6396520488271575 | 0.7604402710042897 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 94691 | 90836 | 3855 | 95179 | 91303 | 3876 | 0.9592886335554593 | 0.9592767312117169 | 0.9592826823466685 | 758 | 618 |
| readlik-nolink | GT | ALL | ALL | Snv | 75017 | 73032 | 1985 | 71671 | 71011 | 660 | 0.973539331084954 | 0.9907912544822871 | 0.9820895344405254 | 91 | 340 |
| readlik-nolink | GT | ALL | ALL | Insertion | 9709 | 8684 | 1025 | 11006 | 9503 | 1503 | 0.8944278504480379 | 0.8634381246592767 | 0.8786598259139013 | 355 | 108 |
| readlik-nolink | GT | ALL | ALL | Deletion | 9965 | 9120 | 845 | 11517 | 9896 | 1621 | 0.9152032112393377 | 0.8592515411999653 | 0.8863452490828814 | 312 | 164 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 985 | 893 | 92 |  | 0.9065989847715736 |  | 0 | 6 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 19674 | 17804 | 1870 | 23508 | 20292 | 3216 | 0.9049506963505134 | 0.8631955079122001 | 0.8835800728339729 | 667 | 278 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 390690 | 375314 | 15376 | 439104 | 375314 | 63790 | 0.960643988840257 | 0.8547268984113103 | 0.9045955984256333 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 200440 | 196019 | 4421 | 193442 | 191322 | 2120 | 0.9779435242466573 | 0.9890406426732561 | 0.9834607801989937 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 95512 | 87413 | 8099 | 132022 | 87380 | 44642 | 0.9152043722254796 | 0.6618593870718517 | 0.7681827716547412 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 99134 | 92353 | 6781 | 124304 | 92124 | 32180 | 0.9315976355236347 | 0.7411185480756854 | 0.8255127723393881 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 9280 | 6398 | 2882 |  | 0.6894396551724138 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 194646 | 179766 | 14880 | 265606 | 185902 | 79704 | 0.923553527943035 | 0.6999164175508084 | 0.7963316826264596 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 94691 | 91492 | 3199 | 93897 | 91905 | 1992 | 0.9662164302837651 | 0.9787852647049426 | 0.9724602368359894 | 369 | 277 |
| readlik | GT | ALL | ALL | Snv | 75017 | 73173 | 1844 | 71681 | 71324 | 357 | 0.9754189050481891 | 0.9950196007310166 | 0.985121765129867 | 78 | 136 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 8933 | 776 | 10496 | 9607 | 889 | 0.9200741579977341 | 0.9153010670731707 | 0.9176814060669704 | 159 | 62 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 9386 | 579 | 10761 | 10082 | 679 | 0.9418966382338184 | 0.9369017749279807 | 0.9393925670555333 | 132 | 77 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 959 | 892 | 67 |  | 0.9301355578727841 |  | 0 | 2 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 18319 | 1355 | 22216 | 20581 | 1635 | 0.9311273762325912 | 0.9264043932301045 | 0.9287598803741659 | 291 | 141 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390682 | 378642 | 12040 | 407410 | 378642 | 28768 | 0.9691820969484133 | 0.9293880857121818 | 0.9488680502999653 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 196264 | 4176 | 192996 | 191474 | 1522 | 0.9791658351626422 | 0.992113826193289 | 0.9855973074188318 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 88902 | 6610 | 116174 | 88713 | 27461 | 0.9307940363514532 | 0.7636218086663109 | 0.8389612592734708 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 93737 | 5397 | 107840 | 93235 | 14605 | 0.9455585369298122 | 0.864567878338279 | 0.9032513212586625 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 8200 | 6358 | 1842 |  | 0.7753658536585366 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 182639 | 12007 | 232214 | 188306 | 43908 | 0.9383136565868294 | 0.8109157931907637 | 0.8699754777048591 |  |  |

</details>

