# Tier 2 results: HG002 chr20 on HPRC v2.1 MC CHM13, 4-haplotype graph

Real reads, real benchmark, run on a 32 GB laptop.

This is the **4-haplotype** graph: CHM13, GRCh38 and 2 recombinants. It is kept as a thin-panel reference rather than the headline configuration -- the caller is tuned on the 34-haplotype graph, whose page is [tier2-chr20-results.md](tier2-chr20-results.md). The two are compared directly in [tier2-chr20-graph-comparison.md](tier2-chr20-graph-comparison.md).

| | |
|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz`, 100,179,277 nodes, **4 haplotypes** (CHM13, GRCh38, 2 recombinants). HG002 itself is **absent** — no circularity |
| chromosome | chr20 component, 2,382,533 nodes |
| reads | 596,017,764 alignments genome-wide (~28.6×); 151 bp paired Illumina |
| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |
| regions | small variants 58.9 Mb; SVs 59.4 Mb |
| engine | `aardvark compare` for small variants; `truvari bench --sizemin 50` for SVs |

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.7`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap down from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every arm on this page was re-run together on one build, so the wall-clock column compares runs made on the same machine in the same session rather than a mixture of vintages.

Two changes since the accuracy results were first produced left the calls untouched. The read path was optimised (vg `44fd008`) — on chr20 `readlik` went **506 s to under 100 s**, so the read-likelihood caller is now near parity with the Poisson caller at matched enumeration rather than 5.9x, and `readlik-support` is *faster* than `poisson`. Then `AD`, `BL`, `GQI` and the explained-share scaling of `GQ` were added — which rescales a quality and does not change a genotype. Both are confirmed by the variant counts below, which are unchanged to the record.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 106,587 | 190 s | 2.9 GB |
| `poisson-z` | panel (`-z`) | yes | 106,686 | 74 s | 2.9 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 104,159 | 126 s | 3.5 GB |
| `readlik-nomismap` | panel (default) | **no** | 106,359 | 104 s | 3.8 GB |
| `readlik-nolink` | panel (default) | **no** | 104,165 | 104 s | 3.8 GB |
| `readlik` | panel (default) | **no** | 104,165 | 105 s | 3.5 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9556 | 0.9917 | 0.9733 | 71,684 | 3,333 | 588 | 0.9659 | 0.9909 | 0.9783 |
| `poisson-z` | 0.9562 | 0.9914 | 0.9735 | 71,734 | 3,283 | 605 | 0.9664 | 0.9906 | 0.9784 |
| `readlik-support` | 0.9573 | 0.9943 | 0.9754 | 71,815 | 3,202 | 403 | 0.9673 | 0.9927 | 0.9798 |
| `readlik-nomismap` | 0.9579 | 0.9950 | **0.9761** | 71,859 | 3,158 | 353 | 0.9676 | 0.9932 | 0.9803 |
| `readlik-nolink` | 0.9576 | 0.9942 | 0.9756 | 71,839 | 3,178 | 408 | 0.9676 | 0.9926 | 0.9799 |
| `readlik` | 0.9570 | 0.9956 | 0.9759 | 71,789 | 3,228 | 312 | 0.9670 | 0.9936 | 0.9801 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7281 | 0.8462 | 0.7827 | 7,069 | 2,640 | 1,462 | 0.7684 | 0.7613 | 0.7648 |
| `poisson-z` | 0.7295 | 0.8497 | 0.7850 | 7,083 | 2,626 | 1,426 | 0.7729 | 0.7694 | 0.7712 |
| `readlik-support` | 0.8086 | 0.8519 | 0.8297 | 7,851 | 1,858 | 1,460 | 0.8624 | 0.7503 | 0.8025 |
| `readlik-nomismap` | 0.8166 | 0.8547 | 0.8352 | 7,928 | 1,781 | 1,429 | 0.8724 | 0.6746 | 0.7608 |
| `readlik-nolink` | 0.8094 | 0.8515 | 0.8299 | 7,858 | 1,851 | 1,464 | 0.8666 | 0.7494 | 0.8038 |
| `readlik` | 0.8175 | 0.8577 | **0.8371** | 7,937 | 1,772 | 1,392 | 0.8728 | 0.7758 | 0.8215 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8251 | 0.8047 | 0.8148 | 8,222 | 1,743 | 2,148 | 0.8750 | 0.6764 | 0.7630 |
| `poisson-z` | 0.8259 | 0.8040 | 0.8148 | 8,230 | 1,735 | 2,159 | 0.8737 | 0.7015 | 0.7782 |
| `readlik-support` | 0.8608 | 0.8901 | 0.8752 | 8,578 | 1,387 | 1,133 | 0.8742 | 0.8375 | 0.8554 |
| `readlik-nomismap` | 0.8633 | 0.8997 | 0.8811 | 8,603 | 1,362 | 1,021 | 0.8796 | 0.8250 | 0.8514 |
| `readlik-nolink` | 0.8607 | 0.8897 | 0.8750 | 8,577 | 1,388 | 1,136 | 0.8737 | 0.8370 | 0.8549 |
| `readlik` | 0.8647 | 0.9010 | **0.8825** | 8,617 | 1,348 | 1,005 | 0.8800 | 0.8556 | 0.8676 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.8746 | — | 0 | 0 | 116 | — | 0.5789 | — |
| `poisson-z` | — | 0.8772 | — | 0 | 0 | 116 | — | 0.5131 | — |
| `readlik-support` | — | 0.8803 | — | 0 | 0 | 116 | — | 0.7585 | — |
| `readlik-nomismap` | — | 0.8812 | — | 0 | 0 | 117 | — | 0.6219 | — |
| `readlik-nolink` | — | 0.8851 | — | 0 | 0 | 111 | — | 0.7600 | — |
| `readlik` | — | 0.8897 | — | 0 | 0 | 106 | — | 0.7618 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7772 | 0.8261 | 0.8009 | 15,291 | 4,383 | 3,726 | 0.8227 | 0.7065 | 0.7602 |
| `poisson-z` | 0.7783 | 0.8275 | 0.8022 | 15,313 | 4,361 | 3,701 | 0.8243 | 0.7181 | 0.7675 |
| `readlik-support` | 0.8351 | 0.8718 | 0.8530 | 16,429 | 3,245 | 2,709 | 0.8684 | 0.7910 | 0.8279 |
| `readlik-nomismap` | 0.8402 | 0.8778 | 0.8586 | 16,531 | 3,143 | 2,567 | 0.8760 | 0.7375 | 0.8008 |
| `readlik-nolink` | 0.8354 | 0.8717 | 0.8531 | 16,435 | 3,239 | 2,711 | 0.8702 | 0.7903 | 0.8283 |
| `readlik` | 0.8414 | 0.8802 | **0.8604** | 16,554 | 3,120 | 2,503 | 0.8765 | 0.8122 | 0.8431 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9185 | 0.9531 | 0.9355 | 86,975 | 7,716 | 4,314 | 0.9040 | 0.8619 | 0.8825 |
| `poisson-z` | 0.9193 | 0.9532 | 0.9359 | 87,047 | 7,644 | 4,306 | 0.9046 | 0.8695 | 0.8867 |
| `readlik-support` | 0.9319 | 0.9660 | 0.9487 | 88,244 | 6,447 | 3,112 | 0.9279 | 0.9137 | 0.9207 |
| `readlik-nomismap` | 0.9335 | 0.9681 | 0.9505 | 88,390 | 6,301 | 2,920 | 0.9319 | 0.8811 | 0.9058 |
| `readlik-nolink` | 0.9322 | 0.9659 | 0.9488 | 88,274 | 6,417 | 3,119 | 0.9290 | 0.9132 | 0.9211 |
| `readlik` | 0.9330 | 0.9692 | **0.9507** | 88,343 | 6,348 | 2,815 | 0.9319 | 0.9269 | 0.9294 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (58.9 Mb vs 59.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

That is exactly where the gap lives. 246 `readlik` calls carry a >=200 bp insertion allele; they contribute **27,951 FP bases and zero TP bases**, which is the whole of the precision difference. The Poisson caller scores better there because it does not emit them — at the two largest sites it emits nothing at all.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7637 | 0.8700 | **0.8134** |
| `sm50-poisson-z` | Deletion | 0.8628 | 0.8094 | **0.8353** |
| `sm50-poisson-z` | ALL | 0.8993 | 0.9385 | **0.9184** |
| `sm50-readlik` | Insertion | 0.8643 | 0.8794 | **0.8718** |
| `sm50-readlik` | Deletion | 0.8724 | 0.9030 | **0.8874** |
| `sm50-readlik` | ALL | 0.9275 | 0.9688 | **0.9477** |

The insertion BASEPAIR precision gap collapses from **-0.006 to -0.009**, and insertion BASEPAIR F1 goes from 0.8134 for `poisson-z` against 0.8718 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it. On chr20, of the 246, only **35 are confirmed true**, **73 are confirmed false**, and **138 fall outside the SV confident region** and cannot be judged at all. See *Known bad output* for the worst of the unjudged ones.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.4889 | 0.5021 | 0.4954 | 374 | 362 | 391 |
| `poisson-z` | 0.4902 | 0.4959 | 0.4930 | 375 | 372 | 390 |
| `readlik-support` | 0.4889 | 0.5188 | **0.5034** | 374 | 333 | 391 |
| `readlik-nomismap` | 0.4876 | 0.4771 | 0.4823 | 373 | 399 | 392 |
| `readlik-nolink` | 0.4915 | 0.5084 | 0.4998 | 376 | 351 | 389 |
| `readlik` | 0.4863 | 0.5180 | 0.5016 | 372 | 335 | 393 |

## Structural variants — aardvark (secondary)

Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

### SV insertion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.3877 | 321 | 507 | 713 | 336 | 377 | 0.4712 | 0.4254 |
| `poisson-z` | 0.4263 | 353 | 475 | 723 | 375 | 348 | 0.5187 | 0.4680 |
| `readlik-nomismap` | 0.4589 | 380 | 448 | 713 | 327 | 386 | 0.4586 | 0.4588 |
| `readlik` | 0.4577 | 379 | 449 | 656 | 307 | 349 | 0.4680 | 0.4628 |

### SV deletion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.5434 | 463 | 389 | 713 | 336 | 377 | 0.4712 | 0.5048 |
| `poisson-z` | 0.5763 | 491 | 361 | 723 | 375 | 348 | 0.5187 | 0.5460 |
| `readlik-nomismap` | 0.5023 | 428 | 424 | 713 | 327 | 386 | 0.4586 | 0.4795 |
| `readlik` | 0.5094 | 434 | 418 | 656 | 307 | 349 | 0.4680 | 0.4878 |

### SV (joint)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 784 | 896 | 713 | 336 | 377 | 0.4712 | 0.4689 |
| `poisson-z` | 0.5024 | 844 | 836 | 723 | 375 | 348 | 0.5187 | 0.5104 |
| `readlik-nomismap` | 0.4810 | 808 | 872 | 713 | 327 | 386 | 0.4586 | 0.4695 |
| `readlik` | 0.4839 | 813 | 867 | 656 | 307 | 349 | 0.4680 | 0.4758 |

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) looks inert on this graph, because it binds only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless.

The two graphs are put side by side in [tier2-chr20-graph-comparison.md](tier2-chr20-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| floor 1e-8, cap 0.1 (original defaults) | 0.9370 | 0.9759 | 0.7783 | 0.8231 | 0.8686 |
| **floor 0.02, cap 0.7 (current defaults)** | 0.9507 | 0.9759 | 0.8371 | 0.8825 | 0.9294 |
| floor 0.05, cap 0.1 | 0.9495 | 0.9745 | 0.8346 | 0.8840 | 0.8954 |
| cap 0.2, floor 1e-8 | 0.9370 | 0.9759 | 0.7783 | 0.8233 | 0.8674 |
| cap 0.4, floor 1e-8 | 0.9370 | 0.9758 | 0.7785 | 0.8234 | 0.8710 |

Sweep rows other than the current one are historical: they were produced at the defaults in force at the time and are kept because the comparison between them is the result. The full grids are in plan §9.20-§9.21.

Raising the floor off 1e-8 changed **1,493 genotypes (1.41%)** on chr20, of which **94% were heterozygous → homozygous** (1/0→1/1: 614, 0/1→1/1: 606, 1/2→1/1: 184), and dropped 1,251 spurious non-reference calls. The failure it corrects is spurious heterozygosity: a few locally misaligned reads, each able to veto the homozygous hypothesis almost without bound, conjuring a second allele that is not there.

The floor was later re-swept at the corrected cap, on both graphs and both benchmarks, and settled at **0.02**. 0.05 wins on small-variant `GT` but costs about 0.01 of SV F1 — which the first sweep never saw, because it was scored on one benchmark only. Plan §9.21 records that as a process rule: a sweep that sets a default has to be scored on every benchmark the project runs.

## Known bad output

Neither benchmark scores these, so they appear in no metric on this page. They are recorded because they are plainly wrong and would mislead anyone reading the VCF.

`readlik` emits a small number of enormous homozygous insertions in and around the chr20 pericentromere, at depths that are physically impossible:

| position | called insertion | GT | DP | GQ |
|---|---|---|---|---|
| chr20:25,849,044 | 61,958 bp | 1/1 | 7,873 | 256 |
| chr20:32,179,077 | 57,716 bp | 1/1 | 5,337 | 256 |
| chr20:1,629,728 | 33,050 bp | 1/1 | 291 | 256 |
| chr20:25,873,453 | 28,685 bp | 1/2 | 5,498 | 256 |
| chr20:25,792,993 | 23,450 bp | 1/1 | 932 | 256 |

Chromosome-median DP is **29**, and the Poisson caller's expected depth (`XD`) never exceeds **167** anywhere on chr20. Median DP rises monotonically with called insertion length — 28 for 1 bp, 28 for 2–15 bp, 35 for 50–199 bp, **330 for >=1 kb** — so these are collapsed-repeat pile-ups, not haplotypes.

The read-likelihood model cannot reject them, and the reason is structural rather than a tuning failure: it computes P(reads | genotype) **conditioned on the reads it is given**, and never asks whether that many reads should be there. The Poisson caller gets this for free, because an observed-vs-expected depth term is the whole of its model. A depth-plausibility guard is the obvious remedy, and the expected depth is already reachable — the read-likelihood caller subclasses `SupportBasedSnarlCaller` and holds a `TraversalSupportFinder` for allele enumeration.

The same blindness has a second consequence, found later and now corrected. Because the model only weighs reads it can see, it had no way to know that a heterozygous deletion produces *no* reads over the deleted interval, and its flat `1/ploidy` mixture asserted that both haplotypes contributed equally everywhere. That cost it 94% of heterozygous deletions above 1 kb and mis-genotyped two thirds of heterozygous insertions above 1 kb. Weighting each haplotype by the reads it is *expected* to contribute at the site is now the default and fixes both, without moving small variants at all — see [tier2-sv-errors.md](tier2-sv-errors.md). It did not remove the need for a depth term: it corrects the *relative* weight between a genotype's haplotypes, while the pile-ups above are a statement about *absolute* depth. That term is now also the default, at `--depth-term 0.1`, and the read arms in the tables on this page carry it — see [tier2-depth-term.md](tier2-depth-term.md). It does not resolve the pile-ups either: it detects them emphatically and still cannot outvote the read evidence at them, which is what the `DR` field and `--depth-quality` are for ([tier2-quality-signals.md](tier2-quality-signals.md)).

Filtering on depth is **not** that remedy, and that has now been tested properly rather than by two spot checks. Sweeping a two-sided cut on DP over a rolling local median, across both chromosomes and both graphs, against the one test a hard filter has to pass — beat lowering the GQ threshold to the same recall:

- a **minimum** fails in all eight dataset-by-benchmark cells. Few reads already means a small likelihood gap, so low depth depresses GQ on its own and a separate cut adds nothing;
- a **maximum** passes in exactly one configuration — 5x the local median, structural calls, 34-haplotype graph, worth about +0.025 precision — and is dominated everywhere else. The two original spot checks (DP 200 moving insertion BASEPAIR precision by 0.0001; DP 58 helping by +0.087 but costing SV insertion recall 0.4976 to 0.4167) were both right and both too narrow to conclude from.

What shipped instead attacks the same blindness from the other side: **GQ is now scaled by the fraction of reads the called genotype explains**, so a pile-up the call does not account for can no longer carry a saturated quality. The giants remain output that no metric charges for — they should be fixed because they are wrong, not because they cost a score — but they no longer look confident. See [tier2-quality-signals.md](tier2-quality-signals.md).

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. `vg call` emits `AD` (per-allele read support, ties split fractionally), `BL` (mean absolute fit), `GQI` (the raw likelihood-ratio quality) and `GQ` (that ratio scaled by the fraction of reads the called genotype explains). The scaling rescales a quality and does not change a genotype, so **the numbers on this page are unaffected by it**; what it changes is how the calls rank. See [tier2-quality-signals.md](tier2-quality-signals.md).

## The genotype mixture

The read-likelihood arms on this page use the **length-weighted mixture**, which became the default after it was found that the flat `1/ploidy` weight breaks heterozygotes whose alleles differ in length. Unlike the `GQ` scaling above, this *does* change genotypes, so these numbers are not comparable with runs made before it. `--flat-mixture` restores the previous model exactly. Derivation and measurements: [tier2-sv-errors.md](tier2-sv-errors.md).

## Raw aardvark summary rows

<details><summary><code>poisson</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson | GT | ALL | ALL | ALL | 94691 | 86975 | 7716 | 91990 | 87676 | 4314 | 0.9185139031164524 | 0.9531035982171975 | 0.9354891215208296 | 1342 | 724 |
| poisson | GT | ALL | ALL | Snv | 75017 | 71684 | 3333 | 70558 | 69970 | 588 | 0.9555700707839556 | 0.991666430454378 | 0.97328368746245 | 199 | 214 |
| poisson | GT | ALL | ALL | Insertion | 9709 | 7069 | 2640 | 9507 | 8045 | 1462 | 0.7280873416417757 | 0.8462185757862628 | 0.7827208504668184 | 713 | 174 |
| poisson | GT | ALL | ALL | Deletion | 9965 | 8222 | 1743 | 11000 | 8852 | 2148 | 0.8250878073256397 | 0.8047272727272727 | 0.8147803626017884 | 430 | 320 |
| poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 925 | 809 | 116 |  | 0.8745945945945945 |  | 0 | 16 |
| poisson | GT | ALL | ALL | JointIndel | 19674 | 15291 | 4383 | 21432 | 17706 | 3726 | 0.7772186642268984 | 0.8261478163493841 | 0.8009366667641129 | 1143 | 510 |
| poisson | BASEPAIR | ALL | ALL | ALL | 390714 | 353211 | 37503 | 409806 | 353211 | 56595 | 0.9040141894070854 | 0.8618980688423302 | 0.8824539049617747 |  |  |
| poisson | BASEPAIR | ALL | ALL | Snv | 200440 | 193613 | 6827 | 190468 | 188743 | 1725 | 0.9659399321492717 | 0.9909433605645043 | 0.9782819098424921 |  |  |
| poisson | BASEPAIR | ALL | ALL | Insertion | 95512 | 73393 | 22119 | 96422 | 73404 | 23018 | 0.7684165340480777 | 0.7612785463898281 | 0.7648308862895322 |  |  |
| poisson | BASEPAIR | ALL | ALL | Deletion | 99134 | 86741 | 12393 | 126980 | 85892 | 41088 | 0.8749873908043658 | 0.6764214836982202 | 0.7629971425745341 |  |  |
| poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11504 | 6660 | 4844 |  | 0.5789290681502086 |  |  |  |
| poisson | BASEPAIR | ALL | ALL | JointIndel | 194646 | 160134 | 34512 | 234906 | 165956 | 68950 | 0.8226935051323941 | 0.7064783360152571 | 0.7601698159969886 |  |  |

</details>

<details><summary><code>poisson-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-z | GT | ALL | ALL | ALL | 94691 | 87047 | 7644 | 92035 | 87729 | 4306 | 0.919274271050047 | 0.9532134514043571 | 0.9359362843204727 | 1349 | 719 |
| poisson-z | GT | ALL | ALL | Snv | 75017 | 71734 | 3283 | 70584 | 69979 | 605 | 0.9562365863737553 | 0.9914286523858098 | 0.9735146793443107 | 200 | 219 |
| poisson-z | GT | ALL | ALL | Insertion | 9709 | 7083 | 2626 | 9488 | 8062 | 1426 | 0.7295293027088269 | 0.8497048903878583 | 0.785044572746258 | 717 | 172 |
| poisson-z | GT | ALL | ALL | Deletion | 9965 | 8230 | 1735 | 11018 | 8859 | 2159 | 0.8258906171600602 | 0.8040479215828644 | 0.8148229131320283 | 432 | 313 |
| poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 945 | 829 | 116 |  | 0.8772486772486773 |  | 0 | 15 |
| poisson-z | GT | ALL | ALL | JointIndel | 19674 | 15313 | 4361 | 21451 | 17750 | 3701 | 0.7783368913286571 | 0.827467250944012 | 0.8021504874990804 | 1149 | 500 |
| poisson-z | BASEPAIR | ALL | ALL | ALL | 390902 | 353624 | 37278 | 406720 | 353624 | 53096 | 0.9046359445589943 | 0.8694531864673486 | 0.8866957029770995 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Snv | 200440 | 193709 | 6731 | 190578 | 188793 | 1785 | 0.9664188784673717 | 0.9906337562572805 | 0.9783765103780948 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Insertion | 95512 | 73826 | 21686 | 95356 | 73369 | 21987 | 0.7729499958120446 | 0.7694219556189438 | 0.7711819406746623 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Deletion | 99134 | 86612 | 12522 | 122114 | 85666 | 36448 | 0.8736861218149172 | 0.7015248046906989 | 0.7781973520547425 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 14032 | 7200 | 6832 |  | 0.5131128848346637 |  |  |  |
| poisson-z | BASEPAIR | ALL | ALL | JointIndel | 194646 | 160438 | 34208 | 231502 | 166235 | 65267 | 0.8242553147765688 | 0.7180715501377958 | 0.7675082436222629 |  |  |

</details>

<details><summary><code>readlik-support</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-support | GT | ALL | ALL | ALL | 94691 | 88244 | 6447 | 91542 | 88430 | 3112 | 0.9319153879460561 | 0.9660046754495204 | 0.9486538861585461 | 559 | 474 |
| readlik-support | GT | ALL | ALL | Snv | 75017 | 71815 | 3202 | 70407 | 70004 | 403 | 0.9573163416292307 | 0.9942761373158919 | 0.9754462620792435 | 82 | 175 |
| readlik-support | GT | ALL | ALL | Insertion | 9709 | 7851 | 1858 | 9860 | 8400 | 1460 | 0.8086311669584921 | 0.8519269776876268 | 0.8297146454667229 | 259 | 159 |
| readlik-support | GT | ALL | ALL | Deletion | 9965 | 8578 | 1387 | 10306 | 9173 | 1133 | 0.8608128449573508 | 0.890064040364836 | 0.8751940986869569 | 218 | 137 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 969 | 853 | 116 |  | 0.8802889576883385 |  | 0 | 3 |
| readlik-support | GT | ALL | ALL | JointIndel | 19674 | 16429 | 3245 | 21135 | 18426 | 2709 | 0.8350615024905967 | 0.8718239886444287 | 0.8530468548076375 | 477 | 299 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 390700 | 362515 | 28185 | 396764 | 362515 | 34249 | 0.9278602508318403 | 0.9136791644403222 | 0.9207151057064196 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 200440 | 193883 | 6557 | 190322 | 188924 | 1398 | 0.9672869686689284 | 0.9926545538613507 | 0.9798065945359036 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 95512 | 82368 | 13144 | 109186 | 81926 | 27260 | 0.8623837842365357 | 0.7503342919421904 | 0.8024665137514548 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 99134 | 86658 | 12476 | 102332 | 85699 | 16633 | 0.8741501402142554 | 0.8374604229371067 | 0.8554120451167406 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 9420 | 7145 | 2275 |  | 0.7584925690021231 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 194646 | 169026 | 25620 | 220938 | 174770 | 46168 | 0.8683764372244999 | 0.7910363993518543 | 0.8279041297538535 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 88390 | 6301 | 91540 | 88620 | 2920 | 0.9334572451447339 | 0.9681013764474546 | 0.9504637234089658 | 560 | 370 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 71859 | 3158 | 70542 | 70189 | 353 | 0.9579028753482544 | 0.9949958889739446 | 0.9760971130919158 | 175 | 83 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 7928 | 1781 | 9832 | 8403 | 1429 | 0.8165619528272736 | 0.8546582587469488 | 0.835175893552639 | 194 | 151 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 8603 | 1362 | 10181 | 9160 | 1021 | 0.8633216256899147 | 0.899715155682153 | 0.8811427635183848 | 191 | 129 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 985 | 868 | 117 |  | 0.8812182741116751 |  | 0 | 7 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 16531 | 3143 | 20998 | 18431 | 2567 | 0.8402460099623869 | 0.8777502619297076 | 0.8585887727423427 | 385 | 287 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390700 | 364103 | 26597 | 413254 | 364103 | 49151 | 0.931924750447914 | 0.8810634621806444 | 0.9057806789940718 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 193954 | 6486 | 190202 | 188910 | 1292 | 0.9676411893833566 | 0.993207221795775 | 0.9802575374245037 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 83321 | 12191 | 122216 | 82442 | 39774 | 0.8723615880727029 | 0.6745597957714211 | 0.7608144290136096 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 87198 | 11936 | 104266 | 86015 | 18251 | 0.8795973127282264 | 0.8249573206989814 | 0.8514015663356435 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 12428 | 7729 | 4699 |  | 0.6219021564209849 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 170519 | 24127 | 238910 | 176186 | 62724 | 0.8760467720888176 | 0.737457620024277 | 0.800800259029337 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 94691 | 88274 | 6417 | 91528 | 88409 | 3119 | 0.9322322079183871 | 0.9659229962415873 | 0.9487786093486892 | 560 | 481 |
| readlik-nolink | GT | ALL | ALL | Snv | 75017 | 71839 | 3178 | 70404 | 69996 | 408 | 0.9576362691123346 | 0.9942048747230271 | 0.9755780074316994 | 77 | 178 |
| readlik-nolink | GT | ALL | ALL | Insertion | 9709 | 7858 | 1851 | 9857 | 8393 | 1464 | 0.8093521474920177 | 0.8514761083493964 | 0.8298799281706507 | 260 | 159 |
| readlik-nolink | GT | ALL | ALL | Deletion | 9965 | 8577 | 1388 | 10301 | 9165 | 1136 | 0.8607124937280481 | 0.8897194447141055 | 0.8749756276268025 | 223 | 135 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 966 | 855 | 111 |  | 0.8850931677018633 |  | 0 | 9 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 19674 | 16435 | 3239 | 21124 | 18413 | 2711 | 0.8353664735183491 | 0.8716625639083507 | 0.8531286418041952 | 483 | 303 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 390700 | 362971 | 27729 | 397452 | 362971 | 34481 | 0.9290273867417456 | 0.9132448698207583 | 0.921068524853074 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 200440 | 193939 | 6501 | 190316 | 188917 | 1399 | 0.9675663540211534 | 0.9926490678660753 | 0.9799472330372575 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 95512 | 82771 | 12741 | 109352 | 81949 | 27403 | 0.8666031493424909 | 0.7494055892896335 | 0.8037546187565543 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 99134 | 86615 | 12519 | 102318 | 85637 | 16681 | 0.8737163838844393 | 0.8369690572528783 | 0.8549480349116425 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 10012 | 7609 | 2403 |  | 0.7599880143827407 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 194646 | 169386 | 25260 | 221682 | 175195 | 46487 | 0.8702259486452328 | 0.7902987161790312 | 0.8283387348212518 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 94691 | 88343 | 6348 | 91284 | 88469 | 2815 | 0.9329608938547486 | 0.9691621751895184 | 0.9507170424144441 | 480 | 398 |
| readlik | GT | ALL | ALL | Snv | 75017 | 71789 | 3228 | 70387 | 70075 | 312 | 0.9569697535225349 | 0.9955673632915168 | 0.9758870610549271 | 121 | 94 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 7937 | 1772 | 9781 | 8389 | 1392 | 0.8174889277989494 | 0.8576832634699928 | 0.8371038811408037 | 180 | 157 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 8617 | 1348 | 10155 | 9150 | 1005 | 0.8647265429001505 | 0.9010339734121122 | 0.8825069828738372 | 179 | 140 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 961 | 855 | 106 |  | 0.8896982310093653 |  | 0 | 7 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 16554 | 3120 | 20897 | 18394 | 2503 | 0.841415065568771 | 0.8802220414413552 | 0.8603811845118405 | 359 | 304 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390700 | 364075 | 26625 | 392776 | 364075 | 28701 | 0.9318530842078321 | 0.9269278163635253 | 0.9293839249702608 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 193823 | 6617 | 190002 | 188790 | 1212 | 0.9669876272201158 | 0.9936211197776865 | 0.9801234748553291 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 83367 | 12145 | 106288 | 82458 | 23830 | 0.8728432029483206 | 0.7757978323046816 | 0.8214642852017403 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 87239 | 11895 | 100796 | 86240 | 14556 | 0.880010894345028 | 0.8555895075201397 | 0.8676283859993165 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 9984 | 7606 | 2378 |  | 0.7618189102564102 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 170606 | 24040 | 217068 | 176304 | 40764 | 0.8764937373488282 | 0.8122063132290342 | 0.843126340568088 |  |  |

</details>

<details><summary><code>poisson</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-poisson | GT | ALL | ALL | ALL | 100207 | 88821 | 11386 | 95410 | 89562 | 5848 | 0.8863752033291088 | 0.938706634524683 | 0.9117906570388599 | 1462 | 884 |
| sv-poisson | GT | ALL | ALL | Snv | 78483 | 72822 | 5661 | 71941 | 70903 | 1038 | 0.9278697297503918 | 0.9855715099873508 | 0.955850591729682 | 270 | 279 |
| sv-poisson | GT | ALL | ALL | Insertion | 9857 | 7066 | 2791 | 10249 | 8423 | 1826 | 0.7168509688546211 | 0.8218362767099229 | 0.7657620259056992 | 701 | 208 |
| sv-poisson | GT | ALL | ALL | Deletion | 10187 | 8149 | 2038 | 11988 | 9326 | 2662 | 0.7999411014037499 | 0.7779446112779446 | 0.7887895354843473 | 430 | 373 |
| sv-poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1232 | 910 | 322 |  | 0.7386363636363636 |  | 0 | 24 |
| sv-poisson | GT | ALL | ALL | SvInsertion | 828 | 321 | 507 | 0 | 0 | 0 | 0.38768115942028986 |  |  | 41 | 0 |
| sv-poisson | GT | ALL | ALL | SvDeletion | 852 | 463 | 389 | 0 | 0 | 0 | 0.5434272300469484 |  |  | 20 | 0 |
| sv-poisson | GT | ALL | ALL | JointIndel | 20044 | 15215 | 4829 | 23469 | 18659 | 4810 | 0.7590800239473159 | 0.7950487877625805 | 0.7766481752437369 | 1131 | 605 |
| sv-poisson | GT | ALL | ALL | JointStructuralVariant | 1680 | 784 | 896 | 0 | 0 | 0 | 0.4666666666666667 |  |  | 61 | 0 |
| sv-poisson | BASEPAIR | ALL | ALL | ALL | 969654 | 766195 | 203459 | 912096 | 766195 | 145901 | 0.7901736083180186 | 0.8400376714731783 | 0.8143430317523582 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Snv | 208362 | 197801 | 10561 | 193834 | 190888 | 2946 | 0.9493141743696067 | 0.9848014280260429 | 0.9667322402100234 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Insertion | 79376 | 60349 | 19027 | 236752 | 190208 | 46544 | 0.760292783712961 | 0.8034060958302358 | 0.7812550933453016 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Deletion | 82226 | 71194 | 11032 | 372090 | 272443 | 99647 | 0.8658331914479604 | 0.7321965115966568 | 0.7934271073872482 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 170138 | 128832 | 41306 |  | 0.7572206091525703 |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 179326 | 88496 | 0 | 0 | 0 | 0.6695715811247769 |  |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 268286 | 106776 | 0 | 0 | 0 | 0.7153110685699964 |  |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | JointIndel | 161602 | 131543 | 30059 | 778980 | 591483 | 187497 | 0.8139936386925904 | 0.7593044750828006 | 0.7856985362615818 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 447612 | 195272 | 0 | 0 | 0 | 0.696256245294641 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | ALL | 1067628 | 864169 | 203459 | 1214080 | 1068179 | 145901 | 0.8094289396681241 | 0.8798258763837639 | 0.8431605693193573 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Snv | 208430 | 197869 | 10561 | 194538 | 191592 | 2946 | 0.9493307105503047 | 0.9848564290781235 | 0.9667673147556781 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Insertion | 104216 | 85189 | 19027 | 266290 | 219746 | 46544 | 0.8174272664466109 | 0.825213113522851 | 0.8213017381630682 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Deletion | 108210 | 97178 | 11032 | 407276 | 307629 | 99647 | 0.8980500877922558 | 0.7553329928598788 | 0.8205319970767486 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 345976 | 304670 | 41306 |  | 0.8806102157375078 |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 181226 | 88496 | 0 | 0 | 0 | 0.6718992147470358 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 270274 | 106776 | 0 | 0 | 0 | 0.7168120938867524 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | JointIndel | 212426 | 182367 | 30059 | 1019542 | 832045 | 187497 | 0.85849660587687 | 0.816096835637963 | 0.8367599515118957 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 451500 | 195272 | 0 | 0 | 0 | 0.6980821680592233 |  |  |  |  |

</details>

<details><summary><code>poisson-z</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-poisson-z | GT | ALL | ALL | ALL | 100207 | 89383 | 10824 | 95469 | 89701 | 5768 | 0.8919835939605018 | 0.9395824822717321 | 0.915164536223702 | 1497 | 896 |
| sv-poisson-z | GT | ALL | ALL | Snv | 78483 | 73280 | 5203 | 71987 | 70921 | 1066 | 0.9337053884280673 | 0.985191770736383 | 0.9587578579480464 | 297 | 290 |
| sv-poisson-z | GT | ALL | ALL | Insertion | 9857 | 7093 | 2764 | 10233 | 8476 | 1757 | 0.7195901389875216 | 0.8283005961106225 | 0.7701279264274401 | 702 | 209 |
| sv-poisson-z | GT | ALL | ALL | Deletion | 10187 | 8166 | 2021 | 12016 | 9354 | 2662 | 0.80160989496417 | 0.778462050599201 | 0.7898664163572933 | 433 | 370 |
| sv-poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1233 | 950 | 283 |  | 0.7704785077047851 |  | 0 | 27 |
| sv-poisson-z | GT | ALL | ALL | SvInsertion | 828 | 353 | 475 | 0 | 0 | 0 | 0.42632850241545894 |  |  | 43 | 0 |
| sv-poisson-z | GT | ALL | ALL | SvDeletion | 852 | 491 | 361 | 0 | 0 | 0 | 0.5762910798122066 |  |  | 22 | 0 |
| sv-poisson-z | GT | ALL | ALL | JointIndel | 20044 | 15259 | 4785 | 23482 | 18780 | 4702 | 0.7612751945719417 | 0.7997615194617154 | 0.7800439296090023 | 1135 | 606 |
| sv-poisson-z | GT | ALL | ALL | JointStructuralVariant | 1680 | 844 | 836 | 0 | 0 | 0 | 0.5023809523809524 |  |  | 65 | 0 |
| sv-poisson-z | BASEPAIR | ALL | ALL | ALL | 969654 | 778806 | 190848 | 916550 | 778806 | 137744 | 0.8031792783817733 | 0.8497146909606678 | 0.8257919079802608 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Snv | 208362 | 198020 | 10342 | 194030 | 190967 | 3063 | 0.9503652297443872 | 0.9842137813740143 | 0.966993388305566 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Insertion | 79376 | 60700 | 18676 | 236464 | 192598 | 43866 | 0.764714775246926 | 0.8144918465390081 | 0.7888188167070147 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Deletion | 82226 | 71138 | 11088 | 373446 | 275117 | 98329 | 0.8651521416583562 | 0.736698210718551 | 0.7957747536320006 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 171994 | 133886 | 38108 |  | 0.7784341314231892 |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 189207 | 78615 | 0 | 0 | 0 | 0.7064654882720613 |  |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 268826 | 106236 | 0 | 0 | 0 | 0.7167508305293525 |  |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | JointIndel | 161602 | 131838 | 29764 | 781904 | 601601 | 180303 | 0.8158191111496145 | 0.769405195522724 | 0.7919326748690676 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 458033 | 184851 | 0 | 0 | 0 | 0.712466012531032 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | ALL | 1067628 | 876780 | 190848 | 1195494 | 1057750 | 137744 | 0.8212411064528095 | 0.8847806848047752 | 0.8518276522382817 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Snv | 208430 | 198088 | 10342 | 194864 | 191801 | 3063 | 0.9503814230197188 | 0.984281344938008 | 0.9670343800965592 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Insertion | 104216 | 85540 | 18676 | 265556 | 221690 | 43866 | 0.8207952713594842 | 0.8348145024025064 | 0.8277455314573465 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Deletion | 108210 | 97122 | 11088 | 408540 | 310211 | 98329 | 0.8975325755475464 | 0.7593161012385569 | 0.8226592392509066 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 326534 | 288426 | 38108 |  | 0.8832954608095941 |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 191107 | 78615 | 0 | 0 | 0 | 0.7085332305114155 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 270814 | 106236 | 0 | 0 | 0 | 0.7182442646863811 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | JointIndel | 212426 | 182662 | 29764 | 1000630 | 820327 | 180303 | 0.8598853247719206 | 0.8198105193727951 | 0.8393698622994008 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 461921 | 184851 | 0 | 0 | 0 | 0.7141944920311949 |  |  |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik-nomismap | GT | ALL | ALL | ALL | 100207 | 90210 | 9997 | 95243 | 90585 | 4658 | 0.9002365104234236 | 0.9510935186837878 | 0.9249664801895188 | 772 | 606 |
| sv-readlik-nomismap | GT | ALL | ALL | Snv | 78483 | 73057 | 5426 | 71942 | 71123 | 819 | 0.9308640087662297 | 0.9886158294181424 | 0.9588711230980156 | 214 | 209 |
| sv-readlik-nomismap | GT | ALL | ALL | Insertion | 9857 | 7855 | 2002 | 10820 | 8829 | 1991 | 0.7968956071827128 | 0.8159889094269871 | 0.8063292451948483 | 256 | 219 |
| sv-readlik-nomismap | GT | ALL | ALL | Deletion | 10187 | 8490 | 1697 | 11150 | 9655 | 1495 | 0.8334151369392363 | 0.8659192825112108 | 0.8493563469935353 | 231 | 166 |
| sv-readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1331 | 978 | 353 |  | 0.734785875281743 |  | 0 | 12 |
| sv-readlik-nomismap | GT | ALL | ALL | SvInsertion | 828 | 380 | 448 | 0 | 0 | 0 | 0.45893719806763283 |  |  | 31 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | SvDeletion | 852 | 428 | 424 | 0 | 0 | 0 | 0.5023474178403756 |  |  | 40 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | JointIndel | 20044 | 16345 | 3699 | 23301 | 19462 | 3839 | 0.8154559968070245 | 0.8352431226127635 | 0.8252309643998686 | 487 | 397 |
| sv-readlik-nomismap | GT | ALL | ALL | JointStructuralVariant | 1680 | 808 | 872 | 0 | 0 | 0 | 0.48095238095238096 |  |  | 71 | 0 |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 969654 | 675125 | 294529 | 861084 | 675125 | 185959 | 0.6962535089836168 | 0.7840408136720691 | 0.7375440942395908 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 208362 | 197911 | 10451 | 193872 | 191220 | 2652 | 0.9498421017268024 | 0.9863208715028472 | 0.9677378428558094 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 79376 | 67875 | 11501 | 332120 | 216021 | 116099 | 0.8551073372303971 | 0.650430567264844 | 0.7388561241353395 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 82226 | 70225 | 12001 | 217466 | 174033 | 43433 | 0.8540485977671296 | 0.8002768248829703 | 0.8262888193085464 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 180808 | 104286 | 76522 |  | 0.5767775762134418 |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 207928 | 59894 | 0 | 0 | 0 | 0.7763663926040429 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 140426 | 234636 | 0 | 0 | 0 | 0.37440743130469095 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 161602 | 138100 | 23502 | 730394 | 494340 | 236054 | 0.8545686315763419 | 0.6768127887140365 | 0.7553741622058465 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 348354 | 294530 | 0 | 0 | 0 | 0.5418613622364221 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | ALL | 1067628 | 773099 | 294529 | 1207748 | 1021789 | 185959 | 0.724127692417209 | 0.8460283105416031 | 0.78034606380216 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Snv | 208430 | 197979 | 10451 | 194870 | 192218 | 2652 | 0.9498584656719282 | 0.9863909272848566 | 0.9677800557631119 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Insertion | 104216 | 92715 | 11501 | 364760 | 248661 | 116099 | 0.8896426652337452 | 0.6817112621998026 | 0.7719195702954588 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Deletion | 108210 | 96209 | 12001 | 254424 | 210991 | 43433 | 0.889095277700767 | 0.829288903562557 | 0.8581513447883985 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 393694 | 317172 | 76522 |  | 0.8056307690744589 |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 209828 | 59894 | 0 | 0 | 0 | 0.7779417325987498 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 142414 | 234636 | 0 | 0 | 0 | 0.37770587455244664 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointIndel | 212426 | 188924 | 23502 | 1012878 | 776824 | 236054 | 0.8893638255204165 | 0.7669472532723586 | 0.8236316859509688 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 352242 | 294530 | 0 | 0 | 0 | 0.544615413159506 |  |  |  |  |

</details>

<details><summary><code>readlik</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik | GT | ALL | ALL | ALL | 100207 | 90120 | 10087 | 94794 | 90306 | 4488 | 0.8993383695749798 | 0.952655231343756 | 0.9252293335125434 | 681 | 717 |
| sv-readlik | GT | ALL | ALL | Snv | 78483 | 72957 | 5526 | 71706 | 70943 | 763 | 0.9295898474828944 | 0.9893593283686163 | 0.9585437683682849 | 151 | 272 |
| sv-readlik | GT | ALL | ALL | Insertion | 9857 | 7857 | 2000 | 10721 | 8795 | 1926 | 0.7970985086740388 | 0.8203525790504617 | 0.8085583821492389 | 250 | 233 |
| sv-readlik | GT | ALL | ALL | Deletion | 10187 | 8493 | 1694 | 11081 | 9610 | 1471 | 0.8337096299204869 | 0.8672502481725476 | 0.8501492513309459 | 228 | 195 |
| sv-readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1286 | 958 | 328 |  | 0.744945567651633 |  | 0 | 17 |
| sv-readlik | GT | ALL | ALL | SvInsertion | 828 | 379 | 449 | 0 | 0 | 0 | 0.4577294685990338 |  |  | 20 | 0 |
| sv-readlik | GT | ALL | ALL | SvDeletion | 852 | 434 | 418 | 0 | 0 | 0 | 0.5093896713615024 |  |  | 32 | 0 |
| sv-readlik | GT | ALL | ALL | JointIndel | 20044 | 16350 | 3694 | 23088 | 19363 | 3725 | 0.8157054480143684 | 0.8386607761607762 | 0.8270238525830534 | 478 | 445 |
| sv-readlik | GT | ALL | ALL | JointStructuralVariant | 1680 | 813 | 867 | 0 | 0 | 0 | 0.48392857142857143 |  |  | 52 | 0 |
| sv-readlik | BASEPAIR | ALL | ALL | ALL | 969654 | 696807 | 272847 | 866978 | 696807 | 170171 | 0.7186140623356373 | 0.8037193562005034 | 0.7587878246703751 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Snv | 208362 | 197759 | 10603 | 193784 | 191175 | 2609 | 0.9491126021059503 | 0.9865365561656277 | 0.9674628006772611 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Insertion | 79376 | 67912 | 11464 | 323230 | 216660 | 106570 | 0.8555734730901028 | 0.6702966927574792 | 0.7516865880981007 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Deletion | 82226 | 70499 | 11727 | 231032 | 194365 | 36667 | 0.8573808770948362 | 0.841290384016067 | 0.8492594227648552 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 172616 | 104252 | 68364 |  | 0.6039532835890068 |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 210803 | 57019 | 0 | 0 | 0 | 0.7871011343354914 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 158861 | 216201 | 0 | 0 | 0 | 0.4235593048615962 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointIndel | 161602 | 138411 | 23191 | 726878 | 515277 | 211601 | 0.8564931127090011 | 0.7088906253869287 | 0.7757330340564366 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 369664 | 273220 | 0 | 0 | 0 | 0.5750088662962525 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | ALL | 1067628 | 794781 | 272847 | 1181670 | 1011499 | 170171 | 0.7444362643167844 | 0.8559910973452826 | 0.7963258190417012 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Snv | 208430 | 197827 | 10603 | 194208 | 191599 | 2609 | 0.9491292040493211 | 0.9865659499093755 | 0.9674855597635371 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Insertion | 104216 | 92752 | 11464 | 354930 | 248360 | 106570 | 0.8899976970906579 | 0.6997436114163356 | 0.7834862180178408 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Deletion | 108210 | 96483 | 11727 | 268324 | 231657 | 36667 | 0.8916273911838093 | 0.8633480419194705 | 0.8772598724520217 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 364208 | 295844 | 68364 |  | 0.8122940737161183 |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 212703 | 57019 | 0 | 0 | 0 | 0.7886008556958646 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 160849 | 216201 | 0 | 0 | 0 | 0.4265985943508818 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointIndel | 212426 | 189235 | 23191 | 987462 | 775861 | 211601 | 0.8908278647623172 | 0.7857122603198908 | 0.8349747968531579 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 373552 | 273220 | 0 | 0 | 0 | 0.5775636545799756 |  |  |  |  |

</details>

