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
| `poisson` | support (Flow) | yes | 106,587 | 216 s | 2.6 GB |
| `poisson-z` | panel (`-z`) | yes | 106,686 | 77 s | 3.0 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 109,535 | 150 s | 3.3 GB |
| `readlik-nomismap` | panel (default) | **no** | 112,614 | 112 s | 3.8 GB |
| `readlik-nolink` | panel (default) | **no** | 109,521 | 114 s | 3.4 GB |
| `readlik` | panel (default) | **no** | 109,476 | 111 s | 3.5 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9556 | 0.9917 | 0.9733 | 71,684 | 3,333 | 588 | 0.9659 | 0.9909 | 0.9783 |
| `poisson-z` | 0.9562 | 0.9914 | 0.9735 | 71,734 | 3,283 | 605 | 0.9664 | 0.9906 | 0.9784 |
| `readlik-support` | 0.9582 | 0.9930 | 0.9753 | 71,885 | 3,132 | 496 | 0.9676 | 0.9919 | 0.9796 |
| `readlik-nomismap` | 0.9585 | 0.9941 | **0.9760** | 71,901 | 3,116 | 416 | 0.9677 | 0.9928 | 0.9801 |
| `readlik-nolink` | 0.9582 | 0.9929 | 0.9752 | 71,881 | 3,136 | 501 | 0.9676 | 0.9919 | 0.9796 |
| `readlik` | 0.9573 | 0.9947 | 0.9756 | 71,813 | 3,204 | 376 | 0.9669 | 0.9932 | 0.9799 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7281 | 0.8462 | 0.7827 | 7,069 | 2,640 | 1,462 | 0.7684 | 0.7613 | 0.7648 |
| `poisson-z` | 0.7295 | 0.8497 | 0.7850 | 7,083 | 2,626 | 1,426 | 0.7729 | 0.7694 | 0.7712 |
| `readlik-support` | 0.8108 | 0.8465 | 0.8282 | 7,872 | 1,837 | 1,528 | 0.8689 | 0.7397 | 0.7992 |
| `readlik-nomismap` | 0.8172 | 0.8502 | 0.8333 | 7,934 | 1,775 | 1,488 | 0.8748 | 0.6679 | 0.7575 |
| `readlik-nolink` | 0.8102 | 0.8459 | 0.8277 | 7,866 | 1,843 | 1,533 | 0.8691 | 0.7382 | 0.7983 |
| `readlik` | 0.8178 | 0.8529 | **0.8350** | 7,940 | 1,769 | 1,453 | 0.8742 | 0.7640 | 0.8154 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8251 | 0.8047 | 0.8148 | 8,222 | 1,743 | 2,148 | 0.8750 | 0.6764 | 0.7630 |
| `poisson-z` | 0.8259 | 0.8040 | 0.8148 | 8,230 | 1,735 | 2,159 | 0.8737 | 0.7015 | 0.7782 |
| `readlik-support` | 0.8624 | 0.8865 | 0.8743 | 8,594 | 1,371 | 1,178 | 0.8768 | 0.8317 | 0.8536 |
| `readlik-nomismap` | 0.8639 | 0.8965 | 0.8799 | 8,609 | 1,356 | 1,061 | 0.8828 | 0.8142 | 0.8471 |
| `readlik-nolink` | 0.8621 | 0.8861 | 0.8739 | 8,591 | 1,374 | 1,183 | 0.8774 | 0.8308 | 0.8535 |
| `readlik` | 0.8650 | 0.8973 | **0.8809** | 8,620 | 1,345 | 1,050 | 0.8849 | 0.8450 | 0.8645 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.8746 | — | 0 | 0 | 116 | — | 0.5789 | — |
| `poisson-z` | — | 0.8772 | — | 0 | 0 | 116 | — | 0.5131 | — |
| `readlik-support` | — | 0.9031 | — | 0 | 0 | 84 | — | 0.8006 | — |
| `readlik-nomismap` | — | 0.8958 | — | 0 | 0 | 92 | — | 0.6130 | — |
| `readlik-nolink` | — | 0.9054 | — | 0 | 0 | 82 | — | 0.8145 | — |
| `readlik` | — | 0.9073 | — | 0 | 0 | 80 | — | 0.7931 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7772 | 0.8261 | 0.8009 | 15,291 | 4,383 | 3,726 | 0.8227 | 0.7065 | 0.7602 |
| `poisson-z` | 0.7783 | 0.8275 | 0.8022 | 15,313 | 4,361 | 3,701 | 0.8243 | 0.7181 | 0.7675 |
| `readlik-support` | 0.8369 | 0.8684 | 0.8524 | 16,466 | 3,208 | 2,790 | 0.8729 | 0.7845 | 0.8264 |
| `readlik-nomismap` | 0.8409 | 0.8746 | 0.8574 | 16,543 | 3,131 | 2,641 | 0.8789 | 0.7309 | 0.7981 |
| `readlik-nolink` | 0.8365 | 0.8680 | 0.8520 | 16,457 | 3,217 | 2,798 | 0.8733 | 0.7837 | 0.8261 |
| `readlik` | 0.8417 | 0.8768 | **0.8589** | 16,560 | 3,114 | 2,583 | 0.8797 | 0.8031 | 0.8396 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9185 | 0.9531 | 0.9355 | 86,975 | 7,716 | 4,314 | 0.9040 | 0.8619 | 0.8825 |
| `poisson-z` | 0.9193 | 0.9532 | 0.9359 | 87,047 | 7,644 | 4,306 | 0.9046 | 0.8695 | 0.8867 |
| `readlik-support` | 0.9330 | 0.9642 | 0.9484 | 88,351 | 6,340 | 3,286 | 0.9303 | 0.9102 | 0.9202 |
| `readlik-nomismap` | 0.9340 | 0.9667 | 0.9501 | 88,444 | 6,247 | 3,057 | 0.9332 | 0.8777 | 0.9046 |
| `readlik-nolink` | 0.9329 | 0.9641 | 0.9482 | 88,338 | 6,353 | 3,299 | 0.9305 | 0.9098 | 0.9200 |
| `readlik` | 0.9333 | 0.9677 | **0.9502** | 88,373 | 6,318 | 2,959 | 0.9332 | 0.9229 | 0.9280 |

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
| `sm50-readlik` | Insertion | 0.8689 | 0.8725 | **0.8707** |
| `sm50-readlik` | Deletion | 0.8767 | 0.8943 | **0.8854** |
| `sm50-readlik` | ALL | 0.9297 | 0.9646 | **0.9468** |

The insertion BASEPAIR precision gap collapses from **0.005 to -0.002**, and insertion BASEPAIR F1 goes from 0.8134 for `poisson-z` against 0.8707 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it. On chr20, of the 246, only **35 are confirmed true**, **73 are confirmed false**, and **138 fall outside the SV confident region** and cannot be judged at all. See *Known bad output* for the worst of the unjudged ones.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.4889 | 0.5021 | 0.4954 | 374 | 362 | 391 |
| `poisson-z` | 0.4902 | 0.4959 | 0.4930 | 375 | 372 | 390 |
| `readlik-support` | 0.5569 | 0.5225 | **0.5391** | 426 | 371 | 339 |
| `readlik-nomismap` | 0.5451 | 0.4711 | 0.5054 | 417 | 458 | 348 |
| `readlik-nolink` | 0.5556 | 0.5056 | 0.5294 | 425 | 399 | 340 |
| `readlik` | 0.5477 | 0.5146 | 0.5306 | 419 | 383 | 346 |

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
| **floor 0.02, cap 0.7 (current defaults)** | 0.9502 | 0.9756 | 0.8350 | 0.8809 | 0.9280 |
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
| readlik-support | GT | ALL | ALL | ALL | 94691 | 88351 | 6340 | 91786 | 88500 | 3286 | 0.9330453791807035 | 0.9641993332316475 | 0.948366572424922 | 540 | 552 |
| readlik-support | GT | ALL | ALL | Snv | 75017 | 71885 | 3132 | 70584 | 70088 | 496 | 0.9582494634549502 | 0.9929729117080358 | 0.9753022228335657 | 66 | 208 |
| readlik-support | GT | ALL | ALL | Insertion | 9709 | 7872 | 1837 | 9952 | 8424 | 1528 | 0.8107941085590689 | 0.8464630225080386 | 0.8282447170050277 | 255 | 179 |
| readlik-support | GT | ALL | ALL | Deletion | 9965 | 8594 | 1371 | 10383 | 9205 | 1178 | 0.8624184646261916 | 0.8865453144563228 | 0.874315475322253 | 219 | 162 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 867 | 783 | 84 |  | 0.903114186851211 |  | 0 | 3 |
| readlik-support | GT | ALL | ALL | JointIndel | 19674 | 16466 | 3208 | 21202 | 18412 | 2790 | 0.8369421571617363 | 0.8684086406942741 | 0.8523850951420779 | 474 | 344 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 390714 | 363489 | 27225 | 399342 | 363489 | 35853 | 0.9303198759194705 | 0.9102198115900657 | 0.9201600899176767 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 200440 | 193944 | 6496 | 190886 | 189344 | 1542 | 0.9675912991418878 | 0.9919218800750186 | 0.9796055375066295 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 95512 | 82994 | 12518 | 112558 | 83264 | 29294 | 0.8689379345003769 | 0.739743065797189 | 0.799152612029216 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 99134 | 86922 | 12212 | 104682 | 87061 | 17621 | 0.8768132023321968 | 0.8316711564547868 | 0.8536458015877726 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 6820 | 5460 | 1360 |  | 0.8005865102639296 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 194646 | 169916 | 24730 | 224060 | 175785 | 48275 | 0.8729488404590898 | 0.7845443184861198 | 0.8263889952306758 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 88444 | 6247 | 91775 | 88718 | 3057 | 0.9340275210949298 | 0.9666902751293925 | 0.9500782526888283 | 530 | 404 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 71901 | 3116 | 70709 | 70293 | 416 | 0.9584627484436861 | 0.9941167319577423 | 0.9759642204067236 | 154 | 95 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 7934 | 1775 | 9930 | 8442 | 1488 | 0.8171799361417241 | 0.8501510574018127 | 0.8333394982623671 | 189 | 157 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 8609 | 1356 | 10253 | 9192 | 1061 | 0.86392373306573 | 0.8965180922656784 | 0.8799191724331195 | 187 | 148 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 883 | 791 | 92 |  | 0.8958097395243488 |  | 0 | 4 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 16543 | 3131 | 21066 | 18425 | 2641 | 0.8408559520178917 | 0.874632108611032 | 0.8574115218055463 | 376 | 309 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390714 | 364607 | 26107 | 415412 | 364607 | 50805 | 0.9331813039717031 | 0.8776997294252453 | 0.9045905974996463 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 193965 | 6475 | 190638 | 189262 | 1376 | 0.9676960686489723 | 0.9927821315792235 | 0.9800786008661931 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 83557 | 11955 | 125272 | 83673 | 41599 | 0.8748324817823938 | 0.6679305830512804 | 0.7575075952341513 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 87516 | 11618 | 107486 | 87517 | 19969 | 0.8828050920975649 | 0.8142176655564446 | 0.847125352900624 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 8994 | 5513 | 3481 |  | 0.6129641983544585 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 171073 | 23573 | 241752 | 176703 | 65049 | 0.8788929646640568 | 0.7309267348356994 | 0.7981097077288735 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 94691 | 88338 | 6353 | 91780 | 88481 | 3299 | 0.9329080905260267 | 0.9640553497494008 | 0.9482260083679108 | 536 | 555 |
| readlik-nolink | GT | ALL | ALL | Snv | 75017 | 71881 | 3136 | 70580 | 70079 | 501 | 0.9581961422077663 | 0.9929016718617172 | 0.9752402413748564 | 61 | 211 |
| readlik-nolink | GT | ALL | ALL | Insertion | 9709 | 7866 | 1843 | 9950 | 8417 | 1533 | 0.8101761252446184 | 0.845929648241206 | 0.8276669468992339 | 256 | 179 |
| readlik-nolink | GT | ALL | ALL | Deletion | 9965 | 8591 | 1374 | 10383 | 9200 | 1183 | 0.862117410938284 | 0.8860637580660695 | 0.8739265776043418 | 219 | 160 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 867 | 785 | 82 |  | 0.9054209919261822 |  | 0 | 5 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 19674 | 16457 | 3217 | 21200 | 18402 | 2798 | 0.8364847006201077 | 0.8680188679245283 | 0.8519600853501408 | 475 | 344 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 390714 | 363542 | 27172 | 399590 | 363542 | 36048 | 0.930455525013181 | 0.909787532220526 | 0.9200054662509616 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 200440 | 193939 | 6501 | 190878 | 189341 | 1537 | 0.9675663540211534 | 0.9919477362503798 | 0.979605361664247 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 95512 | 83013 | 12499 | 112824 | 83287 | 29537 | 0.8691368623837842 | 0.7382028646387293 | 0.7983369175642107 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 99134 | 86978 | 12156 | 104798 | 87066 | 17732 | 0.8773780942966086 | 0.8307982976774366 | 0.8534531100956365 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 6728 | 5480 | 1248 |  | 0.8145065398335315 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 194646 | 169991 | 24655 | 224350 | 175833 | 48517 | 0.8733341553384092 | 0.7837441497659906 | 0.8261173089151967 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 94691 | 88373 | 6318 | 91514 | 88555 | 2959 | 0.9332777138270797 | 0.9676661494416154 | 0.9501608849677771 | 465 | 435 |
| readlik | GT | ALL | ALL | Snv | 75017 | 71813 | 3204 | 70548 | 70172 | 376 | 0.9572896810056387 | 0.9946702954017123 | 0.9756220632591283 | 111 | 105 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 7940 | 1769 | 9875 | 8422 | 1453 | 0.8177979194561746 | 0.8528607594936709 | 0.8349613999409576 | 177 | 165 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 8620 | 1345 | 10228 | 9178 | 1050 | 0.8650275965880582 | 0.8973406335549472 | 0.8808878851633128 | 177 | 161 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 863 | 783 | 80 |  | 0.9073001158748552 |  | 0 | 4 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 16560 | 3114 | 20966 | 18383 | 2583 | 0.8417200365965233 | 0.8768005341982257 | 0.858902232857039 | 354 | 330 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390714 | 364616 | 26098 | 395060 | 364616 | 30444 | 0.9332043387234652 | 0.9229382878550094 | 0.928042923283285 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 193805 | 6635 | 190430 | 189132 | 1298 | 0.9668978247854719 | 0.9931838470829176 | 0.9798645792561853 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 83495 | 12017 | 109518 | 83667 | 25851 | 0.8741833486891699 | 0.7639566098723497 | 0.8153615251017712 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 87728 | 11406 | 103870 | 87766 | 16104 | 0.8849436116771239 | 0.8449600462116107 | 0.8644897553774248 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 6882 | 5458 | 1424 |  | 0.7930834059866317 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 171223 | 23423 | 220270 | 176891 | 43379 | 0.8796635944226956 | 0.8030644209379398 | 0.8396205787586538 |  |  |

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

