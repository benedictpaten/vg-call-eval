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

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.7`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap down from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every arm on this page was re-run together on one build, so the wall-clock column compares runs made on the same machine in the same session rather than a mixture of vintages.

Two changes since the accuracy results were first produced left the calls untouched. The read path was optimised (vg `44fd008`) — on chr20 `readlik` went **506 s to under 100 s**, so the read-likelihood caller is now near parity with the Poisson caller at matched enumeration rather than 5.9x, and `readlik-support` is *faster* than `poisson`. Then `AD`, `BL`, `GQI` and the explained-share scaling of `GQ` were added — which rescales a quality and does not change a genotype. Both are confirmed by the variant counts below, which are unchanged to the record.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 124,445 | 329 s | 3.1 GB |
| `poisson-z` | panel (`-z`) | yes | 124,769 | 112 s | 3.1 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 105,198 | 153 s | 4.2 GB |
| `readlik-nomismap` | panel (default) | **no** | 121,500 | 135 s | 3.9 GB |
| `readlik-nolink` | panel (default) | **no** | 105,251 | 127 s | 3.8 GB |
| `readlik` | panel (default) | **no** | 105,251 | 136 s | 3.6 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9536 | 0.9581 | 0.9558 | 71,533 | 3,484 | 3,050 | 0.9626 | 0.9656 | 0.9641 |
| `poisson-z` | 0.9569 | 0.9583 | 0.9576 | 71,786 | 3,231 | 3,036 | 0.9658 | 0.9655 | 0.9656 |
| `readlik-support` | 0.9549 | 0.9935 | 0.9738 | 71,633 | 3,384 | 455 | 0.9639 | 0.9912 | 0.9774 |
| `readlik-nomismap` | 0.9625 | 0.9682 | 0.9654 | 72,207 | 2,810 | 2,295 | 0.9681 | 0.9720 | 0.9700 |
| `readlik-nolink` | 0.9598 | 0.9936 | 0.9764 | 72,004 | 3,013 | 449 | 0.9662 | 0.9912 | 0.9786 |
| `readlik` | 0.9598 | 0.9971 | **0.9781** | 71,999 | 3,018 | 205 | 0.9660 | 0.9936 | 0.9796 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7516 | 0.8184 | 0.7835 | 7,297 | 2,412 | 1,878 | 0.7817 | 0.5305 | 0.6321 |
| `poisson-z` | 0.7540 | 0.8220 | 0.7866 | 7,321 | 2,388 | 1,836 | 0.7890 | 0.5251 | 0.6306 |
| `readlik-support` | 0.8766 | 0.8733 | 0.8749 | 8,511 | 1,198 | 1,345 | 0.8927 | 0.5956 | 0.7145 |
| `readlik-nomismap` | 0.9012 | 0.9000 | 0.9006 | 8,750 | 959 | 1,028 | 0.9129 | 0.5385 | 0.6775 |
| `readlik-nolink` | 0.8844 | 0.8730 | 0.8787 | 8,587 | 1,122 | 1,351 | 0.9079 | 0.5808 | 0.7084 |
| `readlik` | 0.9013 | 0.9141 | **0.9077** | 8,751 | 958 | 867 | 0.9121 | 0.6604 | 0.7661 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8300 | 0.6849 | 0.7505 | 8,271 | 1,694 | 4,201 | 0.9098 | 0.5445 | 0.6813 |
| `poisson-z` | 0.8341 | 0.6848 | 0.7521 | 8,312 | 1,653 | 4,219 | 0.9096 | 0.4965 | 0.6424 |
| `readlik-support` | 0.8942 | 0.8677 | 0.8808 | 8,911 | 1,054 | 1,467 | 0.9177 | 0.7494 | 0.8251 |
| `readlik-nomismap` | 0.9208 | 0.9256 | 0.9232 | 9,176 | 789 | 778 | 0.9238 | 0.8155 | 0.8663 |
| `readlik-nolink` | 0.9017 | 0.8696 | 0.8853 | 8,985 | 980 | 1,446 | 0.9182 | 0.7735 | 0.8397 |
| `readlik` | 0.9191 | 0.9385 | **0.9287** | 9,159 | 806 | 632 | 0.9241 | 0.8577 | 0.8897 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.6561 | — | 0 | 0 | 489 | — | 0.4840 | — |
| `poisson-z` | — | 0.6549 | — | 0 | 0 | 518 | — | 0.3756 | — |
| `readlik-support` | — | 0.8208 | — | 0 | 0 | 228 | — | 0.6060 | — |
| `readlik-nomismap` | — | 0.7541 | — | 0 | 0 | 376 | — | 0.5230 | — |
| `readlik-nolink` | — | 0.8804 | — | 0 | 0 | 154 | — | 0.7625 | — |
| `readlik` | — | 0.9070 | — | 0 | 0 | 116 | — | 0.7840 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7913 | 0.7382 | 0.7639 | 15,568 | 4,106 | 6,568 | 0.8470 | 0.5347 | 0.6555 |
| `poisson-z` | 0.7946 | 0.7392 | 0.7659 | 15,633 | 4,041 | 6,573 | 0.8504 | 0.4980 | 0.6281 |
| `readlik-support` | 0.8855 | 0.8677 | 0.8765 | 17,422 | 2,252 | 3,040 | 0.9054 | 0.6620 | 0.7648 |
| `readlik-nomismap` | 0.9112 | 0.9020 | 0.9066 | 17,926 | 1,748 | 2,182 | 0.9185 | 0.6415 | 0.7554 |
| `readlik-nolink` | 0.8932 | 0.8718 | 0.8823 | 17,572 | 2,102 | 2,951 | 0.9132 | 0.6716 | 0.7740 |
| `readlik` | 0.9103 | 0.9253 | **0.9177** | 17,910 | 1,764 | 1,615 | 0.9182 | 0.7511 | 0.8263 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9198 | 0.9017 | 0.9107 | 87,101 | 7,590 | 9,618 | 0.9140 | 0.7221 | 0.8068 |
| `poisson-z` | 0.9232 | 0.9019 | 0.9124 | 87,419 | 7,272 | 9,609 | 0.9176 | 0.6876 | 0.7861 |
| `readlik-support` | 0.9405 | 0.9624 | 0.9513 | 89,055 | 5,636 | 3,495 | 0.9443 | 0.8239 | 0.8800 |
| `readlik-nomismap` | 0.9519 | 0.9526 | 0.9522 | 90,133 | 4,558 | 4,477 | 0.9534 | 0.8009 | 0.8706 |
| `readlik-nolink` | 0.9460 | 0.9634 | 0.9546 | 89,576 | 5,115 | 3,400 | 0.9497 | 0.8302 | 0.8859 |
| `readlik` | 0.9495 | 0.9801 | **0.9645** | 89,909 | 4,782 | 1,820 | 0.9522 | 0.8847 | 0.9172 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (58.9 Mb vs 59.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

That is exactly where the gap lives. 246 `readlik` calls carry a >=200 bp insertion allele; they contribute **27,951 FP bases and zero TP bases**, which is the whole of the precision difference. The Poisson caller scores better there because it does not emit them — at the two largest sites it emits nothing at all.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7802 | 0.8282 | **0.8035** |
| `sm50-poisson-z` | Deletion | 0.8949 | 0.6969 | **0.7836** |
| `sm50-poisson-z` | ALL | 0.9109 | 0.8746 | **0.8924** |
| `sm50-readlik` | Insertion | 0.8975 | 0.8901 | **0.8938** |
| `sm50-readlik` | Deletion | 0.9137 | 0.8999 | **0.9068** |
| `sm50-readlik` | ALL | 0.9453 | 0.9700 | **0.9575** |

The insertion BASEPAIR precision gap collapses from **-0.135 to -0.062**, and insertion BASEPAIR F1 goes from 0.8035 for `poisson-z` against 0.8938 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it. On chr20, of the 246, only **35 are confirmed true**, **73 are confirmed false**, and **138 fall outside the SV confident region** and cannot be judged at all. See *Known bad output* for the worst of the unjudged ones.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.4810 | 0.4289 | 0.4535 | 368 | 478 | 397 |
| `poisson-z` | 0.4824 | 0.4029 | 0.4391 | 369 | 535 | 396 |
| `readlik-support` | 0.4745 | 0.4449 | 0.4592 | 363 | 443 | 402 |
| `readlik-nomismap` | 0.4850 | 0.3998 | 0.4383 | 371 | 545 | 394 |
| `readlik-nolink` | 0.4967 | 0.4380 | 0.4655 | 380 | 476 | 385 |
| `readlik` | 0.4902 | 0.4986 | **0.4944** | 375 | 367 | 390 |

## Structural variants — aardvark (secondary)

Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

### SV insertion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.3925 | 325 | 503 | 802 | 314 | 488 | 0.3915 | 0.3920 |
| `poisson-z` | 0.4130 | 342 | 486 | 829 | 333 | 496 | 0.4017 | 0.4073 |
| `readlik-nomismap` | 0.4964 | 411 | 417 | 815 | 324 | 491 | 0.3975 | 0.4415 |
| `readlik` | 0.4940 | 409 | 419 | 708 | 293 | 415 | 0.4138 | 0.4504 |

### SV deletion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.5340 | 455 | 397 | 802 | 314 | 488 | 0.3915 | 0.4518 |
| `poisson-z` | 0.5622 | 479 | 373 | 829 | 333 | 496 | 0.4017 | 0.4686 |
| `readlik-nomismap` | 0.5270 | 449 | 403 | 815 | 324 | 491 | 0.3975 | 0.4532 |
| `readlik` | 0.5282 | 450 | 402 | 708 | 293 | 415 | 0.4138 | 0.4641 |

### SV (joint)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4643 | 780 | 900 | 802 | 314 | 488 | 0.3915 | 0.4248 |
| `poisson-z` | 0.4887 | 821 | 859 | 829 | 333 | 496 | 0.4017 | 0.4409 |
| `readlik-nomismap` | 0.5119 | 860 | 820 | 815 | 324 | 491 | 0.3975 | 0.4475 |
| `readlik` | 0.5113 | 859 | 821 | 708 | 293 | 415 | 0.4138 | 0.4574 |

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and the old cap of 0.1 was telling the model 0.1: overriding the mapper at exactly the sites that matter. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless — see 

The two graphs are put side by side in [tier2-chr20-graph-comparison.md](tier2-chr20-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9645 | 0.9781 | 0.9077 | 0.9287 | 0.9172 |

Only the current row is available here: the preserved old-default arms (`arms.floor-1e-8.json`, `arms.readlik.json`) exist for the 4-haplotype runs alone, so the before-and-after is on [tier2-chr20-4hap-results.md](tier2-chr20-4hap-results.md). Mixing rows from two graphs into one table is exactly what the one-build-per-matrix rule forbids. The full grids are in plan §9.20-§9.21.

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
| readlik-support | GT | ALL | ALL | ALL | 94691 | 89055 | 5636 | 92844 | 89349 | 3495 | 0.940480087864739 | 0.9623562104174744 | 0.9512923987708899 | 814 | 504 |
| readlik-support | GT | ALL | ALL | Snv | 75017 | 71633 | 3384 | 69867 | 69412 | 455 | 0.95489022488236 | 0.9934876264903316 | 0.9738066180631509 | 148 | 288 |
| readlik-support | GT | ALL | ALL | Insertion | 9709 | 8511 | 1198 | 10614 | 9269 | 1345 | 0.8766093315480482 | 0.8732805728283399 | 0.8749417860934054 | 351 | 78 |
| readlik-support | GT | ALL | ALL | Deletion | 9965 | 8911 | 1054 | 11091 | 9624 | 1467 | 0.8942298043151029 | 0.8677305923721936 | 0.8807809293263322 | 315 | 128 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1272 | 1044 | 228 |  | 0.8207547169811321 |  | 0 | 10 |
| readlik-support | GT | ALL | ALL | JointIndel | 19674 | 17422 | 2252 | 22977 | 19937 | 3040 | 0.8855342075836129 | 0.867693780737259 | 0.8765232242114674 | 666 | 216 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 390748 | 368986 | 21762 | 447858 | 368986 | 78872 | 0.9443068166695671 | 0.823890608183844 | 0.8799984736574744 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 200440 | 193209 | 7231 | 188550 | 186887 | 1663 | 0.9639243663939333 | 0.9911800583399629 | 0.9773622294857119 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 95512 | 85264 | 10248 | 141126 | 84053 | 57073 | 0.8927045816232515 | 0.5955883394980372 | 0.7144889717418104 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 99134 | 90976 | 8158 | 118282 | 88646 | 29636 | 0.9177073456130087 | 0.7494462386500059 | 0.8250857327643133 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 17532 | 10625 | 6907 |  | 0.6060346794433037 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 194646 | 176240 | 18406 | 276940 | 183324 | 93616 | 0.9054385910832999 | 0.6619628800462194 | 0.7647903214312664 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 90133 | 4558 | 94411 | 89934 | 4477 | 0.9518644855371683 | 0.9525796782154622 | 0.9522219475850148 | 449 | 210 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 72207 | 2810 | 72144 | 69849 | 2295 | 0.9625418238532599 | 0.968188622754491 | 0.9653569657198386 | 147 | 50 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 8750 | 959 | 10284 | 9256 | 1028 | 0.9012256669069935 | 0.9000388953714508 | 0.9006318901842488 | 171 | 69 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 9176 | 789 | 10454 | 9676 | 778 | 0.920822880080281 | 0.9255787258465659 | 0.9231946780583973 | 131 | 86 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1529 | 1153 | 376 |  | 0.7540876389797253 |  | 0 | 5 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 17926 | 1748 | 22267 | 20085 | 2182 | 0.9111517739148114 | 0.9020074549782189 | 0.9065565556418935 | 302 | 160 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390704 | 372509 | 18195 | 465092 | 372509 | 92583 | 0.9534302182726565 | 0.8009361588674929 | 0.8705555996989938 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 194039 | 6401 | 192594 | 187202 | 5392 | 0.9680652564358412 | 0.9720032815144812 | 0.9700302721984452 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 87195 | 8317 | 156848 | 84470 | 72378 | 0.9129219365105955 | 0.5385468734061002 | 0.6774534198910425 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 91584 | 7550 | 109426 | 89233 | 20193 | 0.9238404583694797 | 0.8154643320600223 | 0.866275935603281 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24286 | 12702 | 11584 |  | 0.5230173762661616 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 178779 | 15867 | 290560 | 186405 | 104155 | 0.918482784131192 | 0.641537031938326 | 0.7554272236138279 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 94691 | 89576 | 5115 | 92892 | 89492 | 3400 | 0.945982194717555 | 0.9633983550790165 | 0.9546108453048061 | 824 | 496 |
| readlik-nolink | GT | ALL | ALL | Snv | 75017 | 72004 | 3013 | 69878 | 69429 | 449 | 0.9598357705586734 | 0.9935745155843041 | 0.9764137800833905 | 147 | 287 |
| readlik-nolink | GT | ALL | ALL | Insertion | 9709 | 8587 | 1122 | 10637 | 9286 | 1351 | 0.8844371201977547 | 0.8729905048415907 | 0.8786765350235057 | 363 | 78 |
| readlik-nolink | GT | ALL | ALL | Deletion | 9965 | 8985 | 980 | 11089 | 9643 | 1446 | 0.9016557952834923 | 0.8696005050049599 | 0.8853380900228667 | 314 | 123 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1288 | 1134 | 154 |  | 0.8804347826086957 |  | 0 | 8 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 19674 | 17572 | 2102 | 23014 | 20063 | 2951 | 0.893158483277422 | 0.8717737029634136 | 0.8823365388993991 | 677 | 209 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 390702 | 371032 | 19670 | 446894 | 371032 | 75862 | 0.9496547240607931 | 0.830246098627415 | 0.8859450140640595 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 200440 | 193674 | 6766 | 188572 | 186910 | 1662 | 0.9662442626222311 | 0.9911863903442717 | 0.9785564167067136 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 95512 | 86719 | 8793 | 144952 | 84182 | 60770 | 0.9079382695368121 | 0.5807577680887466 | 0.7083947154445023 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 99134 | 91029 | 8105 | 114564 | 88616 | 25948 | 0.9182419755078984 | 0.7735065116441465 | 0.839682911159015 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 16310 | 12436 | 3874 |  | 0.7624770079705702 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 194646 | 177748 | 16898 | 275826 | 185234 | 90592 | 0.9131859889234816 | 0.6715610566081515 | 0.7739533565692603 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 94691 | 89909 | 4782 | 91384 | 89564 | 1820 | 0.9494988964104297 | 0.9800840409699728 | 0.9645490714732674 | 399 | 294 |
| readlik | GT | ALL | ALL | Snv | 75017 | 71999 | 3018 | 69773 | 69568 | 205 | 0.9597691189996934 | 0.9970619007352415 | 0.9780601517513223 | 115 | 109 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 8751 | 958 | 10095 | 9228 | 867 | 0.9013286641260686 | 0.9141158989598811 | 0.9076772475667192 | 160 | 77 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 9159 | 806 | 10269 | 9637 | 632 | 0.9191169091821375 | 0.938455545817509 | 0.9286855630907692 | 124 | 100 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1247 | 1131 | 116 |  | 0.9069767441860465 |  | 0 | 8 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 17910 | 1764 | 21611 | 19996 | 1615 | 0.9103385178408051 | 0.9252695386608671 | 0.9177433030372101 | 284 | 185 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390702 | 372028 | 18674 | 420516 | 372028 | 48488 | 0.9522039815511566 | 0.8846940425572392 | 0.9172084445857955 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 193626 | 6814 | 188016 | 186811 | 1205 | 0.966004789463181 | 0.9935909709811931 | 0.9796037081827017 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 87115 | 8397 | 127856 | 84435 | 43421 | 0.9120843454225647 | 0.660391377800025 | 0.7660946730662868 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 91612 | 7522 | 103796 | 89022 | 14774 | 0.9241229043516856 | 0.8576631084049482 | 0.8896535465201579 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 15942 | 12499 | 3443 |  | 0.7840296073265588 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 178727 | 15919 | 247594 | 185956 | 61638 | 0.9182156324815306 | 0.7510521256573262 | 0.8262638503914285 |  |  |

</details>

<details><summary><code>poisson</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-poisson | GT | ALL | ALL | ALL | 100207 | 88731 | 11476 | 101178 | 89926 | 11252 | 0.885477062480665 | 0.8887900531736148 | 0.8871304647451471 | 1788 | 760 |
| sv-poisson | GT | ALL | ALL | Snv | 78483 | 72509 | 5974 | 74051 | 70496 | 3555 | 0.9238816049335525 | 0.9519925456779786 | 0.937726446839757 | 315 | 224 |
| sv-poisson | GT | ALL | ALL | Insertion | 9857 | 7259 | 2598 | 11081 | 8813 | 2268 | 0.7364309627675764 | 0.7953253316487682 | 0.7647459348913181 | 812 | 160 |
| sv-poisson | GT | ALL | ALL | Deletion | 10187 | 8183 | 2004 | 14254 | 9569 | 4685 | 0.8032786885245902 | 0.67132033113512 | 0.7313951901291185 | 594 | 338 |
| sv-poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1792 | 1048 | 744 |  | 0.5848214285714286 |  | 0 | 38 |
| sv-poisson | GT | ALL | ALL | SvInsertion | 828 | 325 | 503 | 0 | 0 | 0 | 0.392512077294686 |  |  | 44 | 0 |
| sv-poisson | GT | ALL | ALL | SvDeletion | 852 | 455 | 397 | 0 | 0 | 0 | 0.534037558685446 |  |  | 23 | 0 |
| sv-poisson | GT | ALL | ALL | JointIndel | 20044 | 15442 | 4602 | 27127 | 19430 | 7697 | 0.7704051087607264 | 0.7162605522173481 | 0.742346854596734 | 1406 | 536 |
| sv-poisson | GT | ALL | ALL | JointStructuralVariant | 1680 | 780 | 900 | 0 | 0 | 0 | 0.4642857142857143 |  |  | 67 | 0 |
| sv-poisson | BASEPAIR | ALL | ALL | ALL | 969654 | 751030 | 218624 | 981332 | 751030 | 230302 | 0.7745340090382755 | 0.7653169365719247 | 0.76989788752969 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Snv | 208362 | 196500 | 11862 | 196980 | 189077 | 7903 | 0.9430702335358655 | 0.9598791755508174 | 0.9514004670122981 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Insertion | 79376 | 61824 | 17552 | 272336 | 186774 | 85562 | 0.7788752267687966 | 0.6858219258562952 | 0.7293927035591942 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Deletion | 82226 | 74175 | 8051 | 387054 | 257463 | 129591 | 0.9020869311410016 | 0.6651862530809655 | 0.7657322689112073 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 191060 | 134021 | 57039 |  | 0.701460274259395 |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 178542 | 89280 | 0 | 0 | 0 | 0.666644263727401 |  |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 249447 | 125615 | 0 | 0 | 0 | 0.6650820397694248 |  |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | JointIndel | 161602 | 135999 | 25603 | 850450 | 578258 | 272192 | 0.8415675548569943 | 0.6799435592921395 | 0.7521712241375851 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 427989 | 214895 | 0 | 0 | 0 | 0.6657328538274401 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | ALL | 1067628 | 849004 | 218624 | 2003506 | 1773204 | 230302 | 0.7952245538708239 | 0.8850505064621718 | 0.8377365239415575 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Snv | 208430 | 196568 | 11862 | 197764 | 189861 | 7903 | 0.9430888067936477 | 0.9600382273821323 | 0.9514880405555497 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Insertion | 104216 | 86664 | 17552 | 306682 | 221120 | 85562 | 0.8315805634451524 | 0.7210074278894751 | 0.7723565639807193 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Deletion | 108210 | 100159 | 8051 | 431678 | 302087 | 129591 | 0.9255983735329452 | 0.6997970709649322 | 0.7970134687911095 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 1067382 | 1010343 | 57039 |  | 0.9465617745099693 |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 180442 | 89280 | 0 | 0 | 0 | 0.6689925182224661 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 251435 | 125615 | 0 | 0 | 0 | 0.6668478981567432 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | JointIndel | 212426 | 186823 | 25603 | 1805742 | 1533550 | 272192 | 0.87947332247465 | 0.8492630730192907 | 0.8641042306162575 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 431877 | 214895 | 0 | 0 | 0 | 0.6677422646620448 |  |  |  |  |

</details>

<details><summary><code>poisson-z</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-poisson-z | GT | ALL | ALL | ALL | 100207 | 89216 | 10991 | 101352 | 90114 | 11238 | 0.8903170437195006 | 0.8891191096376984 | 0.889717673448042 | 1888 | 751 |
| sv-poisson-z | GT | ALL | ALL | Snv | 78483 | 72879 | 5604 | 74061 | 70518 | 3543 | 0.9285960016818929 | 0.9521610564264593 | 0.9402308991934452 | 395 | 221 |
| sv-poisson-z | GT | ALL | ALL | Insertion | 9857 | 7293 | 2564 | 11060 | 8846 | 2214 | 0.7398802881201176 | 0.7998191681735985 | 0.7686830493747834 | 823 | 131 |
| sv-poisson-z | GT | ALL | ALL | Deletion | 10187 | 8223 | 1964 | 14326 | 9616 | 4710 | 0.8072052616079317 | 0.6712271394667039 | 0.7329629387420492 | 597 | 330 |
| sv-poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1905 | 1134 | 771 |  | 0.5952755905511811 |  | 0 | 69 |
| sv-poisson-z | GT | ALL | ALL | SvInsertion | 828 | 342 | 486 | 0 | 0 | 0 | 0.41304347826086957 |  |  | 45 | 0 |
| sv-poisson-z | GT | ALL | ALL | SvDeletion | 852 | 479 | 373 | 0 | 0 | 0 | 0.562206572769953 |  |  | 28 | 0 |
| sv-poisson-z | GT | ALL | ALL | JointIndel | 20044 | 15516 | 4528 | 27291 | 19596 | 7695 | 0.7740969866294153 | 0.7180389139276685 | 0.7450149270539612 | 1420 | 530 |
| sv-poisson-z | GT | ALL | ALL | JointStructuralVariant | 1680 | 821 | 859 | 0 | 0 | 0 | 0.4886904761904762 |  |  | 73 | 0 |
| sv-poisson-z | BASEPAIR | ALL | ALL | ALL | 969654 | 750996 | 218658 | 999698 | 750996 | 248702 | 0.7744989449844997 | 0.7512228693065306 | 0.7626833598056619 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Snv | 208362 | 197615 | 10747 | 197032 | 189068 | 7964 | 0.948421497201985 | 0.9595801697186244 | 0.9539682035170262 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Insertion | 79376 | 61975 | 17401 | 277304 | 187899 | 89405 | 0.780777565007055 | 0.677592101087615 | 0.7255344417193318 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Deletion | 82226 | 73991 | 8235 | 410802 | 259036 | 151766 | 0.8998491961180162 | 0.6305616817834383 | 0.741513838732905 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 183484 | 128293 | 55191 |  | 0.6992053803056397 |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 185653 | 82169 | 0 | 0 | 0 | 0.6931954805803855 |  |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 239343 | 135719 | 0 | 0 | 0 | 0.6381424937743626 |  |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | JointIndel | 161602 | 135966 | 25636 | 871590 | 575228 | 296362 | 0.8413633494634968 | 0.6599754471712617 | 0.7397119877809546 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 424996 | 217888 | 0 | 0 | 0 | 0.661077270549586 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | ALL | 1067628 | 848970 | 218658 | 2415100 | 2166398 | 248702 | 0.7951927075723004 | 0.8970220694795247 | 0.8430435874390028 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Snv | 208430 | 197683 | 10747 | 200736 | 192772 | 7964 | 0.9484383246173775 | 0.9603260003188268 | 0.9543451445838782 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Insertion | 104216 | 86815 | 17401 | 311956 | 222551 | 89405 | 0.8330294772395793 | 0.7134050955904038 | 0.7685905168974295 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Deletion | 108210 | 99975 | 8235 | 459356 | 307590 | 151766 | 0.9238979761574716 | 0.6696113689600223 | 0.7764655921093248 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 1443052 | 1387861 | 55191 |  | 0.9617539769876623 |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 187553 | 82169 | 0 | 0 | 0 | 0.6953567006028429 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 241331 | 135719 | 0 | 0 | 0 | 0.6400503911948018 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | JointIndel | 212426 | 186790 | 25636 | 2214364 | 1918002 | 296362 | 0.8793179742592715 | 0.8661638285304494 | 0.8726913358394996 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 428884 | 217888 | 0 | 0 | 0 | 0.6631146679200708 |  |  |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik-nomismap | GT | ALL | ALL | ALL | 100207 | 90998 | 9209 | 99212 | 91719 | 7493 | 0.9081002325186863 | 0.9244748619118655 | 0.916214391007768 | 1062 | 455 |
| sv-readlik-nomismap | GT | ALL | ALL | Snv | 78483 | 72857 | 5626 | 73519 | 70652 | 2867 | 0.9283156861995592 | 0.9610032780641738 | 0.9443767139275245 | 295 | 182 |
| sv-readlik-nomismap | GT | ALL | ALL | Insertion | 9857 | 8468 | 1389 | 11785 | 9739 | 2046 | 0.8590849142741199 | 0.8263894781501909 | 0.8424200773202462 | 363 | 116 |
| sv-readlik-nomismap | GT | ALL | ALL | Deletion | 10187 | 8813 | 1374 | 11932 | 10132 | 1800 | 0.865122214587219 | 0.8491451558833389 | 0.8570592317371769 | 335 | 140 |
| sv-readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1976 | 1196 | 780 |  | 0.6052631578947368 |  | 0 | 17 |
| sv-readlik-nomismap | GT | ALL | ALL | SvInsertion | 828 | 411 | 417 | 0 | 0 | 0 | 0.4963768115942029 |  |  | 33 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | SvDeletion | 852 | 449 | 403 | 0 | 0 | 0 | 0.5269953051643192 |  |  | 36 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | JointIndel | 20044 | 17281 | 2763 | 25693 | 21067 | 4626 | 0.862153262821792 | 0.8199509594052855 | 0.8405227044364363 | 698 | 273 |
| sv-readlik-nomismap | GT | ALL | ALL | JointStructuralVariant | 1680 | 860 | 820 | 0 | 0 | 0 | 0.5119047619047619 |  |  | 69 | 0 |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 969654 | 685511 | 284143 | 936768 | 685511 | 251257 | 0.7069645461164498 | 0.7317831095852976 | 0.7191597663056762 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 208362 | 197345 | 11017 | 195980 | 189139 | 6841 | 0.9471256755070502 | 0.9650933768751914 | 0.9560251116225983 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 79376 | 70758 | 8618 | 378732 | 211348 | 167384 | 0.8914281394880065 | 0.5580410422145475 | 0.6863940183053092 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 82226 | 73833 | 8393 | 228920 | 180973 | 47947 | 0.8979276627830614 | 0.7905512842914555 | 0.8408252507310016 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 202170 | 117379 | 84791 |  | 0.5805955384082703 |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 197581 | 70241 | 0 | 0 | 0 | 0.7377325238404612 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 155049 | 220013 | 0 | 0 | 0 | 0.4133956519188827 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 161602 | 144591 | 17011 | 809822 | 509700 | 300122 | 0.894735213673098 | 0.6293975713181416 | 0.7389699585287817 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 352630 | 290254 | 0 | 0 | 0 | 0.5485126399163768 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | ALL | 1067628 | 783485 | 284143 | 2215898 | 1964641 | 251257 | 0.7338557999602858 | 0.8866116581178376 | 0.8030338460407556 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Snv | 208430 | 197413 | 11017 | 201292 | 194451 | 6841 | 0.9471429256824834 | 0.9660145460326293 | 0.9564856598667125 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Insertion | 104216 | 95598 | 8618 | 414414 | 247030 | 167384 | 0.9173063637061488 | 0.5960947265295091 | 0.7226127819586154 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Deletion | 108210 | 99817 | 8393 | 267546 | 219599 | 47947 | 0.9224378523241844 | 0.8207896959775142 | 0.8686502058838314 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 1332646 | 1247855 | 84791 |  | 0.9363739507716228 |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 199481 | 70241 | 0 | 0 | 0 | 0.7395800120123683 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 157037 | 220013 | 0 | 0 | 0 | 0.4164885293727622 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointIndel | 212426 | 195415 | 17011 | 2014606 | 1714484 | 300122 | 0.9199203487332059 | 0.8510269501828149 | 0.8841336038319413 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 356518 | 290254 | 0 | 0 | 0 | 0.5512267074022994 |  |  |  |  |

</details>

<details><summary><code>readlik</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik | GT | ALL | ALL | ALL | 100207 | 90764 | 9443 | 95801 | 91040 | 4761 | 0.9057650663127327 | 0.9503032327428732 | 0.9274997811885645 | 950 | 719 |
| sv-readlik | GT | ALL | ALL | Snv | 78483 | 72633 | 5850 | 70946 | 70167 | 779 | 0.9254615649248882 | 0.9890198178896625 | 0.9561856663868713 | 223 | 373 |
| sv-readlik | GT | ALL | ALL | Insertion | 9857 | 8473 | 1384 | 11487 | 9675 | 1812 | 0.8595921680024348 | 0.8422564638286759 | 0.8508360217413324 | 341 | 137 |
| sv-readlik | GT | ALL | ALL | Deletion | 10187 | 8799 | 1388 | 11718 | 10041 | 1677 | 0.8637479140080495 | 0.8568868407578085 | 0.8603036980340681 | 329 | 183 |
| sv-readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1650 | 1157 | 493 |  | 0.7012121212121212 |  | 0 | 26 |
| sv-readlik | GT | ALL | ALL | SvInsertion | 828 | 409 | 419 | 0 | 0 | 0 | 0.4939613526570048 |  |  | 25 | 0 |
| sv-readlik | GT | ALL | ALL | SvDeletion | 852 | 450 | 402 | 0 | 0 | 0 | 0.528169014084507 |  |  | 32 | 0 |
| sv-readlik | GT | ALL | ALL | JointIndel | 20044 | 17272 | 2772 | 24855 | 20873 | 3982 | 0.8617042506485731 | 0.83979078656206 | 0.8506064074360891 | 670 | 346 |
| sv-readlik | GT | ALL | ALL | JointStructuralVariant | 1680 | 859 | 821 | 0 | 0 | 0 | 0.5113095238095238 |  |  | 57 | 0 |
| sv-readlik | BASEPAIR | ALL | ALL | ALL | 969654 | 696727 | 272927 | 917154 | 696727 | 220427 | 0.718531558679694 | 0.759661954262861 | 0.7385245345578353 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Snv | 208362 | 196867 | 11495 | 191440 | 188719 | 2721 | 0.9448315911730546 | 0.98578666945257 | 0.9648747310141292 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Insertion | 79376 | 70621 | 8755 | 365796 | 212450 | 153346 | 0.8897021769804475 | 0.5807881988868112 | 0.7027975611324834 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Deletion | 82226 | 74016 | 8210 | 227368 | 189196 | 38172 | 0.900153236202661 | 0.8321135779881074 | 0.8647971824872497 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 194146 | 118712 | 75434 |  | 0.6114573568345472 |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 199542 | 68280 | 0 | 0 | 0 | 0.7450545511571118 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 164079 | 210983 | 0 | 0 | 0 | 0.43747167135033677 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointIndel | 161602 | 144637 | 16965 | 787310 | 520358 | 266952 | 0.8950198636155493 | 0.660931526336513 | 0.76036674202161 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 363621 | 279263 | 0 | 0 | 0 | 0.5656090367780191 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | ALL | 1067628 | 794701 | 272927 | 2186902 | 1966475 | 220427 | 0.7443613318496705 | 0.8992058171788219 | 0.8144894354670457 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Snv | 208430 | 196935 | 11495 | 196756 | 194035 | 2721 | 0.9448495897903373 | 0.9861706885685824 | 0.9650680326870069 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Insertion | 104216 | 95461 | 8755 | 400804 | 247458 | 153346 | 0.9159917862900131 | 0.6174040179239728 | 0.7376269162686325 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Deletion | 108210 | 100000 | 8210 | 265612 | 227440 | 38172 | 0.924129008409574 | 0.8562866135566164 | 0.8889152502791228 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 1323730 | 1248296 | 75434 |  | 0.943014058758206 |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 201442 | 68280 | 0 | 0 | 0 | 0.74685046084487 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 166067 | 210983 | 0 | 0 | 0 | 0.440437607744331 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointIndel | 212426 | 195461 | 16965 | 1990146 | 1723194 | 266952 | 0.9201368947304003 | 0.8658631075308043 | 0.8921753527618291 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 367509 | 279263 | 0 | 0 | 0 | 0.5682203311213225 |  |  |  |  |

</details>

