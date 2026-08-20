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
| `poisson` | support (Flow) | yes | 124,445 | 329 s | 2.9 GB |
| `poisson-z` | panel (`-z`) | yes | 124,769 | 109 s | 3.3 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 117,324 | 175 s | 3.8 GB |
| `readlik-nomismap` | panel (default) | **no** | 138,405 | 194 s | 3.6 GB |
| `readlik-nolink` | panel (default) | **no** | 117,047 | 139 s | 3.4 GB |
| `readlik` | panel (default) | **no** | 116,945 | 180 s | 3.9 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9536 | 0.9581 | 0.9558 | 71,533 | 3,484 | 3,050 | 0.9626 | 0.9656 | 0.9641 |
| `poisson-z` | 0.9569 | 0.9583 | 0.9576 | 71,786 | 3,231 | 3,036 | 0.9658 | 0.9655 | 0.9656 |
| `readlik-support` | 0.9731 | 0.9914 | 0.9822 | 73,000 | 2,017 | 612 | 0.9780 | 0.9899 | 0.9839 |
| `readlik-nomismap` | 0.9761 | 0.9671 | 0.9716 | 73,221 | 1,796 | 2,427 | 0.9800 | 0.9715 | 0.9757 |
| `readlik-nolink` | 0.9735 | 0.9913 | 0.9823 | 73,029 | 1,988 | 621 | 0.9781 | 0.9898 | 0.9839 |
| `readlik` | 0.9731 | 0.9955 | **0.9842** | 73,002 | 2,015 | 324 | 0.9778 | 0.9927 | 0.9852 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7516 | 0.8184 | 0.7835 | 7,297 | 2,412 | 1,878 | 0.7817 | 0.5305 | 0.6321 |
| `poisson-z` | 0.7540 | 0.8220 | 0.7866 | 7,321 | 2,388 | 1,836 | 0.7890 | 0.5251 | 0.6306 |
| `readlik-support` | 0.8939 | 0.8652 | 0.8793 | 8,679 | 1,030 | 1,476 | 0.9140 | 0.5813 | 0.7107 |
| `readlik-nomismap` | 0.9115 | 0.8935 | 0.9024 | 8,850 | 859 | 1,129 | 0.9184 | 0.5371 | 0.6778 |
| `readlik-nolink` | 0.8951 | 0.8650 | 0.8798 | 8,691 | 1,018 | 1,480 | 0.9159 | 0.5712 | 0.7036 |
| `readlik` | 0.9119 | 0.9087 | **0.9103** | 8,854 | 855 | 948 | 0.9164 | 0.6473 | 0.7587 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8300 | 0.6849 | 0.7505 | 8,271 | 1,694 | 4,201 | 0.9098 | 0.5445 | 0.6813 |
| `poisson-z` | 0.8341 | 0.6848 | 0.7521 | 8,312 | 1,653 | 4,219 | 0.9096 | 0.4965 | 0.6424 |
| `readlik-support` | 0.9135 | 0.8618 | 0.8869 | 9,103 | 862 | 1,581 | 0.9306 | 0.7387 | 0.8236 |
| `readlik-nomismap` | 0.9325 | 0.9177 | 0.9250 | 9,292 | 673 | 887 | 0.9375 | 0.7917 | 0.8584 |
| `readlik-nolink` | 0.9147 | 0.8623 | 0.8877 | 9,115 | 850 | 1,577 | 0.9296 | 0.7598 | 0.8362 |
| `readlik` | 0.9313 | 0.9312 | **0.9312** | 9,280 | 685 | 728 | 0.9382 | 0.8364 | 0.8843 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.6561 | — | 0 | 0 | 489 | — | 0.4840 | — |
| `poisson-z` | — | 0.6549 | — | 0 | 0 | 518 | — | 0.3756 | — |
| `readlik-support` | — | 0.8646 | — | 0 | 0 | 147 | — | 0.5250 | — |
| `readlik-nomismap` | — | 0.7264 | — | 0 | 0 | 356 | — | 0.4244 | — |
| `readlik-nolink` | — | 0.8806 | — | 0 | 0 | 128 | — | 0.6677 | — |
| `readlik` | — | 0.9024 | — | 0 | 0 | 101 | — | 0.6828 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7913 | 0.7382 | 0.7639 | 15,568 | 4,106 | 6,568 | 0.8470 | 0.5347 | 0.6555 |
| `poisson-z` | 0.7946 | 0.7392 | 0.7659 | 15,633 | 4,041 | 6,573 | 0.8504 | 0.4980 | 0.6281 |
| `readlik-support` | 0.9038 | 0.8635 | 0.8832 | 17,782 | 1,892 | 3,204 | 0.9224 | 0.6464 | 0.7601 |
| `readlik-nomismap` | 0.9221 | 0.8954 | 0.9086 | 18,142 | 1,532 | 2,372 | 0.9281 | 0.6305 | 0.7509 |
| `readlik-nolink` | 0.9051 | 0.8644 | 0.8843 | 17,806 | 1,868 | 3,185 | 0.9229 | 0.6550 | 0.7662 |
| `readlik` | 0.9217 | 0.9192 | **0.9205** | 18,134 | 1,540 | 1,777 | 0.9275 | 0.7304 | 0.8172 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9198 | 0.9017 | 0.9107 | 87,101 | 7,590 | 9,618 | 0.9140 | 0.7221 | 0.8068 |
| `poisson-z` | 0.9232 | 0.9019 | 0.9124 | 87,419 | 7,272 | 9,609 | 0.9176 | 0.6876 | 0.7861 |
| `readlik-support` | 0.9587 | 0.9598 | 0.9593 | 90,782 | 3,909 | 3,816 | 0.9600 | 0.8142 | 0.8811 |
| `readlik-nomismap` | 0.9649 | 0.9502 | 0.9575 | 91,363 | 3,328 | 4,799 | 0.9645 | 0.7961 | 0.8722 |
| `readlik-nolink` | 0.9593 | 0.9599 | 0.9596 | 90,835 | 3,856 | 3,806 | 0.9605 | 0.8199 | 0.8846 |
| `readlik` | 0.9625 | 0.9775 | **0.9699** | 91,136 | 3,555 | 2,101 | 0.9630 | 0.8730 | 0.9158 |

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
| `sm50-readlik` | Insertion | 0.9095 | 0.8811 | **0.8951** |
| `sm50-readlik` | Deletion | 0.9271 | 0.8858 | **0.9060** |
| `sm50-readlik` | ALL | 0.9584 | 0.9634 | **0.9609** |

The insertion BASEPAIR precision gap collapses from **-0.122 to -0.053**, and insertion BASEPAIR F1 goes from 0.8035 for `poisson-z` against 0.8951 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it. On chr20, of the 246, only **35 are confirmed true**, **73 are confirmed false**, and **138 fall outside the SV confident region** and cannot be judged at all. See *Known bad output* for the worst of the unjudged ones.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.4810 | 0.4289 | 0.4535 | 368 | 478 | 397 |
| `poisson-z` | 0.4824 | 0.4029 | 0.4391 | 369 | 535 | 396 |
| `readlik-support` | 0.5307 | 0.4410 | 0.4817 | 406 | 502 | 359 |
| `readlik-nomismap` | 0.5529 | 0.3922 | 0.4589 | 423 | 643 | 342 |
| `readlik-nolink` | 0.5490 | 0.4278 | 0.4809 | 420 | 547 | 345 |
| `readlik` | 0.5359 | 0.4926 | **0.5133** | 410 | 410 | 355 |

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
| **floor 0.02, cap 0.7 (current defaults)** | 0.9699 | 0.9842 | 0.9103 | 0.9312 | 0.9158 |

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
| readlik-support | GT | ALL | ALL | ALL | 94691 | 90782 | 3909 | 94983 | 91167 | 3816 | 0.9587183576052635 | 0.9598243896276176 | 0.959271054804979 | 774 | 599 |
| readlik-support | GT | ALL | ALL | Snv | 75017 | 73000 | 2017 | 71504 | 70892 | 612 | 0.9731127611074823 | 0.9914410382635936 | 0.982191402983025 | 102 | 332 |
| readlik-support | GT | ALL | ALL | Insertion | 9709 | 8679 | 1030 | 10952 | 9476 | 1476 | 0.8939128643526625 | 0.8652300949598247 | 0.8793376438398628 | 354 | 105 |
| readlik-support | GT | ALL | ALL | Deletion | 9965 | 9103 | 862 | 11441 | 9860 | 1581 | 0.9134972403411942 | 0.861812778603269 | 0.886902666626015 | 318 | 156 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1086 | 939 | 147 |  | 0.8646408839779005 |  | 0 | 6 |
| readlik-support | GT | ALL | ALL | JointIndel | 19674 | 17782 | 1892 | 23479 | 20275 | 3204 | 0.9038324692487547 | 0.8635376293709273 | 0.8832257017962536 | 672 | 267 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 390682 | 375049 | 15633 | 460642 | 375049 | 85593 | 0.959985358936424 | 0.8141875903630151 | 0.8810957990142413 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 200440 | 196027 | 4413 | 192992 | 191045 | 1947 | 0.9779834364398323 | 0.9899114989222351 | 0.9839113177138715 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 95512 | 87298 | 8214 | 149914 | 87151 | 62763 | 0.9140003350364352 | 0.5813399682484625 | 0.7106675645429072 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 99134 | 92250 | 6884 | 124050 | 91636 | 32414 | 0.9305586378033772 | 0.7387021362353889 | 0.823604872682025 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 13950 | 7324 | 6626 |  | 0.5250179211469534 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 194646 | 179548 | 15098 | 287914 | 186111 | 101803 | 0.9224335460271467 | 0.6464117757385887 | 0.7601411027788858 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 91363 | 3328 | 96459 | 91660 | 4799 | 0.9648541044027416 | 0.9502482920204439 | 0.9574955015147099 | 398 | 258 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 73221 | 1796 | 73782 | 71355 | 2427 | 0.9760587600143967 | 0.9671057981621534 | 0.9715616540913877 | 111 | 64 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 8850 | 859 | 10596 | 9467 | 1129 | 0.911525388814502 | 0.8934503586258966 | 0.9023973720287548 | 162 | 81 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 9292 | 673 | 10780 | 9893 | 887 | 0.9324636226793778 | 0.9177179962894249 | 0.9250320494428458 | 125 | 106 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1301 | 945 | 356 |  | 0.7263643351268255 |  | 0 | 7 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 18142 | 1532 | 22677 | 20305 | 2372 | 0.9221307309138965 | 0.8954006261851215 | 0.908569121143217 | 287 | 194 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390638 | 376751 | 13887 | 473242 | 376751 | 96491 | 0.9644504630885884 | 0.7961064318044467 | 0.8722299393434274 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 196431 | 4009 | 196914 | 191306 | 5608 | 0.9799990021951707 | 0.9715205622759174 | 0.9757413648072043 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 87718 | 7794 | 162624 | 87340 | 75284 | 0.918397688248597 | 0.5370670995670995 | 0.6777782421202233 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 92936 | 6198 | 116758 | 92432 | 24326 | 0.9374785643674218 | 0.7916545333082101 | 0.8584176155767088 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 17544 | 7445 | 10099 |  | 0.42436160510715915 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 180654 | 13992 | 296926 | 187217 | 109709 | 0.9281156561141765 | 0.630517367963735 | 0.7509054814301092 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 94691 | 90835 | 3856 | 95004 | 91198 | 3806 | 0.959278072889715 | 0.9599385289040462 | 0.959608187256204 | 771 | 592 |
| readlik-nolink | GT | ALL | ALL | Snv | 75017 | 73029 | 1988 | 71517 | 70896 | 621 | 0.9734993401495661 | 0.9913167498636688 | 0.9823272587970241 | 99 | 330 |
| readlik-nolink | GT | ALL | ALL | Insertion | 9709 | 8691 | 1018 | 10966 | 9486 | 1480 | 0.8951488309815635 | 0.8650373882910816 | 0.8798355519498191 | 357 | 102 |
| readlik-nolink | GT | ALL | ALL | Deletion | 9965 | 9115 | 850 | 11449 | 9872 | 1577 | 0.9147014550928249 | 0.8622587125513145 | 0.8877062225686344 | 315 | 153 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1072 | 944 | 128 |  | 0.8805970149253731 |  | 0 | 7 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 19674 | 17806 | 1868 | 23487 | 20302 | 3185 | 0.9050523533597642 | 0.8643930685059821 | 0.8842555652881441 | 672 | 262 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 390636 | 375207 | 15429 | 457652 | 375207 | 82445 | 0.9605028722391178 | 0.8198522021099001 | 0.8846217322418802 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 200440 | 196053 | 4387 | 193018 | 191053 | 1965 | 0.9781131510676512 | 0.9898196023168824 | 0.9839315581750652 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 95512 | 87479 | 8033 | 152584 | 87151 | 65433 | 0.915895384873105 | 0.5711673569968018 | 0.7035742763694433 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 99134 | 92159 | 6975 | 120716 | 91715 | 29001 | 0.9296406883612081 | 0.7597584413002418 | 0.8361580729595095 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 10828 | 7230 | 3598 |  | 0.6677133357960843 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 194646 | 179638 | 15008 | 284128 | 186096 | 98032 | 0.92289592388233 | 0.6549724068025679 | 0.7661873335541365 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 94691 | 91136 | 3555 | 93381 | 91280 | 2101 | 0.9624568332787699 | 0.9775007763892012 | 0.9699204736046891 | 341 | 345 |
| readlik | GT | ALL | ALL | Snv | 75017 | 73002 | 2015 | 71375 | 71051 | 324 | 0.9731394217310743 | 0.9954605954465849 | 0.9841734631271587 | 76 | 125 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 8854 | 855 | 10388 | 9440 | 948 | 0.9119373776908023 | 0.908740854832499 | 0.9103363102308537 | 147 | 86 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 9280 | 685 | 10583 | 9855 | 728 | 0.9312594079277471 | 0.9312104318246244 | 0.93123491923224 | 118 | 125 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1035 | 934 | 101 |  | 0.9024154589371981 |  | 0 | 9 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 18134 | 1540 | 22006 | 20229 | 1777 | 0.9217241028768933 | 0.9192492956466418 | 0.9204850358290337 | 265 | 220 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390636 | 376172 | 14464 | 430908 | 376172 | 54736 | 0.9629732026746127 | 0.872975205844403 | 0.9157683581159379 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 195989 | 4451 | 192286 | 190891 | 1395 | 0.977793853522251 | 0.9927451816564908 | 0.9852127965072237 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 87531 | 7981 | 134780 | 87238 | 47542 | 0.916439819080325 | 0.6472622050749369 | 0.7586827272118567 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 93004 | 6130 | 110384 | 92322 | 18062 | 0.938164504609922 | 0.8363712132193072 | 0.8843482574470233 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 10540 | 7197 | 3343 |  | 0.6828273244781784 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 180535 | 14111 | 255704 | 186757 | 68947 | 0.9275042898389898 | 0.7303640146419297 | 0.8172129895884779 |  |  |

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

