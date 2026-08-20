# Tier 2 results: HG002 chr6 on HPRC v2.1 MC CHM13, 34-haplotype graph

Real reads, real benchmark, run on a 32 GB laptop.

This is the **34-haplotype** graph: CHM13, GRCh38 and 32 recombinants from haplotype sampling. It is the primary subject because it is what the caller is tuned for -- both the linkage transition and the panel frequency prior are panel-size effects and have little to work with on a thin panel -- and because it is the better-performing configuration. The 4-haplotype graph has its own page at [tier2-chr6-4hap-results.md](tier2-chr6-4hap-results.md), and the two are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md).

| | |
|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz`, 101,366,693 nodes, **34 haplotypes** (CHM13, GRCh38, 32 recombinants from haplotype sampling; the file is named for the recombinant count, not the total). HG002 itself is **absent** — no circularity |
| chromosome | chr6 component, 5,499,123 nodes |
| reads | 596,017,764 alignments genome-wide (~28.6×); 151 bp paired Illumina |
| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |
| regions | small variants 167.2 Mb; SVs 168.4 Mb |
| engine | `aardvark compare` for small variants; `truvari bench --sizemin 50` for SVs |

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.7`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap down from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every arm on this page was re-run together on one build, so the wall-clock column compares runs made on the same machine in the same session rather than a mixture of vintages.

Two changes since the accuracy results were first produced left the calls untouched. The read path was optimised (vg `44fd008`). Then `AD`, `BL`, `GQI` and the explained-share scaling of `GQ` were added — which rescales a quality and does not change a genotype. Both are confirmed by the variant counts below, which are unchanged to the record.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 294,626 | 603 s | 4.8 GB |
| `poisson-z` | panel (`-z`) | yes | 294,835 | 214 s | 5.2 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 298,098 | 338 s | 9.0 GB |
| `readlik-nomismap` | panel (default) | **no** | 305,480 | 339 s | 8.8 GB |
| `readlik-nolink` | panel (default) | **no** | 297,938 | 272 s | 7.7 GB |
| `readlik` | panel (default) | **no** | 297,484 | 376 s | 5.3 GB |

**Peak RSS in this table is repeatable to about ±0.35 GB, so read it accordingly.** Three
back-to-back runs of one binary on chr6-4hap, identical parameters and a warm cache, gave 7.3, 6.6
and 7.0 GB -- a 0.7 GB spread on a 7 GB measurement. Differences smaller than that are not evidence
of anything, and a single measurement of each of two arms cannot resolve one. Thread count matters
too: the same run at `--threads 6` instead of 5 measured 8.7 GB, because the read and GBWT caches
are per thread. Wall clock is worse still -- a run immediately after a full rebuild took 956 s
against 260 s warm, purely from page cache.

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9573 | 0.9800 | 0.9685 | 214,033 | 9,552 | 4,235 | 0.9639 | 0.9825 | 0.9731 |
| `poisson-z` | 0.9607 | 0.9801 | 0.9703 | 214,796 | 8,789 | 4,230 | 0.9671 | 0.9824 | 0.9747 |
| `readlik-support` | 0.9786 | 0.9933 | 0.9859 | 218,793 | 4,792 | 1,422 | 0.9821 | 0.9914 | 0.9868 |
| `readlik-nomismap` | 0.9791 | 0.9855 | 0.9823 | 218,910 | 4,675 | 3,132 | 0.9825 | 0.9856 | 0.9841 |
| `readlik-nolink` | 0.9787 | 0.9933 | 0.9860 | 218,816 | 4,769 | 1,421 | 0.9821 | 0.9914 | 0.9868 |
| `readlik` | 0.9781 | 0.9961 | **0.9870** | 218,691 | 4,894 | 821 | 0.9817 | 0.9935 | 0.9876 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7744 | 0.8529 | 0.8118 | 21,742 | 6,333 | 4,279 | 0.7996 | 0.6375 | 0.7094 |
| `poisson-z` | 0.7792 | 0.8558 | 0.8157 | 21,876 | 6,199 | 4,192 | 0.8041 | 0.5959 | 0.6845 |
| `readlik-support` | 0.9021 | 0.8871 | 0.8945 | 25,327 | 2,748 | 3,491 | 0.9218 | 0.5414 | 0.6821 |
| `readlik-nomismap` | 0.9213 | 0.9155 | 0.9184 | 25,866 | 2,209 | 2,534 | 0.9327 | 0.4495 | 0.6067 |
| `readlik-nolink` | 0.9030 | 0.8870 | 0.8949 | 25,352 | 2,723 | 3,493 | 0.9230 | 0.5401 | 0.6815 |
| `readlik` | 0.9202 | 0.9226 | **0.9214** | 25,835 | 2,240 | 2,292 | 0.9291 | 0.5903 | 0.7219 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8557 | 0.7206 | 0.7824 | 24,222 | 4,083 | 10,150 | 0.9292 | 0.5149 | 0.6626 |
| `poisson-z` | 0.8593 | 0.7229 | 0.7852 | 24,323 | 3,982 | 10,083 | 0.9306 | 0.4736 | 0.6278 |
| `readlik-support` | 0.9299 | 0.8896 | 0.9093 | 26,322 | 1,983 | 3,515 | 0.9431 | 0.7556 | 0.8390 |
| `readlik-nomismap` | 0.9416 | 0.9376 | 0.9396 | 26,651 | 1,654 | 1,885 | 0.9480 | 0.8398 | 0.8906 |
| `readlik-nolink` | 0.9306 | 0.8899 | 0.9098 | 26,342 | 1,963 | 3,509 | 0.9439 | 0.7584 | 0.8410 |
| `readlik` | 0.9422 | 0.9423 | **0.9422** | 26,668 | 1,637 | 1,730 | 0.9500 | 0.8194 | 0.8799 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.7865 | — | 0 | 0 | 825 | — | 0.3420 | — |
| `poisson-z` | — | 0.7915 | — | 0 | 0 | 836 | — | 0.4509 | — |
| `readlik-support` | — | 0.9002 | — | 0 | 0 | 341 | — | 0.5817 | — |
| `readlik-nomismap` | — | 0.8658 | — | 0 | 0 | 487 | — | 0.4986 | — |
| `readlik-nolink` | — | 0.8994 | — | 0 | 0 | 346 | — | 0.7447 | — |
| `readlik` | — | 0.9144 | — | 0 | 0 | 288 | — | 0.7619 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8153 | 0.7798 | 0.7972 | 45,964 | 10,416 | 15,254 | 0.8656 | 0.5453 | 0.6691 |
| `poisson-z` | 0.8194 | 0.7825 | 0.8005 | 46,199 | 10,181 | 15,111 | 0.8685 | 0.5176 | 0.6486 |
| `readlik-support` | 0.9161 | 0.8890 | 0.9023 | 51,649 | 4,731 | 7,347 | 0.9326 | 0.6314 | 0.7530 |
| `readlik-nomismap` | 0.9315 | 0.9232 | 0.9273 | 52,517 | 3,863 | 4,906 | 0.9405 | 0.5857 | 0.7218 |
| `readlik-nolink` | 0.9169 | 0.8890 | 0.9027 | 51,694 | 4,686 | 7,348 | 0.9337 | 0.6381 | 0.7581 |
| `readlik` | 0.9312 | 0.9315 | **0.9314** | 52,503 | 3,877 | 4,310 | 0.9397 | 0.6924 | 0.7973 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9287 | 0.9308 | 0.9297 | 259,997 | 19,968 | 19,489 | 0.9264 | 0.7386 | 0.8219 |
| `poisson-z` | 0.9322 | 0.9313 | 0.9318 | 260,995 | 18,970 | 19,341 | 0.9296 | 0.7139 | 0.8076 |
| `readlik-support` | 0.9660 | 0.9686 | 0.9673 | 270,442 | 9,523 | 8,769 | 0.9682 | 0.8013 | 0.8769 |
| `readlik-nomismap` | 0.9695 | 0.9713 | 0.9704 | 271,427 | 8,538 | 8,038 | 0.9722 | 0.7630 | 0.8550 |
| `readlik-nolink` | 0.9662 | 0.9686 | 0.9674 | 270,510 | 9,455 | 8,769 | 0.9687 | 0.8063 | 0.8801 |
| `readlik` | 0.9687 | 0.9814 | **0.9750** | 271,194 | 8,771 | 5,131 | 0.9714 | 0.8441 | 0.9033 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (167.2 Mb vs 168.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

The same mechanism was traced site by site on chr20, where 246 `readlik` calls carrying a >=200 bp insertion allele contribute 27,951 FP bases and zero TP bases — the whole of the precision difference there. The size-matched control below is the general test.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.8019 | 0.8645 | **0.8320** |
| `sm50-poisson-z` | Deletion | 0.9239 | 0.7428 | **0.8235** |
| `sm50-poisson-z` | ALL | 0.9271 | 0.9091 | **0.9180** |
| `sm50-readlik` | Insertion | 0.9267 | 0.9024 | **0.9144** |
| `sm50-readlik` | Deletion | 0.9480 | 0.9078 | **0.9275** |
| `sm50-readlik` | ALL | 0.9702 | 0.9725 | **0.9714** |

The insertion BASEPAIR precision gap collapses from **0.006 to -0.038**, and insertion BASEPAIR F1 goes from 0.8320 for `poisson-z` against 0.9144 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5320 | 0.4618 | 0.4944 | 823 | 936 | 724 |
| `poisson-z` | 0.5417 | 0.4442 | 0.4881 | 838 | 1,036 | 709 |
| `readlik-support` | 0.5740 | 0.5057 | 0.5377 | 888 | 860 | 659 |
| `readlik-nomismap` | 0.5992 | 0.4657 | 0.5241 | 927 | 1,066 | 620 |
| `readlik-nolink` | 0.5941 | 0.4981 | 0.5419 | 919 | 928 | 628 |
| `readlik` | 0.5856 | 0.5499 | **0.5672** | 906 | 731 | 641 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and the old cap of 0.1 was telling the model 0.1: overriding the mapper at exactly the sites that matter. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless — see 

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9750 | 0.9870 | 0.9214 | 0.9422 | 0.9033 |

Only the current row is available here: the preserved old-default arms (`arms.floor-1e-8.json`, `arms.readlik.json`) exist for the 4-haplotype runs alone, so the before-and-after is on [tier2-chr6-4hap-results.md](tier2-chr6-4hap-results.md). Mixing rows from two graphs into one table is exactly what the one-build-per-matrix rule forbids. The full grids are in plan §9.20-§9.21.

The floor was later re-swept at the corrected cap, on both graphs and both benchmarks, and settled at **0.02**. 0.05 wins on small-variant `GT` but costs about 0.01 of SV F1 — which the first sweep never saw, because it was scored on one benchmark only. Plan §9.21 records that as a process rule: a sweep that sets a default has to be scored on every benchmark the project runs.

## Known bad output (measured on chr20)

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
| poisson | GT | ALL | ALL | ALL | 279965 | 259997 | 19968 | 281484 | 261995 | 19489 | 0.928676798885575 | 0.9307633826434184 | 0.9297189200272442 | 4157 | 1701 |
| poisson | GT | ALL | ALL | Snv | 223585 | 214033 | 9552 | 212198 | 207963 | 4235 | 0.9572779927097077 | 0.9800422247146533 | 0.9685263646222478 | 670 | 431 |
| poisson | GT | ALL | ALL | Insertion | 28075 | 21742 | 6333 | 29096 | 24817 | 4279 | 0.7744256455921639 | 0.8529351113555128 | 0.8117865954917517 | 2157 | 370 |
| poisson | GT | ALL | ALL | Deletion | 28305 | 24222 | 4083 | 36326 | 26176 | 10150 | 0.8557498675145734 | 0.7205858063095304 | 0.7823729660146876 | 1330 | 819 |
| poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3864 | 3039 | 825 |  | 0.7864906832298136 |  | 0 | 81 |
| poisson | GT | ALL | ALL | JointIndel | 56380 | 45964 | 10416 | 69286 | 54032 | 15254 | 0.8152536360411493 | 0.7798400831336778 | 0.7971537417052307 | 3487 | 1270 |
| poisson | BASEPAIR | ALL | ALL | ALL | 1115456 | 1033389 | 82067 | 1399186 | 1033389 | 365797 | 0.9264273983016811 | 0.7385644224570572 | 0.8218975106595691 |  |  |
| poisson | BASEPAIR | ALL | ALL | Snv | 602392 | 580625 | 21767 | 572210 | 562179 | 10031 | 0.9638657219883399 | 0.9824697226542702 | 0.9730788093742253 |  |  |
| poisson | BASEPAIR | ALL | ALL | Insertion | 257406 | 205833 | 51573 | 322190 | 205394 | 116796 | 0.7996433649565278 | 0.6374934045128651 | 0.7094208177701448 |  |  |
| poisson | BASEPAIR | ALL | ALL | Deletion | 267064 | 248155 | 18909 | 475064 | 244609 | 230455 | 0.9291967468471977 | 0.5148969402017413 | 0.6626170671443703 |  |  |
| poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 74880 | 25611 | 49269 |  | 0.34202724358974357 |  |  |  |
| poisson | BASEPAIR | ALL | ALL | JointIndel | 524470 | 453988 | 70482 | 872134 | 475614 | 396520 | 0.865612904455927 | 0.5453450960517535 | 0.6691308350132357 |  |  |

</details>

<details><summary><code>poisson-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-z | GT | ALL | ALL | ALL | 279965 | 260995 | 18970 | 281690 | 262349 | 19341 | 0.9322415301912739 | 0.9313394156697078 | 0.931790254584566 | 4315 | 1706 |
| poisson-z | GT | ALL | ALL | Snv | 223585 | 214796 | 8789 | 212217 | 207987 | 4230 | 0.9606905651094662 | 0.9800675723434032 | 0.9702823363201399 | 778 | 434 |
| poisson-z | GT | ALL | ALL | Insertion | 28075 | 21876 | 6199 | 29078 | 24886 | 4192 | 0.7791985752448798 | 0.8558360272370864 | 0.8157212233356769 | 2183 | 352 |
| poisson-z | GT | ALL | ALL | Deletion | 28305 | 24323 | 3982 | 36385 | 26302 | 10083 | 0.8593181416710829 | 0.7228803078191562 | 0.7852164979252937 | 1354 | 786 |
| poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 4010 | 3174 | 836 |  | 0.7915211970074812 |  | 0 | 134 |
| poisson-z | GT | ALL | ALL | JointIndel | 56380 | 46199 | 10181 | 69473 | 54362 | 15111 | 0.8194217807733238 | 0.7824910396844817 | 0.8005307067761451 | 3537 | 1272 |
| poisson-z | BASEPAIR | ALL | ALL | ALL | 1115610 | 1037072 | 78538 | 1452654 | 1037072 | 415582 | 0.9296008461738421 | 0.7139153576832473 | 0.807605448661041 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Snv | 602392 | 582558 | 19834 | 572258 | 562194 | 10064 | 0.9670745959441692 | 0.9824135267658993 | 0.9746837165917057 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Insertion | 257406 | 206986 | 50420 | 344826 | 205467 | 139359 | 0.8041226700232318 | 0.5958570409423883 | 0.6844987123195558 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Deletion | 267064 | 248540 | 18524 | 517298 | 244998 | 272300 | 0.9306383488601983 | 0.47361095538741693 | 0.6277525168652335 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 63774 | 28756 | 35018 |  | 0.4509047574246558 |  |  |  |
| poisson-z | BASEPAIR | ALL | ALL | JointIndel | 524470 | 455526 | 68944 | 925898 | 479221 | 446677 | 0.868545388678094 | 0.5175742900405876 | 0.6486261898087171 |  |  |

</details>

<details><summary><code>readlik-support</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-support | GT | ALL | ALL | ALL | 279965 | 270442 | 9523 | 279622 | 270853 | 8769 | 0.9659850338435162 | 0.9686398065960475 | 0.9673105987249317 | 1803 | 1501 |
| readlik-support | GT | ALL | ALL | Snv | 223585 | 218793 | 4792 | 213446 | 212024 | 1422 | 0.9785674352036139 | 0.9933378934250349 | 0.9858973456251097 | 200 | 773 |
| readlik-support | GT | ALL | ALL | Insertion | 28075 | 25327 | 2748 | 30912 | 27421 | 3491 | 0.9021193232413179 | 0.8870665113871635 | 0.8945295960145939 | 926 | 355 |
| readlik-support | GT | ALL | ALL | Deletion | 28305 | 26322 | 1983 | 31847 | 28332 | 3515 | 0.9299417064122947 | 0.8896285364398531 | 0.9093385457361243 | 677 | 335 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3417 | 3076 | 341 |  | 0.9002048580626281 |  | 0 | 38 |
| readlik-support | GT | ALL | ALL | JointIndel | 56380 | 51649 | 4731 | 66176 | 58829 | 7347 | 0.9160872649875843 | 0.8889778771760155 | 0.9023289997839138 | 1603 | 728 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 1115456 | 1079940 | 35516 | 1347758 | 1079940 | 267818 | 0.9681601067186872 | 0.801286284332944 | 0.876854386180007 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 602392 | 591629 | 10763 | 580120 | 575159 | 4961 | 0.9821328968512198 | 0.9914483210370268 | 0.9867686242578648 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 257406 | 237277 | 20129 | 438304 | 237279 | 201025 | 0.9218005796290685 | 0.5413571402496897 | 0.6821182964606843 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 267064 | 251865 | 15199 | 332400 | 251159 | 81241 | 0.9430885480633856 | 0.75559265944645 | 0.838992956387997 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 36456 | 21206 | 15250 |  | 0.5816875137151635 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 524470 | 489142 | 35328 | 807160 | 509644 | 297516 | 0.9326405704806757 | 0.6314039347836861 | 0.7530130043721055 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 271427 | 8538 | 279921 | 271883 | 8038 | 0.9695033307734895 | 0.9712847553416857 | 0.9703932254841712 | 1169 | 744 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 218910 | 4675 | 216080 | 212948 | 3132 | 0.9790907261220565 | 0.9855053683820807 | 0.9822875749632762 | 458 | 164 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 25866 | 2209 | 30002 | 27468 | 2534 | 0.9213178984861977 | 0.9155389640690621 | 0.9184193407263523 | 374 | 256 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 26651 | 1654 | 30209 | 28324 | 1885 | 0.941565094506271 | 0.9376013770730577 | 0.939579055463938 | 337 | 302 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3630 | 3143 | 487 |  | 0.865840220385675 |  | 0 | 22 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 52517 | 3863 | 63841 | 58935 | 4906 | 0.9314827953174885 | 0.9231528328190348 | 0.9272991073546158 | 711 | 580 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115518 | 1084519 | 30999 | 1421342 | 1084519 | 336823 | 0.9722111162706474 | 0.7630246626075919 | 0.8550089480696609 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 591855 | 10537 | 583462 | 575060 | 8402 | 0.9825080678362262 | 0.9855997477127902 | 0.9840514794308006 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 240081 | 17325 | 531434 | 238903 | 292531 | 0.9326938765996131 | 0.44954406379719775 | 0.6066784644111041 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 253165 | 13899 | 299940 | 251904 | 48036 | 0.9479562951202708 | 0.8398479695939188 | 0.8906334831322368 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 44734 | 22306 | 22428 |  | 0.4986363839585103 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 493246 | 31224 | 876108 | 513113 | 362995 | 0.940465612904456 | 0.5856732275016322 | 0.7218288615438032 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 279965 | 270510 | 9455 | 279664 | 270895 | 8769 | 0.9662279213473113 | 0.9686445162766749 | 0.9674347096863941 | 1806 | 1513 |
| readlik-nolink | GT | ALL | ALL | Snv | 223585 | 218816 | 4769 | 213448 | 212027 | 1421 | 0.9786703043585213 | 0.9933426408305536 | 0.9859518894189874 | 200 | 774 |
| readlik-nolink | GT | ALL | ALL | Insertion | 28075 | 25352 | 2723 | 30918 | 27425 | 3493 | 0.9030097951914514 | 0.8870237402160553 | 0.8949453852551508 | 927 | 364 |
| readlik-nolink | GT | ALL | ALL | Deletion | 28305 | 26342 | 1963 | 31859 | 28350 | 3509 | 0.9306482953541777 | 0.8898584387457233 | 0.9097964029610589 | 679 | 334 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3439 | 3093 | 346 |  | 0.8993893573713289 |  | 0 | 41 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 56380 | 51694 | 4686 | 66216 | 58868 | 7348 | 0.9168854203618304 | 0.8890298417300955 | 0.9027428000189495 | 1606 | 739 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 1115480 | 1080589 | 34891 | 1340176 | 1080589 | 259587 | 0.968721088679313 | 0.8063037989040246 | 0.8800817378329864 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 591634 | 10758 | 580124 | 575163 | 4961 | 0.982141197094251 | 0.9914483800015169 | 0.9867728428367101 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 237593 | 19813 | 438878 | 237048 | 201830 | 0.923028212240585 | 0.5401227676028418 | 0.6814724651645153 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 252089 | 14975 | 331464 | 251372 | 80092 | 0.943927298325495 | 0.7583689329761302 | 0.8410347445706214 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 29546 | 22002 | 7544 |  | 0.7446693291816151 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 489682 | 34788 | 799888 | 510422 | 289466 | 0.93367018132591 | 0.63811683635709 | 0.7581060991162865 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 271194 | 8771 | 276063 | 270932 | 5131 | 0.9686710838854857 | 0.9814136628233411 | 0.975000740979651 | 790 | 1019 |
| readlik | GT | ALL | ALL | Snv | 223585 | 218691 | 4894 | 213106 | 212285 | 821 | 0.9781112328644587 | 0.9961474571340084 | 0.9870469582817407 | 173 | 358 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 25835 | 2240 | 29620 | 27328 | 2292 | 0.920213713268032 | 0.9226198514517218 | 0.9214152115449455 | 336 | 300 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 26668 | 1637 | 29971 | 28241 | 1730 | 0.9421656951068715 | 0.9422775349504521 | 0.9422216117098697 | 281 | 323 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3366 | 3078 | 288 |  | 0.9144385026737968 |  | 0 | 38 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 52503 | 3877 | 62957 | 58647 | 4310 | 0.9312344803121675 | 0.9315405753133091 | 0.9313875026636604 | 617 | 661 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115540 | 1083658 | 31882 | 1283790 | 1083658 | 200132 | 0.9714201194040555 | 0.8441084601064037 | 0.9033005047242356 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 591386 | 11006 | 578420 | 574669 | 3751 | 0.9817295050399075 | 0.9935150928391134 | 0.9875871387168419 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 239149 | 18257 | 403284 | 238058 | 165226 | 0.9290731373783051 | 0.5902986480991064 | 0.7219176007105896 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 253709 | 13355 | 308134 | 252481 | 55653 | 0.9499932600425366 | 0.8193870199328863 | 0.8798698110429873 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 28984 | 22083 | 6901 |  | 0.761903118962186 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 492858 | 31612 | 740402 | 512622 | 227780 | 0.9397258184452876 | 0.6923563145426403 | 0.7972945615772612 |  |  |

</details>

