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

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.5`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap down from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every arm on this page was re-run together on one build, so the wall-clock column compares runs made on the same machine in the same session rather than a mixture of vintages.

Two changes since the accuracy results were first produced left the calls untouched. The read path was optimised (vg `44fd008`). Then `AD`, `BL`, `GQI` and the explained-share scaling of `GQ` were added — which rescales a quality and does not change a genotype. Both are confirmed by the variant counts below, which are unchanged to the record.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 294,626 | 605 s | 4.4 GB |
| `poisson-z` | haplotype (`-z`) | yes | 294,835 | 226 s | 4.1 GB |
| `readlik` | support (Flow) | yes | 284,466 | 406 s | 4.7 GB |
| `readlik-nomismap` | support (Flow) | yes | 290,599 | 384 s | 5.4 GB |
| `readlik-z-nolink` | haplotype (`-z`) | **no** | 284,529 | 267 s | 6.5 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 284,529 | 338 s | 4.7 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9573 | 0.9800 | 0.9685 | 214,033 | 9,552 | 4,235 | 0.9639 | 0.9825 | 0.9731 |
| `poisson-z` | 0.9607 | 0.9801 | 0.9703 | 214,796 | 8,789 | 4,230 | 0.9671 | 0.9824 | 0.9747 |
| `readlik` | 0.9600 | 0.9948 | 0.9771 | 214,637 | 8,948 | 1,089 | 0.9659 | 0.9924 | 0.9790 |
| `readlik-nomismap` | 0.9611 | 0.9844 | 0.9726 | 214,897 | 8,688 | 3,298 | 0.9668 | 0.9849 | 0.9758 |
| `readlik-z-nolink` | 0.9644 | 0.9949 | 0.9794 | 215,626 | 7,959 | 1,073 | 0.9682 | 0.9924 | 0.9802 |
| `readlik-z` | 0.9643 | 0.9972 | **0.9804** | 215,594 | 7,991 | 593 | 0.9679 | 0.9940 | 0.9808 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7744 | 0.8529 | 0.8118 | 21,742 | 6,333 | 4,279 | 0.7996 | 0.6375 | 0.7094 |
| `poisson-z` | 0.7792 | 0.8558 | 0.8157 | 21,876 | 6,199 | 4,192 | 0.8041 | 0.5959 | 0.6845 |
| `readlik` | 0.8845 | 0.8926 | 0.8885 | 24,832 | 3,243 | 3,227 | 0.9070 | 0.6072 | 0.7274 |
| `readlik-nomismap` | 0.8847 | 0.8860 | 0.8854 | 24,839 | 3,236 | 3,464 | 0.9092 | 0.4757 | 0.6246 |
| `readlik-z-nolink` | 0.8900 | 0.8926 | 0.8913 | 24,986 | 3,089 | 3,229 | 0.9143 | 0.6056 | 0.7286 |
| `readlik-z` | 0.9075 | 0.9263 | **0.9168** | 25,477 | 2,598 | 2,125 | 0.9202 | 0.6436 | 0.7575 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8557 | 0.7206 | 0.7824 | 24,222 | 4,083 | 10,150 | 0.9292 | 0.5149 | 0.6626 |
| `poisson-z` | 0.8593 | 0.7229 | 0.7852 | 24,323 | 3,982 | 10,083 | 0.9306 | 0.4736 | 0.6278 |
| `readlik` | 0.9130 | 0.8947 | 0.9038 | 25,842 | 2,463 | 3,264 | 0.9271 | 0.7640 | 0.8377 |
| `readlik-nomismap` | 0.9126 | 0.8907 | 0.9015 | 25,832 | 2,473 | 3,413 | 0.9272 | 0.7941 | 0.8555 |
| `readlik-z-nolink` | 0.9181 | 0.8953 | 0.9066 | 25,988 | 2,317 | 3,248 | 0.9307 | 0.7675 | 0.8413 |
| `readlik-z` | 0.9301 | 0.9490 | **0.9394** | 26,326 | 1,979 | 1,488 | 0.9360 | 0.8330 | 0.8815 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.7865 | — | 0 | 0 | 825 | — | 0.3420 | — |
| `poisson-z` | — | 0.7915 | — | 0 | 0 | 836 | — | 0.4509 | — |
| `readlik` | — | 0.8584 | — | 0 | 0 | 546 | — | 0.5846 | — |
| `readlik-nomismap` | — | 0.8200 | — | 0 | 0 | 738 | — | 0.4954 | — |
| `readlik-z-nolink` | — | 0.8953 | — | 0 | 0 | 410 | — | 0.7768 | — |
| `readlik-z` | — | 0.9117 | — | 0 | 0 | 338 | — | 0.8294 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8153 | 0.7798 | 0.7972 | 45,964 | 10,416 | 15,254 | 0.8656 | 0.5453 | 0.6691 |
| `poisson-z` | 0.8194 | 0.7825 | 0.8005 | 46,199 | 10,181 | 15,111 | 0.8685 | 0.5176 | 0.6486 |
| `readlik` | 0.8988 | 0.8916 | 0.8952 | 50,674 | 5,706 | 7,037 | 0.9172 | 0.6726 | 0.7761 |
| `readlik-nomismap` | 0.8987 | 0.8841 | 0.8914 | 50,671 | 5,709 | 7,615 | 0.9184 | 0.5917 | 0.7197 |
| `readlik-z-nolink` | 0.9041 | 0.8940 | 0.8990 | 50,974 | 5,406 | 6,887 | 0.9227 | 0.6846 | 0.7860 |
| `readlik-z` | 0.9188 | 0.9361 | **0.9274** | 51,803 | 4,577 | 3,951 | 0.9283 | 0.7343 | 0.8200 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9287 | 0.9308 | 0.9297 | 259,997 | 19,968 | 19,489 | 0.9264 | 0.7386 | 0.8219 |
| `poisson-z` | 0.9322 | 0.9313 | 0.9318 | 260,995 | 18,970 | 19,341 | 0.9296 | 0.7139 | 0.8076 |
| `readlik` | 0.9477 | 0.9703 | 0.9588 | 265,311 | 14,654 | 8,126 | 0.9521 | 0.8315 | 0.8877 |
| `readlik-nomismap` | 0.9486 | 0.9606 | 0.9546 | 265,568 | 14,397 | 10,913 | 0.9532 | 0.7689 | 0.8512 |
| `readlik-z-nolink` | 0.9523 | 0.9709 | 0.9615 | 266,600 | 13,365 | 7,960 | 0.9561 | 0.8396 | 0.8941 |
| `readlik-z` | 0.9551 | 0.9832 | **0.9689** | 267,397 | 12,568 | 4,544 | 0.9586 | 0.8724 | 0.9135 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (167.2 Mb vs 168.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

The same mechanism was traced site by site on chr20, where 246 `readlik-z` calls carrying a >=200 bp insertion allele contribute 27,951 FP bases and zero TP bases — the whole of the precision difference there. The size-matched control below is the general test.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.8019 | 0.8645 | **0.8320** |
| `sm50-poisson-z` | Deletion | 0.9239 | 0.7428 | **0.8235** |
| `sm50-poisson-z` | ALL | 0.9271 | 0.9091 | **0.9180** |
| `sm50-readlik-z` | Insertion | 0.9168 | 0.9090 | **0.9129** |
| `sm50-readlik-z` | Deletion | 0.9340 | 0.9172 | **0.9255** |
| `sm50-readlik-z` | ALL | 0.9572 | 0.9774 | **0.9672** |

The insertion BASEPAIR precision gap collapses from **-0.048 to -0.044**, and insertion BASEPAIR F1 goes from 0.8320 for `poisson-z` against 0.9129 for `readlik-z`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5320 | 0.4618 | 0.4944 | 823 | 936 | 724 |
| `poisson-z` | 0.5417 | 0.4442 | 0.4881 | 838 | 1,036 | 709 |
| `readlik` | 0.5255 | 0.4766 | 0.4999 | 813 | 885 | 734 |
| `readlik-nomismap` | 0.5307 | 0.4319 | 0.4762 | 821 | 1,080 | 726 |
| `readlik-z-nolink` | 0.5456 | 0.4716 | 0.5059 | 844 | 949 | 703 |
| `readlik-z` | 0.5352 | 0.5187 | **0.5268** | 828 | 760 | 719 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and the old cap of 0.1 was telling the model 0.1: overriding the mapper at exactly the sites that matter. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless — see 

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9689 | 0.9804 | 0.9168 | 0.9394 | 0.9135 |

Only the current row is available here: the preserved old-default arms (`arms.floor-1e-8.json`, `arms.readlik-z.json`) exist for the 4-haplotype runs alone, so the before-and-after is on [tier2-chr6-4hap-results.md](tier2-chr6-4hap-results.md). Mixing rows from two graphs into one table is exactly what the one-build-per-matrix rule forbids. The full grids are in plan §9.20-§9.21.

The floor was later re-swept at the corrected cap, on both graphs and both benchmarks, and settled at **0.02**. 0.05 wins on small-variant `GT` but costs about 0.01 of SV F1 — which the first sweep never saw, because it was scored on one benchmark only. Plan §9.21 records that as a process rule: a sweep that sets a default has to be scored on every benchmark the project runs.

## Known bad output (measured on chr20)

Neither benchmark scores these, so they appear in no metric on this page. They are recorded because they are plainly wrong and would mislead anyone reading the VCF.

`readlik-z` emits a small number of enormous homozygous insertions in and around the chr20 pericentromere, at depths that are physically impossible:

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

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 265311 | 14654 | 273653 | 265527 | 8126 | 0.9476577429321522 | 0.9703054598341696 | 0.9588478868572163 | 1992 | 1239 |
| readlik | GT | ALL | ALL | Snv | 223585 | 214637 | 8948 | 208740 | 207651 | 1089 | 0.9599794261690185 | 0.9947829836159816 | 0.9770713750111547 | 399 | 642 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 24832 | 3243 | 30055 | 26828 | 3227 | 0.8844879786286732 | 0.8926301780069872 | 0.8885404258127988 | 921 | 290 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 25842 | 2463 | 31002 | 27738 | 3264 | 0.9129835718071012 | 0.8947164699051674 | 0.9037577248434675 | 672 | 272 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3856 | 3310 | 546 |  | 0.858402489626556 |  | 0 | 35 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 50674 | 5706 | 64913 | 57876 | 7037 | 0.8987938985455836 | 0.8915933634248918 | 0.8951791515183855 | 1593 | 597 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115502 | 1062070 | 53432 | 1277270 | 1062070 | 215200 | 0.9521004892864379 | 0.8315156544818245 | 0.8877318858629237 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 581848 | 20544 | 566394 | 562093 | 4301 | 0.9658959614337508 | 0.9924063461124235 | 0.9789717125057448 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 233479 | 23927 | 380558 | 231071 | 149487 | 0.9070456788109057 | 0.6071899684148014 | 0.7274284396580225 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 247588 | 19476 | 318858 | 243616 | 75242 | 0.9270736602462332 | 0.7640266200001254 | 0.8376900688892183 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 48116 | 28131 | 19985 |  | 0.5846495968077147 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 481067 | 43403 | 747532 | 502818 | 244714 | 0.917244074970923 | 0.6726374255550264 | 0.7761241237034865 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 265568 | 14397 | 277330 | 266417 | 10913 | 0.9485757148214955 | 0.9606497674250892 | 0.954574562619069 | 2188 | 831 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 214897 | 8688 | 211624 | 208326 | 3298 | 0.9611422948766688 | 0.9844157562469286 | 0.9726398228266656 | 534 | 326 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 24839 | 3236 | 30389 | 26925 | 3464 | 0.8847373107747106 | 0.8860113856987726 | 0.8853738898806476 | 951 | 249 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 25832 | 2473 | 31217 | 27804 | 3413 | 0.9126302773361596 | 0.8906685459845597 | 0.9015156796250134 | 703 | 236 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 4100 | 3362 | 738 |  | 0.82 |  | 0 | 20 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 50671 | 5709 | 65706 | 58091 | 7615 | 0.8987406881873005 | 0.8841049523635589 | 0.8913627464366658 | 1654 | 505 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115496 | 1063277 | 52219 | 1382918 | 1063277 | 319641 | 0.9531876402963345 | 0.7688648206184314 | 0.851161576904388 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 582418 | 19974 | 571232 | 562589 | 8643 | 0.966842189139298 | 0.9848695451235225 | 0.9757726105833051 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 234038 | 23368 | 486828 | 231572 | 255256 | 0.9092173453610249 | 0.47567518712974605 | 0.6245858371672585 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 247624 | 19440 | 306624 | 243487 | 63137 | 0.9272084593955007 | 0.7940898298893759 | 0.8555017016826969 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 57530 | 28498 | 29032 |  | 0.4953589431600904 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 481662 | 42808 | 850982 | 503557 | 347425 | 0.9183785535874311 | 0.591736370452019 | 0.7197306421515836 |  |  |

</details>

<details><summary><code>readlik-z-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z-nolink | GT | ALL | ALL | ALL | 279965 | 266600 | 13365 | 273738 | 265778 | 7960 | 0.952261889879092 | 0.9709210997377054 | 0.9615009766116313 | 1988 | 1253 |
| readlik-z-nolink | GT | ALL | ALL | Snv | 223585 | 215626 | 7959 | 208745 | 207672 | 1073 | 0.9644027998300423 | 0.994859757119931 | 0.9793945500579433 | 372 | 641 |
| readlik-z-nolink | GT | ALL | ALL | Insertion | 28075 | 24986 | 3089 | 30068 | 26839 | 3229 | 0.889973285841496 | 0.8926100838100306 | 0.8912897346495064 | 935 | 300 |
| readlik-z-nolink | GT | ALL | ALL | Deletion | 28305 | 25988 | 2317 | 31009 | 27761 | 3248 | 0.9181416710828475 | 0.8952562159373085 | 0.9065545339293215 | 681 | 274 |
| readlik-z-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3916 | 3506 | 410 |  | 0.8953013278855976 |  | 0 | 38 |
| readlik-z-nolink | GT | ALL | ALL | JointIndel | 56380 | 50974 | 5406 | 64993 | 58106 | 6887 | 0.9041149343738915 | 0.8940347422030065 | 0.8990465841682481 | 1616 | 612 |
| readlik-z-nolink | BASEPAIR | ALL | ALL | ALL | 1115502 | 1066522 | 48980 | 1270302 | 1066522 | 203780 | 0.9560915175409815 | 0.8395814538590036 | 0.8940566785871765 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 583242 | 19150 | 566404 | 562119 | 4285 | 0.9682100691908259 | 0.9924347285683011 | 0.9801727455301588 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 235356 | 22050 | 381508 | 231034 | 150474 | 0.9143376611267803 | 0.6055810100967739 | 0.7285988847665609 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 248560 | 18504 | 318100 | 244146 | 73954 | 0.930713237276458 | 0.7675133605784344 | 0.8412715304061272 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 41190 | 31998 | 9192 |  | 0.7768390386016023 |  |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 483916 | 40554 | 740798 | 507178 | 233620 | 0.9226762255229088 | 0.6846373775307169 | 0.786030342991963 |  |  |

</details>

<details><summary><code>readlik-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z | GT | ALL | ALL | ALL | 279965 | 267397 | 12568 | 270324 | 265780 | 4544 | 0.955108674298573 | 0.983190541720306 | 0.9689461844947999 | 929 | 895 |
| readlik-z | GT | ALL | ALL | Snv | 223585 | 215594 | 7991 | 208482 | 207889 | 593 | 0.9642596775275623 | 0.9971556297426156 | 0.9804317957716053 | 287 | 312 |
| readlik-z | GT | ALL | ALL | Insertion | 28075 | 25477 | 2598 | 28838 | 26713 | 2125 | 0.9074621549421193 | 0.9263125043345586 | 0.9167904432323623 | 347 | 287 |
| readlik-z | GT | ALL | ALL | Deletion | 28305 | 26326 | 1979 | 29178 | 27690 | 1488 | 0.9300830242006712 | 0.9490026732469669 | 0.9394476021045399 | 295 | 260 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3826 | 3488 | 338 |  | 0.9116570831155254 |  | 0 | 36 |
| readlik-z | GT | ALL | ALL | JointIndel | 56380 | 51803 | 4577 | 61842 | 57891 | 3951 | 0.9188187300461157 | 0.9361113806151159 | 0.9273844496619754 | 642 | 583 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 1115544 | 1069394 | 46150 | 1225792 | 1069394 | 156398 | 0.9586300495542982 | 0.872410653683496 | 0.9134904174368822 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 602392 | 583049 | 19343 | 565224 | 561858 | 3366 | 0.9678896798098249 | 0.9940448388603457 | 0.9807929180565034 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 257406 | 236867 | 20539 | 360882 | 232272 | 128610 | 0.9202077651647592 | 0.6436231233477979 | 0.757456577044051 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 267064 | 249980 | 17084 | 294182 | 245061 | 49121 | 0.9360303148309019 | 0.8330251341006588 | 0.8815289300346398 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 38514 | 31943 | 6571 |  | 0.8293867165186686 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 524470 | 486847 | 37623 | 693578 | 509276 | 184302 | 0.9282647243884302 | 0.7342735784583709 | 0.8199513476065446 |  |  |

</details>

