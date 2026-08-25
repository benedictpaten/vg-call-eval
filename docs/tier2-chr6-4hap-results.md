# Tier 2 results: HG002 chr6 on HPRC v2.1 MC CHM13, 4-haplotype graph

Real reads, real benchmark, run on a 32 GB laptop.

This is the **4-haplotype** graph: CHM13, GRCh38 and 2 recombinants. It is kept as a thin-panel reference rather than the headline configuration -- the caller is tuned on the 34-haplotype graph, whose page is [tier2-chr6-results.md](tier2-chr6-results.md). The two are compared directly in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md).

| | |
|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz`, 100,179,277 nodes, **4 haplotypes** (CHM13, GRCh38, 2 recombinants). HG002 itself is **absent** — no circularity |
| chromosome | chr6 component, 5,316,235 nodes |
| reads | 596,017,764 alignments genome-wide (~28.6×); 151 bp paired Illumina |
| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |
| regions | small variants 167.2 Mb; SVs 168.4 Mb |
| engine | `aardvark compare` for small variants; `truvari bench --sizemin 50` for SVs |

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.7`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap down from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every number on this page — accuracy and cost alike — comes from one `vg` build in one pass, which is what the refresh harness exists to guarantee: a table whose rows come from different builds is not a comparison, it is a mixture of vintages.

Build: `vg version v1.4.0-18654-g648296d56`.

The wall column is what the caller costs unaided, and the repeatability note below applies to it harder than to the memory column. It includes snarl decomposition, which is single-threaded — 46 s of a 197 s chr20 run — and which `vg call -r` skips for byte-identical output given `vg snarls -T -P <ref path>`. The whole-genome harness caches one snarl file per contig for exactly that reason; this matrix does not, so these figures include it.

| arm | enumeration | pack? | variants | wall | CPU | peak RSS |
|---|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 288,849 | 329 s | 981 s (3.0x) | 6.1 GB |
| `poisson-z` | panel (`-z`) | yes | 289,002 | 161 s | 366 s (2.3x) | 5.9 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 293,606 | 308 s | 1,004 s (3.3x) | 8.0 GB |
| `readlik-nomismap` | panel (default) | **no** | 296,674 | 267 s | 912 s (3.4x) | 7.3 GB |
| `readlik-nolink` | panel (default) | **no** | 293,633 | 255 s | 875 s (3.4x) | 6.9 GB |
| `readlik` | panel (default) | **no** | 295,204 | 282 s | 892 s (3.2x) | 6.7 GB |

`CPU` is user+sys, with the multiple of wall clock beside it. It is the column that separates work from waiting: this caller has phases that run on one thread and phases that block on a subprocess, so a wall-clock change can come from either doing less or waiting less, and only CPU distinguishes them. A multiple well under `--threads` means the run spent its time parked rather than computing.

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9654 | 0.9926 | 0.9788 | 215,840 | 7,745 | 1,556 | 0.9734 | 0.9919 | 0.9826 |
| `poisson-z` | 0.9660 | 0.9926 | 0.9791 | 215,988 | 7,597 | 1,575 | 0.9738 | 0.9917 | 0.9827 |
| `readlik-support` | 0.9695 | 0.9939 | **0.9815** | 216,760 | 6,825 | 1,289 | 0.9758 | 0.9924 | 0.9840 |
| `readlik-nomismap` | 0.9698 | 0.9916 | 0.9806 | 216,822 | 6,763 | 1,791 | 0.9761 | 0.9907 | 0.9833 |
| `readlik-nolink` | 0.9694 | 0.9938 | 0.9815 | 216,751 | 6,834 | 1,308 | 0.9757 | 0.9924 | 0.9840 |
| `readlik` | 0.9698 | 0.9923 | 0.9809 | 216,834 | 6,751 | 1,636 | 0.9761 | 0.9911 | 0.9836 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7610 | 0.8659 | 0.8101 | 21,365 | 6,710 | 3,719 | 0.7939 | 0.8478 | 0.8199 |
| `poisson-z` | 0.7619 | 0.8671 | 0.8111 | 21,389 | 6,686 | 3,684 | 0.7960 | 0.8476 | 0.8210 |
| `readlik-support` | 0.8352 | 0.8631 | **0.8489** | 23,448 | 4,627 | 3,945 | 0.8859 | 0.7341 | 0.8029 |
| `readlik-nomismap` | 0.8396 | 0.8548 | 0.8471 | 23,572 | 4,503 | 4,253 | 0.8945 | 0.6300 | 0.7393 |
| `readlik-nolink` | 0.8349 | 0.8625 | 0.8485 | 23,440 | 4,635 | 3,965 | 0.8864 | 0.7335 | 0.8027 |
| `readlik` | 0.8400 | 0.8566 | 0.8482 | 23,584 | 4,491 | 4,187 | 0.8954 | 0.8046 | 0.8476 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8506 | 0.8211 | 0.8356 | 24,075 | 4,230 | 5,576 | 0.8989 | 0.7317 | 0.8068 |
| `poisson-z` | 0.8515 | 0.8223 | 0.8367 | 24,102 | 4,203 | 5,539 | 0.8987 | 0.7340 | 0.8080 |
| `readlik-support` | 0.8798 | 0.9024 | **0.8910** | 24,904 | 3,401 | 2,867 | 0.8988 | 0.8230 | 0.8593 |
| `readlik-nomismap` | 0.8785 | 0.9011 | 0.8896 | 24,866 | 3,439 | 2,935 | 0.9044 | 0.8286 | 0.8649 |
| `readlik-nolink` | 0.8798 | 0.9021 | 0.8908 | 24,902 | 3,403 | 2,877 | 0.8983 | 0.8227 | 0.8588 |
| `readlik` | 0.8789 | 0.9021 | 0.8904 | 24,878 | 3,427 | 2,898 | 0.9061 | 0.8621 | 0.8835 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.9081 | — | 0 | 0 | 279 | — | 0.7762 | — |
| `poisson-z` | — | 0.9137 | — | 0 | 0 | 264 | — | 0.7713 | — |
| `readlik-support` | — | 0.9370 | — | 0 | 0 | 183 | — | 0.7958 | — |
| `readlik-nomismap` | — | 0.9303 | — | 0 | 0 | 206 | — | 0.8053 | — |
| `readlik-nolink` | — | 0.9360 | — | 0 | 0 | 186 | — | 0.7965 | — |
| `readlik` | — | 0.9324 | — | 0 | 0 | 199 | — | 0.8210 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8060 | 0.8454 | 0.8252 | 45,440 | 10,940 | 9,574 | 0.8474 | 0.7808 | 0.8127 |
| `poisson-z` | 0.8069 | 0.8469 | 0.8264 | 45,491 | 10,889 | 9,487 | 0.8483 | 0.7819 | 0.8137 |
| `readlik-support` | 0.8576 | 0.8856 | **0.8714** | 48,352 | 8,028 | 6,995 | 0.8925 | 0.7778 | 0.8312 |
| `readlik-nomismap` | 0.8591 | 0.8806 | 0.8697 | 48,438 | 7,942 | 7,394 | 0.8996 | 0.7209 | 0.8004 |
| `readlik-nolink` | 0.8574 | 0.8850 | 0.8710 | 48,342 | 8,038 | 7,028 | 0.8925 | 0.7773 | 0.8309 |
| `readlik` | 0.8596 | 0.8820 | 0.8707 | 48,462 | 7,918 | 7,284 | 0.9008 | 0.8326 | 0.8654 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9333 | 0.9593 | 0.9461 | 261,280 | 18,685 | 11,130 | 0.9229 | 0.9095 | 0.9162 |
| `poisson-z` | 0.9340 | 0.9596 | 0.9466 | 261,479 | 18,486 | 11,062 | 0.9236 | 0.9102 | 0.9169 |
| `readlik-support` | 0.9469 | 0.9697 | **0.9582** | 265,112 | 14,853 | 8,284 | 0.9457 | 0.9053 | 0.9251 |
| `readlik-nomismap` | 0.9475 | 0.9666 | 0.9569 | 265,260 | 14,705 | 9,185 | 0.9493 | 0.8684 | 0.9070 |
| `readlik-nolink` | 0.9469 | 0.9695 | 0.9580 | 265,093 | 14,872 | 8,336 | 0.9457 | 0.9048 | 0.9248 |
| `readlik` | 0.9476 | 0.9675 | 0.9575 | 265,296 | 14,669 | 8,920 | 0.9499 | 0.9385 | 0.9442 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (167.2 Mb vs 168.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

The same mechanism was traced site by site on chr20, where 246 `readlik` calls carrying a >=200 bp insertion allele contribute 27,951 FP bases and zero TP bases — the whole of the precision difference there. The size-matched control below is the general test.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7938 | 0.8825 | **0.8358** |
| `sm50-poisson-z` | Deletion | 0.8945 | 0.8348 | **0.8636** |
| `sm50-poisson-z` | ALL | 0.9219 | 0.9498 | **0.9356** |
| `sm50-readlik` | Insertion | 0.8924 | 0.8845 | **0.8884** |
| `sm50-readlik` | Deletion | 0.9035 | 0.9069 | **0.9052** |
| `sm50-readlik` | ALL | 0.9484 | 0.9704 | **0.9593** |

The insertion BASEPAIR precision gap collapses from **0.043 to -0.002**, and insertion BASEPAIR F1 goes from 0.8358 for `poisson-z` against 0.8884 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5512 | 0.5490 | 846 | 670 | 701 |
| `poisson-z` | 0.5488 | 0.5468 | 0.5478 | 849 | 692 | 698 |
| `readlik-support` | 0.5721 | 0.5712 | 0.5717 | 885 | 650 | 662 |
| `readlik-nomismap` | 0.5947 | 0.5421 | 0.5672 | 920 | 772 | 627 |
| `readlik-nolink` | 0.5779 | 0.5647 | 0.5712 | 894 | 676 | 653 |
| `readlik` | 0.5966 | 0.5691 | **0.5825** | 923 | 689 | 624 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) looks inert on this graph, because it binds only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless.

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9575 | 0.9809 | 0.8482 | 0.8904 | 0.9442 |

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
| poisson | GT | ALL | ALL | ALL | 279965 | 261280 | 18685 | 273564 | 262434 | 11130 | 0.9332595145821799 | 0.9593148221257183 | 0.9461078150154776 | 3409 | 1970 |
| poisson | GT | ALL | ALL | Snv | 223585 | 215840 | 7745 | 211630 | 210074 | 1556 | 0.9653599302278776 | 0.9926475452440581 | 0.978813591900827 | 498 | 542 |
| poisson | GT | ALL | ALL | Insertion | 28075 | 21365 | 6710 | 27729 | 24010 | 3719 | 0.7609973285841496 | 0.8658804861336507 | 0.8100580521287053 | 1948 | 485 |
| poisson | GT | ALL | ALL | Deletion | 28305 | 24075 | 4230 | 31169 | 25593 | 5576 | 0.8505564387917329 | 0.8211043023516956 | 0.8355709195002339 | 963 | 900 |
| poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3036 | 2757 | 279 |  | 0.908102766798419 |  | 0 | 43 |
| poisson | GT | ALL | ALL | JointIndel | 56380 | 45440 | 10940 | 61934 | 52360 | 9574 | 0.8059595601277049 | 0.8454160880937772 | 0.8252164541953574 | 2911 | 1428 |
| poisson | BASEPAIR | ALL | ALL | ALL | 1115470 | 1029518 | 85952 | 1131928 | 1029518 | 102410 | 0.9229454848628829 | 0.9095260475931332 | 0.9161866300495061 |  |  |
| poisson | BASEPAIR | ALL | ALL | Snv | 602392 | 586381 | 16011 | 574920 | 570274 | 4646 | 0.9734209617657605 | 0.9919188756696584 | 0.9825828668978973 |  |  |
| poisson | BASEPAIR | ALL | ALL | Insertion | 257406 | 204345 | 53061 | 241638 | 204862 | 36776 | 0.7938626139250833 | 0.847805394846837 | 0.8199477644159954 |  |  |
| poisson | BASEPAIR | ALL | ALL | Deletion | 267064 | 240074 | 26990 | 327584 | 239704 | 87880 | 0.8989380822574364 | 0.7317329295692098 | 0.8067630952668078 |  |  |
| poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 23016 | 17866 | 5150 |  | 0.7762426138338547 |  |  |  |
| poisson | BASEPAIR | ALL | ALL | JointIndel | 524470 | 444419 | 80051 | 592238 | 462432 | 129806 | 0.8473678189410262 | 0.7808212238998511 | 0.8127345904802796 |  |  |

</details>

<details><summary><code>poisson-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-z | GT | ALL | ALL | ALL | 279965 | 261479 | 18486 | 273626 | 262564 | 11062 | 0.9339703177182862 | 0.9595725552396337 | 0.9465983549566294 | 3428 | 1969 |
| poisson-z | GT | ALL | ALL | Snv | 223585 | 215988 | 7597 | 211674 | 210099 | 1575 | 0.9660218708768478 | 0.9925593129056947 | 0.9791108097522244 | 506 | 548 |
| poisson-z | GT | ALL | ALL | Insertion | 28075 | 21389 | 6686 | 27721 | 24037 | 3684 | 0.7618521816562779 | 0.8671043613145268 | 0.8110779286798757 | 1954 | 482 |
| poisson-z | GT | ALL | ALL | Deletion | 28305 | 24102 | 4203 | 31172 | 25633 | 5539 | 0.851510333863275 | 0.8223084819709996 | 0.8366546765967963 | 968 | 890 |
| poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3059 | 2795 | 264 |  | 0.9136972866949984 |  | 0 | 49 |
| poisson-z | GT | ALL | ALL | JointIndel | 56380 | 45491 | 10889 | 61952 | 52465 | 9487 | 0.8068641362185172 | 0.8468653150826446 | 0.8263809420700016 | 2922 | 1421 |
| poisson-z | BASEPAIR | ALL | ALL | ALL | 1115470 | 1030284 | 85186 | 1131950 | 1030284 | 101666 | 0.9236321909150403 | 0.9101850788462388 | 0.9168593320340658 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Snv | 602392 | 586614 | 15778 | 575096 | 570349 | 4747 | 0.9738077530910105 | 0.991745725930975 | 0.9826948870269363 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Insertion | 257406 | 204887 | 52519 | 241888 | 205019 | 36869 | 0.7959682369486337 | 0.8475782180182564 | 0.8209629096071364 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Deletion | 267064 | 240004 | 27060 | 326620 | 239730 | 86890 | 0.8986759728005272 | 0.7339722001102198 | 0.8080163159300016 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 23904 | 18438 | 5466 |  | 0.7713353413654619 |  |  |  |
| poisson-z | BASEPAIR | ALL | ALL | JointIndel | 524470 | 444891 | 79579 | 592412 | 463187 | 129225 | 0.848267775087231 | 0.781866336265977 | 0.8137146666164571 |  |  |

</details>

<details><summary><code>readlik-support</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-support | GT | ALL | ALL | ALL | 279965 | 265112 | 14853 | 273075 | 264791 | 8284 | 0.9469469397960459 | 0.9696640117183923 | 0.9581708460984673 | 1368 | 1705 |
| readlik-support | GT | ALL | ALL | Snv | 223585 | 216760 | 6825 | 211956 | 210667 | 1289 | 0.9694746964241787 | 0.993918549132839 | 0.9815444622430187 | 182 | 618 |
| readlik-support | GT | ALL | ALL | Insertion | 28075 | 23448 | 4627 | 28826 | 24881 | 3945 | 0.8351914514692788 | 0.8631443835426351 | 0.8489378786657612 | 701 | 588 |
| readlik-support | GT | ALL | ALL | Deletion | 28305 | 24904 | 3401 | 29387 | 26520 | 2867 | 0.8798445504327858 | 0.902439854357369 | 0.8909989739187264 | 485 | 478 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2906 | 2723 | 183 |  | 0.9370268410185822 |  | 0 | 21 |
| readlik-support | GT | ALL | ALL | JointIndel | 56380 | 48352 | 8028 | 61119 | 54124 | 6995 | 0.8576090812344803 | 0.8855511379440109 | 0.8713561604293897 | 1186 | 1087 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 1115436 | 1054864 | 60572 | 1165218 | 1054864 | 110354 | 0.9456965706683306 | 0.9052932584288949 | 0.9250539538220177 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 602392 | 587791 | 14601 | 576678 | 572300 | 4378 | 0.9757616303005352 | 0.9924082416877357 | 0.9840145381909223 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 257406 | 228025 | 29381 | 311208 | 228465 | 82743 | 0.8858573615222645 | 0.7341231587876919 | 0.8028842277086462 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 267064 | 240044 | 27020 | 292394 | 240645 | 51749 | 0.8988257496330467 | 0.8230162041628761 | 0.8592520992255003 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21158 | 16837 | 4321 |  | 0.795774647887324 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 524470 | 468069 | 56401 | 624760 | 485947 | 138813 | 0.8924609605887849 | 0.7778138805301236 | 0.8312027528501652 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 265260 | 14705 | 275075 | 265890 | 9185 | 0.9474755773043059 | 0.9666091066072889 | 0.9569467108829757 | 1621 | 996 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 216822 | 6763 | 213173 | 211382 | 1791 | 0.9697519958852338 | 0.9915983731523222 | 0.9805535172717552 | 421 | 199 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 23572 | 4503 | 29282 | 25029 | 4253 | 0.8396081923419413 | 0.8547571887166178 | 0.847114968391621 | 614 | 459 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 24866 | 3439 | 29663 | 26728 | 2935 | 0.8785020314432079 | 0.9010551865960962 | 0.8896356957145349 | 586 | 331 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2957 | 2751 | 206 |  | 0.9303347987825499 |  | 0 | 7 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 48438 | 7942 | 61902 | 54508 | 7394 | 0.8591344448385952 | 0.8805531323705211 | 0.8697119372935526 | 1200 | 797 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115434 | 1058911 | 56523 | 1219442 | 1058911 | 160531 | 0.9493264505116394 | 0.8683570026290713 | 0.9070383180948367 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 587983 | 14409 | 577622 | 572238 | 5384 | 0.97608035963293 | 0.9906790253833823 | 0.9833255117468193 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 230255 | 27151 | 365884 | 230525 | 135359 | 0.894520718242776 | 0.6300494145685518 | 0.7393457902904469 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 241545 | 25519 | 292186 | 242103 | 50083 | 0.9044461252733427 | 0.8285920612212768 | 0.8648590493203917 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21166 | 17044 | 4122 |  | 0.8052537087782292 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 471800 | 52670 | 679236 | 489672 | 189564 | 0.8995748088546532 | 0.720915852516651 | 0.8003967633841987 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 279965 | 265093 | 14872 | 273066 | 264730 | 8336 | 0.9468790741699855 | 0.9694725817201703 | 0.9580426408596463 | 1367 | 1719 |
| readlik-nolink | GT | ALL | ALL | Snv | 223585 | 216751 | 6834 | 211951 | 210643 | 1308 | 0.9694344432766062 | 0.9938287623082693 | 0.9814800482786757 | 179 | 633 |
| readlik-nolink | GT | ALL | ALL | Insertion | 28075 | 23440 | 4635 | 28826 | 24861 | 3965 | 0.834906500445236 | 0.8624505654617359 | 0.8484550456469928 | 702 | 583 |
| readlik-nolink | GT | ALL | ALL | Deletion | 28305 | 24902 | 3403 | 29383 | 26506 | 2877 | 0.8797738915385974 | 0.9020862403430555 | 0.8907903689746437 | 486 | 480 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2906 | 2720 | 186 |  | 0.9359944941500344 |  | 0 | 23 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 56380 | 48342 | 8038 | 61115 | 54087 | 7028 | 0.8574317133735367 | 0.8850036815838992 | 0.8709995506730477 | 1188 | 1086 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 1115436 | 1054863 | 60573 | 1165854 | 1054863 | 110991 | 0.9456956741579078 | 0.9047985425276235 | 0.9247951816735268 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 587761 | 14631 | 576702 | 572292 | 4410 | 0.9757118288423485 | 0.9923530696963077 | 0.9839620931298112 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 228165 | 29241 | 311624 | 228583 | 83041 | 0.8864012493881261 | 0.7335218083331194 | 0.8027475678684973 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 239912 | 27152 | 292358 | 240513 | 51845 | 0.8983314860857322 | 0.8226660464225367 | 0.8588354115277973 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21272 | 16943 | 4329 |  | 0.7964930424971793 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 468077 | 56393 | 625254 | 486039 | 139215 | 0.8924762140827883 | 0.7773464863879319 | 0.8309424096420399 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 265296 | 14669 | 274642 | 265722 | 8920 | 0.947604164806315 | 0.9675213550731497 | 0.9574591911490689 | 1499 | 1068 |
| readlik | GT | ALL | ALL | Snv | 223585 | 216834 | 6751 | 212892 | 211256 | 1636 | 0.9698056667486639 | 0.9923153523852469 | 0.9809313926718922 | 334 | 250 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 23584 | 4491 | 29200 | 25013 | 4187 | 0.8400356188780054 | 0.8566095890410959 | 0.8482416511222441 | 599 | 473 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 24878 | 3427 | 29607 | 26709 | 2898 | 0.8789259848083377 | 0.9021177424257777 | 0.890370868553556 | 566 | 333 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2943 | 2744 | 199 |  | 0.9323819232076113 |  | 0 | 12 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 48462 | 7918 | 61750 | 54466 | 7284 | 0.8595601277048599 | 0.8820404858299595 | 0.8706552199726805 | 1165 | 818 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115434 | 1059595 | 55839 | 1129012 | 1059595 | 69417 | 0.9499396647403612 | 0.9385152682168125 | 0.9441929099653099 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 588001 | 14391 | 577386 | 572273 | 5113 | 0.9761102405078421 | 0.9911445722618837 | 0.9835699580233636 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 230487 | 26919 | 286806 | 230759 | 56047 | 0.8954220181347754 | 0.804582191446483 | 0.8475750890025944 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 241977 | 25087 | 281294 | 242498 | 38796 | 0.9060637150645539 | 0.8620802434463586 | 0.8835249237495183 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 20774 | 17055 | 3719 |  | 0.820978145759122 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 472464 | 52006 | 588874 | 490312 | 98562 | 0.9008408488569413 | 0.8326263343261886 | 0.8653914202375848 |  |  |

</details>

