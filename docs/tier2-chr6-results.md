# Tier 2 results: HG002 chr6 on HPRC v2.1 MC CHM13

Real reads, real benchmark, run on a 32 GB laptop.

| | |
|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz`, 100,179,277 nodes, **4 haplotypes** (CHM13, GRCh38, 2 recombinants). HG002 itself is **absent** — no circularity |
| chromosome | chr6 component, 5,316,235 nodes |
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
| `poisson` | support (Flow) | yes | 288,849 | 322 s | 5.8 GB |
| `poisson-z` | haplotype (`-z`) | yes | 289,002 | 156 s | 6.6 GB |
| `readlik` | support (Flow) | yes | 286,574 | 288 s | 8.1 GB |
| `readlik-nomismap` | support (Flow) | yes | 287,939 | 288 s | 7.6 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 286,585 | 235 s | 6.6 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9654 | 0.9926 | 0.9788 | 215,840 | 7,745 | 1,556 | 0.9734 | 0.9919 | 0.9826 |
| `poisson-z` | 0.9660 | 0.9926 | 0.9791 | 215,988 | 7,597 | 1,575 | 0.9738 | 0.9917 | 0.9827 |
| `readlik` | 0.9684 | 0.9951 | 0.9816 | 216,524 | 7,061 | 1,026 | 0.9752 | 0.9933 | 0.9841 |
| `readlik-nomismap` | 0.9684 | 0.9945 | 0.9813 | 216,513 | 7,072 | 1,160 | 0.9751 | 0.9928 | 0.9839 |
| `readlik-z` | 0.9690 | 0.9952 | **0.9819** | 216,643 | 6,942 | 1,018 | 0.9755 | 0.9933 | 0.9843 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7610 | 0.8659 | 0.8101 | 21,365 | 6,710 | 3,719 | 0.7939 | 0.8478 | 0.8199 |
| `poisson-z` | 0.7619 | 0.8671 | 0.8111 | 21,389 | 6,686 | 3,684 | 0.7960 | 0.8476 | 0.8210 |
| `readlik` | 0.8339 | 0.8664 | **0.8498** | 23,413 | 4,662 | 3,826 | 0.8858 | 0.7061 | 0.7858 |
| `readlik-nomismap` | 0.8340 | 0.8642 | 0.8488 | 23,414 | 4,661 | 3,904 | 0.8853 | 0.6603 | 0.7564 |
| `readlik-z` | 0.8342 | 0.8656 | 0.8496 | 23,419 | 4,656 | 3,848 | 0.8871 | 0.7060 | 0.7862 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8506 | 0.8211 | 0.8356 | 24,075 | 4,230 | 5,576 | 0.8989 | 0.7317 | 0.8068 |
| `poisson-z` | 0.8515 | 0.8223 | 0.8367 | 24,102 | 4,203 | 5,539 | 0.8987 | 0.7340 | 0.8080 |
| `readlik` | 0.8787 | 0.9064 | 0.8923 | 24,872 | 3,433 | 2,733 | 0.8970 | 0.8698 | 0.8832 |
| `readlik-nomismap` | 0.8782 | 0.9052 | 0.8915 | 24,857 | 3,448 | 2,776 | 0.8957 | 0.8552 | 0.8750 |
| `readlik-z` | 0.8790 | 0.9062 | **0.8924** | 24,880 | 3,425 | 2,738 | 0.8969 | 0.8696 | 0.8831 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.9081 | — | 0 | 0 | 279 | — | 0.7762 | — |
| `poisson-z` | — | 0.9137 | — | 0 | 0 | 264 | — | 0.7713 | — |
| `readlik` | — | 0.9149 | — | 0 | 0 | 267 | — | 0.7718 | — |
| `readlik-nomismap` | — | 0.9134 | — | 0 | 0 | 273 | — | 0.7614 | — |
| `readlik-z` | — | 0.9197 | — | 0 | 0 | 252 | — | 0.7864 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8060 | 0.8454 | 0.8252 | 45,440 | 10,940 | 9,574 | 0.8474 | 0.7808 | 0.8127 |
| `poisson-z` | 0.8069 | 0.8469 | 0.8264 | 45,491 | 10,889 | 9,487 | 0.8483 | 0.7819 | 0.8137 |
| `readlik` | 0.8564 | 0.8880 | 0.8719 | 48,285 | 8,095 | 6,826 | 0.8915 | 0.7811 | 0.8326 |
| `readlik-nomismap` | 0.8562 | 0.8864 | 0.8710 | 48,271 | 8,109 | 6,953 | 0.8906 | 0.7481 | 0.8132 |
| `readlik-z` | 0.8567 | 0.8878 | **0.8720** | 48,299 | 8,081 | 6,838 | 0.8921 | 0.7815 | 0.8331 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9333 | 0.9593 | 0.9461 | 261,280 | 18,685 | 11,130 | 0.9229 | 0.9095 | 0.9162 |
| `poisson-z` | 0.9340 | 0.9596 | 0.9466 | 261,479 | 18,486 | 11,062 | 0.9236 | 0.9102 | 0.9169 |
| `readlik` | 0.9459 | 0.9712 | 0.9584 | 264,809 | 15,156 | 7,852 | 0.9449 | 0.9070 | 0.9256 |
| `readlik-nomismap` | 0.9458 | 0.9703 | 0.9579 | 264,784 | 15,181 | 8,113 | 0.9445 | 0.8873 | 0.9150 |
| `readlik-z` | 0.9463 | 0.9712 | **0.9586** | 264,942 | 15,023 | 7,856 | 0.9454 | 0.9073 | 0.9259 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (167.2 Mb vs 168.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

The same mechanism was traced site by site on chr20, where 246 `readlik-z` calls carrying a >=200 bp insertion allele contribute 27,951 FP bases and zero TP bases — the whole of the precision difference there. The size-matched control below is the general test.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7938 | 0.8825 | **0.8358** |
| `sm50-poisson-z` | Deletion | 0.8945 | 0.8348 | **0.8636** |
| `sm50-poisson-z` | ALL | 0.9219 | 0.9498 | **0.9356** |
| `sm50-readlik-z` | Insertion | 0.8835 | 0.8876 | **0.8855** |
| `sm50-readlik-z` | Deletion | 0.8954 | 0.9071 | **0.9012** |
| `sm50-readlik-z` | ALL | 0.9441 | 0.9714 | **0.9575** |

The insertion BASEPAIR precision gap collapses from **0.142 to -0.005**, and insertion BASEPAIR F1 goes from 0.8358 for `poisson-z` against 0.8855 for `readlik-z`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5512 | 0.5490 | 846 | 670 | 701 |
| `poisson-z` | 0.5488 | 0.5468 | 0.5478 | 849 | 692 | 698 |
| `readlik` | 0.5301 | 0.5638 | 0.5464 | 820 | 622 | 727 |
| `readlik-nomismap` | 0.5339 | 0.5374 | 0.5357 | 826 | 705 | 721 |
| `readlik-z` | 0.5385 | 0.5688 | **0.5532** | 833 | 620 | 714 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) looked inert on this graph, because it binds only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising the default to **0.5** removed 94% of the excess false-positive SNVs. A clamp that is inert on a sparse graph is not thereby harmless — see [tier2-chr20-hap32.md](tier2-chr20-hap32.md) and plan §9.20.

| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.5 (current defaults)** | 0.9586 | 0.9819 | 0.8496 | 0.8924 | 0.9259 |

Sweep rows other than the current one are historical: they were produced at the defaults in force at the time and are kept because the comparison between them is the result. The full grids are in plan §9.20-§9.21.

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

The same blindness has a second consequence, found later and now corrected. Because the model only weighs reads it can see, it had no way to know that a heterozygous deletion produces *no* reads over the deleted interval, and its flat `1/ploidy` mixture asserted that both haplotypes contributed equally everywhere. That cost it 94% of heterozygous deletions above 1 kb and mis-genotyped two thirds of heterozygous insertions above 1 kb. Weighting each haplotype by the reads it is *expected* to contribute at the site is now the default and fixes both, without moving small variants at all — see [tier2-sv-errors.md](tier2-sv-errors.md). It does not remove the need for a depth term: it corrects the *relative* weight between a genotype's haplotypes, while the pile-ups above are a statement about *absolute* depth.

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

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 264809 | 15156 | 272424 | 264572 | 7852 | 0.9458646616541353 | 0.9711772824714415 | 0.9583538581468398 | 1397 | 1453 |
| readlik | GT | ALL | ALL | Snv | 223585 | 216524 | 7061 | 211466 | 210440 | 1026 | 0.9684191694433885 | 0.9951481562047799 | 0.9816017391579184 | 200 | 462 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 23413 | 4662 | 28632 | 24806 | 3826 | 0.8339447907390917 | 0.8663732886281084 | 0.8498498012275432 | 703 | 545 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 24872 | 3433 | 29187 | 26454 | 2733 | 0.8787140081257728 | 0.9063624216260664 | 0.8923240966576904 | 494 | 426 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3139 | 2872 | 267 |  | 0.9149410640331316 |  | 0 | 20 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 48285 | 8095 | 60958 | 54132 | 6826 | 0.8564207165661583 | 0.8880212605400439 | 0.8719347668293203 | 1197 | 991 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115520 | 1054107 | 61413 | 1162162 | 1054107 | 108055 | 0.9449467512908778 | 0.9070224288868506 | 0.925596286048711 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 587434 | 14958 | 575050 | 571187 | 3863 | 0.9751689929481135 | 0.9932823232762369 | 0.9841423203296069 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 228020 | 29386 | 322088 | 227414 | 94674 | 0.8858379369556265 | 0.7060616974243064 | 0.7857985816465864 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 239564 | 27500 | 274938 | 239128 | 35810 | 0.8970284276428122 | 0.869752453280376 | 0.8831798940419862 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24662 | 19033 | 5629 |  | 0.7717541156435002 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467584 | 56886 | 621688 | 485575 | 136113 | 0.8915362175148245 | 0.7810589877880866 | 0.8326490156384717 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 264784 | 15181 | 273023 | 264910 | 8113 | 0.94577536477774 | 0.970284554781099 | 0.9578732057059981 | 1561 | 1278 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 216513 | 7072 | 211822 | 210662 | 1160 | 0.9683699711519109 | 0.9945237038645656 | 0.9812726004256312 | 326 | 347 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 23414 | 4661 | 28758 | 24854 | 3904 | 0.8339804096170971 | 0.864246470547326 | 0.8488437368832656 | 715 | 523 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 24857 | 3448 | 29289 | 26513 | 2776 | 0.8781840664193605 | 0.9052203899074738 | 0.8914972934989441 | 520 | 392 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3154 | 2881 | 273 |  | 0.913443246670894 |  | 0 | 16 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 48271 | 8109 | 61201 | 54248 | 6953 | 0.8561724015608372 | 0.8863907452492606 | 0.8710195604338279 | 1235 | 931 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115520 | 1053632 | 61888 | 1187504 | 1053632 | 133872 | 0.9445209409064831 | 0.8872660639458899 | 0.9149987147333245 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 587407 | 14985 | 575290 | 571160 | 4130 | 0.9751241716357455 | 0.9928210120113334 | 0.983893021985446 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 227882 | 29524 | 344244 | 227289 | 116955 | 0.88530181891642 | 0.6602555164360163 | 0.7563943391554802 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 239219 | 27845 | 279298 | 238847 | 40451 | 0.8957366024623311 | 0.8551690309275397 | 0.8749828519438815 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24898 | 18958 | 5940 |  | 0.7614266206120973 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467101 | 57369 | 648440 | 485094 | 163346 | 0.8906152878143649 | 0.7480938868669422 | 0.813156919675808 |  |  |

</details>

<details><summary><code>readlik-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z | GT | ALL | ALL | ALL | 279965 | 264942 | 15023 | 272405 | 264549 | 7856 | 0.9463397210365582 | 0.9711605880949322 | 0.9585895091049347 | 1396 | 1453 |
| readlik-z | GT | ALL | ALL | Snv | 223585 | 216643 | 6942 | 211458 | 210440 | 1018 | 0.9689514055057361 | 0.9951858052190033 | 0.9818934028040767 | 202 | 460 |
| readlik-z | GT | ALL | ALL | Insertion | 28075 | 23419 | 4656 | 28627 | 24779 | 3848 | 0.8341585040071238 | 0.8655814440912425 | 0.8495795175106301 | 702 | 546 |
| readlik-z | GT | ALL | ALL | Deletion | 28305 | 24880 | 3425 | 29180 | 26442 | 2738 | 0.8789966437025261 | 0.9061686086360521 | 0.8923758341993632 | 492 | 426 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3140 | 2888 | 252 |  | 0.9197452229299363 |  | 0 | 21 |
| readlik-z | GT | ALL | ALL | JointIndel | 56380 | 48299 | 8081 | 60947 | 54109 | 6838 | 0.8566690315714792 | 0.8878041577107979 | 0.8719587468399608 | 1194 | 993 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 1115520 | 1054626 | 60894 | 1162438 | 1054626 | 107812 | 0.9454120051635112 | 0.9072535481462237 | 0.9259398110061731 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 602392 | 587629 | 14763 | 575030 | 571172 | 3858 | 0.975492702426327 | 0.9932907848286177 | 0.9843112950308518 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 257406 | 228339 | 29067 | 322162 | 227431 | 94731 | 0.8870772243071257 | 0.7059522848753111 | 0.7862179448036161 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 267064 | 239537 | 27527 | 274856 | 239016 | 35840 | 0.8969273282808615 | 0.8696044474197397 | 0.8830545868624498 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 25150 | 19778 | 5372 |  | 0.7864015904572564 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467876 | 56594 | 622168 | 486225 | 135943 | 0.8920929700459511 | 0.7815011379563076 | 0.8331430756361053 |  |  |

</details>

