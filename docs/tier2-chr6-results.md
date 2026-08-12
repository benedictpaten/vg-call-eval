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
| `poisson` | support (Flow) | yes | 288,849 | 337 s | 4.8 GB |
| `poisson-z` | haplotype (`-z`) | yes | 289,002 | 169 s | 4.4 GB |
| `readlik` | support (Flow) | yes | 286,465 | 337 s | 5.0 GB |
| `readlik-nomismap` | support (Flow) | yes | 287,943 | 346 s | 4.7 GB |
| `readlik-z-nolink` | haplotype (`-z`) | **no** | 286,474 | 296 s | 4.6 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 286,474 | 286 s | 5.0 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9654 | 0.9926 | 0.9788 | 215,840 | 7,745 | 1,556 | 0.9734 | 0.9919 | 0.9826 |
| `poisson-z` | 0.9660 | 0.9926 | 0.9791 | 215,988 | 7,597 | 1,575 | 0.9738 | 0.9917 | 0.9827 |
| `readlik` | 0.9683 | 0.9950 | 0.9815 | 216,492 | 7,093 | 1,047 | 0.9751 | 0.9932 | 0.9841 |
| `readlik-nomismap` | 0.9684 | 0.9945 | 0.9813 | 216,509 | 7,076 | 1,165 | 0.9751 | 0.9928 | 0.9839 |
| `readlik-z-nolink` | 0.9688 | 0.9951 | 0.9818 | 216,611 | 6,974 | 1,040 | 0.9754 | 0.9932 | 0.9842 |
| `readlik-z` | 0.9682 | 0.9965 | **0.9822** | 216,484 | 7,101 | 735 | 0.9749 | 0.9942 | 0.9845 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7610 | 0.8659 | 0.8101 | 21,365 | 6,710 | 3,719 | 0.7939 | 0.8478 | 0.8199 |
| `poisson-z` | 0.7619 | 0.8671 | 0.8111 | 21,389 | 6,686 | 3,684 | 0.7960 | 0.8476 | 0.8210 |
| `readlik` | 0.8339 | 0.8670 | 0.8502 | 23,413 | 4,662 | 3,806 | 0.8854 | 0.8004 | 0.8408 |
| `readlik-nomismap` | 0.8339 | 0.8644 | 0.8489 | 23,412 | 4,663 | 3,898 | 0.8848 | 0.6617 | 0.7572 |
| `readlik-z-nolink` | 0.8342 | 0.8662 | 0.8499 | 23,419 | 4,656 | 3,829 | 0.8866 | 0.8003 | 0.8412 |
| `readlik-z` | 0.8413 | 0.8721 | **0.8564** | 23,619 | 4,456 | 3,634 | 0.8931 | 0.8467 | 0.8693 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8506 | 0.8211 | 0.8356 | 24,075 | 4,230 | 5,576 | 0.8989 | 0.7317 | 0.8068 |
| `poisson-z` | 0.8515 | 0.8223 | 0.8367 | 24,102 | 4,203 | 5,539 | 0.8987 | 0.7340 | 0.8080 |
| `readlik` | 0.8787 | 0.9059 | 0.8921 | 24,873 | 3,432 | 2,747 | 0.8981 | 0.8274 | 0.8613 |
| `readlik-nomismap` | 0.8782 | 0.9048 | 0.8913 | 24,858 | 3,447 | 2,788 | 0.8962 | 0.8506 | 0.8728 |
| `readlik-z-nolink` | 0.8790 | 0.9057 | 0.8922 | 24,881 | 3,424 | 2,753 | 0.8976 | 0.8275 | 0.8611 |
| `readlik-z` | 0.8799 | 0.9164 | **0.8978** | 24,906 | 3,399 | 2,412 | 0.9013 | 0.8725 | 0.8866 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.9081 | — | 0 | 0 | 279 | — | 0.7762 | — |
| `poisson-z` | — | 0.9137 | — | 0 | 0 | 264 | — | 0.7713 | — |
| `readlik` | — | 0.9146 | — | 0 | 0 | 268 | — | 0.7713 | — |
| `readlik-nomismap` | — | 0.9132 | — | 0 | 0 | 274 | — | 0.7651 | — |
| `readlik-z-nolink` | — | 0.9203 | — | 0 | 0 | 250 | — | 0.7869 | — |
| `readlik-z` | — | 0.9216 | — | 0 | 0 | 245 | — | 0.7873 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8060 | 0.8454 | 0.8252 | 45,440 | 10,940 | 9,574 | 0.8474 | 0.7808 | 0.8127 |
| `poisson-z` | 0.8069 | 0.8469 | 0.8264 | 45,491 | 10,889 | 9,487 | 0.8483 | 0.7819 | 0.8137 |
| `readlik` | 0.8564 | 0.8881 | 0.8720 | 48,286 | 8,094 | 6,821 | 0.8919 | 0.8123 | 0.8502 |
| `readlik-nomismap` | 0.8562 | 0.8863 | 0.8710 | 48,270 | 8,110 | 6,960 | 0.8906 | 0.7474 | 0.8128 |
| `readlik-z-nolink` | 0.8567 | 0.8879 | 0.8720 | 48,300 | 8,080 | 6,832 | 0.8922 | 0.8129 | 0.8507 |
| `readlik-z` | 0.8607 | 0.8958 | **0.8779** | 48,525 | 7,855 | 6,291 | 0.8972 | 0.8565 | 0.8764 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9333 | 0.9593 | 0.9461 | 261,280 | 18,685 | 11,130 | 0.9229 | 0.9095 | 0.9162 |
| `poisson-z` | 0.9340 | 0.9596 | 0.9466 | 261,479 | 18,486 | 11,062 | 0.9236 | 0.9102 | 0.9169 |
| `readlik` | 0.9458 | 0.9711 | 0.9583 | 264,778 | 15,187 | 7,868 | 0.9450 | 0.9262 | 0.9355 |
| `readlik-nomismap` | 0.9458 | 0.9702 | 0.9578 | 264,779 | 15,186 | 8,125 | 0.9445 | 0.8867 | 0.9147 |
| `readlik-z-nolink` | 0.9462 | 0.9711 | 0.9585 | 264,911 | 15,054 | 7,872 | 0.9454 | 0.9266 | 0.9359 |
| `readlik-z` | 0.9466 | 0.9741 | **0.9602** | 265,009 | 14,956 | 7,026 | 0.9476 | 0.9524 | 0.9500 |

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
| `sm50-readlik-z` | Insertion | 0.8895 | 0.8940 | **0.8917** |
| `sm50-readlik-z` | Deletion | 0.8990 | 0.9150 | **0.9070** |
| `sm50-readlik-z` | ALL | 0.9461 | 0.9754 | **0.9605** |

The insertion BASEPAIR precision gap collapses from **0.001 to -0.011**, and insertion BASEPAIR F1 goes from 0.8358 for `poisson-z` against 0.8917 for `readlik-z`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5512 | 0.5490 | 846 | 670 | 701 |
| `poisson-z` | 0.5488 | 0.5468 | 0.5478 | 849 | 692 | 698 |
| `readlik` | 0.5417 | 0.5684 | 0.5547 | 838 | 625 | 709 |
| `readlik-nomismap` | 0.5391 | 0.5383 | 0.5387 | 834 | 711 | 713 |
| `readlik-z-nolink` | 0.5507 | 0.5728 | 0.5616 | 852 | 625 | 695 |
| `readlik-z` | 0.5514 | 0.5881 | **0.5691** | 853 | 587 | 694 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) looked inert on this graph, because it binds only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising the default to **0.5** removed 94% of the excess false-positive SNVs. A clamp that is inert on a sparse graph is not thereby harmless — see [tier2-chr20-hap32.md](tier2-chr20-hap32.md) and plan §9.20.

| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.5 (current defaults)** | 0.9602 | 0.9822 | 0.8564 | 0.8978 | 0.9500 |

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

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 264778 | 15187 | 272385 | 264517 | 7868 | 0.9457539335274052 | 0.9711144152578153 | 0.9582664127322944 | 1406 | 1486 |
| readlik | GT | ALL | ALL | Snv | 223585 | 216492 | 7093 | 211425 | 210378 | 1047 | 0.9682760471409084 | 0.9950478893224548 | 0.9814794381049708 | 204 | 495 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 23413 | 4662 | 28621 | 24815 | 3806 | 0.8339447907390917 | 0.867020719052444 | 0.8501611678243393 | 707 | 541 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 24873 | 3432 | 29199 | 26452 | 2747 | 0.878749337572867 | 0.905921435665605 | 0.8921285353260359 | 495 | 429 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3140 | 2872 | 268 |  | 0.9146496815286624 |  | 0 | 21 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 48286 | 8094 | 60960 | 54139 | 6821 | 0.8564384533522525 | 0.8881069553805774 | 0.8719852673023818 | 1202 | 991 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115520 | 1054213 | 61307 | 1138182 | 1054213 | 83969 | 0.9450417742398164 | 0.9262253312739087 | 0.9355389488051216 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 587375 | 15017 | 575044 | 571132 | 3912 | 0.9750710500803463 | 0.9931970423132839 | 0.9840505840923137 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 227895 | 29511 | 283958 | 227294 | 56664 | 0.8853523227896786 | 0.8004493622296255 | 0.8407628351817676 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 239854 | 27210 | 289250 | 239321 | 49929 | 0.898114309678579 | 0.8273846153846154 | 0.8612998268400726 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24800 | 19128 | 5672 |  | 0.7712903225806451 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467749 | 56721 | 598008 | 485743 | 112265 | 0.8918508208286461 | 0.8122683977471873 | 0.8502013584113115 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 264779 | 15186 | 273028 | 264903 | 8125 | 0.9457575054024611 | 0.9702411474281026 | 0.9578428939653023 | 1564 | 1281 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 216509 | 7076 | 211822 | 210657 | 1165 | 0.9683520808641009 | 0.9945000991398438 | 0.981251925368814 | 326 | 349 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 23412 | 4663 | 28750 | 24852 | 3898 | 0.8339091718610864 | 0.8644173913043478 | 0.8488892614166945 | 717 | 522 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 24858 | 3447 | 29299 | 26511 | 2788 | 0.8782193958664547 | 0.9048431687088296 | 0.8913325160484633 | 521 | 394 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3157 | 2883 | 274 |  | 0.9132087424770352 |  | 0 | 16 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 48270 | 8110 | 61206 | 54246 | 6960 | 0.8561546647747428 | 0.8862856582687971 | 0.870959642766283 | 1238 | 932 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115520 | 1053600 | 61920 | 1188256 | 1053600 | 134656 | 0.9444922547332186 | 0.886677618291008 | 0.9146722598030363 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 587409 | 14983 | 575290 | 571160 | 4130 | 0.975127491732958 | 0.9928210120113334 | 0.983894712021574 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 227751 | 29655 | 343268 | 227154 | 116114 | 0.8847928952705065 | 0.6617395154806157 | 0.7571809264994692 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 239340 | 27724 | 280816 | 238862 | 41954 | 0.8961896773807028 | 0.8505996809298616 | 0.8727997454368863 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24952 | 19091 | 5861 |  | 0.7651090092978519 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467091 | 57379 | 649036 | 485107 | 163929 | 0.8905962209468606 | 0.7474269532044447 | 0.812754826014711 |  |  |

</details>

<details><summary><code>readlik-z-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z-nolink | GT | ALL | ALL | ALL | 279965 | 264911 | 15054 | 272362 | 264490 | 7872 | 0.946228992909828 | 0.9710972896365866 | 0.9585018666409969 | 1400 | 1487 |
| readlik-z-nolink | GT | ALL | ALL | Snv | 223585 | 216611 | 6974 | 211417 | 210377 | 1040 | 0.9688082832032561 | 0.9950808118552434 | 0.9817688131246171 | 202 | 493 |
| readlik-z-nolink | GT | ALL | ALL | Insertion | 28075 | 23419 | 4656 | 28616 | 24787 | 3829 | 0.8341585040071238 | 0.8661937377690803 | 0.849874343357256 | 706 | 543 |
| readlik-z-nolink | GT | ALL | ALL | Deletion | 28305 | 24881 | 3424 | 29193 | 26440 | 2753 | 0.8790319731496202 | 0.9056965710958106 | 0.8921650819472894 | 492 | 429 |
| readlik-z-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3136 | 2886 | 250 |  | 0.920280612244898 |  | 0 | 22 |
| readlik-z-nolink | GT | ALL | ALL | JointIndel | 56380 | 48300 | 8080 | 60945 | 54113 | 6832 | 0.8566867683575736 | 0.8878989252604808 | 0.8720136404788092 | 1198 | 994 |
| readlik-z-nolink | BASEPAIR | ALL | ALL | ALL | 1115520 | 1054613 | 60907 | 1138174 | 1054613 | 83561 | 0.9454003514056225 | 0.9265832816423499 | 0.9358972424827859 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 587564 | 14828 | 575024 | 571117 | 3907 | 0.9753847992669226 | 0.9932055009877848 | 0.9842144890040659 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 228219 | 29187 | 284032 | 227311 | 56721 | 0.8866110347078157 | 0.8003006703470031 | 0.8412478297323296 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 239710 | 27354 | 289008 | 239156 | 49852 | 0.8975751130815085 | 0.8275065050102419 | 0.8611178010601789 |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 25164 | 19802 | 5362 |  | 0.7869178191066603 |  |  |  |
| readlik-z-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467929 | 56541 | 598204 | 486269 | 111935 | 0.8921940244437241 | 0.8128815587993393 | 0.8506931616038326 |  |  |

</details>

<details><summary><code>readlik-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z | GT | ALL | ALL | ALL | 279965 | 265009 | 14956 | 271744 | 264718 | 7026 | 0.9465790366652974 | 0.9741447833254828 | 0.9601641016532888 | 1256 | 1179 |
| readlik-z | GT | ALL | ALL | Snv | 223585 | 216484 | 7101 | 211355 | 210620 | 735 | 0.9682402665652884 | 0.9965224385512527 | 0.982177795851356 | 316 | 237 |
| readlik-z | GT | ALL | ALL | Insertion | 28075 | 23619 | 4456 | 28406 | 24772 | 3634 | 0.8412822796081924 | 0.8720692811377878 | 0.856399176468427 | 489 | 527 |
| readlik-z | GT | ALL | ALL | Deletion | 28305 | 24906 | 3399 | 28857 | 26445 | 2412 | 0.879915209326974 | 0.9164154277991475 | 0.8977944887389234 | 451 | 388 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3126 | 2881 | 245 |  | 0.9216250799744082 |  | 0 | 27 |
| readlik-z | GT | ALL | ALL | JointIndel | 56380 | 48525 | 7855 | 60389 | 54098 | 6291 | 0.8606775452288046 | 0.8958253986653198 | 0.8778998153769754 | 940 | 942 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 1115520 | 1057039 | 58481 | 1109832 | 1057039 | 52793 | 0.9475751219162364 | 0.9524315391879131 | 0.9499971240504873 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 602392 | 587302 | 15090 | 574158 | 570811 | 3347 | 0.974949866532092 | 0.994170594157009 | 0.9844664228863969 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 257406 | 229883 | 27523 | 270338 | 228889 | 41449 | 0.8930755304849148 | 0.8466771227130481 | 0.8692576144390064 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 267064 | 240694 | 26370 | 274986 | 239923 | 35063 | 0.9012596231614893 | 0.8724916904860611 | 0.8866423676518236 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 25268 | 19894 | 5374 |  | 0.7873199303466836 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 524470 | 470577 | 53893 | 570592 | 488706 | 81886 | 0.8972429309588728 | 0.8564894004823062 | 0.8763926469810053 |  |  |

</details>

