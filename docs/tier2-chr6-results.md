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
| `poisson` | support (Flow) | yes | 288,849 | 342 s | 6.9 GB |
| `poisson-z` | haplotype (`-z`) | yes | 289,002 | 168 s | 6.1 GB |
| `readlik` | support (Flow) | yes | 286,543 | 292 s | 9.6 GB |
| `readlik-nomismap` | support (Flow) | yes | 287,902 | 295 s | 10.3 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 286,557 | 242 s | 9.2 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9654 | 0.9926 | 0.9788 | 215,840 | 7,745 | 1,556 | 0.9734 | 0.9919 | 0.9826 |
| `poisson-z` | 0.9660 | 0.9926 | 0.9791 | 215,988 | 7,597 | 1,575 | 0.9738 | 0.9917 | 0.9827 |
| `readlik` | 0.9684 | 0.9951 | 0.9816 | 216,523 | 7,062 | 1,027 | 0.9752 | 0.9933 | 0.9841 |
| `readlik-nomismap` | 0.9683 | 0.9945 | 0.9812 | 216,506 | 7,079 | 1,168 | 0.9751 | 0.9928 | 0.9839 |
| `readlik-z` | 0.9690 | 0.9952 | **0.9819** | 216,645 | 6,940 | 1,018 | 0.9755 | 0.9933 | 0.9843 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7610 | 0.8659 | 0.8101 | 21,365 | 6,710 | 3,719 | 0.7939 | 0.8478 | 0.8199 |
| `poisson-z` | 0.7619 | 0.8671 | 0.8111 | 21,389 | 6,686 | 3,684 | 0.7960 | 0.8476 | 0.8210 |
| `readlik` | 0.8345 | 0.8647 | **0.8493** | 23,428 | 4,647 | 3,883 | 0.8885 | 0.6827 | 0.7721 |
| `readlik-nomismap` | 0.8347 | 0.8620 | 0.8481 | 23,434 | 4,641 | 3,980 | 0.8879 | 0.5923 | 0.7105 |
| `readlik-z` | 0.8348 | 0.8640 | 0.8491 | 23,436 | 4,639 | 3,901 | 0.8895 | 0.6825 | 0.7724 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8506 | 0.8211 | 0.8356 | 24,075 | 4,230 | 5,576 | 0.8989 | 0.7317 | 0.8068 |
| `poisson-z` | 0.8515 | 0.8223 | 0.8367 | 24,102 | 4,203 | 5,539 | 0.8987 | 0.7340 | 0.8080 |
| `readlik` | 0.8782 | 0.9085 | 0.8931 | 24,858 | 3,447 | 2,665 | 0.8941 | 0.8815 | 0.8877 |
| `readlik-nomismap` | 0.8772 | 0.9077 | 0.8922 | 24,829 | 3,476 | 2,697 | 0.8924 | 0.8707 | 0.8814 |
| `readlik-z` | 0.8784 | 0.9083 | **0.8931** | 24,864 | 3,441 | 2,669 | 0.8937 | 0.8824 | 0.8880 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.9081 | — | 0 | 0 | 279 | — | 0.7762 | — |
| `poisson-z` | — | 0.9137 | — | 0 | 0 | 264 | — | 0.7713 | — |
| `readlik` | — | 0.9146 | — | 0 | 0 | 268 | — | 0.7685 | — |
| `readlik-nomismap` | — | 0.9123 | — | 0 | 0 | 277 | — | 0.6892 | — |
| `readlik-z` | — | 0.9191 | — | 0 | 0 | 254 | — | 0.7711 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8060 | 0.8454 | 0.8252 | 45,440 | 10,940 | 9,574 | 0.8474 | 0.7808 | 0.8127 |
| `poisson-z` | 0.8069 | 0.8469 | 0.8264 | 45,491 | 10,889 | 9,487 | 0.8483 | 0.7819 | 0.8137 |
| `readlik` | 0.8564 | 0.8882 | 0.8720 | 48,286 | 8,094 | 6,816 | 0.8913 | 0.7715 | 0.8271 |
| `readlik-nomismap` | 0.8560 | 0.8864 | 0.8709 | 48,263 | 8,117 | 6,954 | 0.8902 | 0.7071 | 0.7882 |
| `readlik-z` | 0.8567 | 0.8880 | **0.8721** | 48,300 | 8,080 | 6,824 | 0.8917 | 0.7718 | 0.8274 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9333 | 0.9593 | 0.9461 | 261,280 | 18,685 | 11,130 | 0.9229 | 0.9095 | 0.9162 |
| `poisson-z` | 0.9340 | 0.9596 | 0.9466 | 261,479 | 18,486 | 11,062 | 0.9236 | 0.9102 | 0.9169 |
| `readlik` | 0.9459 | 0.9712 | 0.9584 | 264,809 | 15,156 | 7,843 | 0.9449 | 0.9010 | 0.9224 |
| `readlik-nomismap` | 0.9457 | 0.9703 | 0.9578 | 264,769 | 15,196 | 8,122 | 0.9443 | 0.8600 | 0.9002 |
| `readlik-z` | 0.9464 | 0.9712 | **0.9586** | 264,945 | 15,020 | 7,842 | 0.9452 | 0.9013 | 0.9227 |

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
| `sm50-readlik-z` | Insertion | 0.8859 | 0.8853 | **0.8856** |
| `sm50-readlik-z` | Deletion | 0.8925 | 0.9098 | **0.9011** |
| `sm50-readlik-z` | ALL | 0.9439 | 0.9715 | **0.9575** |

The insertion BASEPAIR precision gap collapses from **0.165 to -0.003**, and insertion BASEPAIR F1 goes from 0.8358 for `poisson-z` against 0.8856 for `readlik-z`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The per-record breakdown of these errors — what the false positives and false negatives actually are, and how much of each is the metric rather than the caller — is in [tier2-sv-errors.md](tier2-sv-errors.md). The short version: the read-likelihood caller's SV deficit is entirely heterozygous deletions, it is *ahead* on insertions, and after representation is harmonised it leads overall on three of the four datasets.

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5512 | **0.5490** | 846 | 670 | 701 |
| `poisson-z` | 0.5488 | 0.5468 | 0.5478 | 849 | 692 | 698 |
| `readlik` | 0.5036 | 0.5539 | 0.5276 | 779 | 616 | 768 |
| `readlik-nomismap` | 0.5087 | 0.5278 | 0.5181 | 787 | 696 | 760 |
| `readlik-z` | 0.5120 | 0.5599 | 0.5349 | 792 | 613 | 755 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) looked inert on this graph, because it binds only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising the default to **0.5** removed 94% of the excess false-positive SNVs. A clamp that is inert on a sparse graph is not thereby harmless — see [tier2-chr20-hap32.md](tier2-chr20-hap32.md) and plan §9.20.

| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.5 (current defaults)** | 0.9586 | 0.9819 | 0.8491 | 0.8931 | 0.9227 |

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

Filtering on depth is **not** that remedy, and that has now been tested properly rather than by two spot checks. Sweeping a two-sided cut on DP over a rolling local median, across both chromosomes and both graphs, against the one test a hard filter has to pass — beat lowering the GQ threshold to the same recall:

- a **minimum** fails in all eight dataset-by-benchmark cells. Few reads already means a small likelihood gap, so low depth depresses GQ on its own and a separate cut adds nothing;
- a **maximum** passes in exactly one configuration — 5x the local median, structural calls, 34-haplotype graph, worth about +0.025 precision — and is dominated everywhere else. The two original spot checks (DP 200 moving insertion BASEPAIR precision by 0.0001; DP 58 helping by +0.087 but costing SV insertion recall 0.4976 to 0.4167) were both right and both too narrow to conclude from.

What shipped instead attacks the same blindness from the other side: **GQ is now scaled by the fraction of reads the called genotype explains**, so a pile-up the call does not account for can no longer carry a saturated quality. The giants remain output that no metric charges for — they should be fixed because they are wrong, not because they cost a score — but they no longer look confident. See [tier2-quality-signals.md](tier2-quality-signals.md).

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. `vg call` emits `AD` (per-allele read support, ties split fractionally), `BL` (mean absolute fit), `GQI` (the raw likelihood-ratio quality) and `GQ` (that ratio scaled by the fraction of reads the called genotype explains). The scaling rescales a quality and does not change a genotype, so **the numbers on this page are unaffected by it**; what it changes is how the calls rank. See [tier2-quality-signals.md](tier2-quality-signals.md).

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
| readlik | GT | ALL | ALL | ALL | 279965 | 264809 | 15156 | 272423 | 264580 | 7843 | 0.9458646616541353 | 0.9712102135282263 | 0.9583698914488474 | 1407 | 1454 |
| readlik | GT | ALL | ALL | Snv | 223585 | 216523 | 7062 | 211468 | 210441 | 1027 | 0.9684146968714359 | 0.9951434732441787 | 0.9815971634072911 | 209 | 464 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 23428 | 4647 | 28689 | 24806 | 3883 | 0.8344790739091719 | 0.8646519571961379 | 0.8492976130574432 | 692 | 554 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 24858 | 3447 | 29126 | 26461 | 2665 | 0.8782193958664547 | 0.9085009956739682 | 0.8931035872680518 | 506 | 416 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3140 | 2872 | 268 |  | 0.9146496815286624 |  | 0 | 20 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 48286 | 8094 | 60955 | 54139 | 6816 | 0.8564384533522525 | 0.8881798047740136 | 0.8720203801103502 | 1198 | 990 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115520 | 1054002 | 61518 | 1169752 | 1054002 | 115750 | 0.9448526247848537 | 0.901047401500489 | 0.9224302402514887 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 587440 | 14952 | 575054 | 571189 | 3865 | 0.9751789532397509 | 0.9932788920692666 | 0.9841457083285294 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 228705 | 28701 | 334122 | 228089 | 106033 | 0.8884991025850213 | 0.6826518457329958 | 0.7720907440002462 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 238774 | 28290 | 270366 | 238319 | 32047 | 0.8940703352005511 | 0.8814680840046456 | 0.8877244860602951 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24792 | 19052 | 5740 |  | 0.7684737011939335 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467479 | 56991 | 629280 | 485460 | 143820 | 0.8913360154060289 | 0.7714530892448512 | 0.8270729230987488 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 264769 | 15196 | 273020 | 264898 | 8122 | 0.9457217866519029 | 0.9702512636436891 | 0.957829504348 | 1569 | 1281 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 216506 | 7079 | 211823 | 210655 | 1168 | 0.9683386631482434 | 0.994485962336479 | 0.9812381552435295 | 334 | 353 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 23434 | 4641 | 28833 | 24853 | 3980 | 0.834692787177204 | 0.8619637221239552 | 0.8481090871618178 | 697 | 530 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 24829 | 3476 | 29206 | 26509 | 2697 | 0.8771948419007243 | 0.9076559611038828 | 0.8921654694728177 | 538 | 382 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3158 | 2881 | 277 |  | 0.9122862571247625 |  | 0 | 16 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 48263 | 8117 | 61197 | 54243 | 6954 | 0.8560305072720823 | 0.8863669787734693 | 0.870934652449148 | 1235 | 928 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115520 | 1053381 | 62139 | 1224922 | 1053381 | 171541 | 0.9442959337349398 | 0.8599576136276432 | 0.9001556116323328 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 587402 | 14990 | 575292 | 571160 | 4132 | 0.9751158713927144 | 0.9928175604736378 | 0.9838871020114743 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 228543 | 28863 | 384838 | 227922 | 156916 | 0.8878697466259527 | 0.5922544031514559 | 0.7105414325457337 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 238341 | 28723 | 273258 | 237923 | 35335 | 0.8924490009885271 | 0.8706899706504475 | 0.8814352209064914 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 27552 | 18990 | 8562 |  | 0.6892421602787456 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 466884 | 57586 | 685648 | 484835 | 200813 | 0.8902015367895209 | 0.7071193965416657 | 0.7881681888213934 |  |  |

</details>

<details><summary><code>readlik-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z | GT | ALL | ALL | ALL | 279965 | 264945 | 15020 | 272398 | 264556 | 7842 | 0.9463504366617256 | 0.971211242373292 | 0.9586196818172414 | 1394 | 1454 |
| readlik-z | GT | ALL | ALL | Snv | 223585 | 216645 | 6940 | 211458 | 210440 | 1018 | 0.968960350649641 | 0.9951858052190033 | 0.981897995630667 | 200 | 462 |
| readlik-z | GT | ALL | ALL | Insertion | 28075 | 23436 | 4639 | 28683 | 24782 | 3901 | 0.8347640249332146 | 0.8639960952480563 | 0.8491285490253767 | 690 | 555 |
| readlik-z | GT | ALL | ALL | Deletion | 28305 | 24864 | 3441 | 29118 | 26449 | 2669 | 0.8784313725490196 | 0.908338484786043 | 0.8931346347198563 | 504 | 416 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3139 | 2885 | 254 |  | 0.9190825103536158 |  | 0 | 21 |
| readlik-z | GT | ALL | ALL | JointIndel | 56380 | 48300 | 8080 | 60940 | 54116 | 6824 | 0.8566867683575736 | 0.8880210042664917 | 0.8720725113002967 | 1194 | 992 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 1115520 | 1054402 | 61118 | 1169918 | 1054402 | 115516 | 0.9452112019506598 | 0.9012614559310994 | 0.9227132829680788 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 602392 | 587636 | 14756 | 575030 | 571172 | 3858 | 0.9755043227665706 | 0.9932907848286177 | 0.9843172106904154 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 257406 | 228974 | 28432 | 334152 | 228060 | 106092 | 0.8895441442701413 | 0.6825037707390649 | 0.7723902393899625 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 267064 | 238677 | 28387 | 269900 | 238173 | 31727 | 0.8937071263816913 | 0.8824490552056317 | 0.8880424114519687 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 25622 | 19757 | 5865 |  | 0.7710951526032316 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 524470 | 467651 | 56819 | 629674 | 485990 | 143684 | 0.8916639655271036 | 0.7718120805369127 | 0.8274204152222384 |  |  |

</details>

