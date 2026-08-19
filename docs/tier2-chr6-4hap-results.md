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

Every arm on this page was re-run together on one build, so the wall-clock column compares runs made on the same machine in the same session rather than a mixture of vintages.

Two changes since the accuracy results were first produced left the calls untouched. The read path was optimised (vg `44fd008`). Then `AD`, `BL`, `GQI` and the explained-share scaling of `GQ` were added — which rescales a quality and does not change a genotype. Both are confirmed by the variant counts below, which are unchanged to the record.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 288,849 | 360 s | 4.3 GB |
| `poisson-z` | panel (`-z`) | yes | 289,002 | 173 s | 3.8 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 293,247 | 379 s | 5.2 GB |
| `readlik-nomismap` | panel (default) | **no** | 295,242 | 291 s | 6.1 GB |
| `readlik-nolink` | panel (default) | **no** | 293,250 | 300 s | 4.8 GB |
| `readlik` | panel (default) | **no** | 293,250 | 269 s | 4.8 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9654 | 0.9926 | 0.9788 | 215,840 | 7,745 | 1,556 | 0.9734 | 0.9919 | 0.9826 |
| `poisson-z` | 0.9660 | 0.9926 | 0.9791 | 215,988 | 7,597 | 1,575 | 0.9738 | 0.9917 | 0.9827 |
| `readlik-support` | 0.9695 | 0.9941 | 0.9816 | 216,758 | 6,827 | 1,255 | 0.9758 | 0.9925 | 0.9841 |
| `readlik-nomismap` | 0.9693 | 0.9950 | 0.9820 | 216,717 | 6,868 | 1,066 | 0.9756 | 0.9932 | 0.9843 |
| `readlik-nolink` | 0.9694 | 0.9940 | 0.9816 | 216,744 | 6,841 | 1,269 | 0.9757 | 0.9925 | 0.9840 |
| `readlik` | 0.9688 | 0.9958 | **0.9821** | 216,600 | 6,985 | 881 | 0.9752 | 0.9938 | 0.9844 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7610 | 0.8659 | 0.8101 | 21,365 | 6,710 | 3,719 | 0.7939 | 0.8478 | 0.8199 |
| `poisson-z` | 0.7619 | 0.8671 | 0.8111 | 21,389 | 6,686 | 3,684 | 0.7960 | 0.8476 | 0.8210 |
| `readlik-support` | 0.8351 | 0.8637 | 0.8491 | 23,445 | 4,630 | 3,925 | 0.8858 | 0.7347 | 0.8032 |
| `readlik-nomismap` | 0.8417 | 0.8674 | 0.8543 | 23,630 | 4,445 | 3,807 | 0.8914 | 0.6300 | 0.7382 |
| `readlik-nolink` | 0.8348 | 0.8629 | 0.8486 | 23,437 | 4,638 | 3,947 | 0.8863 | 0.7341 | 0.8030 |
| `readlik` | 0.8417 | 0.8699 | **0.8556** | 23,632 | 4,443 | 3,719 | 0.8922 | 0.8070 | 0.8474 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8506 | 0.8211 | 0.8356 | 24,075 | 4,230 | 5,576 | 0.8989 | 0.7317 | 0.8068 |
| `poisson-z` | 0.8515 | 0.8223 | 0.8367 | 24,102 | 4,203 | 5,539 | 0.8987 | 0.7340 | 0.8080 |
| `readlik-support` | 0.8797 | 0.9028 | 0.8911 | 24,899 | 3,406 | 2,854 | 0.8989 | 0.8235 | 0.8596 |
| `readlik-nomismap` | 0.8797 | 0.9112 | 0.8952 | 24,901 | 3,404 | 2,583 | 0.9001 | 0.8584 | 0.8787 |
| `readlik-nolink` | 0.8795 | 0.9024 | 0.8908 | 24,895 | 3,410 | 2,866 | 0.8983 | 0.8230 | 0.8590 |
| `readlik` | 0.8801 | 0.9120 | **0.8958** | 24,911 | 3,394 | 2,553 | 0.9018 | 0.8276 | 0.8631 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.9081 | — | 0 | 0 | 279 | — | 0.7762 | — |
| `poisson-z` | — | 0.9137 | — | 0 | 0 | 264 | — | 0.7713 | — |
| `readlik-support` | — | 0.9307 | — | 0 | 0 | 204 | — | 0.7901 | — |
| `readlik-nomismap` | — | 0.9298 | — | 0 | 0 | 208 | — | 0.8012 | — |
| `readlik-nolink` | — | 0.9303 | — | 0 | 0 | 205 | — | 0.7954 | — |
| `readlik` | — | 0.9319 | — | 0 | 0 | 200 | — | 0.8069 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8060 | 0.8454 | 0.8252 | 45,440 | 10,940 | 9,574 | 0.8474 | 0.7808 | 0.8127 |
| `poisson-z` | 0.8069 | 0.8469 | 0.8264 | 45,491 | 10,889 | 9,487 | 0.8483 | 0.7819 | 0.8137 |
| `readlik-support` | 0.8575 | 0.8857 | 0.8714 | 48,344 | 8,036 | 6,983 | 0.8925 | 0.7782 | 0.8314 |
| `readlik-nomismap` | 0.8608 | 0.8914 | 0.8758 | 48,531 | 7,849 | 6,598 | 0.8958 | 0.7317 | 0.8055 |
| `readlik-nolink` | 0.8573 | 0.8851 | 0.8710 | 48,332 | 8,048 | 7,018 | 0.8924 | 0.7778 | 0.8311 |
| `readlik` | 0.8610 | 0.8931 | **0.8767** | 48,543 | 7,837 | 6,472 | 0.8971 | 0.8170 | 0.8552 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9333 | 0.9593 | 0.9461 | 261,280 | 18,685 | 11,130 | 0.9229 | 0.9095 | 0.9162 |
| `poisson-z` | 0.9340 | 0.9596 | 0.9466 | 261,479 | 18,486 | 11,062 | 0.9236 | 0.9102 | 0.9169 |
| `readlik-support` | 0.9469 | 0.9698 | 0.9582 | 265,102 | 14,863 | 8,238 | 0.9457 | 0.9053 | 0.9251 |
| `readlik-nomismap` | 0.9474 | 0.9719 | 0.9595 | 265,248 | 14,717 | 7,664 | 0.9473 | 0.8761 | 0.9103 |
| `readlik-nolink` | 0.9468 | 0.9696 | 0.9581 | 265,076 | 14,889 | 8,287 | 0.9456 | 0.9049 | 0.9248 |
| `readlik` | 0.9471 | 0.9730 | **0.9599** | 265,143 | 14,822 | 7,353 | 0.9476 | 0.9294 | 0.9384 |

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
| `sm50-readlik` | Insertion | 0.8897 | 0.8906 | **0.8901** |
| `sm50-readlik` | Deletion | 0.8995 | 0.9080 | **0.9037** |
| `sm50-readlik` | ALL | 0.9464 | 0.9727 | **0.9594** |

The insertion BASEPAIR precision gap collapses from **0.041 to -0.008**, and insertion BASEPAIR F1 goes from 0.8358 for `poisson-z` against 0.8901 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5469 | 0.5512 | 0.5490 | 846 | 670 | 701 |
| `poisson-z` | 0.5488 | 0.5468 | 0.5478 | 849 | 692 | 698 |
| `readlik-support` | 0.5714 | 0.5677 | 0.5695 | 884 | 658 | 663 |
| `readlik-nomismap` | 0.5818 | 0.5372 | 0.5586 | 900 | 766 | 647 |
| `readlik-nolink` | 0.5779 | 0.5647 | 0.5712 | 894 | 676 | 653 |
| `readlik` | 0.5779 | 0.5763 | **0.5771** | 894 | 644 | 653 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) looks inert on this graph, because it binds only where `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60. **That reading did not survive the 34-haplotype graph.** There the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless.

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9599 | 0.9821 | 0.8556 | 0.8958 | 0.9384 |

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
| readlik-support | GT | ALL | ALL | ALL | 279965 | 265102 | 14863 | 272960 | 264722 | 8238 | 0.9469112210454879 | 0.969819753810082 | 0.9582285874453902 | 1368 | 1686 |
| readlik-support | GT | ALL | ALL | Snv | 223585 | 216758 | 6827 | 211868 | 210613 | 1255 | 0.9694657512802737 | 0.9940765004625521 | 0.9816168921199427 | 182 | 605 |
| readlik-support | GT | ALL | ALL | Insertion | 28075 | 23445 | 4630 | 28792 | 24867 | 3925 | 0.8350845948352627 | 0.8636774103917755 | 0.8491403717603055 | 701 | 582 |
| readlik-support | GT | ALL | ALL | Deletion | 28305 | 24899 | 3406 | 29355 | 26501 | 2854 | 0.8796679031973149 | 0.9027763583716573 | 0.8910723362826163 | 485 | 474 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2945 | 2741 | 204 |  | 0.930730050933786 |  | 0 | 25 |
| readlik-support | GT | ALL | ALL | JointIndel | 56380 | 48344 | 8036 | 61092 | 54109 | 6983 | 0.8574671869457254 | 0.885696981601519 | 0.8713534995766994 | 1186 | 1081 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 1115468 | 1054906 | 60562 | 1165236 | 1054906 | 110330 | 0.9457070933455733 | 0.9053153180986513 | 0.9250705045459647 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 602392 | 587791 | 14601 | 576420 | 572117 | 4303 | 0.9757616303005352 | 0.9925349571493008 | 0.9840768246954615 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 257406 | 228002 | 29404 | 310882 | 228406 | 82476 | 0.88576800851573 | 0.734703199284616 | 0.8031942641718942 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 267064 | 240077 | 26987 | 292198 | 240636 | 51562 | 0.8989493155198753 | 0.8235374643221377 | 0.859592594290055 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21624 | 17086 | 4538 |  | 0.790140584535701 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 524470 | 468079 | 56391 | 624704 | 486128 | 138576 | 0.8924800274562892 | 0.7781733428951951 | 0.8314162336221081 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 265248 | 14717 | 272990 | 265326 | 7664 | 0.9474327148036362 | 0.9719257115645261 | 0.959522935210648 | 1360 | 1168 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 216717 | 6868 | 212228 | 211162 | 1066 | 0.9692823758302211 | 0.9949771000998926 | 0.9819616800115674 | 386 | 201 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 23630 | 4445 | 28705 | 24898 | 3807 | 0.8416740872662511 | 0.8673750217732102 | 0.8543313072832888 | 499 | 511 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 24901 | 3404 | 29095 | 26512 | 2583 | 0.8797385620915033 | 0.9112218594260182 | 0.8952034882809369 | 475 | 434 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2962 | 2754 | 208 |  | 0.9297771775827144 |  | 0 | 22 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 48531 | 7849 | 60762 | 54164 | 6598 | 0.8607839659453707 | 0.8914123959053356 | 0.8758304881192736 | 974 | 967 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115468 | 1056679 | 58789 | 1206062 | 1056679 | 149383 | 0.9472965607260809 | 0.8761398667730182 | 0.9103298255891588 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 587687 | 14705 | 575808 | 571873 | 3935 | 0.975588985245488 | 0.9931661248193843 | 0.9842990902619778 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 229445 | 27961 | 364580 | 229676 | 134904 | 0.8913739384474333 | 0.629974216907126 | 0.7382170831422498 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 240380 | 26684 | 280742 | 240990 | 39752 | 0.9000838750262109 | 0.8584038013549807 | 0.8787498829117012 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21632 | 17331 | 4301 |  | 0.8011741863905325 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 469825 | 54645 | 666954 | 487997 | 178957 | 0.8958091025225465 | 0.7316801458571356 | 0.8054685896652211 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 279965 | 265076 | 14889 | 272945 | 264658 | 8287 | 0.9468183522940368 | 0.9696385718734544 | 0.9580925961493001 | 1368 | 1701 |
| readlik-nolink | GT | ALL | ALL | Snv | 223585 | 216744 | 6841 | 211858 | 210589 | 1269 | 0.9694031352729388 | 0.9940101388665993 | 0.9815524401327733 | 181 | 619 |
| readlik-nolink | GT | ALL | ALL | Insertion | 28075 | 23437 | 4638 | 28792 | 24845 | 3947 | 0.83479964381122 | 0.8629133092525701 | 0.8486236991995663 | 702 | 582 |
| readlik-nolink | GT | ALL | ALL | Deletion | 28305 | 24895 | 3410 | 29354 | 26488 | 2866 | 0.8795265854089384 | 0.9023642433739865 | 0.8907990646227495 | 485 | 474 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2941 | 2736 | 205 |  | 0.9302958177490649 |  | 0 | 26 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 56380 | 48332 | 8048 | 61087 | 54069 | 7018 | 0.8572543455125932 | 0.8851146725162473 | 0.8709617669280083 | 1187 | 1082 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 1115468 | 1054835 | 60633 | 1165674 | 1054835 | 110839 | 0.9456434429315766 | 0.9049142384577507 | 0.9248306330776426 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 587759 | 14633 | 576426 | 572101 | 4325 | 0.975708508745136 | 0.9924968686353496 | 0.9840310881776533 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 228128 | 29278 | 311306 | 228521 | 82785 | 0.8862575075950055 | 0.7340719420762851 | 0.8030178923328628 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 239897 | 27167 | 292212 | 240504 | 51708 | 0.8982753197735375 | 0.8230462814668802 | 0.8590168869551313 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21620 | 17197 | 4423 |  | 0.7954209065679926 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 468025 | 56445 | 625138 | 486222 | 138916 | 0.8923770663717658 | 0.7777834654108373 | 0.831148999066556 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 265143 | 14822 | 272332 | 264979 | 7353 | 0.947057667922776 | 0.9729998678084103 | 0.9598535132905656 | 1256 | 1292 |
| readlik | GT | ALL | ALL | Snv | 223585 | 216600 | 6985 | 211796 | 210915 | 881 | 0.9687590849117785 | 0.9958403369279873 | 0.9821130585665806 | 314 | 276 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 23632 | 4443 | 28582 | 24863 | 3719 | 0.8417453250222618 | 0.8698831432370023 | 0.8555829523916358 | 488 | 529 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 24911 | 3394 | 29019 | 26466 | 2553 | 0.8800918565624448 | 0.912023157241807 | 0.895773036330979 | 454 | 459 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 2935 | 2735 | 200 |  | 0.9318568994889267 |  | 0 | 28 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 48543 | 7837 | 60536 | 54064 | 6472 | 0.860996807378503 | 0.8930884102021938 | 0.876749045238968 | 942 | 1016 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115468 | 1057071 | 58397 | 1137370 | 1057071 | 80299 | 0.9476479827301186 | 0.9293994038879169 | 0.9384349873359735 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 587463 | 14929 | 575262 | 571692 | 3570 | 0.975217134357694 | 0.9937941320650416 | 0.9844179991663768 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 229645 | 27761 | 284886 | 229889 | 54997 | 0.89215092111295 | 0.8069508505156449 | 0.8474147418378163 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 240831 | 26233 | 291636 | 241353 | 50283 | 0.9017726088128688 | 0.8275830144426614 | 0.8630864397205646 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 21414 | 17278 | 4136 |  | 0.8068553282899038 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 470476 | 53994 | 597936 | 488520 | 109416 | 0.8970503555970789 | 0.8170105161756442 | 0.8551616644791992 |  |  |

</details>

