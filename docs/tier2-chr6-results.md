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
| `poisson` | support (Flow) | yes | 294,626 | 773 s | 3.9 GB |
| `poisson-z` | panel (`-z`) | yes | 294,835 | 221 s | 4.6 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 298,098 | 352 s | 8.1 GB |
| `readlik-nomismap` | panel (default) | **no** | 305,961 | 360 s | 4.8 GB |
| `readlik-nolink` | panel (default) | **no** | 297,938 | 356 s | 4.5 GB |
| `readlik` | panel (default) | **no** | 297,938 | 344 s | 4.8 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9573 | 0.9800 | 0.9685 | 214,033 | 9,552 | 4,235 | 0.9639 | 0.9825 | 0.9731 |
| `poisson-z` | 0.9607 | 0.9801 | 0.9703 | 214,796 | 8,789 | 4,230 | 0.9671 | 0.9824 | 0.9747 |
| `readlik-support` | 0.9786 | 0.9933 | 0.9859 | 218,790 | 4,795 | 1,423 | 0.9821 | 0.9914 | 0.9868 |
| `readlik-nomismap` | 0.9791 | 0.9854 | 0.9822 | 218,918 | 4,667 | 3,157 | 0.9825 | 0.9855 | 0.9840 |
| `readlik-nolink` | 0.9787 | 0.9933 | 0.9859 | 218,814 | 4,771 | 1,423 | 0.9821 | 0.9914 | 0.9868 |
| `readlik` | 0.9782 | 0.9960 | **0.9870** | 218,705 | 4,880 | 852 | 0.9817 | 0.9934 | 0.9875 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7744 | 0.8529 | 0.8118 | 21,742 | 6,333 | 4,279 | 0.7996 | 0.6375 | 0.7094 |
| `poisson-z` | 0.7792 | 0.8558 | 0.8157 | 21,876 | 6,199 | 4,192 | 0.8041 | 0.5959 | 0.6845 |
| `readlik-support` | 0.9020 | 0.8870 | 0.8945 | 25,325 | 2,750 | 3,492 | 0.9218 | 0.5414 | 0.6821 |
| `readlik-nomismap` | 0.9214 | 0.9152 | 0.9183 | 25,869 | 2,206 | 2,544 | 0.9328 | 0.4494 | 0.6065 |
| `readlik-nolink` | 0.9029 | 0.8869 | 0.8949 | 25,350 | 2,725 | 3,496 | 0.9229 | 0.5401 | 0.6814 |
| `readlik` | 0.9205 | 0.9225 | **0.9215** | 25,842 | 2,233 | 2,297 | 0.9289 | 0.5904 | 0.7220 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8557 | 0.7206 | 0.7824 | 24,222 | 4,083 | 10,150 | 0.9292 | 0.5149 | 0.6626 |
| `poisson-z` | 0.8593 | 0.7229 | 0.7852 | 24,323 | 3,982 | 10,083 | 0.9306 | 0.4736 | 0.6278 |
| `readlik-support` | 0.9299 | 0.8896 | 0.9093 | 26,322 | 1,983 | 3,515 | 0.9426 | 0.7553 | 0.8387 |
| `readlik-nomismap` | 0.9422 | 0.9365 | 0.9394 | 26,670 | 1,635 | 1,919 | 0.9482 | 0.8381 | 0.8897 |
| `readlik-nolink` | 0.9306 | 0.8899 | 0.9098 | 26,342 | 1,963 | 3,509 | 0.9439 | 0.7584 | 0.8410 |
| `readlik` | 0.9428 | 0.9417 | **0.9423** | 26,687 | 1,618 | 1,749 | 0.9505 | 0.8180 | 0.8793 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.7865 | — | 0 | 0 | 825 | — | 0.3420 | — |
| `poisson-z` | — | 0.7915 | — | 0 | 0 | 836 | — | 0.4509 | — |
| `readlik-support` | — | 0.8999 | — | 0 | 0 | 342 | — | 0.5786 | — |
| `readlik-nomismap` | — | 0.8657 | — | 0 | 0 | 487 | — | 0.4977 | — |
| `readlik-nolink` | — | 0.8991 | — | 0 | 0 | 347 | — | 0.7442 | — |
| `readlik` | — | 0.9133 | — | 0 | 0 | 292 | — | 0.7595 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8153 | 0.7798 | 0.7972 | 45,964 | 10,416 | 15,254 | 0.8656 | 0.5453 | 0.6691 |
| `poisson-z` | 0.8194 | 0.7825 | 0.8005 | 46,199 | 10,181 | 15,111 | 0.8685 | 0.5176 | 0.6486 |
| `readlik-support` | 0.9161 | 0.8889 | 0.9023 | 51,647 | 4,733 | 7,349 | 0.9324 | 0.6312 | 0.7528 |
| `readlik-nomismap` | 0.9319 | 0.9225 | 0.9272 | 52,539 | 3,841 | 4,950 | 0.9406 | 0.5850 | 0.7214 |
| `readlik-nolink` | 0.9168 | 0.8890 | 0.9027 | 51,692 | 4,688 | 7,352 | 0.9336 | 0.6381 | 0.7581 |
| `readlik` | 0.9317 | 0.9311 | **0.9314** | 52,529 | 3,851 | 4,338 | 0.9399 | 0.6919 | 0.7970 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9287 | 0.9308 | 0.9297 | 259,997 | 19,968 | 19,489 | 0.9264 | 0.7386 | 0.8219 |
| `poisson-z` | 0.9322 | 0.9313 | 0.9318 | 260,995 | 18,970 | 19,341 | 0.9296 | 0.7139 | 0.8076 |
| `readlik-support` | 0.9660 | 0.9686 | 0.9673 | 270,437 | 9,528 | 8,772 | 0.9680 | 0.8012 | 0.8768 |
| `readlik-nomismap` | 0.9696 | 0.9710 | 0.9703 | 271,457 | 8,508 | 8,107 | 0.9723 | 0.7625 | 0.8547 |
| `readlik-nolink` | 0.9662 | 0.9686 | 0.9674 | 270,506 | 9,459 | 8,775 | 0.9687 | 0.8062 | 0.8800 |
| `readlik` | 0.9688 | 0.9812 | **0.9750** | 271,234 | 8,731 | 5,190 | 0.9715 | 0.8438 | 0.9031 |

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
| `sm50-readlik` | Insertion | 0.9266 | 0.9012 | **0.9137** |
| `sm50-readlik` | Deletion | 0.9479 | 0.9070 | **0.9270** |
| `sm50-readlik` | ALL | 0.9701 | 0.9718 | **0.9710** |

The insertion BASEPAIR precision gap collapses from **0.005 to -0.037**, and insertion BASEPAIR F1 goes from 0.8320 for `poisson-z` against 0.9137 for `readlik`.
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
| `readlik-nomismap` | 0.6005 | 0.4662 | 0.5249 | 929 | 1,066 | 618 |
| `readlik-nolink` | 0.5941 | 0.4981 | 0.5419 | 919 | 928 | 628 |
| `readlik` | 0.5837 | 0.5497 | **0.5662** | 903 | 730 | 644 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and the old cap of 0.1 was telling the model 0.1: overriding the mapper at exactly the sites that matter. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless — see 

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9750 | 0.9870 | 0.9215 | 0.9423 | 0.9031 |

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
| readlik-support | GT | ALL | ALL | ALL | 279965 | 270437 | 9528 | 279622 | 270850 | 8772 | 0.9659671744682371 | 0.9686290778264943 | 0.9672962948273991 | 1804 | 1502 |
| readlik-support | GT | ALL | ALL | Snv | 223585 | 218790 | 4795 | 213446 | 212023 | 1423 | 0.9785540174877564 | 0.9933332083993142 | 0.9858882283147561 | 201 | 774 |
| readlik-support | GT | ALL | ALL | Insertion | 28075 | 25325 | 2750 | 30912 | 27420 | 3492 | 0.9020480854853072 | 0.8870341614906833 | 0.8944781252904285 | 926 | 355 |
| readlik-support | GT | ALL | ALL | Deletion | 28305 | 26322 | 1983 | 31847 | 28332 | 3515 | 0.9299417064122947 | 0.8896285364398531 | 0.9093385457361243 | 677 | 335 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3417 | 3075 | 342 |  | 0.8999122036874452 |  | 0 | 38 |
| readlik-support | GT | ALL | ALL | JointIndel | 56380 | 51647 | 4733 | 66176 | 58827 | 7349 | 0.9160517914153955 | 0.8889476547388782 | 0.9022962232293822 | 1603 | 728 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 1115456 | 1079809 | 35647 | 1347668 | 1079809 | 267859 | 0.9680426659590338 | 0.8012425909051785 | 0.8767800565460774 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 602392 | 591623 | 10769 | 580120 | 575156 | 4964 | 0.9821229365595825 | 0.9914431496931669 | 0.986761035661446 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 257406 | 237277 | 20129 | 438304 | 237283 | 201021 | 0.9218005796290685 | 0.541366266335694 | 0.6821255408788638 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 267064 | 251744 | 15320 | 332400 | 251074 | 81326 | 0.942635473145014 | 0.7553369434416366 | 0.8386560230422692 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 36456 | 21095 | 15361 |  | 0.5786427474215493 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 524470 | 489021 | 35449 | 807160 | 509452 | 297708 | 0.9324098613838733 | 0.6311660637296199 | 0.7527686408315898 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 271457 | 8508 | 279971 | 271864 | 8107 | 0.9696104870251638 | 0.9710434294980551 | 0.9703264292326793 | 1165 | 789 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 218918 | 4667 | 216097 | 212940 | 3157 | 0.9791265066976765 | 0.9853908198633021 | 0.9822486756822874 | 458 | 182 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 25869 | 2206 | 30011 | 27467 | 2544 | 0.9214247551202137 | 0.9152310819366233 | 0.9183174752034549 | 371 | 256 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 26670 | 1635 | 30236 | 28317 | 1919 | 0.9422363540010599 | 0.9365326101336156 | 0.9393758240857746 | 336 | 328 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3627 | 3140 | 487 |  | 0.865729252826027 |  | 0 | 23 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 52539 | 3841 | 63874 | 58924 | 4950 | 0.9318730046115644 | 0.922503679118264 | 0.9271646723858462 | 707 | 607 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115518 | 1084576 | 30942 | 1422364 | 1084576 | 337788 | 0.9722622136083864 | 0.7625164866377383 | 0.8547095570243218 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 591831 | 10561 | 583550 | 575104 | 8446 | 0.9824682266696769 | 0.9855265187216177 | 0.9839949963807977 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 240116 | 17290 | 531772 | 238963 | 292809 | 0.9328298485660785 | 0.4493711590681721 | 0.6065497390730433 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 253219 | 13845 | 300606 | 251935 | 48671 | 0.9481584938441722 | 0.8380903907440304 | 0.8897332610806119 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 44870 | 22334 | 22536 |  | 0.49774905281925563 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 493335 | 31135 | 877248 | 513232 | 364016 | 0.9406353080252445 | 0.5850477858028744 | 0.7214035554754228 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 279965 | 270506 | 9459 | 279664 | 270889 | 8775 | 0.966213633847088 | 0.9686230619600663 | 0.9674168476881212 | 1806 | 1516 |
| readlik-nolink | GT | ALL | ALL | Snv | 223585 | 218814 | 4771 | 213448 | 212025 | 1423 | 0.9786613592146164 | 0.993333270866909 | 0.9859427345190287 | 200 | 777 |
| readlik-nolink | GT | ALL | ALL | Insertion | 28075 | 25350 | 2725 | 30918 | 27422 | 3496 | 0.9029385574354408 | 0.8869267093602432 | 0.8948610136833475 | 927 | 364 |
| readlik-nolink | GT | ALL | ALL | Deletion | 28305 | 26342 | 1963 | 31859 | 28350 | 3509 | 0.9306482953541777 | 0.8898584387457233 | 0.9097964029610589 | 679 | 334 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3439 | 3092 | 347 |  | 0.8990985751671998 |  | 0 | 41 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 56380 | 51692 | 4688 | 66216 | 58864 | 7352 | 0.9168499467896417 | 0.8889694333695783 | 0.9026944628434016 | 1606 | 739 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 1115480 | 1080544 | 34936 | 1340222 | 1080544 | 259678 | 0.96868074730161 | 0.8062425478763966 | 0.8800286028190717 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 591629 | 10763 | 580124 | 575163 | 4961 | 0.9821328968512198 | 0.9914483800015169 | 0.9867686534624499 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 237566 | 19840 | 438878 | 237048 | 201830 | 0.9229233195807401 | 0.5401227676028418 | 0.6814438752462909 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 252088 | 14976 | 331464 | 251369 | 80095 | 0.9439235539046821 | 0.7583598822194869 | 0.8410276925365553 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 29546 | 21988 | 7558 |  | 0.7441954917755365 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 489654 | 34816 | 799888 | 510405 | 289483 | 0.9336167940968978 | 0.6380955833816734 | 0.7580735017686909 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 271234 | 8731 | 276105 | 270915 | 5190 | 0.9688139588877182 | 0.9812028032813603 | 0.9749690266881085 | 785 | 1066 |
| readlik | GT | ALL | ALL | Snv | 223585 | 218705 | 4880 | 213125 | 212273 | 852 | 0.9781738488717937 | 0.9960023460410558 | 0.9870075941781121 | 171 | 383 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 25842 | 2233 | 29624 | 27327 | 2297 | 0.9204630454140694 | 0.9224615176883608 | 0.9214611979769031 | 331 | 301 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 26687 | 1618 | 29988 | 28239 | 1749 | 0.9428369546016605 | 0.9416766706682673 | 0.9422564554449738 | 283 | 342 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3368 | 3076 | 292 |  | 0.9133016627078385 |  | 0 | 40 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 52529 | 3851 | 62980 | 58642 | 4338 | 0.9316956367506208 | 0.9311209907907272 | 0.93140822513662 | 614 | 683 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115522 | 1083701 | 31821 | 1284386 | 1083701 | 200685 | 0.9714743411604612 | 0.8437502433069186 | 0.9031187862201384 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 591383 | 11009 | 578540 | 574740 | 3800 | 0.9817245248940889 | 0.9934317419711688 | 0.9875434377141207 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 239104 | 18302 | 403310 | 238131 | 165179 | 0.9288983162785638 | 0.590441595794798 | 0.7219716928862377 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 253834 | 13230 | 308702 | 252528 | 56174 | 0.9504613126441602 | 0.8180316292087515 | 0.8792881189195692 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 29108 | 22107 | 7001 |  | 0.7594819293664972 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 492938 | 31532 | 741120 | 512766 | 228354 | 0.9398783533853223 | 0.6918798575129533 | 0.7970334046755224 |  |  |

</details>

