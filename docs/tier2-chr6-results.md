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

Every number on this page — accuracy and cost alike — comes from one `vg` build in one pass, which is what the refresh harness exists to guarantee: a table whose rows come from different builds is not a comparison, it is a mixture of vintages.

Build: `vg version v1.4.0-18654-g648296d56`.

The wall column is what the caller costs unaided, and the repeatability note below applies to it harder than to the memory column. It includes snarl decomposition, which is single-threaded — 46 s of a 197 s chr20 run — and which `vg call -r` skips for byte-identical output given `vg snarls -T -P <ref path>`. The whole-genome harness caches one snarl file per contig for exactly that reason; this matrix does not, so these figures include it.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 294,626 | 666 s | 6.0 GB |
| `poisson-z` | panel (`-z`) | yes | 294,835 | 198 s | 6.1 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 299,877 | 350 s | 8.2 GB |
| `readlik-nomismap` | panel (default) | **no** | 303,729 | 356 s | 8.5 GB |
| `readlik-nolink` | panel (default) | **no** | 299,880 | 277 s | 8.0 GB |
| `readlik` | panel (default) | **no** | 296,793 | 366 s | 7.6 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9573 | 0.9800 | 0.9685 | 214,033 | 9,552 | 4,235 | 0.9639 | 0.9825 | 0.9731 |
| `poisson-z` | 0.9607 | 0.9801 | 0.9703 | 214,796 | 8,789 | 4,230 | 0.9671 | 0.9824 | 0.9747 |
| `readlik-support` | 0.9788 | 0.9928 | 0.9858 | 218,844 | 4,741 | 1,530 | 0.9822 | 0.9908 | 0.9865 |
| `readlik-nomismap` | 0.9803 | 0.9858 | 0.9830 | 219,183 | 4,402 | 3,085 | 0.9831 | 0.9853 | 0.9842 |
| `readlik-nolink` | 0.9789 | 0.9929 | 0.9859 | 218,878 | 4,707 | 1,519 | 0.9822 | 0.9908 | 0.9865 |
| `readlik` | 0.9800 | 0.9962 | **0.9880** | 219,116 | 4,469 | 819 | 0.9828 | 0.9930 | 0.9879 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7744 | 0.8529 | 0.8118 | 21,742 | 6,333 | 4,279 | 0.7996 | 0.6375 | 0.7094 |
| `poisson-z` | 0.7792 | 0.8558 | 0.8157 | 21,876 | 6,199 | 4,192 | 0.8041 | 0.5959 | 0.6845 |
| `readlik-support` | 0.9025 | 0.8858 | 0.8941 | 25,339 | 2,736 | 3,547 | 0.9218 | 0.5409 | 0.6817 |
| `readlik-nomismap` | 0.9312 | 0.9201 | 0.9256 | 26,144 | 1,931 | 2,427 | 0.9419 | 0.4559 | 0.6144 |
| `readlik-nolink` | 0.9031 | 0.8857 | 0.8943 | 25,355 | 2,720 | 3,550 | 0.9228 | 0.5397 | 0.6811 |
| `readlik` | 0.9307 | 0.9267 | **0.9287** | 26,128 | 1,947 | 2,204 | 0.9407 | 0.5964 | 0.7300 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8557 | 0.7206 | 0.7824 | 24,222 | 4,083 | 10,150 | 0.9292 | 0.5149 | 0.6626 |
| `poisson-z` | 0.8593 | 0.7229 | 0.7852 | 24,323 | 3,982 | 10,083 | 0.9306 | 0.4736 | 0.6278 |
| `readlik-support` | 0.9306 | 0.8880 | 0.9088 | 26,340 | 1,965 | 3,589 | 0.9423 | 0.7530 | 0.8371 |
| `readlik-nomismap` | 0.9498 | 0.9443 | 0.9470 | 26,883 | 1,422 | 1,702 | 0.9524 | 0.8590 | 0.9033 |
| `readlik-nolink` | 0.9315 | 0.8893 | 0.9099 | 26,367 | 1,938 | 3,546 | 0.9433 | 0.7568 | 0.8398 |
| `readlik` | 0.9506 | 0.9491 | **0.9498** | 26,906 | 1,399 | 1,543 | 0.9541 | 0.8259 | 0.8854 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.7865 | — | 0 | 0 | 825 | — | 0.3420 | — |
| `poisson-z` | — | 0.7915 | — | 0 | 0 | 836 | — | 0.4509 | — |
| `readlik-support` | — | 0.9183 | — | 0 | 0 | 262 | — | 0.6035 | — |
| `readlik-nomismap` | — | 0.8879 | — | 0 | 0 | 380 | — | 0.5060 | — |
| `readlik-nolink` | — | 0.9174 | — | 0 | 0 | 266 | — | 0.7647 | — |
| `readlik` | — | 0.9434 | — | 0 | 0 | 178 | — | 0.8104 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8153 | 0.7798 | 0.7972 | 45,964 | 10,416 | 15,254 | 0.8656 | 0.5453 | 0.6691 |
| `poisson-z` | 0.8194 | 0.7825 | 0.8005 | 46,199 | 10,181 | 15,111 | 0.8685 | 0.5176 | 0.6486 |
| `readlik-support` | 0.9166 | 0.8884 | 0.9023 | 51,679 | 4,701 | 7,398 | 0.9323 | 0.6314 | 0.7529 |
| `readlik-nomismap` | 0.9405 | 0.9299 | 0.9352 | 53,027 | 3,353 | 4,509 | 0.9472 | 0.5959 | 0.7316 |
| `readlik-nolink` | 0.9174 | 0.8890 | 0.9030 | 51,722 | 4,658 | 7,362 | 0.9333 | 0.6373 | 0.7574 |
| `readlik` | 0.9407 | 0.9382 | **0.9394** | 53,034 | 3,346 | 3,925 | 0.9475 | 0.6991 | 0.8046 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9287 | 0.9308 | 0.9297 | 259,997 | 19,968 | 19,489 | 0.9264 | 0.7386 | 0.8219 |
| `poisson-z` | 0.9322 | 0.9313 | 0.9318 | 260,995 | 18,970 | 19,341 | 0.9296 | 0.7139 | 0.8076 |
| `readlik-support` | 0.9663 | 0.9681 | 0.9672 | 270,523 | 9,442 | 8,928 | 0.9680 | 0.8026 | 0.8776 |
| `readlik-nomismap` | 0.9723 | 0.9730 | 0.9726 | 272,210 | 7,755 | 7,594 | 0.9758 | 0.7720 | 0.8620 |
| `readlik-nolink` | 0.9665 | 0.9683 | 0.9674 | 270,600 | 9,365 | 8,881 | 0.9685 | 0.8070 | 0.8804 |
| `readlik` | 0.9721 | 0.9829 | **0.9775** | 272,150 | 7,815 | 4,744 | 0.9759 | 0.8506 | 0.9090 |

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
| `sm50-readlik` | Insertion | 0.9387 | 0.9104 | **0.9244** |
| `sm50-readlik` | Deletion | 0.9517 | 0.9184 | **0.9347** |
| `sm50-readlik` | ALL | 0.9747 | 0.9798 | **0.9772** |

The insertion BASEPAIR precision gap collapses from **-0.001 to -0.046**, and insertion BASEPAIR F1 goes from 0.8320 for `poisson-z` against 0.9244 for `readlik`.
There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22). The aardvark block below is kept for continuity with earlier runs.

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5320 | 0.4618 | 0.4944 | 823 | 936 | 724 |
| `poisson-z` | 0.5417 | 0.4442 | 0.4881 | 838 | 1,036 | 709 |
| `readlik-support` | 0.5779 | 0.5092 | 0.5414 | 894 | 854 | 653 |
| `readlik-nomismap` | 0.6167 | 0.4831 | 0.5418 | 954 | 1,027 | 593 |
| `readlik-nolink` | 0.5973 | 0.4928 | 0.5400 | 924 | 953 | 623 |
| `readlik` | 0.6063 | 0.5627 | **0.5837** | 938 | 729 | 609 |


Kept for continuity. These categories are scored against the small-variant truth set and should not be read as the SV result; prefer the truvari table above.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph 23.3% of reads sit at MAPQ 1, meaning p(wrong) = 0.79, and the old cap of 0.1 was telling the model 0.1: overriding the mapper at exactly the sites that matter. Raising it removed 94% of the excess false-positive SNVs, and the default is now **0.7**. A clamp that is inert on a sparse graph is not thereby harmless — see 

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.7 (current defaults)** | 0.9775 | 0.9880 | 0.9287 | 0.9498 | 0.9090 |

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
| readlik-support | GT | ALL | ALL | ALL | 279965 | 270523 | 9442 | 280149 | 271221 | 8928 | 0.9662743557230368 | 0.9681312444449204 | 0.967201908844901 | 1787 | 1593 |
| readlik-support | GT | ALL | ALL | Snv | 223585 | 218844 | 4741 | 213837 | 212307 | 1530 | 0.9787955363731914 | 0.992845017466575 | 0.9857702201489853 | 185 | 828 |
| readlik-support | GT | ALL | ALL | Insertion | 28075 | 25339 | 2736 | 31063 | 27516 | 3547 | 0.902546749777382 | 0.8858127032160448 | 0.8941014345420819 | 926 | 373 |
| readlik-support | GT | ALL | ALL | Deletion | 28305 | 26340 | 1965 | 32041 | 28452 | 3589 | 0.9305776364599894 | 0.8879872663150339 | 0.9087837230698376 | 676 | 358 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3208 | 2946 | 262 |  | 0.9183291770573566 |  | 0 | 34 |
| readlik-support | GT | ALL | ALL | JointIndel | 56380 | 51679 | 4701 | 66312 | 58914 | 7398 | 0.916619368570415 | 0.8884364820846905 | 0.9023079113345289 | 1602 | 765 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 1115428 | 1079775 | 35653 | 1345394 | 1079775 | 265619 | 0.9680364846498385 | 0.8025715886944642 | 0.8775726159795385 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 602392 | 591640 | 10752 | 581228 | 575898 | 5330 | 0.9821511573858882 | 0.9908297604382446 | 0.9864713715123201 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 257406 | 237275 | 20131 | 439730 | 237830 | 201900 | 0.9217928098024133 | 0.5408546153321356 | 0.6817171206051255 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 267064 | 251667 | 15397 | 334482 | 251854 | 82628 | 0.9423471527424138 | 0.7529672747711387 | 0.8370796071493057 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 30910 | 18655 | 12255 |  | 0.6035263668715626 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 524470 | 488942 | 35528 | 805122 | 508339 | 296783 | 0.932259233130589 | 0.6313813310280927 | 0.7528726089220166 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 272210 | 7755 | 280868 | 273274 | 7594 | 0.9723001089421892 | 0.9729623880256918 | 0.9726311357450208 | 1207 | 668 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 219183 | 4402 | 216562 | 213477 | 3085 | 0.9803117382650893 | 0.9857546568650086 | 0.9830256633935279 | 446 | 196 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 26144 | 1931 | 30385 | 27958 | 2427 | 0.9312199465716829 | 0.9201250617080796 | 0.9256392589938967 | 391 | 235 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 26883 | 1422 | 30532 | 28830 | 1702 | 0.9497615262321145 | 0.9442552076509891 | 0.9470003629087007 | 370 | 228 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3389 | 3009 | 380 |  | 0.8878725287695486 |  | 0 | 9 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 53027 | 3353 | 64306 | 59797 | 4509 | 0.9405285562256119 | 0.9298821260846577 | 0.9351750412654318 | 761 | 472 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115432 | 1088445 | 26987 | 1409904 | 1088445 | 321459 | 0.9758057864576236 | 0.771999370169884 | 0.862019945068696 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 592203 | 10189 | 584720 | 576142 | 8578 | 0.9830857647511919 | 0.9853297304692844 | 0.9842064685657081 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 242446 | 14960 | 530632 | 241926 | 288706 | 0.9418816966193484 | 0.4559204872680125 | 0.6144262285773047 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 254347 | 12717 | 295622 | 253946 | 41676 | 0.9523822005212234 | 0.8590226708431713 | 0.9032965677507663 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 38612 | 19539 | 19073 |  | 0.5060343934528126 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 496793 | 27677 | 864866 | 515411 | 349455 | 0.9472286308082445 | 0.5959431865745676 | 0.7316028484966772 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 279965 | 270600 | 9365 | 280179 | 271298 | 8881 | 0.9665493901023342 | 0.9683024066757323 | 0.9674251042540221 | 1795 | 1597 |
| readlik-nolink | GT | ALL | ALL | Snv | 223585 | 218878 | 4707 | 213847 | 212328 | 1519 | 0.9789476038195765 | 0.9928967906961519 | 0.9858728577118383 | 187 | 821 |
| readlik-nolink | GT | ALL | ALL | Insertion | 28075 | 25355 | 2720 | 31072 | 27522 | 3550 | 0.9031166518254675 | 0.885749227600412 | 0.8943486328267493 | 930 | 384 |
| readlik-nolink | GT | ALL | ALL | Deletion | 28305 | 26367 | 1938 | 32040 | 28494 | 3546 | 0.9315315315315316 | 0.8893258426966292 | 0.9099395438689145 | 678 | 358 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3220 | 2954 | 266 |  | 0.9173913043478261 |  | 0 | 34 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 56380 | 51722 | 4658 | 66332 | 58970 | 7362 | 0.9173820503724726 | 0.889012844479286 | 0.9029746800107078 | 1608 | 776 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 1115436 | 1080348 | 35088 | 1338708 | 1080348 | 258360 | 0.9685432422837348 | 0.8070079509497217 | 0.8804275543733374 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 591646 | 10746 | 581262 | 575916 | 5346 | 0.9821611176775256 | 0.9908027705234472 | 0.9864630187252283 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 237531 | 19875 | 440312 | 237634 | 202678 | 0.9227873476142747 | 0.5396945802067625 | 0.6810659615231027 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 251934 | 15130 | 333154 | 252120 | 81034 | 0.9433469130994818 | 0.7567671407217083 | 0.8398188868921951 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 25208 | 19277 | 5931 |  | 0.764717549984132 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 489465 | 35005 | 798674 | 509031 | 289643 | 0.9332564303010659 | 0.6373451495854379 | 0.7574250106316424 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 272150 | 7815 | 277353 | 272609 | 4744 | 0.9720857964388405 | 0.9828954437125252 | 0.9774607352686623 | 841 | 872 |
| readlik | GT | ALL | ALL | Snv | 223585 | 219116 | 4469 | 213832 | 213013 | 819 | 0.9800120759442718 | 0.9961698903812338 | 0.9880249277660768 | 167 | 336 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 26128 | 1947 | 30056 | 27852 | 2204 | 0.9306500445235975 | 0.9266702155975513 | 0.9286558661114266 | 357 | 272 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 26906 | 1399 | 30319 | 28776 | 1543 | 0.9505741035152799 | 0.9491078201787658 | 0.9498403959662836 | 317 | 250 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3146 | 2968 | 178 |  | 0.9434202161474888 |  | 0 | 14 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 53034 | 3346 | 63521 | 59596 | 3925 | 0.9406527137282724 | 0.9382094110609089 | 0.9394294737394558 | 674 | 536 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115436 | 1088504 | 26932 | 1279620 | 1088504 | 191116 | 0.9758551812923377 | 0.8506462856160423 | 0.9089591224589322 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 592060 | 10332 | 580092 | 576055 | 4037 | 0.982848377800502 | 0.9930407590520124 | 0.9879182803532984 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 242130 | 15276 | 405472 | 241815 | 163657 | 0.940654064007832 | 0.5963790348038829 | 0.7299600291119833 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 254810 | 12254 | 307926 | 254322 | 53604 | 0.9541158673576371 | 0.8259192143566961 | 0.8854012324457571 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 23510 | 19052 | 4458 |  | 0.8103785623139089 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 496940 | 27530 | 736908 | 515189 | 221719 | 0.9475089137605582 | 0.6991225498976806 | 0.8045817931443954 |  |  |

</details>

