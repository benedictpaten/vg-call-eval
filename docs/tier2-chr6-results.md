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

**All read-likelihood arms below run at the current clamp defaults, `--mismap-min 0.02` and `--mismap-max 0.95`.** The floor caps how much one read can veto an allele; the cap bounds how far a low-MAPQ read is discounted. Both were set by measurement — the floor from 1e-8, the cap raised from an original 0.1 that was actively wrong on haplotype-rich graphs — and the sweeps are in harness plan §9.20-§9.21. `poisson` and `poisson-z` do not use the read-likelihood model, so neither reaches them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

Every number on this page — accuracy and cost alike — comes from one `vg` build in one pass, which is what the refresh harness exists to guarantee: a table whose rows come from different builds is not a comparison, it is a mixture of vintages.

Build: `vg version v1.4.0-18924-g8acbb43a2` (the code of `0cab3fbd4`, built before that commit, so the string names its parent).

The wall column is what the caller costs unaided, and the repeatability note below applies to it harder than to the memory column. It includes snarl decomposition, which is single-threaded — 46 s of a 197 s chr20 run — and which `vg call -r` skips for byte-identical output given `vg snarls -T -P <ref path>`. The whole-genome harness caches one snarl file per contig for exactly that reason; this matrix does not, so these figures include it.

| arm | enumeration | pack? | variants | wall | CPU | Δ wall | Δ CPU | peak RSS |
|---|---|---|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 294,626 | 600 s | 2,453 s (4.1x) | -67 s | -243 s | 7.0 GB |
| `poisson-z` | panel (`-z`) | yes | 294,835 | 192 s | 531 s (2.8x) | -6 s | +7 s | 6.2 GB |
| `readlik-support` | support (`--enumerate-support`) | yes | 299,876 | 319 s | 1,144 s (3.6x) | -31 s | -88 s | 12.9 GB |
| `readlik-nomismap` | panel (default) | **no** | 303,771 | 297 s | 1,046 s (3.5x) | -59 s | -98 s | 13.2 GB |
| `readlik-nolink` | panel (default) | **no** | 299,824 | 254 s | 912 s (3.6x) | -23 s | -68 s | 12.7 GB |
| `readlik` | panel (default) | **no** | 296,566 | 295 s | 1,029 s (3.5x) | -71 s | -152 s | 12.7 GB |

`CPU` is user+sys, with the multiple of wall clock beside it. It is the column that separates work from waiting: this caller has phases that run on one thread and phases that block on a subprocess, so a wall-clock change can come from either doing less or waiting less, and only CPU distinguishes them. A multiple well under `--threads` means the run spent its time parked rather than computing.

`Δ wall` is against `vg version v1.4.0-18654-g648296d56`. Read it with the repeatability note below: run-to-run variance on this measurement is larger than most of these deltas, and the Poisson arms are the control -- their code is untouched by any read-likelihood change, so a Δ on those rows is the machine and not the caller.

**Peak RSS in this table is repeatable to about ±0.35 GB, so read it accordingly.** Three back-to-back runs of one binary on chr6-4hap, identical parameters and a warm cache, gave 7.3, 6.6 and 7.0 GB -- a 0.7 GB spread on a 7 GB measurement. Differences smaller than that are not evidence of anything, and a single measurement of each of two arms cannot resolve one. Thread count matters too: the same run at `--threads 6` instead of 5 measured 8.7 GB, because the read and GBWT caches are per thread. Wall clock is worse still -- a run immediately after a full rebuild took 956 s against 260 s warm, purely from page cache.

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9573 | 0.9800 | 0.9685 | 214,033 | 9,552 | 4,235 | 0.9639 | 0.9825 | 0.9731 |
| `poisson-z` | 0.9607 | 0.9801 | 0.9703 | 214,796 | 8,789 | 4,230 | 0.9671 | 0.9824 | 0.9747 |
| `readlik-support` | 0.9786 | 0.9926 | 0.9856 | 218,791 | 4,794 | 1,573 | 0.9820 | 0.9906 | 0.9863 |
| `readlik-nomismap` | 0.9803 | 0.9858 | 0.9830 | 219,182 | 4,403 | 3,082 | 0.9831 | 0.9853 | 0.9842 |
| `readlik-nolink` | 0.9787 | 0.9927 | 0.9857 | 218,829 | 4,756 | 1,561 | 0.9820 | 0.9906 | 0.9863 |
| `readlik` | 0.9797 | 0.9963 | **0.9879** | 219,048 | 4,537 | 799 | 0.9826 | 0.9931 | 0.9879 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7744 | 0.8529 | 0.8118 | 21,742 | 6,333 | 4,279 | 0.7996 | 0.6375 | 0.7094 |
| `poisson-z` | 0.7792 | 0.8558 | 0.8157 | 21,876 | 6,199 | 4,192 | 0.8041 | 0.5959 | 0.6845 |
| `readlik-support` | 0.9022 | 0.8866 | 0.8943 | 25,328 | 2,747 | 3,517 | 0.9206 | 0.6139 | 0.7366 |
| `readlik-nomismap` | 0.9311 | 0.9204 | 0.9257 | 26,142 | 1,933 | 2,419 | 0.9415 | 0.4557 | 0.6141 |
| `readlik-nolink` | 0.9027 | 0.8863 | 0.8945 | 25,344 | 2,731 | 3,528 | 0.9221 | 0.6106 | 0.7347 |
| `readlik` | 0.9304 | 0.9270 | **0.9287** | 26,122 | 1,953 | 2,193 | 0.9403 | 0.7300 | 0.8219 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8557 | 0.7206 | 0.7824 | 24,222 | 4,083 | 10,150 | 0.9292 | 0.5149 | 0.6626 |
| `poisson-z` | 0.8593 | 0.7229 | 0.7852 | 24,323 | 3,982 | 10,083 | 0.9306 | 0.4736 | 0.6278 |
| `readlik-support` | 0.9310 | 0.8869 | 0.9084 | 26,352 | 1,953 | 3,629 | 0.9431 | 0.7487 | 0.8347 |
| `readlik-nomismap` | 0.9499 | 0.9445 | 0.9472 | 26,887 | 1,418 | 1,694 | 0.9525 | 0.8617 | 0.9049 |
| `readlik-nolink` | 0.9320 | 0.8883 | 0.9096 | 26,379 | 1,926 | 3,584 | 0.9438 | 0.7501 | 0.8359 |
| `readlik` | 0.9510 | 0.9491 | **0.9500** | 26,917 | 1,388 | 1,543 | 0.9542 | 0.8569 | 0.9030 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8153 | 0.7798 | 0.7972 | 45,964 | 10,416 | 15,254 | 0.8656 | 0.5453 | 0.6691 |
| `poisson-z` | 0.8194 | 0.7825 | 0.8005 | 46,199 | 10,181 | 15,111 | 0.8685 | 0.5176 | 0.6486 |
| `readlik-support` | 0.9166 | 0.8883 | 0.9022 | 51,680 | 4,700 | 7,409 | 0.9321 | 0.6761 | 0.7837 |
| `readlik-nomismap` | 0.9406 | 0.9301 | 0.9353 | 53,029 | 3,351 | 4,497 | 0.9471 | 0.5964 | 0.7319 |
| `readlik-nolink` | 0.9174 | 0.8888 | 0.9029 | 51,723 | 4,657 | 7,375 | 0.9332 | 0.6788 | 0.7859 |
| `readlik` | 0.9407 | 0.9384 | **0.9396** | 53,039 | 3,341 | 3,910 | 0.9474 | 0.7908 | 0.8620 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9287 | 0.9308 | 0.9297 | 259,997 | 19,968 | 19,489 | 0.9264 | 0.7386 | 0.8219 |
| `poisson-z` | 0.9322 | 0.9313 | 0.9318 | 260,995 | 18,970 | 19,341 | 0.9296 | 0.7139 | 0.8076 |
| `readlik-support` | 0.9661 | 0.9679 | 0.9670 | 270,471 | 9,494 | 8,982 | 0.9679 | 0.8356 | 0.8969 |
| `readlik-nomismap` | 0.9723 | 0.9730 | 0.9727 | 272,211 | 7,754 | 7,579 | 0.9758 | 0.7724 | 0.8623 |
| `readlik-nolink` | 0.9664 | 0.9681 | 0.9672 | 270,552 | 9,413 | 8,936 | 0.9684 | 0.8377 | 0.8983 |
| `readlik` | 0.9719 | 0.9830 | **0.9774** | 272,087 | 7,878 | 4,709 | 0.9757 | 0.9115 | 0.9425 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates both callers, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (167.2 Mb vs 168.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

316 `readlik` calls carry an insertion allele of 200 bp or more, against 287 from `poisson-z`. Every base of those alleles inside the small-variant confident region is scored FP, and the size-matched control below measures what that does to each caller's precision.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.8019 | 0.8645 | **0.8320** |
| `sm50-poisson-z` | Deletion | 0.9239 | 0.7428 | **0.8235** |
| `sm50-poisson-z` | ALL | 0.9271 | 0.9091 | **0.9180** |
| `sm50-readlik` | Insertion | 0.9384 | 0.9117 | **0.9248** |
| `sm50-readlik` | Deletion | 0.9521 | 0.9183 | **0.9349** |
| `sm50-readlik` | ALL | 0.9746 | 0.9802 | **0.9774** |

Restricting raises insertion BASEPAIR precision from 0.7300 to 0.9117 for `readlik` and from 0.5959 to 0.8645 for `poisson-z`. `readlik` minus `poisson-z` goes from +0.134 to +0.047, so the unrestricted comparison overstates the difference between them. Insertion BASEPAIR F1 is 0.8320 for `poisson-z` against 0.9248 for `readlik` restricted, and 0.6845 against 0.8219 unrestricted.

Whether those large calls are *correct* is a separate question, and the truvari comparison below is what answers it. Split into single alleles there are 352 insertions of 200 bp or more; truvari confirms **183** and rejects **91**, and the other **78** fall outside the SV confident region or above truvari's 50 kb size cap, so it does not judge them. *Known bad output* lists the largest.

## Structural variants — truvari (GIAB `stvar` benchmark)

The SV metric. Reciprocal-overlap matching, `--sizemin 50`. It replaced aardvark's `Sv*` categories as the primary measure: those are scored against the *small-variant* truth set, which contains no record over 50 bp at all, so they have almost nothing to match (plan §9.22).

**What these errors are made of, per record, is in [tier2-sv-errors.md](tier2-sv-errors.md)** — including the finding that about a quarter of all false positives are the metric rather than the caller, and that harmonising representation with `truvari refine` moves every arm up by roughly 0.05 F1. Read the ranking between arms here; treat the absolute level as benchmark-relative.

| arm | recall | precision | **F1** | TP-base | FP | FN |
|---|---|---|---|---|---|---|
| `poisson` | 0.5320 | 0.4618 | 0.4944 | 823 | 936 | 724 |
| `poisson-z` | 0.5417 | 0.4442 | 0.4881 | 838 | 1,036 | 709 |
| `readlik-support` | 0.5818 | 0.5068 | 0.5417 | 900 | 867 | 647 |
| `readlik-nomismap` | 0.6173 | 0.4834 | 0.5422 | 955 | 1,028 | 592 |
| `readlik-nolink` | 0.5960 | 0.4936 | 0.5400 | 922 | 947 | 625 |
| `readlik` | 0.6063 | 0.5670 | **0.5860** | 938 | 717 | 609 |

## Long reads — ONT, 16-haplotype E821 graph

Same sample, truth and confident regions as above; different reads and a different graph. The ONT reads (45.4x on chr6, mean length 33,449 bp on chr20) are aligned to the E821 graph: 16 haplotypes from haplotype sampling plus CHM13 and GRCh38, 18 in the panel against 34 above. So a figure here is not comparable with one in the short-read tables. The pair of columns is comparable, and it answers what the long-read preset buys on long-read data. The preset sets `--gap-open 1 --gap-extend 1 --mismap-min 0.05 --read-min-mapq 5 --insertion-nats 0.9 --read-phasing --regenotype`.

Build: `vg version v1.4.0-18924-g8acbb43a2` (the code of `0cab3fbd4`, built before that commit, so the string names its parent) — the same build as the short-read arms. Directories: `work/E821-chr6/results-0923-default`, `work/E821-chr6/results-0923-preset`.

| | short-read defaults | `--preset ont` | Δ |
|---|---|---|---|
| `readlik` ALL F1 | 0.9373 | 0.9646 | +0.0273 |
| `readlik` SNV F1 | 0.9871 | 0.9882 | +0.0011 |
| `readlik` indel F1 | 0.7774 | 0.8807 | +0.1033 |
| `readlik` indel recall | 0.8199 | 0.8880 | +0.0681 |
| `readlik` indel precision | 0.7391 | 0.8735 | +0.1344 |
| `readlik-nomismap` ALL F1 | 0.9373 | 0.9644 | +0.0271 |
| `readlik-nolink` ALL F1 | 0.9241 | 0.9345 | +0.0104 |
| `readlik` SV F1 (truvari) | 0.6056 | 0.6047 | -0.0009 |

Small variants are aardvark's GT comparison, with indels from its joint indel row.

| arm | configuration | variants | wall | CPU | peak RSS |
|---|---|---|---|---|---|
| `readlik` | short-read defaults | 312,247 | 295 s | 1,057 s | 6.8 GB |
| `readlik` | `--preset ont` | 302,392 | 374 s | 1,152 s | 7.7 GB |
| `readlik-nomismap` | short-read defaults | 312,427 | 290 s | 1,051 s | 7.3 GB |
| `readlik-nomismap` | `--preset ont` | 302,463 | 373 s | 1,144 s | 8.2 GB |
| `readlik-nolink` | short-read defaults | 315,662 | 278 s | 1,000 s | 7.4 GB |
| `readlik-nolink` | `--preset ont` | 308,390 | 273 s | 989 s | 8.2 GB |

Run serially on the same machine as the short-read arms, `--threads 5`.

*The MAPQ mismapping term.* `readlik` minus `readlik-nomismap` is 0.0000 ALL F1 on ONT at short-read defaults and +0.0002 under the preset, against +0.0047 on short reads. `--no-mismap-term` puts every read on the floor, so the term changes only a read whose `e_r = clamp(10^(−MAPQ/10), --mismap-min, --mismap-max)` is above it (MAPQ below 17 at the 0.02 floor), and it reduces a read on the ceiling to near-silence. At the default ceiling of 0.95 that is MAPQ 0 alone: 4.96% of chr20 short-read alignments against 0.63% of ONT, 8x fewer. Long reads anchor uniquely, so there is almost no ambiguous-placement class for the term to act on. Under the preset, `--read-min-mapq 5` drops those reads before scoring, so no read reaches the ceiling and the term acts only between MAPQ 5 and 13, where 10^(−MAPQ/10) lies above the preset's 0.05 floor.

*The linkage layer.* `readlik` minus `readlik-nolink` is +0.0132 on ONT at short-read defaults and +0.0302 under the preset, against +0.0102 on short reads. `--linkage-weight 0` turns the whole layer off, and read phasing and phase-driven re-genotyping run inside it, so under the preset that arm loses both as well as the transition model. Switch error across depth, and ONT's SV figures against short reads with confidence intervals, are in [coverage.md](coverage.md).

## Calibration: the two mismapping clamps

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. The floor caps that veto; the current default is **0.02**.

The *upper* clamp (`--mismap-max`) **binds hard here**, and looked inert on the 4-haplotype graph. There it reaches only reads whose `e_r` is already large — 6.3% of chr20 reads at MAPQ ≤ 9, against 90% at MAPQ 60 — so it appeared to be a knob with nothing to act on. On this graph the old cap of 0.1 was overriding the mapper at exactly the sites that matter: 23.3% of reads at those sites sit at MAPQ 1, meaning p(wrong) = 0.79, and were being told 0.1. Raising it to 0.5 removed 94% of the excess false-positive SNVs, and the default is now **0.95**. A clamp that is inert on a sparse graph is not thereby harmless.

The two graphs are put side by side in [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md); the grids are in plan §9.20.

| `readlik` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| **floor 0.02, cap 0.95 (current defaults)** | 0.9774 | 0.9879 | 0.9287 | 0.9500 | 0.9425 |

Only the current row is available here: the preserved old-default arms (`arms.floor-1e-8.json`, `arms.readlik.json`) exist for the chr20 4-haplotype run alone, so the before-and-after is on [tier2-chr20-4hap-results.md](tier2-chr20-4hap-results.md). Mixing rows from two datasets into one table is exactly what the one-build-per-matrix rule forbids. The full grids are in plan §9.20-§9.21.

The floor was later re-swept at the corrected cap, on both graphs and both benchmarks, and settled at **0.02**. 0.05 wins on small-variant `GT` but costs about 0.01 of SV F1 — which the first sweep never saw, because it was scored on one benchmark only. Plan §9.21 records that as a process rule: a sweep that sets a default has to be scored on every benchmark the project runs.

## Known bad output

Neither benchmark scores these, so they appear in no metric on this page. They are recorded because they are plainly wrong and would mislead anyone reading the VCF.

`readlik` calls **14 records carrying an insertion of 10 kb or more** on chr6 (15 such alleles), with median DP **2,649** against a median of **29** over all of the contig's records. 3 of them start inside another giant's reference span, so the records overstate the number of independent events. Length is the net insertion, ALT minus REF. The five largest:

| position | net insertion | REF length | GT | DP | GQ |
|---|---|---|---|---|---|
| chr6:32,344,855 | 64,151 bp | 65,217 bp | `2\|1` | 19,914 | 8 |
| chr6:32,362,682 | 38,583 bp | 9,493 bp | `1\|2` | 5,486 | 9 |
| chr6:56,902,800 | 36,712 bp | 8,891 bp | `1\|0` | 8,737 | 6 |
| chr6:31,837,918 | 32,738 bp | 1 bp | `1\|0` | 3,691 | 16 |
| chr6:31,881,914 | 32,735 bp | 1 bp | `0\|1` | 2,935 | 0 |

Median DP by called insertion length: 28 for 1 bp, 27 for 2-15 bp, 25 for 16-49 bp, 34 for 50-199 bp, 68 for 200-999 bp, 419 for >=1 kb. So the giants are collapsed-repeat pile-ups, not haplotypes.

The read-likelihood model cannot reject them, and the reason is structural rather than a tuning failure: it computes P(reads | genotype) **conditioned on the reads it is given**, and never asks whether that many reads should be there. The Poisson caller gets this for free, because an observed-vs-expected depth term is the whole of its model. A depth-plausibility guard is the obvious remedy, and the expected depth is already reachable — the read-likelihood caller subclasses `SupportBasedSnarlCaller` and holds a `TraversalSupportFinder` for allele enumeration.

The same blindness has a second consequence, found later and now corrected. Because the model only weighs reads it can see, it had no way to know that a heterozygous deletion produces *no* reads over the deleted interval, and its flat `1/ploidy` mixture asserted that both haplotypes contributed equally everywhere. That cost it 94% of heterozygous deletions above 1 kb and mis-genotyped two thirds of heterozygous insertions above 1 kb. Weighting each haplotype by the reads it is *expected* to contribute at the site is now the default and fixes both, without moving small variants at all — see [tier2-sv-errors.md](tier2-sv-errors.md). It did not remove the need for a depth term: it corrects the *relative* weight between a genotype's haplotypes, while the pile-ups above are a statement about *absolute* depth. That term is now also the default, at `--depth-term 0.1`, and the read arms in the tables on this page carry it — see [tier2-depth-term.md](tier2-depth-term.md). It does not resolve the pile-ups either: it detects them emphatically and still cannot outvote the read evidence at them, which is what the `DR` field and `--depth-quality` are for ([tier2-quality-signals.md](tier2-quality-signals.md)).

Filtering on depth is **not** that remedy, and that has now been tested properly rather than by two spot checks. Sweeping a two-sided cut on DP over a rolling local median, across both chromosomes and both graphs, against the one test a hard filter has to pass — beat lowering the GQ threshold to the same recall:

- a **minimum** fails in all eight dataset-by-benchmark cells. Few reads already means a small likelihood gap, so low depth depresses GQ on its own and a separate cut adds nothing;
- a **maximum** passes in exactly one configuration — 5x the local median, structural calls, 34-haplotype graph, worth about +0.025 precision — and is dominated everywhere else. The two original spot checks (DP 200 moving insertion BASEPAIR precision by 0.0001; DP 58 helping by +0.087 but costing SV insertion recall 0.4976 to 0.4167) were both right and both too narrow to conclude from.

What shipped instead attacks the same blindness from the other side: **GQ is scaled by the fraction of reads the called genotype explains**, which lowers the quality of a pile-up the call does not account for. It does not reach all of them: 2 of the 14 giants above still carry GQ 256. The giants remain output that no metric charges for, and they should be fixed because they are wrong, not because they cost a score. See [tier2-quality-signals.md](tier2-quality-signals.md).

## Quality fields

Every arm above is scored at **every** GQ, so nothing on this page depends on the quality field. `vg call` emits `AD` (per-allele read support, ties split fractionally), `BL` (mean absolute fit), `GQI` (the raw likelihood-ratio quality) and `GQ` (that ratio scaled by the fraction of reads the called genotype explains). The scaling rescales a quality and does not change a genotype, so **the numbers on this page are unaffected by it**; what it changes is how the calls rank. See [tier2-quality-signals.md](tier2-quality-signals.md).

## The genotype mixture

The read-likelihood arms on this page use the **length-weighted mixture**, which became the default after it was found that the flat `1/ploidy` weight breaks heterozygotes whose alleles differ in length. Unlike the `GQ` scaling above, this *does* change genotypes, so these numbers are not comparable with runs made before it. `--flat-mixture` restores the previous model exactly. Derivation and measurements: [tier2-sv-errors.md](tier2-sv-errors.md).

### At long read length

`w_h ∝ L_h + R − 1` counts the read START positions that yield an overlap, so R is the full read length; the *scoring* window is the read's overlap with the site, but the weight answers a sampling question, not a fit question. At ONT read length the correction all but vanishes. With the measured chr20 mean of 33,449 bp against Illumina's 151:

| | Illumina R=151 | ONT R=33,449 |
|---|---|---|
| 50 bp deletion | 1.33x, w = 0.571 | 1.00x, w = 0.500 |
| 300 bp Alu | 3.00x, w = 0.750 | 1.01x, w = 0.502 |
| 5 kb deletion | 34.3x, w = 0.972 | 1.15x, w = 0.535 |
| 30 kb deletion | 201x, w = 0.995 | 1.90x, w = 0.655 |

That is the model being right rather than distorted: a haplotype carrying a 5 kb deletion really does yield only 13% fewer 33 kb reads over the site, where it yields 34x fewer 151 bp reads. The flat mixture that `--flat-mixture` restores, and that loses large heterozygous deletions on short reads, is very nearly what the length weighting already computes for ONT.

The caveat is the mean, not the term. ONT read lengths are heavily skewed (mean 33,449, median 19,222, p10 1,845, p90 86,190), so one mean R summarises a distribution that spans two orders of magnitude. A per-read R is the obvious refinement and has not been measured.

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
| readlik-support | GT | ALL | ALL | ALL | 279965 | 270471 | 9494 | 280136 | 271154 | 8982 | 0.9660886182201347 | 0.9679370020275866 | 0.9670119268566836 | 1790 | 1615 |
| readlik-support | GT | ALL | ALL | Snv | 223585 | 218791 | 4794 | 213818 | 212245 | 1573 | 0.9785584900597089 | 0.9926432760572075 | 0.9855505631977505 | 187 | 840 |
| readlik-support | GT | ALL | ALL | Insertion | 28075 | 25328 | 2747 | 31020 | 27503 | 3517 | 0.9021549421193232 | 0.8866215344938749 | 0.8943207936717851 | 926 | 372 |
| readlik-support | GT | ALL | ALL | Deletion | 28305 | 26352 | 1953 | 32096 | 28467 | 3629 | 0.9310015898251193 | 0.8869329511465603 | 0.9084331354905436 | 677 | 368 |
| readlik-support | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3202 | 2939 | 263 |  | 0.9178638351030606 |  | 0 | 35 |
| readlik-support | GT | ALL | ALL | JointIndel | 56380 | 51680 | 4700 | 66318 | 58909 | 7409 | 0.9166371053565094 | 0.8882807081033807 | 0.9022361582870019 | 1603 | 775 |
| readlik-support | BASEPAIR | ALL | ALL | ALL | 1115426 | 1079606 | 35820 | 1291956 | 1079606 | 212350 | 0.96788670875522 | 0.835636817352913 | 0.8969129120347332 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Snv | 602392 | 591544 | 10848 | 581268 | 575810 | 5458 | 0.9819917927196908 | 0.9906101832545401 | 0.9862821609109198 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Insertion | 257406 | 236977 | 20429 | 386984 | 237574 | 149410 | 0.9206351056307934 | 0.613911686271267 | 0.7366196366521726 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Deletion | 267064 | 251865 | 15199 | 336660 | 252051 | 84609 | 0.9430885480633856 | 0.748681162003208 | 0.8347148265330043 |  |  |
| readlik-support | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 28240 | 18736 | 9504 |  | 0.6634560906515581 |  |  |  |
| readlik-support | BASEPAIR | ALL | ALL | JointIndel | 524470 | 488842 | 35628 | 751884 | 508361 | 243523 | 0.9320685644555456 | 0.6761162626149778 | 0.7837242382749524 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 279965 | 272211 | 7754 | 280870 | 273291 | 7579 | 0.972303680817245 | 0.9730159860433653 | 0.972659703020187 | 1210 | 662 |
| readlik-nomismap | GT | ALL | ALL | Snv | 223585 | 219182 | 4403 | 216562 | 213480 | 3082 | 0.9803072656931369 | 0.985768509710845 | 0.9830303027486045 | 446 | 193 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 28075 | 26142 | 1933 | 30382 | 27963 | 2419 | 0.9311487088156724 | 0.9203804884471068 | 0.9257332854417091 | 392 | 234 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 28305 | 26887 | 1418 | 30534 | 28840 | 1694 | 0.9499028440204911 | 0.9445208619899129 | 0.9472042080072838 | 372 | 226 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3392 | 3008 | 384 |  | 0.8867924528301887 |  | 0 | 9 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 56380 | 53029 | 3351 | 64308 | 59811 | 4497 | 0.9405640297978006 | 0.9300709087516328 | 0.9352880392702245 | 764 | 469 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 1115432 | 1088404 | 27028 | 1409084 | 1088404 | 320680 | 0.9757690293984752 | 0.7724195292828533 | 0.8622674603765633 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 602392 | 592200 | 10192 | 584712 | 576138 | 8574 | 0.9830807846053733 | 0.985336370726101 | 0.9842072853409028 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 257406 | 242360 | 15046 | 530700 | 241819 | 288881 | 0.9415475940731762 | 0.45566044846429243 | 0.6141189943148412 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 267064 | 254391 | 12673 | 294794 | 254034 | 40760 | 0.9525469550369948 | 0.8617339565934178 | 0.9048676543340668 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 38632 | 19504 | 19128 |  | 0.5048664319734935 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 524470 | 496751 | 27719 | 864126 | 515357 | 348769 | 0.9471485499647263 | 0.5963910355665725 | 0.7319163173316445 |  |  |

</details>

<details><summary><code>readlik-nolink</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nolink | GT | ALL | ALL | ALL | 279965 | 270552 | 9413 | 280169 | 271233 | 8936 | 0.9663779400996553 | 0.9681049652174223 | 0.9672406817507812 | 1800 | 1625 |
| readlik-nolink | GT | ALL | ALL | Snv | 223585 | 218829 | 4756 | 213820 | 212259 | 1561 | 0.9787284477939039 | 0.9926994668412683 | 0.9856644527499139 | 190 | 842 |
| readlik-nolink | GT | ALL | ALL | Insertion | 28075 | 25344 | 2731 | 31040 | 27512 | 3528 | 0.9027248441674087 | 0.886340206185567 | 0.8944574981779606 | 931 | 383 |
| readlik-nolink | GT | ALL | ALL | Deletion | 28305 | 26379 | 1926 | 32094 | 28510 | 3584 | 0.9319554848966614 | 0.8883280363930953 | 0.9096189425671806 | 679 | 365 |
| readlik-nolink | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3215 | 2952 | 263 |  | 0.9181959564541213 |  | 0 | 35 |
| readlik-nolink | GT | ALL | ALL | JointIndel | 56380 | 51723 | 4657 | 66349 | 58974 | 7375 | 0.9173997871585668 | 0.8888453480836185 | 0.9028968629326823 | 1610 | 783 |
| readlik-nolink | BASEPAIR | ALL | ALL | ALL | 1115434 | 1080201 | 35233 | 1289442 | 1080201 | 209241 | 0.9684131916366185 | 0.8377274821201729 | 0.8983423677561754 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Snv | 602392 | 591544 | 10848 | 581302 | 575818 | 5484 | 0.9819917927196908 | 0.9905660052778074 | 0.9862602640261601 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Insertion | 257406 | 237350 | 20056 | 388934 | 237479 | 151455 | 0.9220841783019821 | 0.6105894573372346 | 0.7346833206455544 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Deletion | 267064 | 252062 | 15002 | 336270 | 252223 | 84047 | 0.9438261989635444 | 0.7500609629167039 | 0.83586109340931 |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 24590 | 19268 | 5322 |  | 0.7835705571370476 |  |  |  |
| readlik-nolink | BASEPAIR | ALL | ALL | JointIndel | 524470 | 489412 | 35058 | 749794 | 508970 | 240824 | 0.9331553759032929 | 0.6788131140019792 | 0.7859187205350721 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 279965 | 272087 | 7878 | 277272 | 272563 | 4709 | 0.9718607683103245 | 0.9830166767650538 | 0.9774068907996527 | 847 | 873 |
| readlik | GT | ALL | ALL | Snv | 223585 | 219048 | 4537 | 213772 | 212973 | 799 | 0.9797079410515017 | 0.996262372995528 | 0.9879158115442227 | 166 | 319 |
| readlik | GT | ALL | ALL | Insertion | 28075 | 26122 | 1953 | 30037 | 27844 | 2193 | 0.9304363312555655 | 0.9269900456104139 | 0.9287099912982617 | 361 | 280 |
| readlik | GT | ALL | ALL | Deletion | 28305 | 26917 | 1388 | 30327 | 28784 | 1543 | 0.9509627274333157 | 0.9491212450951297 | 0.9500410939201867 | 320 | 258 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 3136 | 2962 | 174 |  | 0.9445153061224489 |  | 0 | 16 |
| readlik | GT | ALL | ALL | JointIndel | 56380 | 53039 | 3341 | 63500 | 59590 | 3910 | 0.9407413976587442 | 0.9384251968503937 | 0.9395818698169484 | 681 | 554 |
| readlik | BASEPAIR | ALL | ALL | ALL | 1115436 | 1088328 | 27108 | 1193958 | 1088328 | 105630 | 0.9756973954579196 | 0.911529551290749 | 0.9425225838466715 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 602392 | 591938 | 10454 | 579946 | 575957 | 3989 | 0.9826458518705428 | 0.9931217734064895 | 0.9878560399058771 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 257406 | 242033 | 15373 | 331136 | 241721 | 89415 | 0.9402772274150564 | 0.7299749951681485 | 0.8218864854763973 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 267064 | 254844 | 12220 | 296856 | 254379 | 42477 | 0.9542431776652788 | 0.8569104212143261 | 0.9029614316752459 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 23290 | 18937 | 4353 |  | 0.8130957492486045 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 524470 | 496877 | 27593 | 651282 | 515037 | 136245 | 0.947388792495281 | 0.7908049047877878 | 0.8620439770519305 |  |  |

</details>

