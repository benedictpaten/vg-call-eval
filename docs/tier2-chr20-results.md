# Tier 2 results: HG002 chr20 on HPRC v2.1 MC CHM13

Real reads, real benchmark, run on a 32 GB laptop.

| | |
|---|---|
| graph | `hprc-v2.1-mc-chm13-eval.HG002.gbz`, 100,179,277 nodes, **4 haplotypes** (CHM13, GRCh38, 2 recombinants). HG002 itself is **absent** — no circularity |
| chromosome | chr20 component, 2,382,533 nodes, IDs 114,818,865–121,250,404 |
| reads | 596,017,764 alignments genome-wide (~28.6×), 13,279,246 on chr20; 151 bp paired Illumina |
| truth | GIAB HG002 **draft** benchmark, defrabb V0.019-20241113, CHM13v2.0 |
| regions | small variants 58.9 Mb (88.9% of chr20); SVs 59.4 Mb (89.6%) |
| engine | `aardvark compare`; SV runs use `--min-variant-gap 1000` + record-basepair |

**All read-likelihood arms below use `--mismap-min 0.01`**, the current default. That floor caps how much one read can veto an allele; it was raised from 1e-8 after measurement, and the before/after comparison is in the calibration section at the end. `poisson` and `poisson-z` do not use the read-likelihood model, so the change cannot affect them.

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

The read path was optimised after the accuracy results below were first produced (vg `44fd008`); the calls are byte-identical, only the cost changed. `readlik-z` went **506 s to 97 s**, so the read-likelihood caller is now **1.35x** the Poisson caller at matched enumeration rather than 5.9x — and `readlik` is now *faster* than `poisson`.

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 106,587 | 156 s | 2.9 GB |
| `poisson-z` | haplotype (`-z`) | yes | 106,686 | 72 s | 2.9 GB |
| `readlik` | support (Flow) | yes | 104,725 | 118 s | 3.9 GB |
| `readlik-nomismap` | support (Flow) | yes | 106,682 | 118 s | 3.8 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 104,733 | 101 s | 3.8 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9556 | 0.9917 | 0.9733 | 71,684 | 3,333 | 588 | 0.9659 | 0.9909 | 0.9783 |
| `poisson-z` | 0.9562 | 0.9914 | 0.9735 | 71,734 | 3,283 | 605 | 0.9664 | 0.9906 | 0.9784 |
| `readlik` | 0.9580 | 0.9947 | 0.9760 | 71,866 | 3,151 | 373 | 0.9677 | 0.9930 | 0.9802 |
| `readlik-nomismap` | 0.9586 | 0.9942 | 0.9761 | 71,909 | 3,108 | 410 | 0.9681 | 0.9927 | 0.9802 |
| `readlik-z` | 0.9583 | 0.9947 | **0.9761** | 71,887 | 3,130 | 376 | 0.9680 | 0.9930 | 0.9804 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7281 | 0.8462 | 0.7827 | 7,069 | 2,640 | 1,462 | 0.7684 | 0.7613 | 0.7648 |
| `poisson-z` | 0.7295 | 0.8497 | 0.7850 | 7,083 | 2,626 | 1,426 | 0.7729 | 0.7694 | 0.7712 |
| `readlik` | 0.8027 | 0.8454 | 0.8235 | 7,793 | 1,916 | 1,543 | 0.8619 | 0.6388 | 0.7338 |
| `readlik-nomismap` | 0.8019 | 0.8420 | 0.8215 | 7,786 | 1,923 | 1,584 | 0.8620 | 0.5973 | 0.7057 |
| `readlik-z` | 0.8034 | 0.8448 | **0.8236** | 7,800 | 1,909 | 1,548 | 0.8665 | 0.6372 | 0.7344 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8251 | 0.8047 | 0.8148 | 8,222 | 1,743 | 2,148 | 0.8750 | 0.6764 | 0.7630 |
| `poisson-z` | 0.8259 | 0.8040 | 0.8148 | 8,230 | 1,735 | 2,159 | 0.8737 | 0.7015 | 0.7782 |
| `readlik` | 0.8551 | 0.8867 | **0.8706** | 8,521 | 1,444 | 1,176 | 0.8678 | 0.8653 | 0.8665 |
| `readlik-nomismap` | 0.8536 | 0.8848 | 0.8689 | 8,506 | 1,459 | 1,199 | 0.8663 | 0.8356 | 0.8507 |
| `readlik-z` | 0.8552 | 0.8865 | 0.8706 | 8,522 | 1,443 | 1,177 | 0.8674 | 0.8627 | 0.8651 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.8746 | — | 0 | 0 | 116 | — | 0.5789 | — |
| `poisson-z` | — | 0.8772 | — | 0 | 0 | 116 | — | 0.5131 | — |
| `readlik` | — | 0.8772 | — | 0 | 0 | 120 | — | 0.6485 | — |
| `readlik-nomismap` | — | 0.8749 | — | 0 | 0 | 124 | — | 0.6406 | — |
| `readlik-z` | — | 0.8789 | — | 0 | 0 | 118 | — | 0.6564 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7772 | 0.8261 | 0.8009 | 15,291 | 4,383 | 3,726 | 0.8227 | 0.7065 | 0.7602 |
| `poisson-z` | 0.7783 | 0.8275 | 0.8022 | 15,313 | 4,361 | 3,701 | 0.8243 | 0.7181 | 0.7675 |
| `readlik` | 0.8292 | 0.8669 | 0.8476 | 16,314 | 3,360 | 2,839 | 0.8649 | 0.7330 | 0.7935 |
| `readlik-nomismap` | 0.8281 | 0.8643 | 0.8458 | 16,292 | 3,382 | 2,907 | 0.8642 | 0.6961 | 0.7711 |
| `readlik-z` | 0.8296 | 0.8667 | **0.8477** | 16,322 | 3,352 | 2,843 | 0.8670 | 0.7312 | 0.7933 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9185 | 0.9531 | 0.9355 | 86,975 | 7,716 | 4,314 | 0.9040 | 0.8619 | 0.8825 |
| `poisson-z` | 0.9193 | 0.9532 | 0.9359 | 87,047 | 7,644 | 4,306 | 0.9046 | 0.8695 | 0.8867 |
| `readlik` | 0.9312 | 0.9650 | 0.9478 | 88,180 | 6,511 | 3,212 | 0.9264 | 0.8757 | 0.9003 |
| `readlik-nomismap` | 0.9315 | 0.9640 | 0.9474 | 88,201 | 6,490 | 3,317 | 0.9262 | 0.8532 | 0.8882 |
| `readlik-z` | 0.9315 | 0.9649 | **0.9479** | 88,209 | 6,482 | 3,219 | 0.9277 | 0.8745 | 0.9003 |

## Reading the insertion BASEPAIR numbers

The insertion `BASEPAIR` precision above understates the read-likelihood caller, and the reason is a property of the benchmark rather than of either caller.

**The `smvar` truth set contains no record >=50 bp** — that size class lives in the separate `stvar` benchmark. But the two confident regions overlap almost completely (58.9 Mb vs 59.4 Mb). So a >=50 bp insertion called inside the small-variant confident region has every one of its bases scored FP, however right the call is. It cannot be scored correct.

That is exactly where the gap lives. 246 `readlik-z` calls carry a >=200 bp insertion allele; they contribute **27,951 FP bases and zero TP bases**, which is the whole of the precision difference. The Poisson caller scores better there because it does not emit them — at the two largest sites it emits nothing at all.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called allele >=50 bp from REF, applied identically to each) gives the size-matched comparison:

| arm | class | BP recall | BP precision | **BP F1** |
|---|---|---|---|---|
| `sm50-poisson-z` | Insertion | 0.7637 | 0.8700 | **0.8134** |
| `sm50-poisson-z` | Deletion | 0.8628 | 0.8094 | **0.8353** |
| `sm50-poisson-z` | ALL | 0.8993 | 0.9385 | **0.9184** |
| `sm50-readlik-z` | Insertion | 0.8574 | 0.8639 | **0.8606** |
| `sm50-readlik-z` | Deletion | 0.8600 | 0.8874 | **0.8735** |
| `sm50-readlik-z` | ALL | 0.9233 | 0.9605 | **0.9416** |

The insertion BASEPAIR precision gap collapses from **0.139 to 0.008**, and insertion BASEPAIR F1 flips from a 0.047 loss into a 0.047 win. There is no insertion-sequence defect in the likelihood model; what the unrestricted number measures is that one caller emits large insertions and the other does not.

Whether those large calls are *correct* is a separate question, and the `stvar` comparison below is what answers it: they are a net win there (SV insertion recall 0.4976 vs 0.4263), but of the 246, only **35 are confirmed true**, **73 are confirmed false**, and **138 fall outside the SV confident region** and cannot be judged at all. See *Known bad output* for the worst of the unjudged ones.

## Structural variants (GIAB `stvar` benchmark)

Of 176,623 chr20 truth records only **2,052 are >=50 bp** — the rest is the local sequence context an SV-aware haplotype comparison needs to place the SV. The rows below are the SV-specific categories, not the whole benchmark.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

### SV insertion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.3877 | 321 | 507 | 713 | 336 | 377 | 0.4712 | 0.4254 |
| `poisson-z` | 0.4263 | 353 | 475 | 723 | 375 | 348 | 0.5187 | 0.4680 |
| `readlik` | 0.4553 | 377 | 451 | 668 | 315 | 353 | 0.4716 | 0.4633 |
| `readlik-nomismap` | 0.4589 | 380 | 448 | 725 | 334 | 391 | 0.4607 | 0.4598 |
| `readlik-z` | 0.4976 | 412 | 416 | 683 | 347 | 336 | 0.5081 | 0.5028 |

### SV deletion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.5434 | 463 | 389 | 713 | 336 | 377 | 0.4712 | 0.5048 |
| `poisson-z` | 0.5763 | 491 | 361 | 723 | 375 | 348 | 0.5187 | 0.5460 |
| `readlik` | 0.5106 | 435 | 417 | 668 | 315 | 353 | 0.4716 | 0.4903 |
| `readlik-nomismap` | 0.5047 | 430 | 422 | 725 | 334 | 391 | 0.4607 | 0.4817 |
| `readlik-z` | 0.5516 | 470 | 382 | 683 | 347 | 336 | 0.5081 | 0.5290 |

### SV (joint)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 784 | 896 | 713 | 336 | 377 | 0.4712 | 0.4689 |
| `poisson-z` | 0.5024 | 844 | 836 | 723 | 375 | 348 | 0.5187 | 0.5104 |
| `readlik` | 0.4833 | 812 | 868 | 668 | 315 | 353 | 0.4716 | 0.4774 |
| `readlik-nomismap` | 0.4821 | 810 | 870 | 725 | 334 | 391 | 0.4607 | 0.4712 |
| `readlik-z` | 0.5250 | 882 | 798 | 683 | 347 | 336 | 0.5081 | 0.5164 |

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the mismapping floor

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. Raising the floor caps that veto.

The *upper* clamp (`--mismap-max`) is inert here: it binds only where `e_r` is already large, i.e. the 6.3% of reads at MAPQ ≤ 9, while 90% are MAPQ 60.

| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| floor 1e-8 (old default) | 0.9370 | 0.9759 | 0.7783 | 0.8231 | 0.8686 |
| **floor 0.01 (current default)** | 0.9479 | 0.9761 | 0.8236 | 0.8706 | 0.9003 |
| floor 0.05 | 0.9495 | 0.9745 | 0.8346 | 0.8840 | 0.8954 |
| cap 0.2, floor 1e-8 | 0.9370 | 0.9759 | 0.7783 | 0.8233 | 0.8674 |
| cap 0.4, floor 1e-8 | 0.9370 | 0.9758 | 0.7785 | 0.8234 | 0.8710 |

Raising the floor to 0.01 changed **1,493 genotypes (1.41%)**, of which **94% were heterozygous → homozygous** (1/0→1/1: 614, 0/1→1/1: 606, 1/2→1/1: 184), and dropped 1,251 spurious non-reference calls. The failure it corrects is spurious heterozygosity: a few locally misaligned reads, each able to veto the homozygous hypothesis almost without bound, conjuring a second allele that is not there.

Calibrated on one chromosome of one sample. 0.05 is better on indel `GT` but costs SNVs and BASEPAIR, so the optimum lies between and is not worth over-fitting here.

## Known bad output

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

Filtering on depth is **not** that remedy, and the measurement says so plainly: dropping every call above DP 200 removes 195 records including all of the giants above, and moves insertion BASEPAIR precision by 0.0001 (0.6226 → 0.6227). Dropping above DP 58 removes 1,202 records and does help (+0.087), but costs SV insertion recall 0.4976 → 0.4167 — it is a blunt proxy for length that discards real SVs. The giants are bad output that no metric charges for; they should be fixed because they are wrong, not because they cost a score.

## Raw aardvark summary rows

<details><summary><code>poisson</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson | GT | ALL | ALL | ALL | 94691 | 86975 | 7716 | 91990 | 87676 | 4314 | 0.9185139031164524 | 0.9531035982171975 | 0.9354891215208296 | 1342 | 724 |
| poisson | GT | ALL | ALL | Snv | 75017 | 71684 | 3333 | 70558 | 69970 | 588 | 0.9555700707839556 | 0.991666430454378 | 0.97328368746245 | 199 | 214 |
| poisson | GT | ALL | ALL | Insertion | 9709 | 7069 | 2640 | 9507 | 8045 | 1462 | 0.7280873416417757 | 0.8462185757862628 | 0.7827208504668184 | 713 | 174 |
| poisson | GT | ALL | ALL | Deletion | 9965 | 8222 | 1743 | 11000 | 8852 | 2148 | 0.8250878073256397 | 0.8047272727272727 | 0.8147803626017884 | 430 | 320 |
| poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 925 | 809 | 116 |  | 0.8745945945945945 |  | 0 | 16 |
| poisson | GT | ALL | ALL | JointIndel | 19674 | 15291 | 4383 | 21432 | 17706 | 3726 | 0.7772186642268984 | 0.8261478163493841 | 0.8009366667641129 | 1143 | 510 |
| poisson | BASEPAIR | ALL | ALL | ALL | 390714 | 353211 | 37503 | 409806 | 353211 | 56595 | 0.9040141894070854 | 0.8618980688423302 | 0.8824539049617747 |  |  |
| poisson | BASEPAIR | ALL | ALL | Snv | 200440 | 193613 | 6827 | 190468 | 188743 | 1725 | 0.9659399321492717 | 0.9909433605645043 | 0.9782819098424921 |  |  |
| poisson | BASEPAIR | ALL | ALL | Insertion | 95512 | 73393 | 22119 | 96422 | 73404 | 23018 | 0.7684165340480777 | 0.7612785463898281 | 0.7648308862895322 |  |  |
| poisson | BASEPAIR | ALL | ALL | Deletion | 99134 | 86741 | 12393 | 126980 | 85892 | 41088 | 0.8749873908043658 | 0.6764214836982202 | 0.7629971425745341 |  |  |
| poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11504 | 6660 | 4844 |  | 0.5789290681502086 |  |  |  |
| poisson | BASEPAIR | ALL | ALL | JointIndel | 194646 | 160134 | 34512 | 234906 | 165956 | 68950 | 0.8226935051323941 | 0.7064783360152571 | 0.7601698159969886 |  |  |

</details>

<details><summary><code>poisson-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| poisson-z | GT | ALL | ALL | ALL | 94691 | 87047 | 7644 | 92035 | 87729 | 4306 | 0.919274271050047 | 0.9532134514043571 | 0.9359362843204727 | 1349 | 719 |
| poisson-z | GT | ALL | ALL | Snv | 75017 | 71734 | 3283 | 70584 | 69979 | 605 | 0.9562365863737553 | 0.9914286523858098 | 0.9735146793443107 | 200 | 219 |
| poisson-z | GT | ALL | ALL | Insertion | 9709 | 7083 | 2626 | 9488 | 8062 | 1426 | 0.7295293027088269 | 0.8497048903878583 | 0.785044572746258 | 717 | 172 |
| poisson-z | GT | ALL | ALL | Deletion | 9965 | 8230 | 1735 | 11018 | 8859 | 2159 | 0.8258906171600602 | 0.8040479215828644 | 0.8148229131320283 | 432 | 313 |
| poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 945 | 829 | 116 |  | 0.8772486772486773 |  | 0 | 15 |
| poisson-z | GT | ALL | ALL | JointIndel | 19674 | 15313 | 4361 | 21451 | 17750 | 3701 | 0.7783368913286571 | 0.827467250944012 | 0.8021504874990804 | 1149 | 500 |
| poisson-z | BASEPAIR | ALL | ALL | ALL | 390902 | 353624 | 37278 | 406720 | 353624 | 53096 | 0.9046359445589943 | 0.8694531864673486 | 0.8866957029770995 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Snv | 200440 | 193709 | 6731 | 190578 | 188793 | 1785 | 0.9664188784673717 | 0.9906337562572805 | 0.9783765103780948 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Insertion | 95512 | 73826 | 21686 | 95356 | 73369 | 21987 | 0.7729499958120446 | 0.7694219556189438 | 0.7711819406746623 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Deletion | 99134 | 86612 | 12522 | 122114 | 85666 | 36448 | 0.8736861218149172 | 0.7015248046906989 | 0.7781973520547425 |  |  |
| poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 14032 | 7200 | 6832 |  | 0.5131128848346637 |  |  |  |
| poisson-z | BASEPAIR | ALL | ALL | JointIndel | 194646 | 160438 | 34208 | 231502 | 166235 | 65267 | 0.8242553147765688 | 0.7180715501377958 | 0.7675082436222629 |  |  |

</details>

<details><summary><code>readlik</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik | GT | ALL | ALL | ALL | 94691 | 88180 | 6511 | 91808 | 88596 | 3212 | 0.9312395053384165 | 0.9650139421401185 | 0.9478259431176715 | 699 | 400 |
| readlik | GT | ALL | ALL | Snv | 75017 | 71866 | 3151 | 70475 | 70102 | 373 | 0.9579961875308264 | 0.9947073430294431 | 0.9760066773246536 | 90 | 130 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 7793 | 1916 | 9978 | 8435 | 1543 | 0.8026573282521372 | 0.8453597915413911 | 0.823455319172178 | 331 | 144 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 8521 | 1444 | 10378 | 9202 | 1176 | 0.8550928248871049 | 0.8866833686644826 | 0.8706016183924501 | 278 | 123 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 977 | 857 | 120 |  | 0.8771750255885363 |  | 0 | 3 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 16314 | 3360 | 21333 | 18494 | 2839 | 0.8292162244586764 | 0.8669197956218065 | 0.8476489518805115 | 609 | 270 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390680 | 361921 | 28759 | 413298 | 361921 | 51377 | 0.9263873246646872 | 0.8756901799669972 | 0.9003256308008428 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 193971 | 6469 | 190350 | 189020 | 1330 | 0.9677260027938536 | 0.9930128710270554 | 0.9802063795769204 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 82320 | 13192 | 128162 | 81868 | 46294 | 0.8618812295837173 | 0.6387852873706715 | 0.7337500273363614 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 86027 | 13107 | 98330 | 85081 | 13249 | 0.8677850182581153 | 0.865259839316587 | 0.8665205890978996 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11062 | 7174 | 3888 |  | 0.6485264870728621 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 168347 | 26299 | 237554 | 174123 | 63431 | 0.8648880531837284 | 0.7329828165385555 | 0.7934910050939721 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 88201 | 6490 | 92046 | 88729 | 3317 | 0.9314612793190483 | 0.9639636703387436 | 0.9474338023807497 | 773 | 353 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 71909 | 3108 | 70622 | 70212 | 410 | 0.958569390938054 | 0.9941944436577836 | 0.976056956245589 | 147 | 104 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 7786 | 1923 | 10027 | 8443 | 1584 | 0.8019363477186116 | 0.8420265283733919 | 0.8214926123528226 | 341 | 137 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 8506 | 1459 | 10406 | 9207 | 1199 | 0.8535875564475665 | 0.8847780126849895 | 0.8689029686928008 | 285 | 110 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 991 | 867 | 124 |  | 0.8748738647830474 |  | 0 | 2 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 16292 | 3382 | 21424 | 18517 | 2907 | 0.8280979973569178 | 0.8643110530246453 | 0.8458170936178713 | 626 | 249 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390680 | 361861 | 28819 | 424144 | 361861 | 62283 | 0.9262337462885226 | 0.8531560036214116 | 0.8881942603556106 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 194047 | 6393 | 190466 | 189069 | 1397 | 0.9681051686290162 | 0.9926653575966314 | 0.9802314452962927 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 82336 | 13176 | 137092 | 81884 | 55208 | 0.8620487478013233 | 0.5972923292387594 | 0.7056542333968251 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 85878 | 13256 | 101576 | 84879 | 16697 | 0.8662820021385196 | 0.8356206190438686 | 0.8506751137037882 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11272 | 7221 | 4051 |  | 0.6406139105748758 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 168214 | 26432 | 249940 | 173984 | 75956 | 0.8642047614644021 | 0.6961030647355365 | 0.7710985908204083 |  |  |

</details>

<details><summary><code>readlik-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z | GT | ALL | ALL | ALL | 94691 | 88209 | 6482 | 91794 | 88575 | 3219 | 0.9315457646450033 | 0.964932348519511 | 0.9479451792168735 | 700 | 408 |
| readlik-z | GT | ALL | ALL | Snv | 75017 | 71887 | 3130 | 70470 | 70094 | 376 | 0.9582761240785422 | 0.9946643961969632 | 0.9761312568926417 | 85 | 133 |
| readlik-z | GT | ALL | ALL | Insertion | 9709 | 7800 | 1909 | 9976 | 8428 | 1548 | 0.8033783087856627 | 0.8448275862068966 | 0.8235817618228264 | 332 | 144 |
| readlik-z | GT | ALL | ALL | Deletion | 9965 | 8522 | 1443 | 10374 | 9197 | 1177 | 0.8551931761164074 | 0.8865432812801234 | 0.8705860881110364 | 283 | 121 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 974 | 856 | 118 |  | 0.8788501026694046 |  | 0 | 10 |
| readlik-z | GT | ALL | ALL | JointIndel | 19674 | 16322 | 3352 | 21324 | 18481 | 2843 | 0.8296228524956796 | 0.8666760457700244 | 0.8477447624549234 | 615 | 275 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 390680 | 362428 | 28252 | 414422 | 362428 | 51994 | 0.9276850619432784 | 0.874538513881985 | 0.9003281571775998 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 200440 | 194029 | 6411 | 190338 | 189011 | 1327 | 0.9680153661943723 | 0.9930281919532621 | 0.9803622615939924 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 95512 | 82765 | 12747 | 128514 | 81892 | 46622 | 0.8665403300108887 | 0.6372224037848017 | 0.7343963241757471 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 99134 | 85988 | 13146 | 98534 | 85010 | 13524 | 0.8673916113543285 | 0.8627478839791342 | 0.8650635157403482 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11722 | 7694 | 4028 |  | 0.6563726326565432 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 194646 | 168753 | 25893 | 238770 | 174596 | 64174 | 0.8669738910637773 | 0.7312308916530552 | 0.7933377477757283 |  |  |

</details>

<details><summary><code>poisson</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-poisson | GT | ALL | ALL | ALL | 100207 | 88821 | 11386 | 95410 | 89562 | 5848 | 0.8863752033291088 | 0.938706634524683 | 0.9117906570388599 | 1462 | 884 |
| sv-poisson | GT | ALL | ALL | Snv | 78483 | 72822 | 5661 | 71941 | 70903 | 1038 | 0.9278697297503918 | 0.9855715099873508 | 0.955850591729682 | 270 | 279 |
| sv-poisson | GT | ALL | ALL | Insertion | 9857 | 7066 | 2791 | 10249 | 8423 | 1826 | 0.7168509688546211 | 0.8218362767099229 | 0.7657620259056992 | 701 | 208 |
| sv-poisson | GT | ALL | ALL | Deletion | 10187 | 8149 | 2038 | 11988 | 9326 | 2662 | 0.7999411014037499 | 0.7779446112779446 | 0.7887895354843473 | 430 | 373 |
| sv-poisson | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1232 | 910 | 322 |  | 0.7386363636363636 |  | 0 | 24 |
| sv-poisson | GT | ALL | ALL | SvInsertion | 828 | 321 | 507 | 0 | 0 | 0 | 0.38768115942028986 |  |  | 41 | 0 |
| sv-poisson | GT | ALL | ALL | SvDeletion | 852 | 463 | 389 | 0 | 0 | 0 | 0.5434272300469484 |  |  | 20 | 0 |
| sv-poisson | GT | ALL | ALL | JointIndel | 20044 | 15215 | 4829 | 23469 | 18659 | 4810 | 0.7590800239473159 | 0.7950487877625805 | 0.7766481752437369 | 1131 | 605 |
| sv-poisson | GT | ALL | ALL | JointStructuralVariant | 1680 | 784 | 896 | 0 | 0 | 0 | 0.4666666666666667 |  |  | 61 | 0 |
| sv-poisson | BASEPAIR | ALL | ALL | ALL | 969654 | 766195 | 203459 | 912096 | 766195 | 145901 | 0.7901736083180186 | 0.8400376714731783 | 0.8143430317523582 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Snv | 208362 | 197801 | 10561 | 193834 | 190888 | 2946 | 0.9493141743696067 | 0.9848014280260429 | 0.9667322402100234 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Insertion | 79376 | 60349 | 19027 | 236752 | 190208 | 46544 | 0.760292783712961 | 0.8034060958302358 | 0.7812550933453016 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Deletion | 82226 | 71194 | 11032 | 372090 | 272443 | 99647 | 0.8658331914479604 | 0.7321965115966568 | 0.7934271073872482 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 170138 | 128832 | 41306 |  | 0.7572206091525703 |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 179326 | 88496 | 0 | 0 | 0 | 0.6695715811247769 |  |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 268286 | 106776 | 0 | 0 | 0 | 0.7153110685699964 |  |  |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | JointIndel | 161602 | 131543 | 30059 | 778980 | 591483 | 187497 | 0.8139936386925904 | 0.7593044750828006 | 0.7856985362615818 |  |  |
| sv-poisson | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 447612 | 195272 | 0 | 0 | 0 | 0.696256245294641 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | ALL | 1067628 | 864169 | 203459 | 1214080 | 1068179 | 145901 | 0.8094289396681241 | 0.8798258763837639 | 0.8431605693193573 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Snv | 208430 | 197869 | 10561 | 194538 | 191592 | 2946 | 0.9493307105503047 | 0.9848564290781235 | 0.9667673147556781 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Insertion | 104216 | 85189 | 19027 | 266290 | 219746 | 46544 | 0.8174272664466109 | 0.825213113522851 | 0.8213017381630682 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Deletion | 108210 | 97178 | 11032 | 407276 | 307629 | 99647 | 0.8980500877922558 | 0.7553329928598788 | 0.8205319970767486 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 345976 | 304670 | 41306 |  | 0.8806102157375078 |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 181226 | 88496 | 0 | 0 | 0 | 0.6718992147470358 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 270274 | 106776 | 0 | 0 | 0 | 0.7168120938867524 |  |  |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | JointIndel | 212426 | 182367 | 30059 | 1019542 | 832045 | 187497 | 0.85849660587687 | 0.816096835637963 | 0.8367599515118957 |  |  |
| sv-poisson | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 451500 | 195272 | 0 | 0 | 0 | 0.6980821680592233 |  |  |  |  |

</details>

<details><summary><code>poisson-z</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-poisson-z | GT | ALL | ALL | ALL | 100207 | 89383 | 10824 | 95469 | 89701 | 5768 | 0.8919835939605018 | 0.9395824822717321 | 0.915164536223702 | 1497 | 896 |
| sv-poisson-z | GT | ALL | ALL | Snv | 78483 | 73280 | 5203 | 71987 | 70921 | 1066 | 0.9337053884280673 | 0.985191770736383 | 0.9587578579480464 | 297 | 290 |
| sv-poisson-z | GT | ALL | ALL | Insertion | 9857 | 7093 | 2764 | 10233 | 8476 | 1757 | 0.7195901389875216 | 0.8283005961106225 | 0.7701279264274401 | 702 | 209 |
| sv-poisson-z | GT | ALL | ALL | Deletion | 10187 | 8166 | 2021 | 12016 | 9354 | 2662 | 0.80160989496417 | 0.778462050599201 | 0.7898664163572933 | 433 | 370 |
| sv-poisson-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1233 | 950 | 283 |  | 0.7704785077047851 |  | 0 | 27 |
| sv-poisson-z | GT | ALL | ALL | SvInsertion | 828 | 353 | 475 | 0 | 0 | 0 | 0.42632850241545894 |  |  | 43 | 0 |
| sv-poisson-z | GT | ALL | ALL | SvDeletion | 852 | 491 | 361 | 0 | 0 | 0 | 0.5762910798122066 |  |  | 22 | 0 |
| sv-poisson-z | GT | ALL | ALL | JointIndel | 20044 | 15259 | 4785 | 23482 | 18780 | 4702 | 0.7612751945719417 | 0.7997615194617154 | 0.7800439296090023 | 1135 | 606 |
| sv-poisson-z | GT | ALL | ALL | JointStructuralVariant | 1680 | 844 | 836 | 0 | 0 | 0 | 0.5023809523809524 |  |  | 65 | 0 |
| sv-poisson-z | BASEPAIR | ALL | ALL | ALL | 969654 | 778806 | 190848 | 916550 | 778806 | 137744 | 0.8031792783817733 | 0.8497146909606678 | 0.8257919079802608 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Snv | 208362 | 198020 | 10342 | 194030 | 190967 | 3063 | 0.9503652297443872 | 0.9842137813740143 | 0.966993388305566 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Insertion | 79376 | 60700 | 18676 | 236464 | 192598 | 43866 | 0.764714775246926 | 0.8144918465390081 | 0.7888188167070147 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Deletion | 82226 | 71138 | 11088 | 373446 | 275117 | 98329 | 0.8651521416583562 | 0.736698210718551 | 0.7957747536320006 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 171994 | 133886 | 38108 |  | 0.7784341314231892 |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 189207 | 78615 | 0 | 0 | 0 | 0.7064654882720613 |  |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 268826 | 106236 | 0 | 0 | 0 | 0.7167508305293525 |  |  |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | JointIndel | 161602 | 131838 | 29764 | 781904 | 601601 | 180303 | 0.8158191111496145 | 0.769405195522724 | 0.7919326748690676 |  |  |
| sv-poisson-z | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 458033 | 184851 | 0 | 0 | 0 | 0.712466012531032 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | ALL | 1067628 | 876780 | 190848 | 1195494 | 1057750 | 137744 | 0.8212411064528095 | 0.8847806848047752 | 0.8518276522382817 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Snv | 208430 | 198088 | 10342 | 194864 | 191801 | 3063 | 0.9503814230197188 | 0.984281344938008 | 0.9670343800965592 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Insertion | 104216 | 85540 | 18676 | 265556 | 221690 | 43866 | 0.8207952713594842 | 0.8348145024025064 | 0.8277455314573465 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Deletion | 108210 | 97122 | 11088 | 408540 | 310211 | 98329 | 0.8975325755475464 | 0.7593161012385569 | 0.8226592392509066 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 326534 | 288426 | 38108 |  | 0.8832954608095941 |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 191107 | 78615 | 0 | 0 | 0 | 0.7085332305114155 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 270814 | 106236 | 0 | 0 | 0 | 0.7182442646863811 |  |  |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | JointIndel | 212426 | 182662 | 29764 | 1000630 | 820327 | 180303 | 0.8598853247719206 | 0.8198105193727951 | 0.8393698622994008 |  |  |
| sv-poisson-z | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 461921 | 184851 | 0 | 0 | 0 | 0.7141944920311949 |  |  |  |  |

</details>

<details><summary><code>readlik</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik | GT | ALL | ALL | ALL | 100207 | 90028 | 10179 | 95084 | 90460 | 4624 | 0.8984202700410151 | 0.9513693155525641 | 0.9241369764909731 | 841 | 648 |
| sv-readlik | GT | ALL | ALL | Snv | 78483 | 72987 | 5496 | 71762 | 71025 | 737 | 0.929972095867895 | 0.989729940637106 | 0.9589209259924636 | 168 | 229 |
| sv-readlik | GT | ALL | ALL | Insertion | 9857 | 7784 | 2073 | 10804 | 8826 | 1978 | 0.7896926042406411 | 0.8169196593854128 | 0.803075425080407 | 327 | 215 |
| sv-readlik | GT | ALL | ALL | Deletion | 10187 | 8445 | 1742 | 11217 | 9646 | 1571 | 0.8289977422204771 | 0.8599447267540341 | 0.8441877091839035 | 285 | 187 |
| sv-readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1301 | 963 | 338 |  | 0.7401998462720983 |  | 0 | 17 |
| sv-readlik | GT | ALL | ALL | SvInsertion | 828 | 377 | 451 | 0 | 0 | 0 | 0.4553140096618358 |  |  | 23 | 0 |
| sv-readlik | GT | ALL | ALL | SvDeletion | 852 | 435 | 417 | 0 | 0 | 0 | 0.5105633802816901 |  |  | 38 | 0 |
| sv-readlik | GT | ALL | ALL | JointIndel | 20044 | 16229 | 3815 | 23322 | 19435 | 3887 | 0.8096687287966474 | 0.8333333333333334 | 0.8213306072046722 | 612 | 419 |
| sv-readlik | GT | ALL | ALL | JointStructuralVariant | 1680 | 812 | 868 | 0 | 0 | 0 | 0.48333333333333334 |  |  | 61 | 0 |
| sv-readlik | BASEPAIR | ALL | ALL | ALL | 969654 | 698490 | 271164 | 868684 | 698490 | 170194 | 0.7203497329975435 | 0.8040783530029332 | 0.7599146620480021 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Snv | 208362 | 197831 | 10531 | 193780 | 191232 | 2548 | 0.9494581545579328 | 0.9868510682216947 | 0.9677935559406592 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Insertion | 79376 | 67387 | 11989 | 321062 | 215574 | 105488 | 0.8489593831888732 | 0.671440407148775 | 0.7498365068499772 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Deletion | 82226 | 70351 | 11875 | 237828 | 199704 | 38124 | 0.8555809597937392 | 0.8396992784701549 | 0.8475657279498962 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 171640 | 102371 | 69269 |  | 0.5964285714285714 |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 207680 | 60142 | 0 | 0 | 0 | 0.7754404044477302 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 164125 | 210937 | 0 | 0 | 0 | 0.43759431773946705 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointIndel | 161602 | 137738 | 23864 | 730530 | 517649 | 212881 | 0.8523285602900954 | 0.7085937606942905 | 0.7738433767828098 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 371805 | 271079 | 0 | 0 | 0 | 0.5783391716079417 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | ALL | 1067628 | 796464 | 271164 | 1187354 | 1017160 | 170194 | 0.7460126560936956 | 0.8566611137032427 | 0.7975173049558452 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Snv | 208430 | 197899 | 10531 | 194218 | 191670 | 2548 | 0.9494746437652929 | 0.9868807216632856 | 0.9678163816099847 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Insertion | 104216 | 92227 | 11989 | 352738 | 247250 | 105488 | 0.8849600829047363 | 0.7009451774404799 | 0.7822768709454622 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Deletion | 108210 | 96335 | 11875 | 275324 | 237200 | 38124 | 0.890259680251363 | 0.861530415074603 | 0.8756594684461012 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 365074 | 295805 | 69269 |  | 0.8102603855656662 |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 209580 | 60142 | 0 | 0 | 0 | 0.77702226737159 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 166113 | 210937 | 0 | 0 | 0 | 0.44055960747911416 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointIndel | 212426 | 188562 | 23864 | 993136 | 780255 | 212881 | 0.8876597026729308 | 0.785647685714746 | 0.8335441484892968 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 375693 | 271079 | 0 | 0 | 0 | 0.5808739401210937 |  |  |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik-nomismap | GT | ALL | ALL | ALL | 100207 | 90088 | 10119 | 95554 | 90707 | 4847 | 0.8990190306066442 | 0.9492747556355569 | 0.923463658151373 | 912 | 513 |
| sv-readlik-nomismap | GT | ALL | ALL | Snv | 78483 | 73065 | 5418 | 72006 | 71181 | 825 | 0.9309659416688965 | 0.9885426214482127 | 0.9588907599994149 | 217 | 160 |
| sv-readlik-nomismap | GT | ALL | ALL | Insertion | 9857 | 7781 | 2076 | 10901 | 8859 | 2042 | 0.7893882520036523 | 0.8126777359875241 | 0.8008637125588992 | 333 | 191 |
| sv-readlik-nomismap | GT | ALL | ALL | Deletion | 10187 | 8432 | 1755 | 11302 | 9683 | 1619 | 0.8277216059683911 | 0.8567510175190232 | 0.8419861721085401 | 287 | 152 |
| sv-readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1345 | 984 | 361 |  | 0.7315985130111524 |  | 0 | 10 |
| sv-readlik-nomismap | GT | ALL | ALL | SvInsertion | 828 | 380 | 448 | 0 | 0 | 0 | 0.45893719806763283 |  |  | 30 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | SvDeletion | 852 | 430 | 422 | 0 | 0 | 0 | 0.5046948356807511 |  |  | 45 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | JointIndel | 20044 | 16213 | 3831 | 23548 | 19526 | 4022 | 0.8088704849331471 | 0.8291999320536776 | 0.8189090581106386 | 620 | 353 |
| sv-readlik-nomismap | GT | ALL | ALL | JointStructuralVariant | 1680 | 810 | 870 | 0 | 0 | 0 | 0.48214285714285715 |  |  | 75 | 0 |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 969654 | 673323 | 296331 | 844444 | 673323 | 171121 | 0.6943951141334951 | 0.7973566038718968 | 0.7423226308611773 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 208362 | 197930 | 10432 | 193884 | 191247 | 2637 | 0.9499332891794089 | 0.9863990839883642 | 0.9678228172818298 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 79376 | 67407 | 11969 | 327590 | 214248 | 113342 | 0.8492113485184438 | 0.654012637748405 | 0.7389383872588566 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 82226 | 70063 | 12163 | 220152 | 174311 | 45841 | 0.8520784180186316 | 0.7917756822558959 | 0.8208209848423379 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 167020 | 103853 | 63167 |  | 0.6217997844569513 |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 206730 | 61092 | 0 | 0 | 0 | 0.7718932723973385 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 140448 | 234614 | 0 | 0 | 0 | 0.37446608827340544 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 161602 | 137470 | 24132 | 714762 | 492412 | 222350 | 0.8506701649732058 | 0.6889174298577708 | 0.7612967338495156 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 347178 | 295706 | 0 | 0 | 0 | 0.5400321053253775 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | ALL | 1067628 | 771297 | 296331 | 1193444 | 1022323 | 171121 | 0.7224398385954658 | 0.8566158110476906 | 0.7838272050911815 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Snv | 208430 | 197998 | 10432 | 195090 | 192453 | 2637 | 0.9499496233747541 | 0.9864831616177149 | 0.9678717641086927 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Insertion | 104216 | 92247 | 11969 | 360168 | 246826 | 113342 | 0.8851519920165809 | 0.6853079673929944 | 0.7725147130918243 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Deletion | 108210 | 96047 | 12163 | 257344 | 211503 | 45841 | 0.8875981887071435 | 0.821868782641134 | 0.8534698301329081 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 380842 | 317675 | 63167 |  | 0.8341385666496868 |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 208630 | 61092 | 0 | 0 | 0 | 0.7735001223481955 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 142436 | 234614 | 0 | 0 | 0 | 0.37776422225169076 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointIndel | 212426 | 188294 | 24132 | 998354 | 776004 | 222350 | 0.88639808686319 | 0.7772834084903751 | 0.8282625348188362 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 351066 | 295706 | 0 | 0 | 0 | 0.5427971526287471 |  |  |  |  |

</details>

<details><summary><code>readlik-z</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik-z | GT | ALL | ALL | ALL | 100207 | 90632 | 9575 | 95070 | 90574 | 4496 | 0.9044477930683485 | 0.9527085305564321 | 0.9279511013023877 | 877 | 658 |
| sv-readlik-z | GT | ALL | ALL | Snv | 78483 | 73468 | 5015 | 71755 | 71051 | 704 | 0.9361008116407374 | 0.9901888370148422 | 0.9623854591692552 | 192 | 233 |
| sv-readlik-z | GT | ALL | ALL | Insertion | 9857 | 7812 | 2045 | 10799 | 8853 | 1946 | 0.7925332251192047 | 0.8197981294564312 | 0.8059351492991361 | 329 | 213 |
| sv-readlik-z | GT | ALL | ALL | Deletion | 10187 | 8470 | 1717 | 11215 | 9674 | 1541 | 0.8314518503975655 | 0.8625947391885868 | 0.8467370336216856 | 292 | 185 |
| sv-readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1301 | 996 | 305 |  | 0.765564950038432 |  | 0 | 27 |
| sv-readlik-z | GT | ALL | ALL | SvInsertion | 828 | 412 | 416 | 0 | 0 | 0 | 0.4975845410628019 |  |  | 25 | 0 |
| sv-readlik-z | GT | ALL | ALL | SvDeletion | 852 | 470 | 382 | 0 | 0 | 0 | 0.5516431924882629 |  |  | 39 | 0 |
| sv-readlik-z | GT | ALL | ALL | JointIndel | 20044 | 16282 | 3762 | 23315 | 19523 | 3792 | 0.8123129115944921 | 0.8373579240832082 | 0.8246453033514645 | 621 | 425 |
| sv-readlik-z | GT | ALL | ALL | JointStructuralVariant | 1680 | 882 | 798 | 0 | 0 | 0 | 0.525 |  |  | 64 | 0 |
| sv-readlik-z | BASEPAIR | ALL | ALL | ALL | 969654 | 731162 | 238492 | 906176 | 731162 | 175014 | 0.7540442260847684 | 0.8068653330037432 | 0.77956104764291 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Snv | 208362 | 198690 | 9672 | 193762 | 191246 | 2516 | 0.9535807872836698 | 0.9870149977807826 | 0.9700098762332708 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Insertion | 79376 | 68250 | 11126 | 325092 | 218828 | 106264 | 0.8598316871598468 | 0.6731263765334121 | 0.755109225378431 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Deletion | 82226 | 70241 | 11985 | 239678 | 203132 | 36546 | 0.8542431834213022 | 0.8475204232345063 | 0.8508685243083177 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 205694 | 128917 | 76777 |  | 0.6267416648030569 |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 229193 | 38629 | 0 | 0 | 0 | 0.8557661431846525 |  |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 172796 | 202266 | 0 | 0 | 0 | 0.46071316209053437 |  |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | JointIndel | 161602 | 138491 | 23111 | 770464 | 550877 | 219587 | 0.8569881560871772 | 0.7149938219047224 | 0.7795779412567909 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 401989 | 240895 | 0 | 0 | 0 | 0.62529009899142 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | ALL | 1067628 | 829136 | 238492 | 1261934 | 1086920 | 175014 | 0.7766150756630587 | 0.8613128737319068 | 0.8167740990682276 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Snv | 208430 | 198758 | 9672 | 194116 | 191600 | 2516 | 0.9535959314877896 | 0.9870386779039337 | 0.970029147079162 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Insertion | 104216 | 93090 | 11126 | 356626 | 250362 | 106264 | 0.8932409610808321 | 0.7020295772041298 | 0.7861758356336824 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Deletion | 108210 | 96225 | 11985 | 273324 | 236778 | 36546 | 0.8892431383421125 | 0.8662905562628967 | 0.8776168014715069 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 437868 | 361091 | 76777 |  | 0.8246572026272758 |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 231093 | 38629 | 0 | 0 | 0 | 0.8567821683066268 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 174784 | 202266 | 0 | 0 | 0 | 0.4635565574857446 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | JointIndel | 212426 | 189315 | 23111 | 1067818 | 848231 | 219587 | 0.8912044664965683 | 0.7943591510912908 | 0.8399996488629476 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 405877 | 240895 | 0 | 0 | 0 | 0.627542627077239 |  |  |  |  |

</details>

