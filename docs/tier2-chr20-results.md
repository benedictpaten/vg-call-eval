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
| `readlik` | support (Flow) | yes | 105,930 | 115 s | 3.8 GB |
| `readlik-nomismap` | support (Flow) | yes | 106,682 | 115 s | 3.5 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 105,936 | 97 s | 3.5 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9556 | 0.9917 | 0.9733 | 71,684 | 3,333 | 588 | 0.9659 | 0.9909 | 0.9783 |
| `poisson-z` | 0.9562 | 0.9914 | 0.9735 | 71,734 | 3,283 | 605 | 0.9664 | 0.9906 | 0.9784 |
| `readlik` | 0.9588 | 0.9947 | 0.9764 | 71,927 | 3,090 | 371 | 0.9683 | 0.9930 | 0.9805 |
| `readlik-nomismap` | 0.9586 | 0.9942 | 0.9761 | 71,909 | 3,108 | 410 | 0.9681 | 0.9927 | 0.9802 |
| `readlik-z` | 0.9591 | 0.9947 | **0.9766** | 71,948 | 3,069 | 375 | 0.9686 | 0.9930 | 0.9807 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7281 | 0.8462 | 0.7827 | 7,069 | 2,640 | 1,462 | 0.7684 | 0.7613 | 0.7648 |
| `poisson-z` | 0.7295 | 0.8497 | 0.7850 | 7,083 | 2,626 | 1,426 | 0.7729 | 0.7694 | 0.7712 |
| `readlik` | 0.8028 | 0.8443 | 0.8230 | 7,794 | 1,915 | 1,557 | 0.8625 | 0.6242 | 0.7243 |
| `readlik-nomismap` | 0.8019 | 0.8420 | 0.8215 | 7,786 | 1,923 | 1,584 | 0.8620 | 0.5973 | 0.7057 |
| `readlik-z` | 0.8035 | 0.8436 | **0.8231** | 7,801 | 1,908 | 1,563 | 0.8669 | 0.6226 | 0.7247 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8251 | 0.8047 | 0.8148 | 8,222 | 1,743 | 2,148 | 0.8750 | 0.6764 | 0.7630 |
| `poisson-z` | 0.8259 | 0.8040 | 0.8148 | 8,230 | 1,735 | 2,159 | 0.8737 | 0.7015 | 0.7782 |
| `readlik` | 0.8550 | 0.8866 | 0.8705 | 8,520 | 1,445 | 1,178 | 0.8677 | 0.8549 | 0.8612 |
| `readlik-nomismap` | 0.8536 | 0.8848 | 0.8689 | 8,506 | 1,459 | 1,199 | 0.8663 | 0.8356 | 0.8507 |
| `readlik-z` | 0.8551 | 0.8866 | **0.8706** | 8,521 | 1,444 | 1,177 | 0.8676 | 0.8564 | 0.8620 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.8746 | — | 0 | 0 | 116 | — | 0.5789 | — |
| `poisson-z` | — | 0.8772 | — | 0 | 0 | 116 | — | 0.5131 | — |
| `readlik` | — | 0.8736 | — | 0 | 0 | 125 | — | 0.6340 | — |
| `readlik-nomismap` | — | 0.8749 | — | 0 | 0 | 124 | — | 0.6406 | — |
| `readlik-z` | — | 0.8753 | — | 0 | 0 | 123 | — | 0.6439 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7772 | 0.8261 | 0.8009 | 15,291 | 4,383 | 3,726 | 0.8227 | 0.7065 | 0.7602 |
| `poisson-z` | 0.7783 | 0.8275 | 0.8022 | 15,313 | 4,361 | 3,701 | 0.8243 | 0.7181 | 0.7675 |
| `readlik` | 0.8292 | 0.8662 | 0.8473 | 16,314 | 3,360 | 2,860 | 0.8652 | 0.7195 | 0.7856 |
| `readlik-nomismap` | 0.8281 | 0.8643 | 0.8458 | 16,292 | 3,382 | 2,907 | 0.8642 | 0.6961 | 0.7711 |
| `readlik-z` | 0.8296 | 0.8660 | **0.8474** | 16,322 | 3,352 | 2,863 | 0.8673 | 0.7192 | 0.7863 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9185 | 0.9531 | 0.9355 | 86,975 | 7,716 | 4,314 | 0.9040 | 0.8619 | 0.8825 |
| `poisson-z` | 0.9193 | 0.9532 | 0.9359 | 87,047 | 7,644 | 4,306 | 0.9046 | 0.8695 | 0.8867 |
| `readlik` | 0.9319 | 0.9649 | 0.9481 | 88,241 | 6,450 | 3,231 | 0.9268 | 0.8688 | 0.8969 |
| `readlik-nomismap` | 0.9315 | 0.9640 | 0.9474 | 88,201 | 6,490 | 3,317 | 0.9262 | 0.8532 | 0.8882 |
| `readlik-z` | 0.9322 | 0.9648 | **0.9482** | 88,270 | 6,421 | 3,238 | 0.9281 | 0.8684 | 0.8973 |

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
| `sm50-readlik-z` | Insertion | 0.8578 | 0.8624 | **0.8601** |
| `sm50-readlik-z` | Deletion | 0.8603 | 0.8865 | **0.8732** |
| `sm50-readlik-z` | ALL | 0.9238 | 0.9596 | **0.9413** |

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
| `readlik` | 0.4553 | 377 | 451 | 700 | 327 | 373 | 0.4671 | 0.4612 |
| `readlik-nomismap` | 0.4589 | 380 | 448 | 725 | 334 | 391 | 0.4607 | 0.4598 |
| `readlik-z` | 0.4976 | 412 | 416 | 710 | 357 | 353 | 0.5028 | 0.5002 |

### SV deletion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.5434 | 463 | 389 | 713 | 336 | 377 | 0.4712 | 0.5048 |
| `poisson-z` | 0.5763 | 491 | 361 | 723 | 375 | 348 | 0.5187 | 0.5460 |
| `readlik` | 0.5059 | 431 | 421 | 700 | 327 | 373 | 0.4671 | 0.4857 |
| `readlik-nomismap` | 0.5047 | 430 | 422 | 725 | 334 | 391 | 0.4607 | 0.4817 |
| `readlik-z` | 0.5446 | 464 | 388 | 710 | 357 | 353 | 0.5028 | 0.5229 |

### SV (joint)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 784 | 896 | 713 | 336 | 377 | 0.4712 | 0.4689 |
| `poisson-z` | 0.5024 | 844 | 836 | 723 | 375 | 348 | 0.5187 | 0.5104 |
| `readlik` | 0.4810 | 808 | 872 | 700 | 327 | 373 | 0.4671 | 0.4739 |
| `readlik-nomismap` | 0.4821 | 810 | 870 | 725 | 334 | 391 | 0.4607 | 0.4712 |
| `readlik-z` | 0.5214 | 876 | 804 | 710 | 357 | 353 | 0.5028 | 0.5120 |

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

## Calibration: the mismapping floor

MAPQ measures confidence that a read is in the right *place*, not that its path through a given site is right. A locally misaligned read is still MAPQ 60, so the mismapping term cannot discount it, yet it vetoes any allele it does not match by `ln(e_r)` — **−13.8 nats from one read** at the old floor of 1e-8. Raising the floor caps that veto.

The *upper* clamp (`--mismap-max`) is inert here: it binds only where `e_r` is already large, i.e. the 6.3% of reads at MAPQ ≤ 9, while 90% are MAPQ 60.

| `readlik-z` variant | ALL GT F1 | SNV GT F1 | Insertion GT F1 | Deletion GT F1 | ALL BP F1 |
|---|---|---|---|---|---|
| floor 1e-8 (old default) | 0.9370 | 0.9759 | 0.7783 | 0.8231 | 0.8686 |
| **floor 0.01 (current default)** | 0.9482 | 0.9766 | 0.8231 | 0.8706 | 0.8973 |
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
| readlik | GT | ALL | ALL | ALL | 94691 | 88241 | 6450 | 91934 | 88703 | 3231 | 0.931883705948823 | 0.9648552222246394 | 0.9480828877768925 | 721 | 368 |
| readlik | GT | ALL | ALL | Snv | 75017 | 71927 | 3090 | 70564 | 70193 | 371 | 0.9588093365503819 | 0.9947423615441302 | 0.9764453785799374 | 110 | 107 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 7794 | 1915 | 9997 | 8440 | 1557 | 0.8027603254712122 | 0.8442532759827949 | 0.8229841380906304 | 332 | 141 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 8520 | 1445 | 10384 | 9206 | 1178 | 0.8549924736578023 | 0.8865562403697997 | 0.8704883267233403 | 279 | 118 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 989 | 864 | 125 |  | 0.8736097067745198 |  | 0 | 2 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 16314 | 3360 | 21370 | 18510 | 2860 | 0.8292162244586764 | 0.8661675245671502 | 0.8472891932378063 | 611 | 261 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390680 | 362099 | 28581 | 416770 | 362099 | 54671 | 0.9268429405139756 | 0.8688221321112364 | 0.8968951637872313 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 194083 | 6357 | 190438 | 189112 | 1326 | 0.9682847734983038 | 0.9930371039393399 | 0.9805047486846131 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 82381 | 13131 | 131238 | 81918 | 49320 | 0.8625198927883407 | 0.6241942120422439 | 0.7242548155028488 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 86019 | 13115 | 99474 | 85039 | 14435 | 0.8677043194060564 | 0.8548867040633733 | 0.861247824483117 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11358 | 7201 | 4157 |  | 0.6340024652227505 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 168400 | 26246 | 242070 | 174158 | 67912 | 0.8651603423651141 | 0.7194530507704383 | 0.7856077077431753 |  |  |

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
| readlik-z | GT | ALL | ALL | ALL | 94691 | 88270 | 6421 | 91920 | 88682 | 3238 | 0.9321899652554096 | 0.9647737162750217 | 0.9482019986045973 | 722 | 376 |
| readlik-z | GT | ALL | ALL | Snv | 75017 | 71948 | 3069 | 70560 | 70185 | 375 | 0.9590892730980978 | 0.9946853741496599 | 0.9765630583838066 | 105 | 110 |
| readlik-z | GT | ALL | ALL | Insertion | 9709 | 7801 | 1908 | 9996 | 8433 | 1563 | 0.8034813060047379 | 0.8436374549819928 | 0.8230698844294271 | 333 | 141 |
| readlik-z | GT | ALL | ALL | Deletion | 9965 | 8521 | 1444 | 10378 | 9201 | 1177 | 0.8550928248871049 | 0.8865870109847754 | 0.8705551687709313 | 284 | 116 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 986 | 863 | 123 |  | 0.8752535496957403 |  | 0 | 9 |
| readlik-z | GT | ALL | ALL | JointIndel | 19674 | 16322 | 3352 | 21360 | 18497 | 2863 | 0.8296228524956796 | 0.8659644194756554 | 0.8474041811010998 | 617 | 266 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 390680 | 362595 | 28085 | 417522 | 362595 | 54927 | 0.9281125217569366 | 0.8684452555793467 | 0.8972880542240678 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 200440 | 194142 | 6298 | 190430 | 189103 | 1327 | 0.9685791259229695 | 0.9930315601533372 | 0.9806529372768994 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 95512 | 82802 | 12710 | 131586 | 81922 | 49664 | 0.8669277158891029 | 0.6225738300427097 | 0.7247075505566006 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 99134 | 86010 | 13124 | 99204 | 84963 | 14241 | 0.8676135331974902 | 0.8564473206725536 | 0.861994266870871 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 12074 | 7775 | 4299 |  | 0.6439456683783336 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 194646 | 168812 | 25834 | 242864 | 174660 | 68204 | 0.8672770054355086 | 0.7191679293761117 | 0.7863088021629087 |  |  |

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
| sv-readlik | GT | ALL | ALL | ALL | 100207 | 90113 | 10094 | 95325 | 90662 | 4663 | 0.8992685141756563 | 0.951083136637818 | 0.9244503538187512 | 879 | 578 |
| sv-readlik | GT | ALL | ALL | Snv | 78483 | 73071 | 5412 | 71900 | 71158 | 742 | 0.9310423913458966 | 0.9896801112656467 | 0.9594661760950861 | 196 | 189 |
| sv-readlik | GT | ALL | ALL | Insertion | 9857 | 7787 | 2070 | 10855 | 8846 | 2009 | 0.7899969564776301 | 0.8149239981575311 | 0.8022668984982884 | 327 | 208 |
| sv-readlik | GT | ALL | ALL | Deletion | 10187 | 8447 | 1740 | 11243 | 9686 | 1557 | 0.8291940708746441 | 0.8615138308280708 | 0.8450450367916328 | 287 | 166 |
| sv-readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1327 | 972 | 355 |  | 0.7324792765636775 |  | 0 | 15 |
| sv-readlik | GT | ALL | ALL | SvInsertion | 828 | 377 | 451 | 0 | 0 | 0 | 0.4553140096618358 |  |  | 27 | 0 |
| sv-readlik | GT | ALL | ALL | SvDeletion | 852 | 431 | 421 | 0 | 0 | 0 | 0.505868544600939 |  |  | 42 | 0 |
| sv-readlik | GT | ALL | ALL | JointIndel | 20044 | 16234 | 3810 | 23425 | 19504 | 3921 | 0.8099181800039912 | 0.8326147278548559 | 0.821109643285971 | 614 | 389 |
| sv-readlik | GT | ALL | ALL | JointStructuralVariant | 1680 | 808 | 872 | 0 | 0 | 0 | 0.48095238095238096 |  |  | 69 | 0 |
| sv-readlik | BASEPAIR | ALL | ALL | ALL | 969654 | 686964 | 282690 | 856030 | 686964 | 169066 | 0.7084630187675192 | 0.8024999123862482 | 0.7525552067060893 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Snv | 208362 | 197987 | 10375 | 193854 | 191329 | 2525 | 0.9502068515372285 | 0.9869747335623716 | 0.9682418636834994 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Insertion | 79376 | 67512 | 11864 | 323482 | 216389 | 107093 | 0.8505341664986897 | 0.6689367569138314 | 0.7488837834477954 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Deletion | 82226 | 70210 | 12016 | 226714 | 187970 | 38744 | 0.8538661737163428 | 0.8291062748661309 | 0.8413040904150735 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 170708 | 101433 | 69275 |  | 0.5941900789652506 |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 206490 | 61332 | 0 | 0 | 0 | 0.7709971548267133 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 153890 | 221172 | 0 | 0 | 0 | 0.4103054961579686 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointIndel | 161602 | 137722 | 23880 | 720904 | 505792 | 215112 | 0.8522295516144602 | 0.7016079810904087 | 0.7696184994230988 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 360380 | 282504 | 0 | 0 | 0 | 0.5605676918386521 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | ALL | 1067628 | 784938 | 282690 | 1176474 | 1007408 | 169066 | 0.735216760894244 | 0.8562943167464815 | 0.7911499238368707 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Snv | 208430 | 198055 | 10375 | 199574 | 197049 | 2525 | 0.9502230964832318 | 0.9873480513493742 | 0.9684299063901441 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Insertion | 104216 | 92352 | 11864 | 355896 | 248803 | 107093 | 0.8861595148537652 | 0.699089059725313 | 0.7815864741215764 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Deletion | 108210 | 96194 | 12016 | 264072 | 225328 | 38744 | 0.8889566583495055 | 0.8532824381229361 | 0.8707543142131323 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 356932 | 287657 | 69275 |  | 0.8059154124595161 |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 208390 | 61332 | 0 | 0 | 0 | 0.7726103172896538 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 155878 | 221172 | 0 | 0 | 0 | 0.4134146664898555 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointIndel | 212426 | 188546 | 23880 | 976900 | 761788 | 215112 | 0.8875843823260806 | 0.7798014126317945 | 0.830209250025766 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 364268 | 282504 | 0 | 0 | 0 | 0.5632092916823858 |  |  |  |  |

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
| sv-readlik-z | GT | ALL | ALL | ALL | 100207 | 90714 | 9493 | 95312 | 90757 | 4555 | 0.9052660991747083 | 0.9522095853617593 | 0.9281446471826293 | 916 | 585 |
| sv-readlik-z | GT | ALL | ALL | Snv | 78483 | 73548 | 4935 | 71891 | 71174 | 717 | 0.9371201406674057 | 0.990026567998776 | 0.9628471278241326 | 221 | 191 |
| sv-readlik-z | GT | ALL | ALL | Insertion | 9857 | 7819 | 2038 | 10852 | 8873 | 1979 | 0.7932433803388454 | 0.8176373018798379 | 0.8052556398416698 | 329 | 204 |
| sv-readlik-z | GT | ALL | ALL | Deletion | 10187 | 8471 | 1716 | 11236 | 9705 | 1531 | 0.8315500147246491 | 0.8637415450338198 | 0.8473401408233207 | 294 | 164 |
| sv-readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1333 | 1005 | 328 |  | 0.7539384846211553 |  | 0 | 26 |
| sv-readlik-z | GT | ALL | ALL | SvInsertion | 828 | 412 | 416 | 0 | 0 | 0 | 0.4975845410628019 |  |  | 29 | 0 |
| sv-readlik-z | GT | ALL | ALL | SvDeletion | 852 | 464 | 388 | 0 | 0 | 0 | 0.5446009389671361 |  |  | 43 | 0 |
| sv-readlik-z | GT | ALL | ALL | JointIndel | 20044 | 16290 | 3754 | 23421 | 19583 | 3838 | 0.8127120335262422 | 0.8361299688313907 | 0.8242547027423466 | 623 | 394 |
| sv-readlik-z | GT | ALL | ALL | JointStructuralVariant | 1680 | 876 | 804 | 0 | 0 | 0 | 0.5214285714285715 |  |  | 72 | 0 |
| sv-readlik-z | BASEPAIR | ALL | ALL | ALL | 969654 | 718545 | 251109 | 892598 | 718545 | 174053 | 0.741032368246818 | 0.8050040443738391 | 0.771694700824593 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Snv | 208362 | 198650 | 9712 | 193834 | 191340 | 2494 | 0.953388813699235 | 0.9871333202637308 | 0.9699676686987858 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Insertion | 79376 | 68318 | 11058 | 328898 | 219757 | 109141 | 0.8606883692803871 | 0.6681615576865776 | 0.7523025921085452 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Deletion | 82226 | 70055 | 12171 | 225134 | 189475 | 35659 | 0.8519811251915452 | 0.841609885668091 | 0.846763749649183 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 203666 | 129120 | 74546 |  | 0.6339791619612503 |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 230922 | 36900 | 0 | 0 | 0 | 0.8622219235163654 |  |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 157728 | 217334 | 0 | 0 | 0 | 0.42053847097279917 |  |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | JointIndel | 161602 | 138373 | 23229 | 757698 | 538352 | 219346 | 0.8562579671043675 | 0.7105099921076735 | 0.7766049054965654 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 388650 | 254234 | 0 | 0 | 0 | 0.604541410269971 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | ALL | 1067628 | 816519 | 251109 | 1256598 | 1082545 | 174053 | 0.7647972889433399 | 0.8614887179511665 | 0.8102685913192695 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Snv | 208430 | 198718 | 9712 | 194476 | 191982 | 2494 | 0.953404020534472 | 0.9871757954709064 | 0.9699960440830117 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Insertion | 104216 | 93158 | 11058 | 360236 | 251095 | 109141 | 0.8938934520611038 | 0.6970291697664864 | 0.7832811001632389 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Deletion | 108210 | 96039 | 12171 | 258526 | 222867 | 35659 | 0.8875242583864708 | 0.8620680318420585 | 0.8746109535489213 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 443360 | 368814 | 74546 |  | 0.831861241429087 |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 232822 | 36900 | 0 | 0 | 0 | 0.8631924722492047 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 159716 | 217334 | 0 | 0 | 0 | 0.42359368783980905 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | JointIndel | 212426 | 189197 | 23229 | 1062122 | 842776 | 219346 | 0.890648978938548 | 0.7934832345060172 | 0.8392631255145847 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 392538 | 254234 | 0 | 0 | 0 | 0.6069186668563264 |  |  |  |  |

</details>

