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

**Read the caveats before the numbers.** The benchmark is a *draft*: its own README reports known errors in highly homozygous regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute values are benchmark-relative; the arm-to-arm comparison is what this table is for.

## Cost

| arm | enumeration | pack? | variants | wall | peak RSS |
|---|---|---|---|---|---|
| `poisson` | support (Flow) | yes | 106,587 | 152 s | 2.9 GB |
| `poisson-z` | haplotype (`-z`) | yes | 106,686 | 74 s | 2.9 GB |
| `readlik` | support (Flow) | yes | 107,121 | 679 s | 4.2 GB |
| `readlik-nomismap` | support (Flow) | yes | 108,500 | 517 s | 4.3 GB |
| `readlik-z` | haplotype (`-z`) | **no** | 107,123 | 434 s | 3.3 GB |

## Small variants (GIAB `smvar` benchmark)

`GT` is the genotype-aware comparison — the one that matters for a genotyper. `BASEPAIR` weights by bases, so it penalises a call that finds the right locus with the wrong sequence. Bold marks the best GT F1 in each class.

### SNV

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9556 | 0.9917 | 0.9733 | 71,684 | 3,333 | 588 | 0.9659 | 0.9909 | 0.9783 |
| `poisson-z` | 0.9562 | 0.9914 | 0.9735 | 71,734 | 3,283 | 605 | 0.9664 | 0.9906 | 0.9784 |
| `readlik` | 0.9582 | 0.9943 | 0.9759 | 71,883 | 3,134 | 404 | 0.9681 | 0.9930 | 0.9804 |
| `readlik-nomismap` | 0.9572 | 0.9928 | 0.9747 | 71,808 | 3,209 | 510 | 0.9674 | 0.9920 | 0.9796 |
| `readlik-z` | 0.9583 | 0.9942 | **0.9759** | 71,892 | 3,125 | 413 | 0.9683 | 0.9930 | 0.9805 |

### Insertion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7281 | 0.8462 | 0.7827 | 7,069 | 2,640 | 1,462 | 0.7684 | 0.7613 | 0.7648 |
| `poisson-z` | 0.7295 | 0.8497 | **0.7850** | 7,083 | 2,626 | 1,426 | 0.7729 | 0.7694 | 0.7712 |
| `readlik` | 0.7519 | 0.8066 | 0.7783 | 7,300 | 2,409 | 2,033 | 0.8244 | 0.5757 | 0.6780 |
| `readlik-nomismap` | 0.7462 | 0.7887 | 0.7669 | 7,245 | 2,464 | 2,269 | 0.8168 | 0.5714 | 0.6724 |
| `readlik-z` | 0.7522 | 0.8062 | 0.7783 | 7,303 | 2,406 | 2,037 | 0.8263 | 0.5742 | 0.6775 |

### Deletion (<50 bp)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.8251 | 0.8047 | 0.8148 | 8,222 | 1,743 | 2,148 | 0.8750 | 0.6764 | 0.7630 |
| `poisson-z` | 0.8259 | 0.8040 | 0.8148 | 8,230 | 1,735 | 2,159 | 0.8737 | 0.7015 | 0.7782 |
| `readlik` | 0.8077 | 0.8385 | 0.8228 | 8,049 | 1,916 | 1,780 | 0.8422 | 0.7713 | 0.8052 |
| `readlik-nomismap` | 0.7882 | 0.8339 | 0.8104 | 7,854 | 2,111 | 1,839 | 0.8259 | 0.7251 | 0.7722 |
| `readlik-z` | 0.8081 | 0.8386 | **0.8231** | 8,053 | 1,912 | 1,778 | 0.8417 | 0.7711 | 0.8049 |

### Indel

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | — | 0.8746 | — | 0 | 0 | 116 | — | 0.5789 | — |
| `poisson-z` | — | 0.8772 | — | 0 | 0 | 116 | — | 0.5131 | — |
| `readlik` | — | 0.8574 | — | 0 | 0 | 145 | — | 0.6282 | — |
| `readlik-nomismap` | — | 0.8521 | — | 0 | 0 | 152 | — | 0.6130 | — |
| `readlik-z` | — | 0.8687 | — | 0 | 0 | 133 | — | 0.6406 | — |

### Indel (joint)

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.7772 | 0.8261 | 0.8009 | 15,291 | 4,383 | 3,726 | 0.8227 | 0.7065 | 0.7602 |
| `poisson-z` | 0.7783 | 0.8275 | **0.8022** | 15,313 | 4,361 | 3,701 | 0.8243 | 0.7181 | 0.7675 |
| `readlik` | 0.7802 | 0.8245 | 0.8017 | 15,349 | 4,325 | 3,958 | 0.8334 | 0.6602 | 0.7368 |
| `readlik-nomismap` | 0.7675 | 0.8135 | 0.7898 | 15,099 | 4,575 | 4,260 | 0.8214 | 0.6393 | 0.7190 |
| `readlik-z` | 0.7805 | 0.8248 | 0.8021 | 15,356 | 4,318 | 3,948 | 0.8341 | 0.6596 | 0.7367 |

### ALL

| arm | GT recall | GT precision | **GT F1** | TP | FN | FP | BP recall | BP precision | BP F1 |
|---|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.9185 | 0.9531 | 0.9355 | 86,975 | 7,716 | 4,314 | 0.9040 | 0.8619 | 0.8825 |
| `poisson-z` | 0.9193 | 0.9532 | 0.9359 | 87,047 | 7,644 | 4,306 | 0.9046 | 0.8695 | 0.8867 |
| `readlik` | 0.9212 | 0.9532 | 0.9369 | 87,232 | 7,459 | 4,362 | 0.9108 | 0.8303 | 0.8687 |
| `readlik-nomismap` | 0.9178 | 0.9490 | 0.9331 | 86,907 | 7,784 | 4,770 | 0.9044 | 0.8166 | 0.8583 |
| `readlik-z` | 0.9214 | 0.9532 | **0.9370** | 87,248 | 7,443 | 4,361 | 0.9112 | 0.8297 | 0.8686 |

## Structural variants (GIAB `stvar` benchmark)

Of 176,623 chr20 truth records only **2,052 are >=50 bp** — the rest is the local sequence context an SV-aware haplotype comparison needs to place the SV. The rows below are the SV-specific categories, not the whole benchmark.

**Precision here is recomputed, not read from aardvark.** Its summary leaves `query_total`/`query_tp`/`query_fp` at zero for the `Sv*` categories, so its own precision and F1 come out as 0/0. The per-variant `BD` decisions *are* in its annotated query VCF, so precision is counted from those over query variants of >=50 bp; recall is the published summary value; F1 is derived from the two.

### SV insertion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.3877 | 321 | 507 | 713 | 336 | 377 | 0.4712 | 0.4254 |
| `poisson-z` | 0.4263 | 353 | 475 | 723 | 375 | 348 | 0.5187 | 0.4680 |
| `readlik` | 0.4203 | 348 | 480 | 750 | 348 | 402 | 0.4640 | 0.4411 |
| `readlik-nomismap` | 0.4106 | 340 | 488 | 799 | 352 | 447 | 0.4406 | 0.4251 |
| `readlik-z` | 0.4614 | 382 | 446 | 755 | 384 | 371 | 0.5086 | 0.4838 |

### SV deletion (>=50 bp)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.5434 | 463 | 389 | 713 | 336 | 377 | 0.4712 | 0.5048 |
| `poisson-z` | 0.5763 | 491 | 361 | 723 | 375 | 348 | 0.5187 | 0.5460 |
| `readlik` | 0.4859 | 414 | 438 | 750 | 348 | 402 | 0.4640 | 0.4747 |
| `readlik-nomismap` | 0.4777 | 407 | 445 | 799 | 352 | 447 | 0.4406 | 0.4584 |
| `readlik-z` | 0.5176 | 441 | 411 | 755 | 384 | 371 | 0.5086 | 0.5131 |

### SV (joint)

| arm | recall | truth TP | truth FN | SV calls | TP | FP | precision\* | F1\* |
|---|---|---|---|---|---|---|---|---|
| `poisson` | 0.4667 | 784 | 896 | 713 | 336 | 377 | 0.4712 | 0.4689 |
| `poisson-z` | 0.5024 | 844 | 836 | 723 | 375 | 348 | 0.5187 | 0.5104 |
| `readlik` | 0.4536 | 762 | 918 | 750 | 348 | 402 | 0.4640 | 0.4587 |
| `readlik-nomismap` | 0.4446 | 747 | 933 | 799 | 352 | 447 | 0.4406 | 0.4426 |
| `readlik-z` | 0.4899 | 823 | 857 | 755 | 384 | 371 | 0.5086 | 0.4991 |

\* recomputed as described above. The per-variant counts are shared across the three SV rows because they are counted over all >=50 bp query variants, not split by insertion/deletion; only recall is category-specific.

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
| readlik | GT | ALL | ALL | ALL | 94691 | 87232 | 7459 | 93168 | 88806 | 4362 | 0.9212279942127551 | 0.9531813498196806 | 0.9369323150367264 | 1719 | 288 |
| readlik | GT | ALL | ALL | Snv | 75017 | 71883 | 3134 | 70615 | 70211 | 404 | 0.9582228028313582 | 0.9942788359413722 | 0.9759179035265743 | 171 | 104 |
| readlik | GT | ALL | ALL | Insertion | 9709 | 7300 | 2409 | 10512 | 8479 | 2033 | 0.7518796992481203 | 0.8066019786910198 | 0.7782801193442401 | 816 | 102 |
| readlik | GT | ALL | ALL | Deletion | 9965 | 8049 | 1916 | 11024 | 9244 | 1780 | 0.8077270446562971 | 0.8385341074020319 | 0.822842324340229 | 732 | 79 |
| readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1017 | 872 | 145 |  | 0.8574237954768928 |  | 0 | 3 |
| readlik | GT | ALL | ALL | JointIndel | 19674 | 15349 | 4325 | 22553 | 18595 | 3958 | 0.7801667174951713 | 0.8245022835099544 | 0.8017220245300664 | 1548 | 184 |
| readlik | BASEPAIR | ALL | ALL | ALL | 390680 | 355832 | 34848 | 428576 | 355832 | 72744 | 0.9108016791235795 | 0.8302658104980214 | 0.8686710869374165 |  |  |
| readlik | BASEPAIR | ALL | ALL | Snv | 200440 | 194045 | 6395 | 190426 | 189091 | 1335 | 0.9680951905807225 | 0.9929894027076135 | 0.9803842917830978 |  |  |
| readlik | BASEPAIR | ALL | ALL | Insertion | 95512 | 78739 | 16773 | 135942 | 78265 | 57677 | 0.8243885585057374 | 0.5757234703035118 | 0.6779740935231903 |  |  |
| readlik | BASEPAIR | ALL | ALL | Deletion | 99134 | 83488 | 15646 | 106676 | 82281 | 24395 | 0.8421732200859443 | 0.7713168847725824 | 0.8051892262612188 |  |  |
| readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11872 | 7458 | 4414 |  | 0.6282008086253369 |  |  |  |
| readlik | BASEPAIR | ALL | ALL | JointIndel | 194646 | 162227 | 32419 | 254490 | 168004 | 86486 | 0.8334463590312671 | 0.660159534755786 | 0.7367506554582858 |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-nomismap | GT | ALL | ALL | ALL | 94691 | 86907 | 7784 | 93576 | 88806 | 4770 | 0.9177957778458354 | 0.9490253911259298 | 0.9331493680496662 | 2036 | 281 |
| readlik-nomismap | GT | ALL | ALL | Snv | 75017 | 71808 | 3209 | 70740 | 70230 | 510 | 0.9572230294466587 | 0.9927905004240882 | 0.9746823966752685 | 271 | 103 |
| readlik-nomismap | GT | ALL | ALL | Insertion | 9709 | 7245 | 2464 | 10737 | 8468 | 2269 | 0.7462148521989906 | 0.7886746763527988 | 0.7668574787958532 | 861 | 101 |
| readlik-nomismap | GT | ALL | ALL | Deletion | 9965 | 7854 | 2111 | 11071 | 9232 | 1839 | 0.788158554942298 | 0.8338903441423539 | 0.8103797721393763 | 904 | 74 |
| readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1028 | 876 | 152 |  | 0.8521400778210116 |  | 0 | 3 |
| readlik-nomismap | GT | ALL | ALL | JointIndel | 19674 | 15099 | 4575 | 22836 | 18576 | 4260 | 0.7674595913388228 | 0.8134524435102469 | 0.7897869914432555 | 1765 | 178 |
| readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 390680 | 353322 | 37358 | 432658 | 353322 | 79336 | 0.9043769837206921 | 0.8166311497764979 | 0.8582672001049387 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 200440 | 193904 | 6536 | 190486 | 188969 | 1517 | 0.9673917381760128 | 0.9920361601377529 | 0.9795589683243778 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 95512 | 78018 | 17494 | 135998 | 77713 | 58285 | 0.8168397688248598 | 0.5714275209929558 | 0.6724421551548665 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 99134 | 81871 | 17263 | 111320 | 80719 | 30601 | 0.8258619646135533 | 0.725107797340999 | 0.772212282609532 |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 11948 | 7324 | 4624 |  | 0.6129896216940074 |  |  |  |
| readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 194646 | 159889 | 34757 | 259266 | 165756 | 93510 | 0.8214348098599509 | 0.6393279489019 | 0.719030149138684 |  |  |

</details>

<details><summary><code>readlik-z</code> — small variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| readlik-z | GT | ALL | ALL | ALL | 94691 | 87248 | 7443 | 93154 | 88793 | 4361 | 0.921396964864665 | 0.9531850484144535 | 0.9370214846211546 | 1750 | 293 |
| readlik-z | GT | ALL | ALL | Snv | 75017 | 71892 | 3125 | 70614 | 70201 | 413 | 0.9583427756375221 | 0.9941513014416404 | 0.9759186763347181 | 190 | 107 |
| readlik-z | GT | ALL | ALL | Insertion | 9709 | 7303 | 2406 | 10511 | 8474 | 2037 | 0.7521886909053456 | 0.8062030254019599 | 0.7782597814597823 | 822 | 101 |
| readlik-z | GT | ALL | ALL | Deletion | 9965 | 8053 | 1912 | 11016 | 9238 | 1778 | 0.8081284495735073 | 0.8385984023238925 | 0.8230815279461438 | 738 | 79 |
| readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1013 | 880 | 133 |  | 0.8687068114511353 |  | 0 | 6 |
| readlik-z | GT | ALL | ALL | JointIndel | 19674 | 15356 | 4318 | 22540 | 18592 | 3948 | 0.780522517027549 | 0.8248447204968944 | 0.8020717781582589 | 1560 | 186 |
| readlik-z | BASEPAIR | ALL | ALL | ALL | 390680 | 356003 | 34677 | 429054 | 356003 | 73051 | 0.9112393774956487 | 0.8297393801246463 | 0.8685817594487968 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Snv | 200440 | 194085 | 6355 | 190424 | 189087 | 1337 | 0.9682947515465975 | 0.9929788261983784 | 0.9804814552290002 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Insertion | 95512 | 78917 | 16595 | 136346 | 78286 | 58060 | 0.826252198676606 | 0.5741715928593432 | 0.6775242521370762 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Deletion | 99134 | 83443 | 15691 | 106574 | 82180 | 24394 | 0.8417192890431133 | 0.7711073995533619 | 0.8048675988773047 |  |  |
| readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 12046 | 7717 | 4329 |  | 0.6406275942221484 |  |  |  |
| readlik-z | BASEPAIR | ALL | ALL | JointIndel | 194646 | 162360 | 32286 | 254966 | 168183 | 86783 | 0.8341296507505934 | 0.6596291270208576 | 0.7366868352970594 |  |  |

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
| sv-readlik | GT | ALL | ALL | ALL | 100207 | 89047 | 11160 | 96644 | 90795 | 5849 | 0.8886305347929785 | 0.9394789122966765 | 0.9133475565043434 | 1908 | 455 |
| sv-readlik | GT | ALL | ALL | Snv | 78483 | 73004 | 5479 | 71955 | 71175 | 780 | 0.930188703286062 | 0.989159891598916 | 0.9587683648098505 | 263 | 179 |
| sv-readlik | GT | ALL | ALL | Insertion | 9857 | 7306 | 2551 | 11388 | 8897 | 2491 | 0.7411991478137364 | 0.7812609764664559 | 0.7607029711216176 | 802 | 145 |
| sv-readlik | GT | ALL | ALL | Deletion | 10187 | 7975 | 2212 | 11928 | 9737 | 2191 | 0.7828605084912142 | 0.8163145539906104 | 0.7992376092134456 | 734 | 117 |
| sv-readlik | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1373 | 986 | 387 |  | 0.718135469774217 |  | 0 | 14 |
| sv-readlik | GT | ALL | ALL | SvInsertion | 828 | 348 | 480 | 0 | 0 | 0 | 0.42028985507246375 |  |  | 47 | 0 |
| sv-readlik | GT | ALL | ALL | SvDeletion | 852 | 414 | 438 | 0 | 0 | 0 | 0.4859154929577465 |  |  | 62 | 0 |
| sv-readlik | GT | ALL | ALL | JointIndel | 20044 | 15281 | 4763 | 24689 | 19620 | 5069 | 0.7623727798842547 | 0.794685892502734 | 0.7781940446384581 | 1536 | 276 |
| sv-readlik | GT | ALL | ALL | JointStructuralVariant | 1680 | 762 | 918 | 0 | 0 | 0 | 0.45357142857142857 |  |  | 109 | 0 |
| sv-readlik | BASEPAIR | ALL | ALL | ALL | 969654 | 685833 | 283821 | 853480 | 685833 | 167647 | 0.7072966233316214 | 0.8035724328631016 | 0.7523670777902227 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Snv | 208362 | 197902 | 10460 | 193830 | 191284 | 2546 | 0.9497989076703046 | 0.9868647784140742 | 0.967977140575296 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Insertion | 79376 | 64336 | 15040 | 311866 | 211667 | 100199 | 0.8105220721628704 | 0.6787113696266986 | 0.7387834979709825 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Deletion | 82226 | 68075 | 14151 | 245410 | 194105 | 51305 | 0.8279011504876803 | 0.7909416894177091 | 0.808999513103892 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 169858 | 99700 | 70158 |  | 0.586960873199967 |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 202028 | 65794 | 0 | 0 | 0 | 0.7543368356595052 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 162445 | 212617 | 0 | 0 | 0 | 0.43311505831035935 |  |  |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointIndel | 161602 | 132411 | 29191 | 727134 | 505472 | 221662 | 0.8193648593458002 | 0.695156601121664 | 0.7521674741083577 |  |  |
| sv-readlik | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 364473 | 278411 | 0 | 0 | 0 | 0.566934314744184 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | ALL | 1067628 | 783807 | 283821 | 1223424 | 1055777 | 167647 | 0.7341574031404197 | 0.862969011560996 | 0.7933687436216117 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Snv | 208430 | 197970 | 10460 | 194582 | 192036 | 2546 | 0.9498152857074318 | 0.9869155420336927 | 0.9680100653110828 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Insertion | 104216 | 89176 | 15040 | 346714 | 246515 | 100199 | 0.855684347892838 | 0.7110038821622432 | 0.7766636419243227 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Deletion | 108210 | 94059 | 14151 | 283324 | 232019 | 51305 | 0.8692265040199612 | 0.818917564343296 | 0.843322397507004 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 398804 | 328646 | 70158 |  | 0.8240789961986339 |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 203928 | 65794 | 0 | 0 | 0 | 0.7560673582429316 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 164433 | 212617 | 0 | 0 | 0 | 0.43610396499138043 |  |  |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointIndel | 212426 | 183235 | 29191 | 1028842 | 807180 | 221662 | 0.8625827346934932 | 0.7845519525835842 | 0.8217190421596994 |  |  |
| sv-readlik | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 368361 | 278411 | 0 | 0 | 0 | 0.5695376423221784 |  |  |  |  |

</details>

<details><summary><code>readlik-nomismap</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik-nomismap | GT | ALL | ALL | ALL | 100207 | 88803 | 11404 | 97236 | 90840 | 6396 | 0.88619557515942 | 0.9342218931260027 | 0.9095752181339906 | 2237 | 389 |
| sv-readlik-nomismap | GT | ALL | ALL | Snv | 78483 | 73006 | 5477 | 72151 | 71206 | 945 | 0.9302141865117286 | 0.9869024684342559 | 0.9577202039036629 | 361 | 152 |
| sv-readlik-nomismap | GT | ALL | ALL | Insertion | 9857 | 7255 | 2602 | 11642 | 8893 | 2749 | 0.7360251597849244 | 0.7638721869094657 | 0.7496901700164892 | 841 | 128 |
| sv-readlik-nomismap | GT | ALL | ALL | Deletion | 10187 | 7795 | 2392 | 12034 | 9740 | 2294 | 0.7651909296161775 | 0.8093734419145754 | 0.7866623018062657 | 902 | 102 |
| sv-readlik-nomismap | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1409 | 1001 | 408 |  | 0.7104329311568488 |  | 0 | 7 |
| sv-readlik-nomismap | GT | ALL | ALL | SvInsertion | 828 | 340 | 488 | 0 | 0 | 0 | 0.4106280193236715 |  |  | 58 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | SvDeletion | 852 | 407 | 445 | 0 | 0 | 0 | 0.47769953051643194 |  |  | 75 | 0 |
| sv-readlik-nomismap | GT | ALL | ALL | JointIndel | 20044 | 15050 | 4994 | 25085 | 19634 | 5451 | 0.7508481341049691 | 0.7826988239984054 | 0.766442720856991 | 1743 | 237 |
| sv-readlik-nomismap | GT | ALL | ALL | JointStructuralVariant | 1680 | 747 | 933 | 0 | 0 | 0 | 0.4446428571428571 |  |  | 133 | 0 |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | ALL | 969654 | 699035 | 270619 | 852700 | 699035 | 153665 | 0.7209117891536568 | 0.8197900785739416 | 0.7671780565137178 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Snv | 208362 | 197811 | 10551 | 193902 | 191154 | 2748 | 0.9493621677657155 | 0.9858278924405112 | 0.9672514594370061 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Insertion | 79376 | 63781 | 15595 | 295178 | 204667 | 90511 | 0.8035300342672849 | 0.6933680694360691 | 0.7443954497843586 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Deletion | 82226 | 66886 | 15340 | 265012 | 200944 | 64068 | 0.8134410040619755 | 0.758244909664468 | 0.784873740046385 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 167108 | 114253 | 52855 |  | 0.6837075424276515 |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 197525 | 70297 | 0 | 0 | 0 | 0.7375234297406487 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 182340 | 192722 | 0 | 0 | 0 | 0.48615962160922727 |  |  |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointIndel | 161602 | 130667 | 30935 | 727298 | 519864 | 207434 | 0.8085729137015631 | 0.7147881611114014 | 0.7587936381795614 |  |  |
| sv-readlik-nomismap | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 379865 | 263019 | 0 | 0 | 0 | 0.590876425607108 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | ALL | 1067628 | 797009 | 270619 | 1225158 | 1071493 | 153665 | 0.7465231335259098 | 0.8745753608922278 | 0.8054917590336325 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Snv | 208430 | 197879 | 10551 | 194764 | 192016 | 2748 | 0.949378688288634 | 0.9858906163356678 | 0.9672902245659644 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Insertion | 104216 | 88621 | 15595 | 330550 | 240039 | 90511 | 0.8503588700391495 | 0.7261806080774467 | 0.7833792047716938 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Deletion | 108210 | 92870 | 15340 | 302368 | 238300 | 64068 | 0.8582386101099714 | 0.7881124986771086 | 0.8216820480939278 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 397476 | 344621 | 52855 |  | 0.8670234177661041 |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 199425 | 70297 | 0 | 0 | 0 | 0.7393723908320419 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 184328 | 192722 | 0 | 0 | 0 | 0.48886885028510807 |  |  |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointIndel | 212426 | 181491 | 30935 | 1030394 | 822960 | 207434 | 0.8543728168868218 | 0.7986847749501647 | 0.8255907892724897 |  |  |
| sv-readlik-nomismap | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 383753 | 263019 | 0 | 0 | 0 | 0.5933358277723835 |  |  |  |  |

</details>

<details><summary><code>readlik-z</code> — structural variants</summary>

| compare_label | comparison | region_label | filter | variant_type | truth_total | truth_tp | truth_fn | query_total | query_tp | query_fp | metric_recall | metric_precision | metric_f1 | truth_fn_gt | query_fp_gt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sv-readlik-z | GT | ALL | ALL | ALL | 100207 | 89626 | 10581 | 96621 | 90887 | 5734 | 0.8944085742512998 | 0.9406547230933234 | 0.9169489150178886 | 2001 | 501 |
| sv-readlik-z | GT | ALL | ALL | Snv | 78483 | 73475 | 5008 | 71952 | 71194 | 758 | 0.936190002930571 | 0.9894651990215699 | 0.9620906449219292 | 321 | 189 |
| sv-readlik-z | GT | ALL | ALL | Insertion | 9857 | 7331 | 2526 | 11384 | 8919 | 2465 | 0.7437354164553109 | 0.7834680252986648 | 0.7630848676004899 | 811 | 164 |
| sv-readlik-z | GT | ALL | ALL | Deletion | 10187 | 7997 | 2190 | 11919 | 9749 | 2170 | 0.7850201236870521 | 0.8179377464552395 | 0.8011409443139051 | 747 | 126 |
| sv-readlik-z | GT | ALL | ALL | Indel | 0 | 0 | 0 | 1366 | 1025 | 341 |  | 0.7503660322108345 |  | 0 | 22 |
| sv-readlik-z | GT | ALL | ALL | SvInsertion | 828 | 382 | 446 | 0 | 0 | 0 | 0.4613526570048309 |  |  | 56 | 0 |
| sv-readlik-z | GT | ALL | ALL | SvDeletion | 852 | 441 | 411 | 0 | 0 | 0 | 0.5176056338028169 |  |  | 66 | 0 |
| sv-readlik-z | GT | ALL | ALL | JointIndel | 20044 | 15328 | 4716 | 24669 | 19693 | 4976 | 0.7647176212332868 | 0.7982893510073371 | 0.7811429435699463 | 1558 | 312 |
| sv-readlik-z | GT | ALL | ALL | JointStructuralVariant | 1680 | 823 | 857 | 0 | 0 | 0 | 0.48988095238095236 |  |  | 122 | 0 |
| sv-readlik-z | BASEPAIR | ALL | ALL | ALL | 969654 | 714442 | 255212 | 887022 | 714442 | 172580 | 0.7368009619926283 | 0.8054388729930035 | 0.7695925406479105 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Snv | 208362 | 198557 | 9805 | 193824 | 191302 | 2522 | 0.9529424751154241 | 0.9869881954763084 | 0.9696665846517162 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Insertion | 79376 | 65198 | 14178 | 314550 | 215119 | 99431 | 0.8213817778673654 | 0.6838944523923065 | 0.7463592792968194 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Deletion | 82226 | 67856 | 14370 | 243496 | 195659 | 47837 | 0.8252377593461923 | 0.8035409205900713 | 0.8142448289864349 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | Indel | 0 | 0 | 0 | 201754 | 123747 | 78007 |  | 0.613355869028619 |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | SvInsertion | 267822 | 224249 | 43573 | 0 | 0 | 0 | 0.837306121229772 |  |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | SvDeletion | 375062 | 165805 | 209257 | 0 | 0 | 0 | 0.4420735771685748 |  |  |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | JointIndel | 161602 | 133054 | 28548 | 759800 | 534525 | 225275 | 0.8233437704978899 | 0.7035075019742038 | 0.7587229086316648 |  |  |
| sv-readlik-z | BASEPAIR | ALL | ALL | JointStructuralVariant | 642884 | 390054 | 252830 | 0 | 0 | 0 | 0.606725319031116 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | ALL | 1067628 | 812416 | 255212 | 1271190 | 1098610 | 172580 | 0.7609541900362299 | 0.8642374468018156 | 0.8093139189536587 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Snv | 208430 | 198625 | 9805 | 194488 | 191966 | 2522 | 0.9529578275680084 | 0.987032618979063 | 0.9696959714365555 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Insertion | 104216 | 90038 | 14178 | 346426 | 246995 | 99431 | 0.8639556306133416 | 0.7129805499587213 | 0.7812409509574124 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Deletion | 108210 | 93840 | 14370 | 277542 | 229705 | 47837 | 0.8672026614915442 | 0.8276405012574674 | 0.8469598382006276 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | Indel | 0 | 0 | 0 | 452734 | 374727 | 78007 |  | 0.8276979418378121 |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | SvInsertion | 269722 | 226149 | 43573 | 0 | 0 | 0 | 0.8384521841006666 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | SvDeletion | 377050 | 167793 | 209257 | 0 | 0 | 0 | 0.4450152499668479 |  |  |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | JointIndel | 212426 | 183878 | 28548 | 1076702 | 851427 | 225275 | 0.8656096711325356 | 0.7907731201390914 | 0.8265008113716849 |  |  |
| sv-readlik-z | RECORD_BP | ALL | ALL | JointStructuralVariant | 646772 | 393942 | 252830 | 0 | 0 | 0 | 0.609089447285906 |  |  |  |  |

</details>

