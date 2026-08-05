# Tier 2: HG002 chr20 on HPRC v2.1 MC CHM13

Truth: GIAB HG002 draft benchmark, defrabb V0.019-20241113, CHM13v2.0 small
variants, restricted to the benchmark BED (58.9 Mb, 88.9% of chr20).

**The benchmark is a draft** — its README notes known errors in homozygous
regions, homopolymers and tandem repeats, and excludes VDJ and TSPY2. Absolute
numbers should be read with that in mind; the arm-to-arm comparison is the point.

**Recall is bounded by the graph, not only by the caller.** This graph carries 4
haplotypes (CHM13, GRCh38, 2 recombinants), so `-z` arms can only propose alleles
present in those walks. Compare down a column for the caller effect, across a row
for what the sampled graph costs.

## Cost

| arm | variants | wall | peak RSS |
|---|---|---|---|
| `poisson` | 106,587 | 152 s | 2.9 GB |
| `poisson-z` | 106,686 | 74 s | 2.9 GB |
| `readlik` | 107,121 | 679 s | 4.2 GB |
| `readlik-nomismap` | 108,500 | 517 s | 4.3 GB |
| `readlik-z` | 107,123 | 434 s | 3.3 GB |

## ALL

| arm | comparison | recall | precision | F1 | TP | FN | FP |
|---|---|---|---|---|---|---|---|
| `poisson` | GT | 0.9185 | 0.9531 | 0.9355 | 86,975 | 7,716 | 4,314 |
| `poisson` | BASEPAIR | 0.9040 | 0.8619 | 0.8825 | 353,211 | 37,503 | 56,595 |
| `poisson-z` | GT | 0.9193 | 0.9532 | 0.9359 | 87,047 | 7,644 | 4,306 |
| `poisson-z` | BASEPAIR | 0.9046 | 0.8695 | 0.8867 | 353,624 | 37,278 | 53,096 |
| `readlik` | GT | 0.9212 | 0.9532 | 0.9369 | 87,232 | 7,459 | 4,362 |
| `readlik` | BASEPAIR | 0.9108 | 0.8303 | 0.8687 | 355,832 | 34,848 | 72,744 |
| `readlik-nomismap` | GT | 0.9178 | 0.9490 | 0.9331 | 86,907 | 7,784 | 4,770 |
| `readlik-nomismap` | BASEPAIR | 0.9044 | 0.8166 | 0.8583 | 353,322 | 37,358 | 79,336 |
| `readlik-z` | GT | 0.9214 | 0.9532 | 0.9370 | 87,248 | 7,443 | 4,361 |
| `readlik-z` | BASEPAIR | 0.9112 | 0.8297 | 0.8686 | 356,003 | 34,677 | 73,051 |

## Snv

| arm | comparison | recall | precision | F1 | TP | FN | FP |
|---|---|---|---|---|---|---|---|
| `poisson` | GT | 0.9556 | 0.9917 | 0.9733 | 71,684 | 3,333 | 588 |
| `poisson` | BASEPAIR | 0.9659 | 0.9909 | 0.9783 | 193,613 | 6,827 | 1,725 |
| `poisson-z` | GT | 0.9562 | 0.9914 | 0.9735 | 71,734 | 3,283 | 605 |
| `poisson-z` | BASEPAIR | 0.9664 | 0.9906 | 0.9784 | 193,709 | 6,731 | 1,785 |
| `readlik` | GT | 0.9582 | 0.9943 | 0.9759 | 71,883 | 3,134 | 404 |
| `readlik` | BASEPAIR | 0.9681 | 0.9930 | 0.9804 | 194,045 | 6,395 | 1,335 |
| `readlik-nomismap` | GT | 0.9572 | 0.9928 | 0.9747 | 71,808 | 3,209 | 510 |
| `readlik-nomismap` | BASEPAIR | 0.9674 | 0.9920 | 0.9796 | 193,904 | 6,536 | 1,517 |
| `readlik-z` | GT | 0.9583 | 0.9942 | 0.9759 | 71,892 | 3,125 | 413 |
| `readlik-z` | BASEPAIR | 0.9683 | 0.9930 | 0.9805 | 194,085 | 6,355 | 1,337 |

## Indel

| arm | comparison | recall | precision | F1 | TP | FN | FP |
|---|---|---|---|---|---|---|---|
| `poisson` | GT | — | 0.8746 | — | 0 | 0 | 116 |
| `poisson` | BASEPAIR | — | 0.5789 | — | 0 | 0 | 4,844 |
| `poisson-z` | GT | — | 0.8772 | — | 0 | 0 | 116 |
| `poisson-z` | BASEPAIR | — | 0.5131 | — | 0 | 0 | 6,832 |
| `readlik` | GT | — | 0.8574 | — | 0 | 0 | 145 |
| `readlik` | BASEPAIR | — | 0.6282 | — | 0 | 0 | 4,414 |
| `readlik-nomismap` | GT | — | 0.8521 | — | 0 | 0 | 152 |
| `readlik-nomismap` | BASEPAIR | — | 0.6130 | — | 0 | 0 | 4,624 |
| `readlik-z` | GT | — | 0.8687 | — | 0 | 0 | 133 |
| `readlik-z` | BASEPAIR | — | 0.6406 | — | 0 | 0 | 4,329 |

## Insertion

| arm | comparison | recall | precision | F1 | TP | FN | FP |
|---|---|---|---|---|---|---|---|
| `poisson` | GT | 0.7281 | 0.8462 | 0.7827 | 7,069 | 2,640 | 1,462 |
| `poisson` | BASEPAIR | 0.7684 | 0.7613 | 0.7648 | 73,393 | 22,119 | 23,018 |
| `poisson-z` | GT | 0.7295 | 0.8497 | 0.7850 | 7,083 | 2,626 | 1,426 |
| `poisson-z` | BASEPAIR | 0.7729 | 0.7694 | 0.7712 | 73,826 | 21,686 | 21,987 |
| `readlik` | GT | 0.7519 | 0.8066 | 0.7783 | 7,300 | 2,409 | 2,033 |
| `readlik` | BASEPAIR | 0.8244 | 0.5757 | 0.6780 | 78,739 | 16,773 | 57,677 |
| `readlik-nomismap` | GT | 0.7462 | 0.7887 | 0.7669 | 7,245 | 2,464 | 2,269 |
| `readlik-nomismap` | BASEPAIR | 0.8168 | 0.5714 | 0.6724 | 78,018 | 17,494 | 58,285 |
| `readlik-z` | GT | 0.7522 | 0.8062 | 0.7783 | 7,303 | 2,406 | 2,037 |
| `readlik-z` | BASEPAIR | 0.8263 | 0.5742 | 0.6775 | 78,917 | 16,595 | 58,060 |

## Deletion

| arm | comparison | recall | precision | F1 | TP | FN | FP |
|---|---|---|---|---|---|---|---|
| `poisson` | GT | 0.8251 | 0.8047 | 0.8148 | 8,222 | 1,743 | 2,148 |
| `poisson` | BASEPAIR | 0.8750 | 0.6764 | 0.7630 | 86,741 | 12,393 | 41,088 |
| `poisson-z` | GT | 0.8259 | 0.8040 | 0.8148 | 8,230 | 1,735 | 2,159 |
| `poisson-z` | BASEPAIR | 0.8737 | 0.7015 | 0.7782 | 86,612 | 12,522 | 36,448 |
| `readlik` | GT | 0.8077 | 0.8385 | 0.8228 | 8,049 | 1,916 | 1,780 |
| `readlik` | BASEPAIR | 0.8422 | 0.7713 | 0.8052 | 83,488 | 15,646 | 24,395 |
| `readlik-nomismap` | GT | 0.7882 | 0.8339 | 0.8104 | 7,854 | 2,111 | 1,839 |
| `readlik-nomismap` | BASEPAIR | 0.8259 | 0.7251 | 0.7722 | 81,871 | 17,263 | 30,601 |
| `readlik-z` | GT | 0.8081 | 0.8386 | 0.8231 | 8,053 | 1,912 | 1,778 |
| `readlik-z` | BASEPAIR | 0.8417 | 0.7711 | 0.8049 | 83,443 | 15,691 | 24,394 |

## Raw aardvark summary rows

<details><summary><code>poisson</code></summary>

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

<details><summary><code>poisson-z</code></summary>

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

<details><summary><code>readlik</code></summary>

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

<details><summary><code>readlik-nomismap</code></summary>

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

<details><summary><code>readlik-z</code></summary>

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

