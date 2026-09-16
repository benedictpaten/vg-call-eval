# Splitting homozygous anchors: how far to push the threshold

> **The run-length tables in this file use a superseded definition of "phased".** They count an
> anchor as phased only when it carries two slots, which throws every *haploid* anchor onto the
> unphased side -- and a haploid anchor's single slot names a haplotype (`src/anchor.cpp:415`). The
> corrected measurement is [phased-run-lengths.md](phased-run-lengths.md), and it changes the
> magnitudes by two orders of magnitude, though not the conclusion. What remains valid here is the
> split/unsplit site accounting, the candidate-pool identity, the gap statistics, and the
> characterisation of what is left haploid.

The goal is long-range phasing, so the measure is how far one haplotype can be followed before the
chain breaks. `--anchors-hom-split` buys that by splitting a homozygous site's single slot when the
reads crossing it carry a confident cross-site phase; `--split-min-q` is the confidence the partition
must reach and `--split-min-side` is the per-strand read count it must reach on **both** sides.

Everything below is chr20 ONT at the shipped preset, vg `a59c58f51`. chr6 is the hold-out and is
never used to choose a value.

## What is being counted

- **phased anchor**: a snarl whose anchor record carries two slots, so a read pinned to it names a
  haplotype. Genuinely heterozygous, or homozygous and split.
- **haploid anchor**: one slot. The read is placed but names no haplotype, so the chain breaks.
- **run**: a maximal consecutive stretch of one kind. Reported in *anchors* over all snarls in node
  order -- which counts the nested and off-reference sites that have no reference position -- and in
  *reference bases* over the positioned subset, first to last POS of the run.
- **candidate pool**: 89,407 homozygous sites are eligible for splitting on chr20, and that number is
  identical in every arm. The threshold does not change what is offered, only where the split/unsplit
  boundary falls inside it, so `split + unsplit = 89,407` on every row and the two columns are one
  number reported twice.
- **held-out agreement**: the run's own leave-one-out self-check at heterozygous sites, where the
  cross-site phase can be compared against the partition the site's own alleles make.

## `--split-min-q`

Site counts. `no-split` is `--read-phasing` without `--anchors-hom-split`; 2.0 is the shipped default.

| arm | hom split | hom unsplit | phased anchors | haploid anchors | % haploid | held-out agree |
|---|---|---|---|---|---|---|
| no-split | 0 | (89,407) | 76,077 | 96,008 | 55.8% | -- |
| 3.0 | 78,795 | 10,612 | 154,453 | 17,630 | 10.2% | 95.203% |
| **2.0** | 80,630 | 8,777 | 156,276 | 15,807 | 9.2% | 95.054% |
| 1.5 | 82,062 | 7,345 | 157,691 | 14,392 | 8.4% | 94.967% |
| 1.0 | 84,142 | 5,265 | 159,751 | 12,331 | 7.2% | 94.869% |
| 0.5 | 85,678 | 3,729 | 161,255 | 10,827 | 6.3% | 94.710% |
| 0.25 | 86,838 | 2,569 | 162,404 | 9,677 | 5.6% | 94.587% |
| 0.1 | 87,137 | 2,270 | 162,695 | 9,386 | 5.5% | 94.473% |
| 0.0 | 87,509 | 1,898 | 163,035 | 9,046 | 5.3% | 94.367% |

The no-split `unsplit` figure is bracketed because the splitting code never runs on that arm and so
reports no counter; it is the candidate pool, which every other arm measures at 89,407.

Runs of **phased** anchors -- the stretch a haplotype survives, which is the thing being bought:

| arm | runs | mean | median | p90 | N50 | max | | runs (bp) | mean bp | p90 bp | N50 bp | max bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| no-split | 24,283 | 3.1 | 2 | 7 | 5 | 329 | | 11,416 | 3,005 | 7,470 | 18,716 | 176,375 |
| 3.0 | 2,003 | 77.1 | 3 | 261 | 430 | 2,214 | | 1,287 | 45,131 | 151,914 | 322,024 | 1,057,364 |
| **2.0** | 1,826 | 85.6 | 4 | 281 | 449 | 2,214 | | 1,188 | 49,612 | 161,544 | 339,123 | 1,167,187 |
| 1.5 | 1,700 | 92.8 | 6 | 308 | 461 | 2,214 | | 1,124 | 53,073 | 172,346 | 364,530 | 1,167,187 |
| 1.0 | 1,556 | 102.7 | 8 | 326 | 474 | 2,214 | | 1,048 | 57,714 | 194,620 | 378,652 | 1,167,187 |
| 0.5 | 1,479 | 109.0 | 10 | 341 | 476 | 2,214 | | 1,002 | 61,193 | 201,861 | 384,195 | 1,236,925 |
| 0.25 | 1,437 | 113.0 | 10 | 349 | 476 | 2,214 | | 968 | 64,044 | 215,242 | 396,240 | 1,236,925 |
| 0.1 | 1,426 | 114.1 | 11 | 354 | 480 | 2,214 | | 961 | 64,671 | 223,352 | 396,240 | 1,236,925 |
| 0.0 | 1,437 | 113.5 | 11 | 349 | 480 | 2,214 | | 962 | 64,680 | 223,352 | 396,240 | 1,236,925 |

Runs of **haploid** anchors -- the breaks:

| arm | runs | total | mean | median | p90 | N50 | max | | total bp | mean bp | p90 bp | max bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| no-split | 24,284 | 96,008 | 4.0 | 1 | 7 | 10 | 826 | | 16,309,549 | 1,429 | 2,786 | 2,011,359 |
| 3.0 | 2,003 | 17,630 | 8.8 | 2 | 20 | 34 | 826 | | 5,689,964 | 4,418 | 9,400 | 2,003,255 |
| **2.0** | 1,826 | 15,807 | 8.7 | 2 | 19 | 40 | 826 | | 5,301,319 | 4,459 | 9,029 | 2,003,255 |
| 1.5 | 1,700 | 14,392 | 8.5 | 1 | 18 | 44 | 826 | | 4,854,316 | 4,315 | 7,465 | 2,003,255 |
| 1.0 | 1,556 | 12,331 | 7.9 | 1 | 15 | 50 | 826 | | 4,134,107 | 3,941 | 5,000 | 2,003,255 |
| 0.5 | 1,479 | 10,827 | 7.3 | 1 | 12 | 51 | 826 | | 3,548,593 | 3,538 | 3,074 | 2,003,255 |
| 0.25 | 1,437 | 9,677 | 6.7 | 1 | 10 | 50 | 826 | | 3,111,561 | 3,211 | 1,939 | 2,003,255 |
| 0.1 | 1,426 | 9,386 | 6.6 | 1 | 9 | 51 | 826 | | 2,987,698 | 3,106 | 1,353 | 2,003,255 |
| 0.0 | 1,437 | 9,046 | 6.3 | 1 | 9 | 50 | 826 | | 2,929,299 | 3,042 | 1,273 | 2,003,255 |

**The haploid N50 in bases is not a usable statistic here and is left out of the table.** One run --
the centromere, 2,003,255 bp with no reads and therefore no anchors -- is a single 2 Mb block that
crosses half the total haploid mass as soon as the rest falls below about 3.9 Mb, so the N50 snaps
from 46 kb at q=2.0 to 2,003,255 at q<=0.5. That is the N50 correctly computed on a distribution one
artefact dominates, not a regression; **p90 is the statistic that tracks the real change**, and it
falls monotonically 9,029 -> 1,273.

Gaps between consecutive phased anchors, which is the read-bridging view of the same thing (33 kb is
the ONT read N50, beyond which no single read spans):

| arm | gaps >33kb | bases in them | gaps >10kb |
|---|---|---|---|
| no-split | 70 | 5.46M | 883 |
| **2.0** | 36 | 3.92M | 475 |
| 0.5 | 23 | 3.30M | 419 |
| 0.0 | 18 | 3.02M | 398 |

### Reading it

Hom-splitting itself is the large effect and it is not in question: 55.8% of anchors haploid becomes
9.2%, the phased-run N50 goes 5 anchors to 449, and in bases 18.7 kb to 339 kb -- **18x**. Everything
the threshold does after that is a second-order trim.

Within that, range gains monotonically as the threshold falls and saturates below 0.25: from 0.25 to
0.0 the phased-run N50 moves 396,240 -> 396,240 bp and p90 moves 215,242 -> 223,352, while held-out
agreement keeps falling 94.587% -> 94.367%. So the bottom of the range is all cost and no gain.

Agreement falls smoothly and slowly across the whole sweep -- 95.20% to 94.37%, 0.84 points for a
2.3x reduction in haploid anchors. Per unit of agreement given up, the range bought is largest around
0.5: relative to the 2.0 default it removes 4,980 haploid anchors and 13 of the 36 unbridgeable gaps
for 0.34 points of agreement, where going on to 0.25 costs a further 0.12 points for 1,150 anchors
and 4 gaps. **0.5 is the knee**, and it is what chr6 is being asked to confirm.

## `--split-min-side`

At q=0.5, varying the per-strand minimum:

| arm | hom split | hom unsplit | phased anchors | haploid anchors | phased N50 (bp) | haploid p90 (bp) | gaps >33kb | held-out agree |
|---|---|---|---|---|---|---|---|---|
| 1 | 86,916 | 2,491 | 162,208 | 9,869 | 394,596 | 2,316 | 21 | 94.7102% |
| **2** | 85,678 | 3,729 | 161,255 | 10,827 | 384,195 | 3,074 | 23 | 94.7102% |
| 3 | 84,668 | 4,739 | 160,348 | 11,735 | 394,596 | 3,769 | 28 | 94.7102% |
| 4 | 83,211 | 6,196 | 158,934 | 13,149 | 384,195 | 4,480 | 30 | 94.7102% |

The parameter moves range monotonically -- side 1 splits 3,705 more sites than side 4 and removes 9
of the unbridgeable gaps -- and the held-out agreement is **identical to four decimal places on every
arm**, 2913190/3075899.

That identity is not evidence that lowering it is free. It is evidence that **the check cannot see
this parameter**: agreement is evaluated per read at heterozygous sites, and `--split-min-side` is a
per-site gate on homozygous ones, so no arm of this sweep changes a single term in the statistic.
Lowering it would buy measured range against unmeasured risk.

**Declined, and it stays at 2** -- not because 1 is worse but because nothing here can say it is not.
Making this decidable needs an accuracy check that reads the split sites themselves, most directly
switch error measured against the T2T-Q100 truth on a phasing whose blocks are carried by split
homozygous anchors. That is worth building; it is not a reason to move the default first.

## What is left haploid, and why

Characterised at the shipped default. Of the 15,807 survivors:

| | count | share |
|---|---|---|
| no VCF line at all (nested, or off the reference) | 8,356 | 52.9% |
| homozygous, stayed collapsed for want of a confident partition | 3,912 | 24.7% |
| half-missing or nested haploid | 2,442 | 15.4% |
| heterozygous, one slot lost to `min-reads` | 1,097 | 6.9% |

So a **majority of what remains haploid is not a homozygous site that the threshold could rescue** --
it has no VCF line, and lowering `--split-min-q` cannot reach it. That is the ceiling the 0.0 arm runs
into: 9,046 haploid anchors at zero threshold against 8,356 with no line.

They are also not isolated. 9,247 (58.5%) have every one of their reads also placed at a phaseable
snarl, and only 18 of 15,807 have fewer than half; the median distance to the nearest phaseable site
is 2,081 bp and p90 is 14,300 bp. The chain routes around them.

## Status

chr20 chooses `--split-min-q 0.5`. chr6 confirmation is running; the default does not move until it
reports.
