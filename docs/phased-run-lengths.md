# The length distribution of phased runs

The goal is long-range phasing, so the measure is how far a single haplotype can be followed before
the chain breaks. This file is chr20 ONT at the shipped preset, vg `a59c58f51`. chr6 is the hold-out
and never chooses a value.

## What counts as phased

An anchor is **phased** when a read pinned to it is thereby assigned to a specific haplotype. The
anchor file says which, and `src/anchor.cpp:415` states the rule outright:

```cpp
// Where the one slot lands. A homozygote's is slot 0 and names no haplotype; a nested haploid
// site's names the strand `nested_strand` gave, so slot 0 would be wrong for every `.|a` site
const int base_slot = (haploid && haploid_slot == 1) ? 1 : 0;
```

So a one-slot record is **not** automatically unphased. A haploid call -- a hemizygous site, a
half-missing genotype, a nested snarl reachable under only one of the parent's two settled strands --
emits one slot, and that slot **is** a haplotype. Only a collapsed homozygote's slot names nothing.
Per snarl:

| anchor record | phased? | why |
|---|---|---|
| >= 2 slots | **yes** | heterozygous, or a homozygote split by read phase |
| 1 slot, `slot == 1` | **yes** | only a haplotype index is ever 1 |
| 1 slot, `slot == 0`, GT haploid or half-missing | **yes** | haploid call sitting on haplotype 0 |
| 1 slot, `slot == 0`, GT heterozygous | **yes** | het whose other slot fell below `min-reads` |
| 1 slot, `slot == 0`, GT homozygous | no | the collapsed homozygote: reads could be either haplotype |
| 1 slot, `slot == 0`, no VCF line | ambiguous | nested or off-reference; see the bounds below |

The data corroborate the code. If slot 0 on a collapsed homozygote were a real haplotype index it
would be roughly 50/50 with slot 1; it is not, and the classes that *are* indexed are balanced:

| one-slot snarls, no-split arm | slot 0 | slot 1 | reading |
|---|---|---|---|
| GT homozygous | 35,827 | **3** | slot 0 is a placeholder, not an index |
| GT half-missing (`./X`) | 1,027 | 1,415 | balanced -- a genuine haplotype index |
| GT heterozygous | 734 | 363 | both occur -- a genuine index |
| no VCF line | 54,658 | 1,981 | mostly placeholders (nested collapsed homozygotes) |

**The ambiguous class has no reference position**, because it has no VCF line, so it cannot enter the
measurement in bases at all. The measurement in **bases is exact**; only the measurement in
**anchors** carries bounds, reported lo (every ambiguous snarl unphased) and hi (every one phased).

> An earlier version of this measurement defined phased as "two slots", which threw every haploid
> anchor onto the wrong side and understated the phased-run N50 by two orders of magnitude -- 339 kb
> instead of 2.1 Mb at the default. The tables below supersede it.

## Snarl classification

| arm | snarls | phased | unphased | ambiguous | % phased (lo-hi) |
|---|---|---|---|---|---|
| no-split | 172,085 | 81,600 | 35,827 | 54,658 | 47.4 - 79.2 |
| 3.0 | 172,083 | 160,310 | 4,380 | 7,393 | 93.2 - 97.5 |
| **2.0** | 172,083 | 162,135 | 3,703 | 6,245 | 94.2 - 97.8 |
| 1.5 | 172,083 | 163,561 | 3,038 | 5,484 | 95.0 - 98.2 |
| 1.0 | 172,082 | 165,629 | 2,085 | 4,368 | 96.3 - 98.8 |
| 0.5 | 172,082 | 167,159 | 1,406 | 3,517 | 97.1 - 99.2 |
| 0.25 | 172,081 | 168,316 | 822 | 2,943 | 97.8 - 99.5 |
| 0.1 | 172,081 | 168,612 | 619 | 2,850 | 98.0 - 99.6 |
| 0.0 | 172,081 | 168,957 | 540 | 2,584 | 98.2 - 99.7 |

Splitting is what converts the mass: 35,827 unphased collapsed homozygotes become 3,703.

## Phased runs, in reference bases

Exact -- every positioned snarl has a VCF line, so nothing is ambiguous.

| arm | runs | total | mean | median | p90 | **N50** | max |
|---|---|---|---|---|---|---|---|
| no-split | 11,000 | 35.4M | 3,220 | 273 | 8,078 | **19,237** | 176,375 |
| 3.0 | 577 | 59.4M | 103,013 | 38 | 177,678 | **1,324,573** | 3,644,510 |
| **2.0** | 475 | 60.3M | 126,944 | 50 | 188,482 | **2,118,890** | 4,560,785 |
| 1.5 | 402 | 61.0M | 151,817 | 112 | 227,035 | **2,402,560** | 5,797,010 |
| 1.0 | 318 | 62.0M | 194,977 | 337 | 294,600 | **3,208,292** | 5,797,010 |
| 0.5 | 250 | 62.9M | 251,501 | 907 | 442,147 | **3,575,641** | 8,587,683 |
| 0.25 | 212 | 63.6M | 299,804 | 1,636 | 683,844 | **3,575,641** | 8,969,836 |
| 0.1 | 199 | 63.7M | 320,242 | 1,848 | 918,832 | **3,575,641** | 8,987,586 |
| 0.0 | 199 | 63.8M | 320,614 | 2,105 | 918,832 | **3,575,641** | 8,989,146 |

### Read-walkable: the same runs, broken wherever no read can bridge

A run contiguous in anchor order can still cross a stretch that holds no anchors at all -- chr20's
centromere is 2 Mb of exactly that -- where the phase chain is carried by the panel rather than by a
read. Breaking each run wherever consecutive phased anchors are further apart than the ONT read N50
gives the runs reads could actually walk:

| arm | runs | mean | median | p90 | **N50** | max |
|---|---|---|---|---|---|---|
| no-split | 11,005 | 3,187 | 273 | 8,067 | **18,932** | 132,424 |
| 3.0 | 584 | 101,062 | 46 | 184,741 | **1,311,284** | 3,644,510 |
| **2.0** | 482 | 124,233 | 61 | 202,073 | **1,879,205** | 4,560,785 |
| 1.5 | 409 | 148,196 | 136 | 273,297 | **2,118,890** | 5,797,010 |
| 1.0 | 326 | 188,486 | 385 | 393,916 | **2,313,753** | 5,797,010 |
| 0.5 | 258 | 241,545 | 1,347 | 583,584 | **2,798,104** | 6,447,419 |
| 0.25 | 220 | 286,373 | 1,819 | 918,832 | **2,798,104** | 6,447,419 |
| 0.1 | 207 | 305,177 | 2,245 | 942,264 | **2,798,104** | 6,447,419 |
| 0.0 | 207 | 305,535 | 2,507 | 942,264 | **2,798,104** | 6,447,419 |

At a 10 kb cut the ordering is the same and 0.5 is the outright maximum: N50 449,404 (q=2.0) ->
473,828 (0.5) -> 462,636 (0.25 and below, i.e. it turns over).

## Phased runs, in anchors

Bounded, because the no-VCF-line class cannot be resolved from the files.

| arm | runs (lo) | N50 (lo) | max (lo) | | runs (hi) | N50 (hi) | max (hi) |
|---|---|---|---|---|---|---|---|
| no-split | 24,453 | 6 | 829 | | 15,250 | 33 | 1,034 |
| 3.0 | 1,586 | 1,199 | 5,249 | | 1,354 | 3,619 | 8,901 |
| **2.0** | 1,393 | 1,289 | 5,249 | | 1,170 | 4,120 | 9,857 |
| 1.5 | 1,242 | 1,289 | 5,249 | | 990 | 4,876 | 13,232 |
| 1.0 | 1,085 | 1,318 | 5,249 | | 725 | 6,070 | 13,232 |
| 0.5 | 970 | 1,517 | 5,249 | | 503 | 6,176 | 16,288 |
| 0.25 | 915 | 1,517 | 5,249 | | 345 | 6,176 | 16,918 |
| 0.1 | 894 | 1,517 | 5,249 | | 303 | 6,176 | 16,941 |
| 0.0 | 900 | 1,517 | 5,249 | | 265 | 6,176 | 16,944 |

The no-split arm is near its **lo** bound -- 54,658 of its 54,658+1,981 no-VCF-line one-slot snarls
sit at slot 0, which the 35,827:3 asymmetry says are overwhelmingly collapsed nested homozygotes.

## Unphased runs -- the breaks

| arm | runs | total bp | mean | median | p90 | max |
|---|---|---|---|---|---|---|
| no-split | 11,000 | 13,171,166 | 1,197 | 1 | 2,762 | 96,578 |
| **2.0** | 475 | 2,280,740 | 4,802 | 295 | 16,409 | 86,187 |
| 0.5 | 250 | 719,519 | 2,878 | 1 | 10,709 | 70,358 |
| 0.0 | 199 | 117,580 | 591 | 1 | 1,589 | 36,286 |

## What this says

**Hom-splitting is worth roughly 100x, not 18x.** Read-walkable phased-run N50 goes 18,932 bp to
1,879,205 bp at the shipped default -- from a run that a single ONT read spans to one that 57 reads
must be chained across. The earlier two-slot-only measurement put this at 18x because it counted
every haploid anchor as a break.

**The threshold then saturates exactly at 0.5.** Read-walkable N50 is 2,798,104 at 0.5, 0.25, 0.1 and
0.0 -- four arms, bit-identical -- while held-out agreement keeps falling, 94.710% -> 94.367%. Below
0.5 there is nothing left to buy and the price is still charged. At the 10 kb cut, 0.5 is not merely
the saturation point but the maximum, with 0.25 turning over.

That is a sharper case for 0.5 than the gap statistics made, and it is the same answer.

## `--split-min-side`

At q=0.5, the corrected metric separates the arms where held-out agreement could not:

| arm | phased | unphased | N50 bp | N50 bp (read-walkable) | mean bp | held-out agree |
|---|---|---|---|---|---|---|
| 1 | 168,305 | 1,009 | 3,575,641 | 2,798,104 | 293,041 | 94.7102% |
| **2** | 167,159 | 1,406 | 3,575,641 | 2,798,104 | 251,501 | 94.7102% |
| 3 | 166,185 | 1,807 | 3,256,172 | 2,711,536 | 214,393 | 94.7102% |
| 4 | 164,734 | 2,672 | 3,256,172 | 2,711,536 | 183,535 | 94.7102% |

Side 1 and side 2 are **tied on N50 under both cuts**, differing only in mean and in run count, and
the held-out agreement is identical to four decimals across all four arms -- because it is scored per
read at *heterozygous* sites while this gate applies to *homozygous* ones, so no arm of the sweep
changes a term in the statistic. Nothing here separates 1 from 2 on the metric that matters.
**Stays at 2.** Deciding it needs switch error measured against T2T-Q100 on a phasing whose blocks
are carried by split homozygous anchors.

## How much of this rests on the classification

`phase_haploid_slot` (`src/graph_caller.cpp:1326`) returns 0 for a haplotype-0 call **and** for three
refusals: a genuinely haploid locus with no second haplotype, a missing `render_phase`, and a phase
naming a different allele than the site settled on. Those refusals are indistinguishable from a real
slot 0 in the file, so the rule above over-counts by however many there are.

Bounded by re-running with **every** slot-0 haploid charged as unphased, which is the worst case:

| N50, reference bases | as classified | strict bound |
|---|---|---|
| chr20 q=2.0 | 2,118,890 | 592,712 |
| chr20 q=0.5 | 3,575,641 | 782,641 |
| chr20 gain | **+68.7%** | **+32.0%** |
| chr6 q=2.0 | 2,327,638 | 973,156 |
| chr6 q=0.5 | 3,541,082 | 1,226,834 |
| chr6 gain | **+52.1%** | **+26.1%** |

The absolute figures move by about 4x; **every comparison survives**, on both contigs. The truth lies
between the two columns, and the decision does not depend on where.

The empirical check says it is nearer the lenient end: refusals can only add to slot 0, yet the
half-missing class runs 1,027 slot-0 against 1,415 slot-1 -- a slot-0 *deficit*, not the excess a
large refused population would leave.

## chr6, held out

Two arms, the default and the value chr20 chose. Nothing else varied.

| chr6 | q=2.0 | q=0.5 | change |
|---|---|---|---|
| homozygous sites split | 208,019 | 220,738 | +12,719 |
| left collapsed | 15,350 | 2,631 | -82.9% |
| phased snarls | 407,998 | 420,712 | +12,714 |
| unphased snarls | 9,629 | 1,784 | -81.5% |
| phased runs | 613 | 131 | -78.6% |
| **run N50, read-walkable** | **2,130,682** | **2,883,125** | **+35.3%** |
| run N50, raw | 2,327,638 | 3,541,082 | +52.1% |
| mean run | 263,867 | 1,281,927 | 4.9x |
| longest run | 6,455,869 | 10,305,202 | +59.6% |
| unphased bases | 4,765,471 | 489,694 | -89.7% |
| held-out agreement | 96.1577% | 95.8719% | -0.286 pts |

The hold-out moves further than chr20 did and pays the same price for it, -0.286 points against
chr20's -0.344.

**The two chr6 VCFs are byte-identical.** `--split-min-q` moves anchors and nothing else, so none of
this is bought by re-genotyping.

## Adopted

`--split-min-q` defaults to **0.5**, was 2.0 (vg `src/anchor.hpp`). It is reachable only through
`--anchors-hom-split`, which is opt-in, so nothing changes for a run that does not ask for it.

`--split-min-side` stays at 2: nothing separates it from 1 on N50 under either cut, and the held-out
agreement cannot see the parameter at all.
