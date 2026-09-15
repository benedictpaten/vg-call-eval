# Re-fitting `--phase-min-q` under `--realign`

`--phase-min-q` decides which heterozygous sites are allowed to carry a read-backed phase link
(`read_phasing.cpp:103`). Its default, 9.5, was fitted against the greedy walk. The exact walk
lowers every per-read score by about half a phred, which moved the distribution the threshold sits
in, and the threshold did not move with it. See `anchors-regression-audit.md` for the diagnosis;
this is the fit.

**Result: 8.5, on a plateau spanning 8.0-9.0.** chr20 switch error 0.5794% -> 0.3826%, a 34%
relative reduction, confirmed on chr6 held out.

## The parameter's own header asked for this

`read_phasing.hpp:98-102`, written when the value was fitted:

> "Fitted on chr20; the distribution is tight (median 10.09, 25th percentile 9.95) so this is
> **sensitive and wants re-fitting whenever the per-read scores move**."

They moved. Measured over chr20 het sites: greedy median 10.06 (reproducing the recorded fitted
distribution), realign median **8.98**, interquartile range 8.82-9.01. The gate at 9.5 went from
sitting comfortably below the distribution to sitting above almost all of it.

## Method

- **Sweep on chr20 ONT, confirm on chr6 ONT.** chr6 sees exactly two arms -- the shipped value and
  the one chr20 chose -- so it cannot become a second fitting set.
- **Two objectives, not one.** Switch error against the T2T-Q100 phased truth
  (`scripts/tier2/phasing_benchmark.py`, whatshap 2.8), *and* small-variant F1 from aardvark. The
  gate reaches both: it decides which sites carry a phase link, and `--regenotype` is computed FROM
  the phasing chain, so the genotypes move with it. A value that wins on one and loses the other
  has not won.
- **One binary** (vg `2a9657e35`), pinned before the sweep, one flag varied.

Three controls, all reproducing:

| control | expected | got |
|---|---|---|
| `pq95` F1 against the standing `bd-c20` arm | 0.86659 indel | 0.86659 |
| `pq95` switch error against `p3b`, a separately-run 9.5 arm | -- | 0.5794% both |
| `pq85` switch error against `an-q85`, separately run | -- | 0.3826% both |

The `p3b`/`an-q85` pair were produced by a different calling path with a different sample name, so
their agreement says the harness is not contributing anything.

**Short reads are not exposed.** `read_phasing` defaults to `false` and is turned on only by
`--preset ont` (`call_main.cpp:555`, `:1291`), and `--phase-min-q` gates nothing else. No
short-read arm is needed, and there is no short-read regression to check for.

**One trap in the measurement.** The truth must be subset to the contig. HG002 is male, so the
whole-genome T2T truth mixes haploid chrX/chrY with diploid autosomes and whatshap refuses the file
outright -- "Inconsistent ploidy (2 and 1)", reported while reading the TRUTH, which reads like a
fault in the calls and is not. All 115,916 chr20 call records are uniformly diploid.

## chr20, the sweep

| `--phase-min-q` | reliable het sites | chain breaks | switch error | switches/pairs | indel F1 | SNV F1 | ALL F1 |
|---|---|---|---|---|---|---|---|
| 9.5 (**shipped**) | 4,521 (5.9%) | 487 | 0.5794% | 341/58,849 | 0.86659 | 0.98563 | 0.95863 |
| 9.0 | 21,210 (27.9%) | 624 | 0.4133% | 243/58,794 | **0.86691** | 0.98528 | 0.95843 |
| **8.5** | 60,203 (79.1%) | 1,084 | **0.3826%** | 225/58,803 | 0.86667 | 0.98541 | 0.95848 |
| 8.0 | 62,595 (82.2%) | 1,464 | 0.3901% | 229/58,708 | 0.86564 | 0.98524 | 0.95809 |
| 7.0 | 65,794 (86.4%) | 2,482 | 0.4673% | 274/58,638 | 0.86317 | 0.98547 | 0.95763 |
| 6.0 | 72,103 (94.7%) | 6,848 | 0.8068% | 473/58,627 | 0.85858 | 0.98506 | 0.95613 |

A U-curve, and the two ends fail by **different mechanisms**. Above the plateau, sites that should
carry links are excluded and the chain falls back on panel-inherited phase. Below it, marginal
sites are admitted that break chains rather than extend them -- 6.0 admits 11,900 more sites than
8.5 and takes the chain breaks from 1,084 to 6,848.

### Most of these differences are not real, and saying which is the point

Switches are rare per-site events, so Poisson is the right first approximation. Against the best arm:

| arm | switch error | vs 8.5 | sigma | verdict |
|---|---|---|---|---|
| 9.5 | 0.5794% | +116 | **4.88** | clearly worse |
| 9.0 | 0.4133% | +18 | 0.83 | indistinguishable |
| 8.5 | 0.3826% | -- | -- | best |
| 8.0 | 0.3901% | +4 | 0.19 | indistinguishable |
| 7.0 | 0.4673% | +49 | 2.19 | marginally worse |
| 6.0 | 0.8068% | +248 | **9.39** | clearly worse |

Adjacent pairs share sites, so this if anything understates the uncertainty. Indel F1 says the same
thing: 9.0, 8.5 and 8.0 span 0.00127, well inside the 0.003 this project requires before moving a
fitted default, while 6.0 is 0.008 below the best.

**So the sweep establishes two claims of different strength.** That 9.5 is wrong is a 4.9-sigma
result. That 8.5 specifically is right is not -- 8.0 through 9.0 are one plateau and the data
cannot separate them.

### Why 8.5 within the plateau

Not because it had the lowest number; that would be reading noise. Three reasons:

1. **It is the midpoint of the plateau (8.0-9.0)**, so it is the furthest from both failure modes.
2. **It ties the best F1.** 0.86667 against 9.0's 0.86691 is a difference of 0.00024.
3. **It is off the cliff.** The reliability distribution has its interquartile range inside
   8.82-9.01, so 9.0 sits ON the median. A gate on a mode that tight is exactly what got stranded
   this time: any future change that nudges the scores swings it wildly. 8.5 sits below the 25th
   percentile, which is the same relationship 9.5 had to the greedy distribution it was fitted
   against (median 10.09, gate 0.6 below).

Reason 3 is the one that matters for the next person. The failure being fixed here was not a bad
number, it was a number placed where a distribution could move across it.

## chr6, held out

Two arms only, the shipped value and the one chr20 chose. Nothing else was run here and no second
candidate was tried.

| chr6 ONT | reliable het sites | chain breaks | switch error | switches/pairs | indel F1 | SNV F1 | ALL F1 |
|---|---|---|---|---|---|---|---|
| 9.5 (**shipped**) | 11,405 (5.8%) | 1,128 | 0.4547% | 728/160,112 | 0.88351 | 0.98816 | 0.96518 |
| **8.5** | 162,622 (82.3%) | 1,464 | **0.2866%** | 459/160,138 | 0.88360 | 0.98813 | 0.96518 |

**-37.0% relative on switch error, at 7.8 sigma**, in a single 172 Mb block, with indel F1 moving
+0.00009 and ALL F1 identical to five decimal places.

The two contigs agree on everything that matters:

| | chr20 (fitted) | chr6 (held out) |
|---|---|---|
| reliable het sites, 9.5 -> 8.5 | 5.9% -> 79.1% | 5.8% -> 82.3% |
| switch error | 0.5794% -> 0.3826% | 0.4547% -> 0.2866% |
| relative change | **-34.0%** | **-37.0%** |
| significance | 4.9 sigma | 7.8 sigma |
| indel F1 | +0.00008 | +0.00009 |

The fitted contig understates the gain rather than overstating it, which is the right direction for
a hold-out to disagree in.

## Adopted

`--phase-min-q` now defaults to **8.5 when the exact walk is in use** and stays at 9.5 otherwise.

Keyed on `--realign` rather than on `--preset ont`, because the value thresholds a quantity the
WALK produces: `--read-phasing --realign` without the preset wants 8.5 too, and `--read-phasing`
alone still wants 9.5 -- greedy's distribution median is 10.06, and greedy at 9.5 scores 0.3545% on
chr20, the best number in the whole experiment. Keying it on nothing is how it came to be stranded.

The struct default in `read_phasing.hpp` stays 9.5 and is now documented as the greedy walk's value.

**What is NOT claimed.** That 8.5 is the optimum. 8.0 through 9.0 are one plateau on chr20 and the
measurement cannot separate them; 8.5 is its midpoint and the point furthest from both failure
modes. What is claimed, on two contigs, is that 9.5 is wrong.
