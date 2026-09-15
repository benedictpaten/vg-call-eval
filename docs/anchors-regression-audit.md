# Auditing the anchors for regressions since the exact walk replaced greedy

Asked 2026-09-15: a collaborator suspects a bug in recent versions of the anchors code, introduced
when the exact anchor-constrained walk (`--realign`, on under `--preset ont`) replaced the greedy
read/allele scoring for long reads.

**There is one, it is real, and it is not in the anchor code.** It is a threshold that was fitted
against greedy and is now on the wrong side of the score distribution. Under `--preset ont`, the
fraction of heterozygous sites reliable enough to carry a read-backed phase link falls from
**77.4% to 5.9%**.

## Method

Two arms, same binary (vg `aca9fc6c6`), same graph, same reads, same everything, one flag apart:

| arm | command | file |
|---|---|---|
| realign (shipped) | `--preset ont` | `work/ont-preset/p3b.*` |
| greedy | `--preset ont --no-realign` | `work/ont-preset/an-greedy.*` |

Comparisons are restricted to snarls both arms called **identically** wherever the measure could
otherwise be confounded by the genotype having moved. That matters: the uncontrolled version of the
switch-rate measure below reads 0.43 percentage points and the controlled version reads 0.13.

## What is NOT wrong

Every structural invariant holds in both arms, at chr20 scale.

| check | realign | greedy |
|---|---|---|
| `score` within `[0, phred(--mismap-min)]` | 12,855,782 clean | 12,831,054 clean |
| `reliability` == dedup mean of its own scores | clean | clean |
| slot always in {0, 1} | clean | clean |
| half-missing genotype has exactly one slot | 2,653 clean | 2,733 clean |
| half-missing slot **is the strand its GT names** | clean | clean |
| `explained` within [0, 1] | 474,950 clean | 472,118 clean |
| **`offset` inside the real read** (from the GAF's own query lengths) | 12,855,782 clean | 12,831,054 clean |

That offset check had never been run: `check_anchors.py` can do it but wants FASTA/FASTQ and prints
"bounds check SKIPPED (no --reads)", which is what every anchors run to date has reported. The
lengths here come from column 2 of `chr20.ont.gaf.gz`.

Structure is stable across the walk change: at the **76,189** snarls where both arms called the same
genotype, only **122 (0.16%)** differ in their (pin node, slot) set, and those are slots crossing the
`min-reads=2` boundary as reads move between them.

`vg call` is deterministic across thread counts: `-t 2` and `-t 6` produce **byte-identical** anchor
files and VCFs. The exact walk's hoisted `ReadScratch` is a local in `compute`, and the calculator's
only mutable member is mutex-guarded, so there is no shared scratch to race on -- confirmed by
reading and then by running.

Two things that look like bugs and are not:

- **46 homozygous snarls carry two anchor slots** with `--anchors-hom-split` off. All 46 are
  symbolic-allele collapse: 17 of them are homozygous on every block, and **17/17 have two slots
  carrying different traversal indices with `AT` listing more traversals than `ALT`+1**. The snarl is
  heterozygous in traversal space and two distinct traversals spell the same ALT, so the VCF prints
  `1|1` while the anchors correctly keep them apart. Present in both walks (46 realign, 36 greedy),
  so it predates the replacement. The anchor header's "a homozygous site collapses to one slot" is
  true of *traversal* homozygosity and reads as false against the VCF's GT -- a documentation trap,
  and a plausible thing for a consumer to report as a bug.
- **19 sites where `reliability` appeared to disagree with its own scores**, realign only. Every one
  was off by *exactly* 0.050, which is the giveaway: `anchor.cpp:817` sets `std::fixed`, then `score`
  gets `setprecision(1)` and `reliability` `setprecision(2)`. A true value of 6.85 prints as `6.9` in
  the score column, so a mean recomputed from the file carries up to 0.05 of rounding. My checker's
  tolerance sat exactly on that boundary. Not a defect -- but note the `score` column is published to
  **one decimal place**, which is coarse for a value a consumer thresholds on.

## The regression

`--preset ont` turns on `--realign`. vg's own counters, from the two logs:

| chr20 ONT, read phasing | realign (shipped) | greedy |
|---|---|---|
| het sites | 76,135 | 76,514 |
| **reliable** | **4,521 (5.9%)** | **59,212 (77.4%)** |
| sites hung off the chain | 71,614 (94.1%) | 17,302 (22.6%) |
| chain breaks | 487 | 787 |

The read-phasing chain is now built almost entirely from *hung* sites rather than reliable ones.

### The code predicted this, in writing

`read_phasing.hpp:98-102`, describing the default, before any of this was measured:

> "A site below this is not allowed to carry a link. **Fitted on chr20; the distribution is tight
> (median 10.09, 25th percentile 9.95) so this is sensitive and wants re-fitting whenever the
> per-read scores move.** `--mismap-min` moves them directly..."

The exact walk moves them, and harder than `--mismap-min` does. Measured over het sites:

| | 25th pct | median | 75th pct | above 9.5 |
|---|---|---|---|---|
| fitted, as recorded in the header | 9.95 | **10.09** | -- | -- |
| greedy, measured here | 9.86 | **10.06** | 10.14 | 82.7% |
| realign (**shipped**) | 8.82 | **8.98** | 9.01 | **8.2%** |

Greedy reproduces the recorded fitted distribution to two decimal places, which is the check that
this is the same quantity the parameter was fitted against. Under realign the median falls *below*
the gate, and because the distribution is as tight as the comment says -- the 25th to 75th
percentile spans 8.82 to 9.01, a fifth of a phred -- essentially the whole distribution crosses at
once. That is why a 0.5 phred shift in the mean costs 92% of the reliable sites -- 77.4% of het
sites down to 5.9% -- rather than a few percent.

### Mechanism, established rather than guessed

`read_phasing.cpp:103` gates a phase link on `sites[...].reliability >= params.reliability`, whose
default is `--phase-min-q 9.5`. `graph_caller.cpp:6064-6084` builds that reliability as the mean over
the site's reads of

    win  = max(r0, r1) / (r0 + r1 + e)
    score = -10 * log10(1 - win)

which is the same quantity the anchor file publishes in its `score` column, and `r0`/`r1` are
`rel_at`, the walk's output.

**The exact walk finds the best correspondence against EVERY allele, the wrong one included.** So the
losing allele's `rel` rises, the winner's share `win` falls, `1 - win` rises, and the score drops.
Measured, paired on (read, snarl) over 3,083,685 observations at same-genotype sites:

| | |
|---|---|
| mean `score` | realign 10.677, greedy 11.180 (**-0.503**) |
| observations that fall under realign | 1,263,797 (40.98%) |
| observations that rise | **732 (0.02%)** |
| of those that change, the share that falls | **99.94%** |

It is essentially one-directional, which is what a systematic mechanism looks like and a bug usually
does not.

### Why a 0.5 phred shift costs 92% of the reliable sites

Because the shift moves a *mode*, not a tail. The second-largest spike in the score distribution:

| | greedy | realign |
|---|---|---|
| mode value | **10.2** | **9.0** |
| reads at it | 4,658,061 | 4,461,235 |

`--phase-min-q 9.5` sits between them. Under greedy the mode clears the gate; under realign it does
not. That also explains the otherwise odd shape of the loss -- at thresholds of 5.0, 11.0, 12.0 and
13.0 the two arms differ by 0.01-0.13%, and **only at 9.5** do they differ by 43.6%.

### The blast radius is exactly one parameter

Every other threshold on this quantity defaults to zero, i.e. off: `--min-confidence`
(`read_likelihood_caller.hpp:340`), `--anchors-min-gqn` and `--anchors-min-q`
(`anchor.hpp:236-237`). `--phase-min-q 9.5` is the **only non-zero default that thresholds the
per-read score**, which is why a shift this large has exactly one visible consequence rather than
many.

`--phase-break 10` sums `phase_link`, which is built from `q0`/`p` and so is also walk-dependent,
but it moved the other way and by little: chain breaks are 487 under realign against 787 under
greedy.

### This is the same class of defect as Part 1, and here it bites

`docs/next-realign-work.md` Part 1 asked whether the ONT preset's parameters were still right after
the walk replacement, checked the gap parameters and `--mismap-min`, and found nothing moved.
**`--phase-min-q` was never on that list**, and it is the one that moved -- because it thresholds a
quantity the walk changes directly, rather than feeding the walk as the gap parameters do.

## Confirmed by moving the gate, and sized

One more arm, `--preset ont --phase-min-q 8.5`, same binary and everything else
(`work/ont-preset/an-q85.*`):

| chr20 ONT | het sites | reliable | hung off the chain | chain breaks |
|---|---|---|---|---|
| realign, `--phase-min-q 9.5` (**shipped**) | 76,135 | **4,521 (5.9%)** | 71,614 (94.1%) | 487 |
| greedy, `--phase-min-q 9.5` | 76,514 | 59,212 (77.4%) | 17,302 | 787 |
| realign, `--phase-min-q 8.5` | 76,135 | **60,203 (79.1%)** | 15,932 | 1,084 |

Moving the gate below the mode restores the population outright -- to slightly *more* than greedy
had. That is the diagnosis confirmed: the code is doing what it should, and one fitted number is on
the wrong side of a distribution the walk moved.

**It changes a lot of output.** Between the shipped gate and 8.5, on the 115,099 snarls both emit:

| | |
|---|---|
| identical GT string | 79,438 |
| **same alleles, phase FLIPPED** | **35,081 (30.48%)** |
| genotype actually differs | 580 (0.50%) |

Both runs report a single phase block for the contig, so a block-orientation difference would flip
~100% of sites, not 30%. This is genuine re-phasing of nearly a third of chr20.

## What is NOT established

**Which phasing is better.** That needs `whatshap compare` against the phased truth via
`scripts/tier2/phasing_benchmark.py`, and whatshap is not installed on this machine. The
anchor-implied switch rate I measured (realign 6.88% against greedy 6.75%, controlled and paired,
99.94% of per-read changes downward) is a measure of per-read slot assignment and is **unmoved by
the gate** -- 7.06% at `--phase-min-q 8.5` -- so it is the wrong instrument for this question and
must not be quoted as though it answered it.

What can be said without whatshap: read phasing was measured to take chr20 from 3.79% to 0.518%
switch error, and that gain is produced by reliable sites carrying phase links. With 94.1% of sites
now hung rather than reliable, the mechanism that produced the gain is mostly switched off, so the
expectation is a regression back toward the panel's rate. **That is an inference from the
mechanism, not a measurement**, and the whatshap run is the thing that would settle it.

## Status

- The anchor *file* is correct: every invariant above passes, including the offset check that had
  never been run.
- The regression is a fitted threshold left on the wrong side of a distribution the walk moved.
- **No default has been changed.** `--phase-min-q` is a fitted parameter and this project fits on
  chr20 and confirms on chr6, which has not been done. 8.5 restores the population on chr20 (79.1% of het
  sites reliable, against greedy's 77.4%) but it was chosen to sit below the mode, not fitted, and
  it re-phases 30% of the contig -- so it is the start of a sweep, not a value to ship.
- **A second, unrelated defect** was found by code review and verified by hand: anchors for
  `reported_inline` and `no_reference` records are collected from the PRE-linkage genotype. See
  below.

## A separate defect: anchors off the render path use the pre-linkage genotype

`graph_caller.cpp:6670-6691` (the render path) rebuilds the genotype from
`linkage_collector->settled_traversals` before calling `collect_anchors_for_record`.
`hand_off_deferred_records` (`graph_caller.cpp:7303`, `:7315`) passes `pr.genotype` unchanged. The
barrier only writes back `pr.genotype` on the branch where the *ploidy* changed, so for a record
whose genotype linkage moved at unchanged ploidy, the anchor rows name the alleles the caller
abandoned.

Walk-independent -- it predates the exact walk -- and small on chr20: `inline_unrendered` is 97 and
`no_ref_unrendered` is 0, so at most 97 chains are exposed and only the ~5% linkage moves. But
`no_reference` is the off-reference population, and the code's own comment says those are "where an
assembler most needs help", so on a gRef run this is the wrong genotype on exactly the sites the
anchors exist for.


## Also found by code review, verified by hand

A six-way adversarial code review of the anchor pipeline ran alongside the measurements (39 agents:
one reader per surface, then three refutation lenses per finding). Most findings were refuted. These
survived, and I checked each against the source myself rather than taking the report.

### Fixed here

- **`--anchors-hom-split` was silently inert without `--read-phasing`.** The split is decided by
  each read's phase across the OTHER sites it crosses; with no phasing chain
  `read_strand_log_odds` returns 0.0 for every read, nothing clears `phase_min` on either side, and
  the flag ran, split nothing, and reported the reads as having no opinion -- blaming the data for a
  switch that was never armed. `call_main.cpp:1172` set the flag with no validation. Now refused,
  with the message naming the flag that fixes it. Not silently disarmed like `--regenotype`,
  because no preset turns this on, so asking for it is always deliberate.
- **Anchors off the render path used the pre-linkage genotype** (above). Both call sites now go
  through one `settled_genotype_for()`, so they cannot drift apart again.

### Real, not fixed here

- **The band can disconnect the DP.** The centre for read row `i` is projected back along the
  diagonal from the next shared visit, so between consecutive shared visits it moves by
  `(allele gap) - (read gap) + 1`. An allele run of more than `2 * band` = 128 visits that the read
  does not share therefore moves the centre further than the band is wide, and the out-of-band reset
  at `allele_likelihood.cpp:1013-1014` leaves every predecessor of the new band `NONE` -- the
  interior dies and the score falls back to `P[0]`, which is allele-independent. Reachable in
  principle; the author's own banded-against-unbanded measurement bounds the VCF impact at **one
  chr20 record in 115,255**, but that measurement was of records, and the anchors have never been
  compared against an unbanded run. A cheap guard would be to widen the band for any row whose
  centre jumps further than it.
- **`--anchors-min-gqn` is a no-op for a negative threshold.** `anchor.cpp:368` gates the filter on
  `params.min_gqn > 0.0`, a guard from when `gqn` was non-negative. v7 made `gqn` signed in [-1, 1],
  so a threshold like `-0.5` -- "drop the sites where linkage overrode the reads hardest", which is
  the population with a 37.8% false-positive rate -- is discarded without a word, while the file
  header still prints `min-gqn=-0.5` as though it had been applied.
- **Two counters are reported as something they do not count.** `single_pin` (`anchor.cpp:604`) is
  incremented inside the per-slot loop and before the `min_reads` test, so it counts up to twice per
  site and counts slots that emitted nothing, while the run report prints it as a count of sites.
  `hom_split_no_opinion` (`anchor.cpp:503`) fires once per READ but is reported as read placements,
  understating the withheld rows by up to 2x.

### Considered and rejected

- **The DP maximises the integer score alone, ignoring the `nats` term** (`better()`,
  `allele_likelihood.cpp:932`), while the value it emits is `log_base * score + nat_adjust`. The
  mechanics are real, but the framing as a walk regression does not survive: the greedy walk chooses
  its correspondence with the same nats-blind, purely structural test
  (`allele_likelihood.cpp:775-793`), so the convention "choose the correspondence, then charge nats
  onto it" predates the exact walk and is unchanged by it. `--insertion-nats 0.9` was also swept and
  fitted under exactly this behaviour. Worth documenting in the header; not a regression.

## The fix for the pre-linkage genotype, gated

`settled_genotype_for()` now serves both call sites. chr20, before against after:

| | |
|---|---|
| VCF | **byte-identical** -- the render path was only factored, and that is the proof |
| snarls in the anchor file | 172,340 before, 172,340 after, none gained or lost |
| snarls whose rows changed | **4** |
| of those, with a VCF line | **0** |

Four of the 97 `reported_inline` chains, all of them off the render path, which is the population
the fix targets and about the 5% of them the linkage layer moves. Reads moved between slots at those
four; no site gained or lost anchors.

A TAP fixture for this needs a `reported_inline` or `no_reference` record whose genotype linkage
moves at unchanged ploidy, and none of the existing fixtures produces one -- so the guard here is
the chr20 gate above, not a unit test. That is a gap, and it is the same gap as the outstanding
coherent-haploid fixture.
