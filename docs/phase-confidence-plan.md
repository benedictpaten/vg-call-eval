# Plan: longer-range phasing — the split thresholds, and decoupling mapQ for phase confidence

Written 2026-09-15. Nothing here is started.

**The goal, stated as a number.** After `--anchors-hom-split`, chr20 has **36 gaps between
consecutive phaseable anchors that exceed the ONT read N50 of 33 kb**, spanning 3.92 Mb. A gap
longer than a read cannot be bridged by one read, so that is the quantity limiting phasing range.
Splitting already halved it (70 -> 36). The question is what halves it again.

Three of those gaps are over 100 kb and one is the 2.01 Mb centromere at chr20:27,047,383-29,059,127.
**Nothing here will touch those** — there are no confidently called variants to anchor on. The
addressable population is the 33 gaps between 33 kb and 100 kb.

## Where the evidence says the constraint is

Of the 15,807 anchors left unphaseable after splitting:

| | |
|---|---|
| no VCF line (nested / off-reference) | 8,356 (52.9%) |
| homozygous, stayed collapsed | 3,912 (24.7%) |
| half-missing — genuinely haploid, nothing to split | 2,442 (15.4%) |
| heterozygous, one slot lost to `min-reads` | 1,097 (6.9%) |

And **the reads are already there**: 99.9% of them have at least half their read placements also
occurring at a phaseable anchor, 58.5% have all of them, and only 18 of 15,807 have none. Median
distance to the nearest phaseable site is 2,081 bp.

So these sites are not short of spanning reads. They are declining on **confidence**, which points
at the split's own thresholds before it points anywhere else.

## Two independent levers

### A. `phase_min` / `phase_min_side` — currently not reachable at all

`AnchorParams::phase_min = 2.0` (nats) and `phase_min_side = 2` (reads per side) decide whether a
homozygous site splits. **Neither is a command-line flag.** They cannot be swept today.

### B. Raw mapQ for phase confidence, decoupled from calling

`allele_likelihood.cpp:333` stores `min(max(mapq, min_mismap), max_mismap)` once, and `:1514`
copies that floored value into `PhaseReadEvidence::mismap`. Phasing therefore inherits the calling
floor, and the fitted `--mismap-min 0.05` caps every read's per-site contribution at phred(0.05).

## Gate 0 for B — answered already, at no cost, and it demotes B

How much information does the 0.05 floor actually destroy on this data?

| chr20 ONT, 85,373 alignments | |
|---|---|
| MAPQ 60 | 80,823 (**94.67%**) |
| at or below the 0.05 floor (MAPQ >= 13) | 82,468 (**96.60%**) |
| above the floor, own `e_r` already used | 2,905 (3.40%) |

**The floor is not erasing a distribution here; it is rescaling a near-constant.** 94.67% of reads
share one MAPQ, so lifting the floor moves almost every read from `e = 0.05` to `e ≈ 1e-6` by
almost exactly the same factor. That matters because **`--regeno-temper` is fitted per run**
(`temper = -1` means fit it, via `fit_calibration`): a near-uniform inflation of every `Lambda` is
absorbed by a correspondingly smaller tau, and the calibrated log-odds the split actually
thresholds comes out close to where it started.

**Prediction to falsify: lifting the floor for phase confidence is close to a no-op on ONT for the
split decision.** The mechanism by which it could help — restoring distinctions among reads the
floor flattens — needs a spread of MAPQ that this dataset does not have.

One effect is **not** absorbed: the per-read `score`, and so the site `reliability`, has its ceiling
raised from phred(0.05) = 13.01 to phred(1e-6) = 60. That moves the distribution `--phase-min-q`
sits in — the exact coupling that stranded it at 9.5 — so B cannot be measured without re-fitting
that gate in the same experiment. But a rescaled gate is not new information.

**So B is demoted from first bet to falsification test**, run after A and cheap to run.

## Lever A, first: expose the thresholds and sweep them

**Change.** Two options on `AnchorParams`, defaulted to today's values so nothing moves until
asked: `--anchors-phase-min` (2.0) and `--anchors-phase-min-side` (2). Note that the run's
self-check at `graph_caller.cpp:1562` hardcodes its own "confident" cut at 2.0 to match
`phase_min`; it must follow the flag, or the check stops measuring the thing being swept.

**The accuracy measure is free, and it is held out.** The run already reports, at heterozygous
sites — where the cross-site phase can be compared against the partition the site's own alleles
make, leave-one-out so a site never judges its own reads — `phase_checked`, `phase_agree`,
`phase_confident`, `phase_confident_agree`. That is ground truth for exactly the inference the
split makes blind at homozygous sites. Relaxing `phase_min` must lower `phase_confident_agree /
phase_confident`; the sweep is asking how much range that buys per point of accuracy.

**Sweep** `phase_min` at 0.5, 1.0, 1.5, 2.0 (current), 3.0 on chr20 ONT. Then
`phase_min_side` at 1, 2, 3 at whichever `phase_min` won.

| measure | role |
|---|---|
| gaps > 33 kb between phaseable anchors, and their total span | **primary** — the goal |
| gaps > 10 kb, p99 gap, count of phaseable anchors | secondary range |
| `phase_confident_agree / phase_confident` | **accuracy of the split itself**, held out, free |
| whatshap switch error | the chain must not degrade |
| indel / SNV / ALL F1 | guard: splitting must not touch genotypes at all |

**Adoption.** Range must improve materially — call it gaps > 33 kb falling by 10% or more, since
splitting itself bought 49% — while the held-out split agreement falls by no more than ~1 point and
switch error does not regress beyond noise (about 0.026 percentage points, from 225 switches).
chr20 fits, chr6 confirms, exactly two arms on chr6.

**Cost.** Six chr20 arms at about 6 minutes, plus scoring; two chr6 arms. Half a day.

## Lever B, second: `--phase-mismap-min`

**Change, and it is small.** `PhaseReadEvidence` gains one `vector<float>` holding the
phase-specific mismap; `AlleleLikelihoodParams` gains `phase_mismap_min`, defaulted to
`min_mismap`; `allele_likelihood.cpp` computes and stores both; `graph_caller.cpp:6065` reads the
new one. One flag. **No `CallInfo` field**, so it does not cross the ploidy swap in
`run_deferred_descent` — the hazard that has silently dropped a field before. Memory is about 27 MB
on chr20, one float per retained read placement.

**Gate 1, byte identity.** With the default equal to the calling floor, the VCF and the anchor file
must be byte-identical to today's. If they are not, the change is not the change it claims to be.

**Gate 2.** `--phase-mismap-min 1e-6`, with `--phase-min-q` re-fitted in the same experiment
because its distribution moves. If the split decisions and the range measures come back within
noise of lever A's best, Gate 0's prediction is confirmed and B stops there, recorded as a negative.

**What would make B matter, and is worth saying so the negative is not over-read.** A dataset with a
real MAPQ spread — short reads, or repeat-rich regions where ONT MAPQ actually varies — has the
distribution this one lacks. The 3.4% of reads below MAPQ 13 are where the floor is already inert,
so they are not it either. Neither scheme captures confidence varying *within* a read, which is
what a chimera actually looks like.

## Order

1. **A**, because it is two flags against a new field, because the evidence points at it, and
   because its accuracy measure is already in the run.
2. **B** as falsification, only if A leaves range on the table.
3. If both stall, the remaining population is the 52.9% of survivors with no VCF line — nested and
   off-reference sites, where the limit may be that the anchor file is the only place they appear
   at all.
