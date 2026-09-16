# Using raw mapQ for phase confidence: measured and rejected

The proposal was that the mismapping clamp should be used less aggressively when judging PHASE than
when judging a genotype. The argument is good: the floor exists because a read that fits one allele
perfectly must not claim certainty about a genotype the aligner never tried, which is an argument
about LOCAL alignment. Long-range phase is a different claim -- forty well-phased heterozygous sites
really should make a read's haplotype close to certain -- and the floor damps exactly that, because
the mismap term never leaves the denominator.

`--phase-mismap-min` decouples the two. Measured on chr20 ONT, vg `b661b0086`.

## Result: the default is the best value, and lowering the floor costs range

| `--phase-mismap-min` | held-out agreement | sites split | left collapsed | fitted temper | **phased-run N50 (read-walkable)** |
|---|---|---|---|---|---|
| **0.05 = `--mismap-min`, the default** | 94.7102% | 85,678 | 3,729 | 0.08 | **2,798,104 bp** |
| 0.01 | 94.7166% | 85,076 | 5,838 | 0.03 | 2,313,753 (-17%) |
| 0.001 | 94.6139% | 84,283 | 6,829 | 0.02 | 1,526,812 (-45%) |
| 0.0 (no floor at all) | 94.7602% | 85,245 | 5,914 | 0.02 | 1,707,287 (-39%) |
| 0.5 (reachability probe) | 79.4631% | 71,357 | 17,893 | 0.04 | 805,834 (-71%) |

**The hypothesis is falsified.** Trusting raw mapQ more does not lengthen phase blocks; it shortens
them by 17% to 45%.

## Why, and it is not "the temper absorbs it"

The temper column shows the compensation happening: τ falls 0.08 -> 0.03 -> 0.02 as the floor drops,
because the calibration fit sees a larger raw Λ and tempers harder to restore agreement. And
agreement IS restored -- it stays flat at 94.6-94.8% across every arm from 0.0 to 0.05.

But the compensation is not neutral, because **agreement is a property of the sign and range is a
property of the ranking**. τ is one scalar over the whole run. A read's raw Λ grows with the number
of sites it crosses, so a smaller floor inflates the well-covered reads more than the thinly-covered
ones, and a single τ chosen to keep the population calibrated then pushes the thin ones below
`--split-min-q`. Fewer reads clear the gate, fewer sites have two confident reads on each side,
fewer splits: 85,678 -> 84,283. Collapsed sites nearly double, 3,729 -> 6,829.

So the argument for the proposal was right about the raw log-odds and wrong about what limits phase
range. The clamp is not what stops forty SNPs from phasing a read confidently; the clamp is holding
the *spread* of confidences in a place where a single temper can calibrate them all at once.

## The null result that was not one

Every arm of the first sweep came back **byte-identical**, which looked exactly like the predicted
"the temper absorbs it" outcome and would have been reported as such.

It was a dead wire. `PhaseReadEvidence` is built only when anchors are NOT armed -- "the anchor
evidence already carries everything it wants" (`src/allele_likelihood.cpp:1328`) -- so with
`--anchors-out` the phasing converts its evidence from `AnchorSiteEvidence`, which carried only the
genotype-clamped probability. The new floor never reached the code under the flag every arm was run
with.

What caught it was refusing a null without proving the knob was connected: the floor was **raised**
to 0.5, where it must damp Λ hard. Under `--anchors-out` it still changed nothing; without
`--anchors-out` the same value moved the temper 0.08 -> 0.04, the reliable het sites 60,203 -> 0, and
the VCF. `AnchorRead` now carries both floors.

**A sweep that reports no effect has not measured anything until one arm is shown to move something.**

## Status

`--phase-mismap-min` ships, defaulting to `--mismap-min`, which is byte-identical to not having it.
It is kept because the question recurs and the answer should be re-measurable, not because any value
other than the default is recommended. **Do not re-tune it without re-reading this file.**
