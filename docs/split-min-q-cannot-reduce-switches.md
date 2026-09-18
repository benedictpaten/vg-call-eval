# --split-min-q cannot reduce hom-split switches, and why

Swept on chr20 ONT, vg `df32ee7e5`, with `--no-anchors-phase-hets` so the run's own held-out check
stays a measurement rather than a mirror (see the last section).

## The result: identical at every threshold

| `--split-min-q` | sites split | left collapsed | agreement over ALL placed reads |
|---|---|---|---|
| 0.5 (default) | 87,456 | 5,479 | **94.2999%** |
| 1.0 | 85,653 | 7,282 | **94.2999%** |
| 2.0 | 82,735 | 10,200 | **94.2999%** |
| 4.0 | 79,567 | 13,368 | **94.2999%** |
| 8.0 | 74,176 | 18,759 | **94.2999%** |

2,939,707/3,117,403 on every row. Raising the threshold to 8 discards 13,280 split sites, 15% of
them, and changes the accuracy of the output by nothing at all.

## Why: the gate and the placement use different rules

```
GATE       -- which sites may split      (anchor.cpp)
    if (|lo| < params.phase_min) continue;        // --split-min-q
    lo > 0 ? ++side0 : ++side1;
    split if side0 >= phase_min_side && side1 >= phase_min_side;

PLACEMENT  -- which slot each read joins (anchor.cpp)
    if (isnan(lo)) continue;                      // phase break: dropped
    if (lo == 0.0) coin = hash(read.name) & 1;    // no opinion: deterministic coin
    best_slot = coin >= 0 ? coin : (lo > 0.0 ? 0 : 1);   // NO |lo| threshold
```

A site qualifies on its confident reads and then places every read with any nonzero opinion by
sign alone. `--split-min-q` never reaches the placement, so the reads that land in slots, and the
slots they land in, are invariant under the sweep. **The rising "confident" column in the run's own
report is not an improvement in the output** -- it is the accuracy of a progressively smaller
reported SUBSET of the same unchanged placements, and it is easy to quote as a gain.

## What the sweep does measure: the confidence gradient

**"correct" below means AGREES WITH THE ALLELE MATCH, not agrees with truth.** The check compares
the read's cross-site strand against the slot the site's own sequence put it in; neither side is
ground truth and both carry error. From the 6.44% het-to-het flip rate, solving `flip = 2p(1-p)`
puts the allele match's own error at about **3.3% per site**, so at |lo| >= 8 the 4.49%
disagreement could be mostly the ALLELE side and the strand's own error nearer 1.2%. Independence
is doubtful -- both read the same `rel_at` values and mismapping is correlated along a read -- so
treat that as a bound rather than a decomposition.

Differencing the confident bands:

| `\|lo\|` | reads | share | agrees with allele match |
|---|---|---|---|
| < 0.5 | 50,546 | 1.6% | **71.7%** |
| 0.5 - 1.0 | 47,121 | 1.5% | 84.8% |
| 1.0 - 2.0 | 89,202 | 2.9% | 89.1% |
| 2.0 - 4.0 | 185,166 | 5.9% | 91.3% |
| 4.0 - 8.0 | 394,123 | 12.6% | 93.7% |
| >= 8.0 | 2,351,245 | 75.4% | **95.5%** |

So a PLACEMENT-side threshold -- the change `--split-min-q` cannot make -- is the only lever, and
this bounds it. Excluding the worst 1.6% moves the output 94.30% -> 94.67%; excluding 24% of all
reads reaches only 95.51%.

## Two reasons not to add one

**The apparent ceiling may not be the strand's.** Reads at |lo| >= 8 nats agree with the allele
match only 95.5% of the time, which LOOKS like the correlated-mismapping floor
[[read-strand-confidence-has-a-floor]] identified. But with the allele match itself erring ~3.3% per
site, most of that 4.5% may be the reference rather than the estimator, and the strand may have more
headroom than this table shows. **This argument is weaker than it first appears and should not be
used to close the question**; a measurement against real truth, not against the allele partition,
is what would settle it.

**Even the worst reads beat a coin.** The sub-0.5 band is 71.7% correct against 50% for the
deterministic coin the `lo == 0.0` branch uses, so coining them would make the output worse.
Dropping them costs coverage at split homozygotes, which is what produced the deletion-shaped holes.

The untested lever is upstream: `--regeno-ceiling`, the per-read error floor in
`calibrated_log_odds`, is pinned at 1 (the plain logistic) because fitting it cost ALL F1 -0.0014
and indel -0.0047. It is the only term that models this floor directly.

## The self-check is no longer held out under the new default

`--anchors-phase-hets` is on by default, so a read's slot at a het site partly follows its strand --
and the self-check asks whether the strand agrees with the slot. Measured: 94.30% with
`--no-anchors-phase-hets`, 98.77% with the default tilt, and **100.0000%** under
`--anchors-strict-hets`, where it is an identity (3,118,918/3,118,918).

vg `df32ee7e5` suppresses the report under both flags rather than relabelling it, because a number
that looks like an accuracy gets quoted as one. **94.30% is the split's blind accuracy; the other
two figures measure the placement rule against itself.**
