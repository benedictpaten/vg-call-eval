# Hom splitting does not bridge switch events

`--anchors-hom-split` splits 87,370 of chr20's 92,603 homozygous sites (94.3%) and 222,021 of
chr6's 224,478 (98.9%). The gate is a crash detector, not a quality filter. The question was
whether that aggression extends haplotypes **across the phasing's switch events**, where the
labels are wrong on one side. Measured on both contigs: **it does not, and no setting of
`--split-min-q` or `--split-min-side` changes a single junction.**

## What makes this one run per contig

A hom site is never a phase site -- `phase_sites` skips `trav_first == trav_second`
(`graph_caller.cpp:6178`) -- so `read_strand_log_odds`' leave-one-out never fires there and `lo`
is the read's global calibrated lambda, **the same number at every hom site it crosses**. The
gate is then a pure function of (pinned read set, one number per read), and both survive a single
run: an unsplit site emits one slot holding every pinned read and a split one emits the same reads
across two, so the read set is in the file whatever the gate did, and `reliability` is computed
over the site's one distinct allele so it is invariant to splitting. Every combination was swept
offline from one chr20 run and one chr6 run; a third chr20 run supplied the untilted het control.

## The result

For each true switch, take the two consecutive assessed het sites bracketing it and ask what
connects them: a read pinned at **both** (direct), a chain through split homozygotes
(split-bridged), or nothing.

| | chr20 | chr6 |
|---|---|---|
| true switches (whatshap `switch/flip`) | 34/75 | 27/193 |
| switch junctions resolved | 28 | 26 |
| **direct** | **28** | **26** |
| split-bridged | 0 | 0 |
| unbridged | 0 | 0 |

Identical at `--split-min-q` 0.5 through 8 and at a coverage-scaled `--split-min-side` of
`max(2, 0.33 x n)`, which collapses 30% of chr20's splits and 27% of chr6's. Also identical in
every het-gap stratum, including the 15-40 kb junctions.

**Why.** The library is ultralong: reads span a median 21 kb of positioned sites, 89 kb at the
90th percentile, 535 kb at most, while consecutive assessed het sites sit a median **260 bp**
apart on chr20 and 273 bp on chr6. The het chain bridges every switch directly and would do so
with hom splitting switched off entirely. Splits are never the only path across a switch.

## What the splits are worth, since not this

Connectivity over the haplotype-carrying sites, hets alone against hets plus splits:

| min reads per link | hets only, N50 | + splits (default), N50 | + splits (q=8, frac=0.33) |
|---|---|---|---|
| 2 | 37.2 Mb | 37.2 Mb | 37.2 Mb |
| 5 | 13.8 Mb | **21.7 Mb** | 21.7 Mb |
| 10 | 2.31 Mb | **4.67 Mb** | 4.37 Mb |
| 20 | 171 kb | 70 kb | 85 kb |

At a permissive link rule the chromosome is one chain either way. Where a link must be
well-supported the splits nearly double the N50 -- and the strict gate keeps almost all of that
(4.67 -> 4.37 Mb) while dropping 30% of the splits. At 20 reads per link the splits *fragment* the
graph and the strict gate is the better of the two.

## Is the gate selective? Yes, but what it selects has no hom analogue

Run the hom gate on **het** sites, where the site's own alleles are held-out ground truth for the
inference the split makes blind. Requires `--no-anchors-phase-hets`, or the slot partly follows
the strand and the check measures itself. 76,796 sites, 3,130,329 read observations, baseline
agreement 95.0054%.

| `--split-min-q` | side rule | sites kept | agreement kept | agreement dropped | disagreements removed |
|---|---|---|---|---|---|
| 0.5 | >= 2 (default) | 99.5% | 95.08% | **71.62%** | 1.7% |
| 8 | >= 2 | 96.9% | 95.47% | 76.67% | 11.7% |
| 0.5 | max(2, 0.33n) | 91.7% | 95.43% | 90.04% | 15.6% |
| 8 | max(2, 0.33n) | 87.4% | 95.82% | 89.02% | 26.3% |

So the gate does pick out bad sites -- but weakly. The worst 10% of sites carry **64.3%** of all
disagreements and the worst 1% sit at 43.7% agreement, worse than a coin; the best gate setting
spends 12.6% of sites to reach 26.3% of the disagreements.

**And that population is not the split's problem.** Inside the worst 1% of sites, agreement is
flat in `|lo|`:

| `\|lo\|` band | worst 1% of sites | all sites |
|---|---|---|
| 0 - 0.5 | 40.2% | 72.2% |
| 0.5 - 2 | 43.8% | 88.2% |
| 2 - 8 | 46.2% | 93.8% |
| >= 8 | **44.9%** | **95.9%** |

A read with overwhelming strand confidence disagrees with those sites half the time. That is not
the strand failing -- it is the site's own allele partition being random, the *reference* of the
check rather than the inference. A homozygous site has no allele partition, so it cannot fail this
way. The proxy's worst decile is an artefact of the measuring instrument.

## The floor of 2

It is inert at this depth and correctly so. The smaller confident side has a median of **14** reads
and sits at a mean fraction of **0.407** of the confident reads -- near the 0.5 a balanced diploid
site gives. It scales with depth, which is what makes the fixed floor safe:

| site depth | sites | median smaller side | smaller side <= 2 |
|---|---|---|---|
| < 10 | 387 | 3 | **45.5%** |
| 10-19 | 3,226 | 4 | 22.2% |
| 20-29 | 12,390 | 9 | 1.9% |
| 30-44 | 46,685 | 15 | 1.0% |
| 45+ | 24,279 | 20 | 1.1% |

So a coverage-scaled `--split-min-side` is a **low-coverage** guard, not a switch lever: it would
change 4.2% of chr20's splits at 44x and much more below 20x. Nothing measured here argues for it
at this depth, and it is a flag with no use on the shipped presets.

## Conclusion

Leave the gate alone. The premise that splitting extends haplotypes across switches is false for
this data at both contigs and every setting; tightening costs real linkage at strict link rules and
buys no measured accuracy; and the sites that look worst in the only available ground truth fail on
an axis homozygous sites do not have.

## Corrections to the earlier sweep

`split-min-q-cannot-reduce-switches.md` reported agreement identical to four decimals at every
threshold. That number is the het self-check, which does not depend on `--split-min-q` at all -- a
constant, not a null result. The parameter does move the output (87,456 -> 74,176 sites split); it
was never measured. This document measures it. The doc's other finding stands: `--split-min-q`
never reaches placement, so both levers can only collapse whole sites.

The switch/flip figures here (34/75 chr20, 27/193 chr6) are from runs with `--anchors-out`, which
implies off-reference nesting and so calls a different VCF from the plain `vg call` arm's 27/76.
Each arm is compared against its own phasing, which is the only correct pairing.
