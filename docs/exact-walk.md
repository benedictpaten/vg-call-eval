# The read-to-allele walk: greedy, exact, and what is actually load-bearing

Measured 2026-09-13 on chr20 and chr6, ONT and short reads. The question: the caller walks each
read against each allele with a **greedy single-pass** correspondence, and both node paths are
known, so an exact intersection ought to be possible. Is greedy costing accuracy?

**Answer: yes, but only in one of the two freedoms "exact" contains.** Optimising the
correspondence BETWEEN anchors gains accuracy; letting the optimiser RE-CHOOSE anchors destroys
it, by an order of magnitude more than the first gains.

## What the walk does

For each read node it searches ahead in the allele for a match (an *anchor*). Allele nodes skipped
between two anchors are a deletion; read nodes the allele lacks are an insertion; where neither
applies the two nodes are scored as substitutes. It is O(m*n) worst case already, because the
anchor search rescans the allele per read step -- so an exact DP is **not** asymptotically worse.

Sizes, from 234,001 chr20 allele traversals: nodes per traversal median **3**, p90 3, p99 12,
p99.9 126, max 3,640; 0.219% exceed 50 nodes. So the median DP is nine cells.

## The greedy walk is not an approximation to a standard alignment

This took four attempts to discover, each caught by one check -- *an optimal DP can never score
lower than greedy* -- rather than by reading the code:

| attempt | defect | cells where greedy won |
|---|---|---|
| 1 | charged gap_open per skipped node; greedy sums them into ONE affine gap | 49 |
| 2 | affine fixed, but forbade I<->D adjacency | 7,887 |
| 3 | allowed I<->D -- the diagnosis was simply wrong, count unmoved | 7,887 |
| 4 | pre-anchor state | **0** |

The real cause: the free-flank rule is gated on `have_anchor`, which is set **only by a true
match**. A leading *substitution* does not close the flank, so allele nodes skipped after it stay
free. That is path-dependent and has no Needleman-Wunsch formulation; it had to be reverse-
engineered into a fourth DP state. Every "disagreement rate" measured before that point --
41.78%, 13.17%, 15.52% -- was substantially measuring these bugs.

## Verified disagreement, once the DP dominates greedy everywhere

| | ONT chr20 | short reads chr20 |
|---|---|---|
| (read, allele) cells | 15,604,944 | 24,650,799 |
| correspondence differs | **41.83%** | **16.04%** |
| exact ever worse | 0 | 0 |
| reads changing preferred allele | **0.295%** | **0.340%** |
| \|delta\| >= 64 score units | 8% | 45% |

It is **not** ONT-specific, and short reads flip slightly MORE reads. For scale, base-level WFA
realignment moved 0.044% of reads and was rejected as not worth it; this is 7-8x that.

## Two freedoms, opposite signs

**Free substitution** -- any read node may pair with any allele node -- is destructive. It lets an
insertion allele explain reads that do not carry the insertion: insertion FP 1,418 -> 3,695.

**Anchor-constrained** -- a node the two paths share may pair only with itself, and the DP
optimises only where they genuinely diverge -- is the improvement.

chr20 ONT indel F1 by gap scale:

| gap_open | greedy | exact, free substitution | exact, anchor-constrained |
|---|---|---|---|
| 1 | 0.85452 | 0.79922 | **0.86758** |
| 2 | 0.81411 | 0.75580 | 0.82313 |
| 3 | 0.79297 | 0.73564 | -- |
| 4 | 0.78663 | 0.73154 | -- |
| 6 | 0.78490 | 0.72965 | -- |

**The gap scale is not the explanation.** Greedy and free-substitution run parallel at an offset of
-0.055 to -0.058 across a sixfold range, and all three walks share the same optimum at gap_open 1.
A mis-set parameter would show as curves crossing; none do. So the anchored walk needs no
re-tuning -- it wins at the shipped values.

## The anchored walk, validated

| | chr20 ONT | chr6 ONT (held out) | chr20 short reads |
|---|---|---|---|
| SNV F1 | +0.00015 | -0.00005 | -0.00002 |
| **Indel F1** | **+0.00521** | **+0.00362** | **+0.00055** |
| ALL F1 | +0.00128 | +0.00073 | +0.00012 |
| SV F1 (truvari) | **+0.0033** | +0.0000 | -- |

Precision AND recall up on indels on every arm; no class regresses anywhere. chr6 carries 69% of
chr20's indel magnitude. The effect scales with read length, as the mechanism predicts: long reads
cross more nodes per site, so there is more correspondence for greedy to get wrong.

Cost: computing BOTH walks ran 226s against a ~221s baseline, so the exact walk alone is
**essentially free**.

## Recommendation

Worth implementing, as a **follow-up and not in the current PR** -- it is a change to the scoring
core, and the branch already carries seven commits of measured work.

Implementation notes for whoever picks it up. The DP must reproduce the walk's own semantics or it
loses to greedy: affine gaps coalesced across consecutive skipped nodes, I<->D adjacency permitted,
and a pre-anchor state in which allele deletions are free until the first true match. Keep the
anchor constraint -- shared nodes pair only with themselves. Cap m*n (200,000 was used here;
28,330 ONT cells and 18,105 short-read cells exceeded it and fell back to greedy) or use an
O(ND) diff, since traversals reach 3,640 nodes. Gate it on the same invariant used here: greedy
must never beat it.
