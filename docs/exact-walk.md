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

## The full design space, chr20 ONT, gap_open 1

| variant | SNV F1 | SNV FP | Indel F1 | Indel FP | ALL F1 |
|---|---|---|---|---|---|
| greedy (shipped) | 0.98582 | 391 | 0.86237 | 3396 | 0.95787 |
| **anchored + substitution** | **0.98597** | **377** | **0.86758** | **3279** | **0.95915** |
| anchored, no substitution | 0.98569 | 392 | 0.86646 | 3304 | 0.95867 |
| free substitution | 0.98410 | 528 | 0.79922 | 5132 | 0.94186 |

**Keep substitution between mutually-unique nodes.** Forbidding it -- which would make the problem
a pure LCS and let Myers apply directly -- costs **+15 SNV false positives** and -0.00112 indel F1.
The mechanism is exactly as predicted: a SNP bubble then scores as insertion + deletion rather than
a base-level mismatch, which at gap_open 1 is -2 against -4, so the disagreement gets CHEAPER and
discriminates less. It still beats greedy (+0.0041 indel), capturing about 80% of the gain, so it
is a defensible simplification -- just not a free one.

## Affine across nodes is not a variable at the ONT preset

`score_gap(L) = -(gap_open + (L-1)*gap_extend)`. At `gap_open == gap_extend == 1` this is exactly
`-L`, so two adjacent skipped nodes cost `-L1 - L2` separately and `-(L1+L2)` coalesced -- the same
number. **The difference between affine and linear is `(gap_open - gap_extend)` per extra node**,
which is zero under `--preset ont`.

Measured: affine and linear give **byte-identical VCFs** in both anchored rows above, and the path
is reachable (10.79% of allele pairs have a leftover gap run of >= 2 nodes). Not a data null, an
arithmetic identity.

This also explains the insertion-coalescing result in `indel-uncertainty.md`, which had been filed
as a lucky null: byte-identical on ONT (1/1) and moving 53 genotypes on short reads (6/1). Not a
coincidence -- `gap_open - gap_extend` is 0 and 5 respectively. **To test affine against linear at
all, run at `gap_open != gap_extend`.**

## How to implement the anchoring

Not as a per-cell predicate. The set-membership test used in the measurement above -- "does this
node appear anywhere in the other sequence" -- answers a global question locally and is wrong
under repeats: it forbids a substitution on account of an occurrence already consumed elsewhere.
Make the anchors a PARTITION and the predicate disappears.

**Phase 1, anchor selection.** Two-pointer scan matching shared nodes in order. Exact whenever both
sequences are node-simple, which is what the data are (0 repeats in 234,001 traversals; 0 in 19,999
reads averaging 936 nodes). **Guard it**: if any node occurs more than once in either sequence, fall
back to LCS on node symbols -- O(MN) DP or Myers O(ND). Common case O(m+n); the pathological case
correct rather than silently wrong.

**Phase 2, score each divergence block.** The anchors partition both sequences into blocks of `p`
read-only against `q` allele-only nodes. **Inside a block there are no shared nodes by
construction**, so substitution needs no test at all. A small `p x q` DP picks how many to pair;
the leftovers on each side are one contiguous run, so each is a single `score_gap(total_bases)`.
Median divergence is 1 node against 1, so the typical block is one comparison.

**Phase 3, flanks.** Before the first anchor and after the last, allele nodes are free and read
nodes are charged -- reproducing the `have_anchor` semantics as the first and last block rather
than as a DP state.

What the partition formulation removes, relative to the four-state DP measured here: the Gotoh
I/D states (a leftover run is contiguous), the pre-anchor state (the flank is just a block), the
substitution predicate (blocks hold no shared nodes), and the O(MN) floor. Three of the four
things that each cost a failed attempt were artefacts of forcing this into one global alignment.

**Ship the invariant that caught all four**: in debug builds, assert the new walk never scores
lower than the greedy walk.

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
