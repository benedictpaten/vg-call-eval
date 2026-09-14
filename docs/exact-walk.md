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

## Implemented and measured (2026-09-14)

The walk was actually replaced, not just simulated, and the numbers moved. Both earlier
headline figures were artefacts and are superseded.

**+0.0052 indel F1 was too high.** The instrumented arm used the exact walk's SCORE with the
greedy walk's `--insertion-nats` bookkeeping -- `score_read_against_allele` ran first and
accumulated `nat_adjust`, then the exact score overwrote `score` alone. The two walks disagree
on 13-42% of cells, so they disagree on how many insertion events exist. With nats carried
along the DP's own path the gain is **+0.0042**.

**~2% runtime was far too low.** That came from an instrumented run that capped the DP at
`m*n <= 200000` and used the greedy score above it -- 28,330 ONT cells took that fallback.
Uncapped, the walk is 1.55x greedy on chr20 ONT and **3.73x on short reads**.

| chr20 ONT / chr6 ONT / chr20 short reads | greedy | DP | delta |
|---|---|---|---|
| indel F1 | 0.86237 / 0.88005 / 0.92858 | **0.86659 / 0.88351 / 0.92918** | +0.0042 / +0.0035 / +0.0006 |
| SNV F1 | 0.98582 / 0.98810 / 0.98518 | 0.98563 / 0.98816 / 0.98516 | -0.0002 / +0.0001 / -0.00002 |
| wall | 196.7s / 409.6s / 143.3s | 305.2s / 575.8s / 534.1s | 1.55x / 1.41x / **3.73x** |

Indel precision AND recall improve on every arm and nothing regresses; chr6 carries 82% of
chr20's gain. The cost profile is the problem, and its shape is diagnostic: short reads have
the SMALLEST accuracy change and the LARGEST slowdown, which is impossible if the cost were
per-cell DP work and follows directly if it is per-call overhead -- the function allocated ten
containers per (read, allele) pair, ~250M times on the short-read arm.

## Optional anchoring is load-bearing, and that kills the simple designs

A shared node visit guarantees identical BASES, so it is the mapper's own alignment. But the
read's alignment INSIDE that node can still be bad -- an ONT homopolymer run -- and the walk is
sometimes right to GAP it rather than pay for those edits. So a shared visit must be protected
from being SUBSTITUTED away, and must remain free to be gapped.

Isolated to one predicate in otherwise identical code:

| chr20 ONT indel F1 | |
|---|---|
| anchors optional (shipped design) | **0.86659** |
| greedy | 0.86237 |
| partition, anchors forced by construction | 0.86164 |
| same DP with anchors forced | 0.85948 |

**Optional anchoring is worth +0.0071**, more than the whole gain over greedy. That retires two
simpler designs that looked strictly better on paper:

- the **partition** (two-pointer anchor chain, tiny DP per divergence block, no band, no
  membership predicate, no pre-anchor state) forces every shared visit to match, and the same
  effect reproduces inside the global DP, so it was not an implementation bug.
- the **fast path** (read walks the allele exactly -> sum the per-node scores) forces matches by
  construction. Its byte gate FAILED against the DP, and since forcing is what costs, it cannot
  be repaired.

Ordering is NOT the reason either design fails. Both sequences are DAG paths, so shared nodes
appear in the same topological order, and repeats are absent in practice -- 0 in 234,001
traversals and in 19,999 reads averaging 936 nodes.

## The band is load-bearing, and the ONT measurement hid it

The band was very nearly deleted on the strength of an ONT-only measurement: it appeared to
recover 42s of 151s and to reach at most 0.1% of traversals. Both figures are true of chr20 ONT
and both are misleading, because **ONT calling is I/O-bound** -- 76% of it is the `gbz-base`
subprocess and 63% is blocked in `waitpid`, so DP cost is largely hidden behind the fetch.

Short reads are the honest measure of DP cost, and there the band is worth 1.79x the CPU:

| chr20 short reads, alone, `-t 6` | real | user CPU |
|---|---|---|
| banded (the shipped reference) | 484.5s | 2303.7s |
| unbanded, perfect-match cuts | 613.6s | 2083.4s |
| **unbanded, exhaustive** | **1343.4s** | **4123.9s** |

And it costs nothing. Against the unbanded walk over the same chr20 short-read call it moves a
**single record of 115,255**; on ONT, `dp-c20` against `bd-c20` leaves indel F1 identical at
0.86659 with SNV F1 marginally *better* banded.

## Perfect-match cuts: cheaper than the band, and much coarser

The proposal was to bound the DP by forcing every node the read matches PERFECTLY -- no edits
at all against that node -- to be paired, cutting the rectangle into independent blocks. A
perfect match earns `match` per base, which no gap or substitution can beat, so in a standard
global alignment the pairing is provably on the optimal path.

**This walk is not a standard global alignment.** Allele sequence before the first match and
after the last is FREE -- it is outside the read's window, and charging it would penalise a read
for being short. Taking a match closes that flank and makes every later allele node chargeable.
So for a read `[X, Y]` against an allele `[X, A, B, C, Y]`:

- force the `X` pairing: `match(X) - gap(A+B+C) + match(Y)`
- decline it: `-gap(len X)` to insert X in the flank, then A, B, C consumed **free**, `match(Y)`

When the intervening stretch is long the second wins, so forcing perfect matches is not exact
here and no smarter choice of cuts repairs it.

Measured, both arms unbanded so the cuts are the only variable: cuts move **43 records** where
the band moves **1**. Cuts are the cheaper bound and the worse approximation, so the band stays
and the cuts are dropped.

**The cut-bound cost four separate off-by-one errors**, each of which showed up only as an
allele becoming unreachable -- a relative likelihood of exactly 0, indistinguishable from
"scored, and hopeless". The subtlest: a row must range PAST its own forced pairing, because a
deletion is taken in the row before the match that follows it, so a row stopping at its own
column leaves the next row's diagonal never computed. That class of bug is now pinned by a test
over 11 node layouts asserting every read reaches every allele.

## Shipped: band + hoist, byte-gated on three arms

Work that depends only on the READ (its own per-node edit scores, and its sorted node keys) or
only on the ALLELE (its sorted node keys) was being recomputed inside the (read, allele) loop --
15.6M times on chr20 ONT, 24.6M on the short-read arm. Computing each once is a pure reordering
of identical arithmetic, and the byte gates prove it:

| chr20/chr6, alone, `-t 6` | before | after | bytes |
|---|---|---|---|
| chr20 short reads | 484.5s / 2303.7s user | **444.5s / 2139.9s user** | identical |
| chr20 ONT | 247.3s / 1356.9s user | **230.5s / 1356.2s user** | identical |
| chr6 ONT (hold-out) | -- | 402.5s / 2056.1s user | identical to `bd-c6` |

Short reads are the honest measure: **-7.1% CPU**. ONT user CPU is flat at 1356s because that arm
is I/O-bound, so its 17s of wall improvement is fetch variance, not the hoist.

Accuracy is unchanged by construction -- byte-identical output on all three arms -- so the walk's
gains stand as previously measured: chr20 ONT indel F1 0.86237 -> 0.86659, chr6 hold-out
0.88005 -> 0.88351, short reads +0.0006 with SNV F1 unchanged.

## The walk did not move the --insertion-nats optimum

`--insertion-nats` was fitted to 0.9 against the GREEDY walk, and the walk changes what a gap
costs in context, so it was re-swept on chr20 ONT with the shipped binary. 0.9 doubles as a
control and reproduces `bd-c20` to every digit, so the harness and binary agree with the
stored arm.

| `--insertion-nats` | indel F1 | SNV F1 | ALL F1 |
|---|---|---|---|
| 0.0 | 0.85668 | 0.98590 | 0.95637 |
| 0.3 | 0.86055 | 0.98576 | 0.95720 |
| 0.6 | 0.86402 | 0.98561 | 0.95795 |
| **0.9 (shipped)** | **0.86659** | 0.98563 | 0.95863 |
| 1.2 | 0.86704 | 0.98561 | 0.95875 |
| 1.5 | 0.86462 | 0.98544 | 0.95811 |

Smooth and unimodal with a peak at 1.2, but **0.9 is within 0.00045 indel F1 of it**, against
0.0026-0.0039 between neighbouring grid points elsewhere. That is a tenth of what the walk
itself bought (+0.0042) and the same magnitude that previously flipped between contigs.

**0.9 stands.** Moving it cannot be justified on this evidence: chr6 is the hold-out and may
confirm a chosen value but never choose between two candidates, and refining the chr20 grid
further would be fitting the fifth decimal. The useful result is the negative one -- the new
walk did not shift the optimum, so the preset needs no re-tuning.
