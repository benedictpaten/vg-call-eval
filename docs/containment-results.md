# Containment: bounding switch damage without fixing the decisions

**The idea** (user's): convert a modest number of relink errors into *minor* switch errors that do
not affect a large number of bases. The chain will not become globally correct; it will become
locally inconsistent instead of globally wrong.

**It works.** Measured on chr20, held out on chr6.

| chr20 | blocks | block N50 | mis-phased sequence |
|---|---|---|---|
| today (one block, always relink) | 1 | 66.21 Mb | **24.72 Mb** |
| refuse the 259 relinks whose spanning links disagree | 260 | 3.00 Mb | **3.30 Mb** |

| chr6 (hold-out) | blocks | block N50 | mis-phased sequence |
|---|---|---|---|
| today | 1 | 172.09 Mb | **52.39 Mb** |
| refuse the 110 disagreeing relinks | 111 | 5.51 Mb | **23.82 Mb** |

**7.5x less mis-phased sequence on chr20, 2.2x on chr6** -- and the number of parity-changing
events does not move (13 and 19). The errors are not fixed. They are contained: each now damages
a block instead of the whole contig. That is exactly the trade the idea asked for.

## Why it works, and the control that nearly killed it

The chain has no redundancy. Stage 1 is `o[m+1] = o[m] ^ (d[m] < 0)` and stage 2 orients each block
against the previous one, so one wrong sign flips everything downstream to the contig end. Cutting
anywhere bounds that -- so **cutting is not evidence that the rule is good.**

At EQUAL CUT COUNT the rule is no better than random, and on chr20 it is worse: 259 random relink
cuts give 2.28 Mb mis-phased against the rule's 3.30 Mb. The rule's cuts cluster, which is why it
reaches a far higher N50 for the same count.

So the comparison has to be at **equal N50**, and there the rule wins:

| chr20, random-cut frontier | N50 | mis-phased |
|---|---|---|
| 25 cuts | 6.58 Mb | 10.95 Mb |
| 50 cuts | 3.27 Mb | 7.18 Mb |
| 100 cuts | 2.05 Mb | 3.83 Mb |
| 259 cuts | 0.85 Mb | 2.28 Mb |
| **259 RULE cuts** | **3.00 Mb** | **3.30 Mb** |

At a matched N50 of ~3.0-3.3 Mb the rule gives 3.30 Mb against random's 7.18 Mb -- **2.2x better**.
chr6 holds out in the same direction but weaker: 23.82 Mb at N50 5.51 against roughly 30 Mb for
random interpolated to that N50, about **1.3x**.

Worth noting where the N50 lands: 3.00 Mb (chr20) and 5.51 Mb (chr6) against mean correct-phase runs
of 2.36 Mb and 5.93 Mb. Blocks cannot usefully exceed the correct-run length, so the rule sits
almost exactly at the ceiling.

## The signal it cuts on

A spanning link `phase_link(rel[m+1-L], rel[m+L])` for L in {2, 5, 10, 20, 50}, judged against the
SETTLED orientation `o[a] ^ o[b]`. It disagrees at 769 of chr20's 61,012 junctions (1.26%) and 1,078
of chr6's 162,518 (0.66%); of those, 259 and 110 are at relinks.

Per-event recall is the replicated part:

| events with >= 1 disagreeing spanning link | chr20 | chr6 |
|---|---|---|
| **relink events** | **2 of 4** | **2 of 4** |
| cascade events | 1 of 9 | 4 of 15 |

The split is mechanistically right: a span helps where stage 1 had no evidence and not where the
reads are uniformly wrong.

## What does NOT work

**Widening `--phase-relink` is inert.** 10, 20 and 30 give byte-identical output: 27 parity changes,
24,766,989 mis-phased bases. Stage 2's vote is a SUM, so more pairs add more of the same correlated
evidence and average the disagreement away. The signal has to be read as a disagreement, not
summed.

**Nothing reaches the cascade errors.** At those junctions every span agrees at full strength, and
read count, per-read evidence, magnitude and span all sit on the background. See
`lsweep-and-nested-strand-results.md`.
