# Bounding switch damage: running plan

**The problem.** The phasing chain has no redundancy anywhere. Stage 1 is
`o[m+1] = o[m] ^ (d[m] < 0)` and stage 2 orients each block against the previous one, so a single
wrong sign -- cascade or relink -- flips everything downstream to the end of the contig. That is
why **13 parity-changing events produce 24.72 Mb of mis-phased sequence on chr20** (19 events,
52.39 Mb on chr6).

**The goal**, in the user's words: convert a modest number of relink and hang errors into *minor*
switch errors that do not affect a large number of bases, by ensuring parity with the cascade
chain. The chain will not become globally correct; it will become locally inconsistent instead of
globally wrong.

**The metric is mis-phased BASES, not switch count.** A good result may well raise the switch count
while cutting mis-phased bases sharply. Report both, plus block N50 and the whatshap figures, so
nothing hides.

## What is already established

| | chr20 | chr6 (hold-out) |
|---|---|---|
| parity changes, all sites | 27 | 28 |
| parity changes, reliable backbone only | **13** | **19** |
| mis-phased sequence, all sites | 24.77 Mb | 52.66 Mb |
| mis-phased sequence, backbone only | **24.72 Mb** | **52.39 Mb** |
| mean correct-phase run, backbone | 4.73 Mb | 8.60 Mb |

The rehung sites add only 50 kb (chr20) and 270 kb (chr6) of mis-phased sequence despite
contributing 14 and 9 extra switches and 71 and 177 extra flips. **Containment already holds for
stage 3.** The entire prize is in the cascade and the relinks.

Four facts that constrain the design:

- **`|d|` is a read counter, not a confidence.** Correlation with the site's read count is
  r = 0.903 and the per-read contribution is a near-constant 0.74 regardless of depth. Raising
  `--phase-break` is a coverage filter: catching all switch junctions needs a threshold of 39,
  which breaks 92.2% of cascade links.
- **Unanimity is the norm** -- 56.9% of linked junctions have zero minority votes -- so the 46/0
  votes at switch junctions are unremarkable.
- **Single-site bypass does not help.** Triangle closure is 99.98%, and at all switch junctions both
  bypasses close at full strength. A bad *site* corrupts both its adjacent links and cancels (a
  flip); a bad *link* does not, and these do not.
- **The failing regions are unassessable, not empty.** Inside switch intervals the ratio of
  whatshap-assessed to truth variants is 0.028 (chr20) and 0.007 (chr6) against a control of 0.50 --
  an 18x to 70x collapse, with truth variant density up to 15x the chromosome average.

## Step 1 -- the L-sweep (RUNNING)

The bypass test skipped ONE site, a few hundred bases, so it stayed inside the mis-specified region
and used the same mis-aligned reads. Reads span a median 21 kb and these regions are ~2 kb, so a
link from a clean site ~10 kb before to a clean site ~10 kb after is carried by reads whose alleles
at BOTH endpoints are fine. Nothing in the algorithm computes it.

Instrumented: at each junction, `phase_link(rel[m+1-L], rel[m+L])` for L in {2, 5, 10, 20, 50},
dumped with the shared-read count, the span in bases, and the parity the cascade implies over the
same stretch.

**Decisive question.** At the 13 + 19 backbone switch junctions, does a spanning link disagree with
the cascade? If yes at L = 10-20, the information exists and is unused. If it agrees at every L, the
reads are wrong over a span longer than a read and no link-based scheme can recover it -- which
closes the question rather than leaving it open.

## Step 2 -- global orientation solve (only if step 1 is positive)

Replace the chain with a solve. Every link, adjacent and spanning, is a signed constraint between
two sites; choosing `o[]` to maximise agreement is weighted correlation clustering. A greedy
maximum-weight spanning tree, or a few rounds of iterative local flipping, is enough to test the
idea. The point is not that more decisions come out right -- it is that a wrong constraint is
outvoted locally instead of propagating.

Prototype offline from the same dump before touching the caller.

## Step 3 -- nested haploid strand asymmetry

chr20's 3,117 half-missing records split 2,105 on strand 0 against 1,012 on strand 1 (67.5/32.5),
while top-level phased hets are balanced at 49.9/50.1. That is ~19 sigma from parity.

It may be legitimate -- a nested chain sits inside one parent traversal, so if nested structure is
preferentially carried by one kind of allele the bias is real. Condition on the parent's genotype
and on which of the parent's traversals carries the child; if the bias survives, it is a bug.

This corner is invisible to every gate in use: whatshap never assesses a half-missing record and the
genotype does not move, so neither switch error nor F1 can see it. It has produced a silent strand
bug before.

Context established: a nested DIPLOID site is a first-class phase site and its hang reaches
`phase_flips` normally (stage 3 sets `decided[t] = 1`, and the flip set is `decided[t] && o[t]`),
with the swap cascading down the nesting tree -- 4,369 nested strands carried on chr20. A nested
HAPLOID site is excluded from the chain entirely (`phase_sites` skips `ploidy != 2`), so its parity
with the parent is guaranteed by construction but never verified against its own reads.

## Gates

- chr20 to develop, chr6 held out, as always.
- Every arm is phasing-and-anchors only, so the VCF genotypes must not move: gate on the VCF body,
  not on F1.
- All instrumentation is env-gated and reverted before the PR is touched.
