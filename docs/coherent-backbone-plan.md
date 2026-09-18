# A coherent backbone for long-range phasing

Target: a phasing whose backbone is internally consistent, so that long-range phase is as accurate
as the reads allow. Everything below rests on measurements already in this repo; each design choice
names the one that motivates it.

## What the measurements say

1. **Coherence predicts switches; reliability does not.** Ranked worst-first against chr20 switch
   positions, coherence enriches 9.3x in its worst 0.1% and 5.1x in its worst 1%; reliability gives
   2.2x and 1.1x and sits BELOW the base rate at 5%. r = 0.437, so it is new information.
   [phase-coherence.md](phase-coherence.md)

2. **Stage 1 is the weak part, and the sites it fails on are perfectly phaseable.** Demoting 2.6% of
   sites from the cascade to stage 3's hanging rule cut true switches 51 -> 33 on chr20 and 62 -> 31
   on chr6 **with flips flat** (77 -> 76, 198 -> 194). If demotion merely converted a propagating
   error into a local one, flips would have risen as switches fell. They did not, so stage 3 phases
   these sites CORRECTLY and the cascade could not.

3. **The difference is one neighbour against six, not weighting.** Stage 1 chains adjacent sites and
   reads only the sign; stage 3 aggregates up to six neighbours weighted by |dv| plus a panel prior.
   A lone link has no competing evidence to weight against, which is why `--phase-break` -- the only
   thing a single link's magnitude can drive -- was measured at 20/40/80 and does nothing.

4. **The cascade is already at a local optimum of read likelihood.** `--phase-cp` aggregates, for
   every junction, the log10 gain of flipping everything downstream over reads SPANNING it, using
   each read's full span. It fires on **3 junctions in all of chr20**. So the problem is not finding
   a better assignment given the sites; it is which sites are in the backbone at all.

5. **Over-demotion hurts.** Threshold 0.80 demotes 2,489 sites and scores 39 switches against 33 at
   0.70 with 1,592. There is a real trade between removing a bad site and fragmenting the chain
   around it.

6. **Residual switches are regional, not per-site.** After the gate, 20 of chr20's 33 sit in two
   megabases, and 6 of the 7 "isolated" ones lie inside a single 4.5 Mb window. No per-site gate
   will reach those.

## The algorithm

### Step 1 -- candidate backbone
Sites passing `--phase-min-q` may carry links, as now. This is a necessary condition, not a
sufficient one: reliability says the reads can separate the site's alleles, nothing more.

### Step 2 -- phase the backbone by aggregation, not by cascade
Replace the single-link sign-only chain. Build a local graph in which each backbone site is joined
to its k nearest backbone neighbours on each side (k ~ 3), with edge weight `d_ij = phase_link(i,j)`.
Choose orientations maximising

    sum over edges of |d_ij| * agree(o_i, o_j, sign d_ij)

by local search from the current cascade as a starting point:

- **single-site moves**: flip site i if the net weighted evidence from its neighbours says so. This
  is stage 3's rule applied to the backbone itself, iterated to a fixed point.
- **suffix moves**: single-site flips cannot cross a switch, which requires inverting a whole
  suffix. `--phase-cp`'s S(j) already computes exactly that gain, so alternate the two move types.

Both moves increase the same objective, so the search terminates. Finding 3 junctions with positive
S(j) says the suffix moves will be rare -- which is the point: they are a safety net, and the
single-site moves are where the gain should be.

### Step 3 -- prune to a coherent fixed point
Measure per-site coherence against the orientation just derived; demote sites below the threshold;
re-derive; repeat until no site is demoted. The invariant is **every site in the backbone is
coherent with the backbone it belongs to**, which the one-shot pass does not give.

Fragmentation is the risk and must be measured, not assumed away: each round removes sites, the
surviving links span further, `phase_link` falls off with distance, and more links drop under
`--phase-break`. Report chain breaks per round and stop when switch error stops improving.

### Step 4 -- hang everything else
Stage 3 unchanged. Demoted sites are NOT dropped: they are phased by hanging and appear in the
output like any other site.

## What this cannot fix, and should not be expected to

The residual switches are three regions -- 26-35 Mb pericentromeric at 0.64x depth, 50-52 Mb, and
65.3-65.5 Mb at 6.7x het density and normal depth. 65.50-65.60 Mb is independently the single most
fragmented 100 kb in the anchor graph. Those are a mapping failure and a probable graph collapse of
paralogous sequence. A phasing algorithm reading those reads cannot recover the truth from them, and
the honest goal for this work is the DISTRIBUTED background, which the coherence gate has already
largely removed: isolated switches went from one per 2.3 Mb to one per 5.1 Mb.

## Gates

chr20 selects, chr6 validates, and nothing ships on chr20 alone. Required of any change:
switch error down on BOTH contigs; F1 not down on either; the whatshap assessed denominator not
shrinking (a re-genotyping arm can flatter itself by dropping sites); and chain breaks reported
beside the switch count, because a chain cut into fragments trivially has fewer switches in it.
