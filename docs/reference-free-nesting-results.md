# Results: phasing the chains the reference does not cross

Implemented on vg branch `snarl-tree-order`. chr20 ONT, `--preset ont`. Off-reference admission is
still behind `VG_CALL_NO_REF_NESTED` (and a gref selection), NOT a default -- see the cost below.

## What was in the way, and what each change did

**1. Every off-reference child of one parent shared its parent's reference start as `anchor_position`.**
The linkage layer sorts each group on `(position, record key)`, so they all tied and fell back to the
record key -- a hash of the snarl. Their order was arbitrary. Fixed by adding the child's base offset
along the settled traversal that reaches it, composed down the tree.

**2. `site_gap` refused any pair with an unpositioned member.** Right for a MIXED pair -- an anchor
differenced against a real coordinate is not a distance -- and wrong when BOTH are anchored. A run is
grouped by `(parent, chain, ploidy)` before it reaches the model, so two unpositioned sites in one run
are in the same chain by construction, and the difference of their offsets is a real distance on the
haplotype the chain sits on. Before this the model got `SIZE_MAX` at every step inside such a chain
and forgot at every one.

**3. `records_for_render()` excluded `no_reference` records, and fed three consumers.** The exclusion
is right for exactly one of them: such a chain has no REF or POS to write. It is genotyped, it is
anchored, and its strand is meaningful, so excluding it from the phase-site construction and the
nested-strand cascade let the VCF's constraint decide what gets INFERRED. Split with a `for_phasing`
flag.

Before (3) the phase chain was **identical to the digit** with off-reference chains admitted --
76,135 het sites either way -- because the sites were genotyped and anchored and then dropped on the
way to the phaser.

## What it buys

| | reference-gated | off-reference phased |
|---|---|---|
| het sites in the phase chain | 76,135 | **77,374** (+1,239) |
| reliable sites | 60,203 | **60,695** (+492) |
| anchors | 636,484 | 657,898 |
| chain breaks | 1,084 | 1,269 |

**+1,239 heterozygous sites phased that previously were not**, in nested chains inside non-reference
alleles -- the population that matters for assembling a complex locus, and the one a reference
coordinate cannot order because the sequence is not on the reference.

## What it costs, measured

| | reference-gated | off-reference phased |
|---|---|---|
| whatshap true switches | 55 | **68** |
| whatshap flips | 85 | 79 |
| raw switch rate | 0.3826% | 0.3845% |
| held-out cross-site agreement | 94.7102% | **94.4171%** |
| re-genotyping round 1, genotypes moved | -- | 7,194 |
| convergence at 12 passes | period-3 limit cycle from round 7 | **no convergence, ~1,200/round** |

Three things to read carefully.

**Switch error here cannot see the gain.** whatshap assesses only sites it can place on the
reference, so all 1,239 new hets are outside its denominator. What this table shows is whether
admitting them DEGRADED the sites whatshap can see, and it did: +13 true switches. The gain is real
and this metric is structurally blind to it.

**Held-out agreement does cover the new sites**, and it falls 0.29 points. That is the honest measure
of the new population's phase quality, and it says they are harder than average rather than free.

**Convergence is worse, and this was a known and accepted trade.** The baseline enters a clean
period-3 limit cycle at round 7 with small numbers (403/258/0 ...). With off-reference chains in the
chain the iteration decays 7,194 -> 2,947 -> 1,778 -> 1,418 -> 1,090 and then plateaus around 1,200 a
round without converging or entering a cycle the detector recognises. The default `--regeno-passes 2`
runs one round, so the default is affected by the 7,194 but not by the plateau.

## Status and recommendation

Off-reference admission stays behind the env var. The coverage gain is exactly where it was wanted,
but it costs +13 true switches on reference-placeable sites, 0.29 points of held-out agreement, and
convergence -- and none of those was in the original bargain. Turning it on by default should wait
on chr6 as a hold-out and on a decision about whether the 1,239 sites are worth the 13 switches.

What is NOT in question is the plumbing: with off-reference chains skipped the default path is
**byte-identical**, VCF and anchors, so every change above is inert until asked for.
