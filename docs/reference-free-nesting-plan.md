# Plan: genotype, order and phase the snarl tree without the reference

The caller currently gates descent, ordering and phasing on reference expressibility. Only the VCF
needs that -- the anchors and the mosaic do not -- and the sites it excludes are nested chains inside
non-reference alleles, which is precisely the population that matters for assembling complex loci.

## Sizing, measured first

`VG_CALL_NO_REF_NESTED=1` admits exactly the excluded population. chr20 ONT, vg `42e2f2ced`:

| | baseline | off-reference admitted |
|---|---|---|
| child calls | 18,561 | **29,657** (+60%) |
| skipped for no reference path | 9,122 | 0 |
| anchors | 636,484 | 657,942 |
| VCF records | 115,979 | 115,979 |
| **het sites in the phase chain** | **76,135** | **76,135** |
| reliable / chain breaks | 60,203 / 1,084 | 60,203 / 1,084 |

**8,329 new snarls**, and the split of them is the point:

| | count | share |
|---|---|---|
| heterozygous, two slots | **1,009** | 12.1% |
| split homozygous, two slots | 1,889 | 22.7% |
| one slot, no haplotype at all | 5,431 | 65.2% |

Two things follow. **Admitting them is purely additive** -- 0 of the 172,082 existing snarls changed a
slot. And **they are genotyped and anchored but never phased by reads**: the chain's numbers are
identical to the digit, so the 2,898 two-slot anchors carry the PANEL's frame, at the panel's 3.79%
switch rate rather than the reads' 0.38%, and the 5,431 one-slot ones carry no haplotype at all.

## Why they are excluded, in four places

1. **Descent.** A child the reference does not thread is skipped outright
   (`src/graph_caller.cpp:7985-7997`), because `get_ref_position` reaches `get_ref_interval`, whose
   `assert(start_steps.size() > 0 && end_steps.size() > 0)` fires. `no_ref_nested` admits it.
2. **The linkage key.** `LinkageCollector::record` takes `(contig, position)`, a reference coordinate.
   An off-reference chain is given **its parent's reference start** as an anchor -- "a place in the
   contig, not a coordinate for this snarl".
3. **The transition distance.** `site_gap` (`src/linkage_model.cpp:32-46`) returns `SIZE_MAX` when
   either site is `unpositioned`, so `rho = 1` and the HMM transition goes uniform: the chain
   forgets. Correct as a refusal, but it means such a site neither gives linkage nor receives it.
4. **The phase order.** `PhaseSite` sorts on `(phase_set, position, record_key)` and **carries no
   `unpositioned` flag at all**, so every off-reference child of one parent shares a position and
   ties break on `record_key`, a hash. Their order in the chain would be arbitrary. Today it does not
   arise because they never reach `phase_sites`; under this plan they would.

## The design

Exactly the shape proposed, mapped onto what exists:

- **(a) genotype the parent chain with haplotype linkage** -- already done, unchanged.
- **(b) recursively genotype each child chain crossed by one or both parent traversals, the crossing
  count setting haploid or diploid** -- already done by `child_ploidy`; only the reference gate is
  removed.
- **(c) order by snarl tree order** -- NEW, and the core of the work.
- **(d) phase over that ordering** -- NEW: the phaser currently orders by reference position.
- **(e) read assignments feed back into (a)** -- the re-genotyping loop already does this; it simply
  now covers the new sites.

So the genuinely new component is a **tree order and a tree distance**, replacing reference position
in both roles.

### The tree key

Every record gains a `tree_key`: the path from the root, as a sequence of offsets.

- Top level: `(contig rank, reference offset)` on a linear reference; `(chain id, offset)` on a gref.
- Each descent appends the child's **start offset along an alignment of the parent's two settled
  traversals**. The two traversals share endpoints, so a node-level alignment is well defined, and a
  child crossed by only one of them takes its offset in that alignment.

Sorting lexicographically on this key gives exactly the stated order: parent chain first, then
children by position along the parent's alignment, recursively. It is total, deterministic, and
reference-free.

### The tree distance

`site_gap` needs a length, not a coordinate. Between two siblings it is the **graph distance along
the parent traversal** between their anchors, in bases -- summing node lengths, which the traversal
already carries. Between a parent and its own child it is the offset into the parent. Across a
subtree boundary it is the parent-level distance plus the residual offsets.

This is the one quantity with real modelling risk: `switch_probability` is calibrated on reference
distances, and a traversal distance is not interchangeable with one. It must be checked rather than
assumed.

### DECIDED: the unstable, genotype-dependent order

The sibling offset comes from an alignment of the parent's two SETTLED traversals, recomputed each
re-genotyping round. It is therefore inside the loop -- order feeds the linkage HMM, which settles the
genotype, which picks the traversals, which define the order.

I raised the alternative (freeze the order on round 1, as the temper already is) because the loop
already fails to converge: a period-3 limit cycle at round 7, whose mechanism is a parent flipping and
its child leaving the chain. A genotype-dependent ordering adds a second discrete state change to the
same loop. **The call is to take the unstable order anyway** -- convergence is not the priority, and a
meaningful biological order is. Recorded here so a later convergence surprise is attributed to a known
decision rather than re-debugged from scratch.

### Computing the sibling offset without a general aligner

The two settled traversals share endpoints, and in a DAG the child chains visited by BOTH appear in
the SAME relative order in each. So no alignment algorithm is needed:

1. Walk each traversal and list the child snarls it enters, in order.
2. The children on both traversals are a common subsequence, in the same order in both -- these are
   the anchors.
3. Merge: emit the anchors in order, and between consecutive anchors emit the private children of
   each traversal.
4. The merged index is the offset.

**Ties are real and must break deterministically.** Between one pair of anchors, `trav_first` and
`trav_second` may each contribute private children -- a heterozygous insertion carrying its own
sub-variation on each allele, which is the complex-locus case this whole change exists for. So the
offset is the triple `(merged index, slot, start node id)`, never a bare integer.

## Phases, each with its own gate

**Phase 1 -- tree key, computed but unused.** Derive and store it; sort by it ONLY where it provably
agrees with reference order. *Gate: chr20 VCF and anchors byte-identical.*

**Phase 2 -- linkage layer ordered by tree key**, with `site_gap` taking a tree distance. Off-reference
chains still skipped, so the tree order must reduce exactly to reference order.
*Gate: byte-identical on chr20 and chr6. A difference here is a bug in the key, not a result.*

**Phase 3 -- `PhaseCall` carries the tree key and `unpositioned`; the phaser sorts on it.**
*Gate: byte-identical again -- still no off-reference sites in play.*

**Phase 4 -- admit off-reference chains.** Remove the descent gate; keep `VG_CALL_NO_REF_NESTED` as a
kill switch. VCF emission is untouched: those records still have no REF or POS and are still not
emitted. Anchors and mosaic gain them. *This is the first phase that may move anything.*

**Phase 5 -- re-measure the re-genotyping loop.** More nested sites in the chain changes the coupling
that already produces a period-3 limit cycle at round 7, and the default `--regeno-passes 2` only
avoids it by stopping early. Convergence must be re-measured, not assumed.

## Testing

**Byte-identity is the gate for phases 1-3**, on chr20 and chr6. Nothing should move until phase 4,
and anything that does is a defect in the key or the distance.

**Phase 4 is measured on the anchors, not the VCF**, because the VCF cannot express these sites:

| metric | baseline | target |
|---|---|---|
| snarls anchored | 172,082 | ~180,411 |
| het sites entering the phase chain | 76,135 | ~77,144 (+1,009) |
| two-slot anchors with READ-derived phase | 0 of the 2,898 new | all of them |
| phased-run N50, read-walkable | 2,798,104 bp | no worse |
| whatshap switch error, chr20 then chr6 held out | 0.3826% / 55 true switches | no worse |
| held-out cross-site agreement | 94.7102% | no worse |

**Switch error is the gate that matters**, and it must be read with its denominator: whatshap assesses
only sites it can place on the reference, so the new sites are invisible to it. A flat switch error
with 1,009 more hets phased is a real gain; a flat switch error is not by itself evidence the new
sites are phased CORRECTLY. The honest check for those is the run's own held-out cross-site
agreement, which does cover them.

**Fixtures.** `nest.gbz` has a nested chain under a deletion but the reference threads it. A new
fixture is needed with a chain inside a non-reference allele -- four requirements, or it silently
tests nothing: the parent must be heterozygous, the child must sit on the ALT traversal only, the
reference must not thread the child, and reads must span both.

**Regression risks to gate explicitly:** the tree distance feeding a switch probability calibrated on
reference distances; the limit cycle's period changing; and per-site memory, since `tree_key` is a
variable-length path on every record and chr20 retains 180,411 of them.

## What this does not fix

The phase chain would still be a sign-only pairwise cascade, still with no independent witness at a
mismapped locus. The 52 chr20 switches live in pericentromeric repeats and are not in this
population. This plan is about COVERAGE -- phasing sites that are currently unphased -- not about the
accuracy of sites already phased.
