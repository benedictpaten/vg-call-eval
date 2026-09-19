# Audit of nested calling: what the strand asymmetry is and is not

Status: **the defect is NOT fixed.** Nothing was changed in vg; all instrumentation is reverted and
`src/` is at `b9a9c3653`.

## The anomaly being audited

A nested haploid site's strand is inherited from its parent and never checked. Measured against the
site's own reads (leave-one-out, restricted to reads that support the called allele):

| | sites | placements | agree |
|---|---|---|---|
| slot 0, parent het | 5,697 | 107,816 | **70.77%** |
| slot 1, parent het | 2,596 | 59,496 | **92.12%** |

Het sites, where the strand IS checkable, reach ~95%. The gap is 21 points and no mechanism has
been found for it.

## Exonerated, with data

| component | verdict |
|---|---|
| `phase_haploid_slot`'s three `return 0` fallbacks | **clean** -- 26 of 10,883 records, and the stamping is faithful: `real_strand 0 -> slot 0` and `1 -> slot 1` with zero exceptions |
| `relate_to_parent` | **clean** -- symmetric, and the 64-bit crossing-mask overflow is unreachable here (max traversal index 1 and 16) |
| `order_arbitrary` (the sorted-pair coin flip) | **not involved** -- 0 for every child and every parent |
| `frame_flipped` defaulting for an unfound parent | **clean** -- the lookup hits 100% of the time |
| last-writer-wins on `linkage_phased` | **consistent** -- `render_phases`, `phase_index`, the flip and the nested cascade all take the last entry |
| pass 2 failing to re-derive the strand | **it does re-derive** -- 2,466 of 5,586 strand-carrying keys change strand between passes |
| reference-vs-alt sequence attracting mismapped reads | **no** -- at the SAME carrying traversal (1), slot 0 scores 72.56% and slot 1 scores 95.72% |
| a sign error on a subpopulation | **no** -- the per-site distribution is a broad spread (slot 0: 2.4% fully inverted, 49.7% clean) rather than a spike at zero |

## What the audit did establish

**The resolution runs twice.** 344,417 PhaseCall appends over 173,036 keys is ~1.99 each: once
during re-genotyping, once at render. `linkage_phased` is never cleared, it accumulates, and every
consumer takes the last entry. That part is coherent by design.

**33,232 keys settle on a DIFFERENT pair between the two passes**, 1,576 of them going heterozygous
-> homozygous, and 1,147 changing ploidy. 6.2% of nested children have a parent whose pair changed
this way.

**The concrete defect: 1,442 children hold strand 1 under a parent whose final settled pair is
homozygous.** `nested_strand_of` cannot produce that -- with `trav_first == trav_second` its first
test fires and returns 0. So those strands index a pair that is not the one the record ships with.
2,358 children in total sit under a now-homozygous parent, where the chain is on both copies and no
strand is meaningful, and every one of them holds a real strand rather than the fallback.

This is distinct from the swap case the read-phasing cascade already handles: the cascade carries a
later **swap** of the parent's pair, and nothing carries a later **change** of it.

## Where to look next

The strand is re-derived in pass 2, so the stale values must arise *within* a pass -- a parent
revised after its own children's generation was resolved. `resolve_generation` builds `pinned_phase`
from the accumulated `phasing_out` at entry, so a revision that lands later in the same pass is not
seen by children already grouped. That is the ordering to check, and it is the only remaining path
by which a child can hold a strand its final parent cannot produce.

Note that none of this is visible to F1 or to switch error: whatshap never assesses a half-missing
record and the genotype does not move. 3,117 chr20 records carry these strands.

## A caution recorded

Two of these dumps came from runs that differ in `--anchors-out`, which implies off-reference
nesting and so descends into more nested chains. Counts from such runs are NOT comparable -- an
apparent 6,753/4,104 against 2,910/3,901 discrepancy was exactly this and nothing more.

## Round 2 of the chase: two more hypotheses killed, root cause still not found

**Hypothesis: the group-wide stamp.** `resolve_generation` derives the strand from
`kv.second.front()` -- one child -- and writes it to every member of the group. If members sat on
different parent strands they would all inherit the front one's. **Refuted: 0 of 20,068 children
are stamped with a strand differing from their own relation.** Groups are homogeneous in
`carrying` by construction (the group key includes the chain), so the front child is representative.
This also explains the `g_nest_both` counter reading 0 -- it is only reached when the front child
has `carrying == -2`, so it under-reports, but nothing is mis-stamped.

**Hypothesis: two sources for the parent's pair.** `carrying` comes from `relate(child,
entries[pidx])`, which reads the parent's **Entry** (the order the barrier settled), while
`nested_strand_of` compares it against **`pin->second.trav_first/trav_second`** -- the
**PhaseCall**, which `apply_read_phasing` rewrites in place when it swaps a pair. At derivation
**48.9% of diploid parents have those two orders disagreeing**, and the strand splits accordingly:
74.5% on strand 0 where the pair is swapped against 29.6% where it is not.

**Refuted, by building the fix and measuring it.** Deriving `carrying` from the PhaseCall's own pair
instead produced **byte-identical output**. `carrying` is a traversal VALUE, not an index into the
pair, so `relate_to_parent(mask, ta, tb)` and `relate_to_parent(mask, tb, ta)` return the same
thing -- swapping the order cannot change which traversal crosses the child. The 74.5/29.6 split is
therefore a *consequence*, not a cause: where the pair is swapped, `trav_first` is the other
traversal and the same child correctly gets the other index. The field is doing its job.

### Where the imbalance enters, precisely

Three stages in one run, keyed by record_key:

| stage | strand 0 | strand 1 | % on 0 |
|---|---|---|---|
| A. derived in `resolve_generation` | 5,427 | 6,314 | 46.2% |
| B. after the read-phasing cascade | 6,753 | 4,104 | **62.2%** |
| C. stamped on the anchors | 6,779 | 4,104 | 62.3% |

B -> C is exact (0 mismatches). The cascade is also exact: **`B == A XOR parent_flipped` holds for
all 11,709 children, 0 violations.** So the shift is not an error in the cascade -- it is the
correlation it acts on:

| at derivation | strand 0 | strand 1 |
|---|---|---|
| children whose parent will flip | 1,889 | 3,732 (33.6% on 0) |
| children whose parent will not | 3,515 | 2,573 (57.7% on 0) |

and that correlation is itself the correct consequence of the swap state, per the refutation above.

### What is still unexplained

Every step now reconciles, yet the output is still wrong: **slot 0 agrees with its own reads 72.55%
and slot 1 92.80%** (supporting reads only). Nothing found so far accounts for a 20-point gap.

**The sharpest untested hypothesis.** `agree` compares `lo > 0` -- the read's strand in the CHAIN's
frame -- against `slot == 0` -- the child's position in its PARENT's pair. Those are commensurable
only if the parent's pair order is the chain frame, which holds only where the parent is a **phase
site** and so received an `o[]`. A parent excluded from `phase_sites` (homozygous, ploidy != 2, or
otherwise) keeps the panel's order, which is related to the chain frame only by accident. Proxies
for this (parent ploidy, parent homozygosity in the FINAL call) do not separate it cleanly, so the
test has to be the real one: condition the agreement on whether the parent's record_key is in
`phase_sites`. That is the next thing to run.

**Status: NOT FIXED.** `src/` is at `b9a9c3653`; the one fix attempted was measured, found inert,
and reverted.
