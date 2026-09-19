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
