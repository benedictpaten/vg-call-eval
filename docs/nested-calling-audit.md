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

## Round 3: the error is proven real and quantified. The cause is still not found.

### Proven real

Restricting the check to **confident** reads -- |lo| >= 8, which agree with the allele partition at
het sites 95.9% of the time -- separates a wrong label from a noisy one:

| slot | all supporting reads | confident reads only |
|---|---|---|
| 1 | 92.80% | **95.05%** -- at the ceiling |
| 0 | 72.55% | **71.45%** -- moves the wrong way |

At slot-0 sites **28.5% of confident reads say the opposite strand**. Those reads are right ~96% of
the time. This is not a measurement artefact and not weak evidence: the label is wrong.

### Quantified

Treating slot 0 as a mixture of correct assignments (at slot 1's 95.05% ceiling) and random ones:

| population | observed | implied random fraction |
|---|---|---|
| slot 1 | 95.05% | **0%** |
| slot 0, parent flipped | 79.21% | 35.2% |
| slot 0, parent not flipped | 63.13% | **70.9%** |
| slot 0 overall | 71.45% | **52.4%** |

chr20 emits 6,779 slot-0 records, so roughly **3,500 nested haploid records carry a haplotype label
with no information in it** -- written identically to the 4,104 slot-1 records that are ~95%
correct. Invisible to F1 and to whatshap, which never assesses a half-missing record.

### Sixteen hypotheses eliminated with data

`phase_haploid_slot`'s fallbacks (26 sites) · `order_arbitrary` (0 everywhere) · homozygous parents
(score *better*) · `frame_flipped` missing a parent (hits 100%) · the 64-bit crossing-mask overflow
(max traversal index 16, properly tested) · reads from the other haplotype (restricting to
supporting reads leaves the gap) · reference-vs-alt (survives conditioning: 79.92/90.73 on
reference, 71.92/92.95 off it) · the parent being a phase site (79.38 vs 82.82, gap in both) ·
the group-wide stamp (0 of 20,068 mis-stamped) · two sources for the parent's pair (fix built,
byte-identical, because `carrying` is a VALUE and so order-independent) · a global sign bias in
`lo` (50.39/49.61) · my tie-break admitting flat reads (moves slot 0 by 1.2 points) · the called
allele · weak reads · the parent's flip status (slot 1 is ~95% either way) · whether the parent's
order was ever `decided` (all 78,176 phase sites are decided).

The cascade is exact: **`B == A XOR parent_flipped` holds for all 11,709 children, 0 violations**,
and derivation is balanced at 46.2% on strand 0. Every step reconciles individually.

### The lead I would follow next

`nested_strand_of` returns 0 on its FIRST match and 1 on the second, so slot 0 is the bucket that
collects any case where `carrying` wrongly equals `trav_first`. `carrying` comes from
`child.parent_crossing`, a bitmask where **bit i means "the parent's candidate traversal i crosses
this child"**, built by `child_crossing_mask(pr.travs, ...)` during the parent's descent. It is
indexed by `traversal_of(trav_arena, parent.trav_offset, parent.num_alleles, parent.final_i)`.

**Those two index spaces have never been verified to be the same list.** The mask indexes
`pr.travs` as it stood at descent; `trav_arena` holds what the linkage layer recorded. If the
parent's candidate set is re-scored between descent and settlement -- and re-genotyping does
re-score sites -- the bit positions and the arena positions need not correspond, and a mis-indexed
bit would make `first` true for a child that does not cross `trav_first`, which lands it in slot 0
with a meaningless strand. That matches the signature exactly: slot 1 clean, slot 0 a 50/50 mixture.

The test: at render time, independently recompute `crossings_of_child(travs[trav_first], child)`
for a sample of slot-0 children and compare it against the mask bit. A disagreement confirms it.

**Status: NOT FIXED.** `src/` is at `b9a9c3653`.

## Round 4: it is a PLOIDY bug, not a strand bug

### The discriminator

A wrong **label** gives a clean read set with the wrong sign. A wrong **ploidy** gives a 50/50 read
set where no strand exists at all. Confident reads only (|lo| >= 8):

| slot | sites | near-50/50 (0.35-0.65) | clean (>=0.9 or <=0.1) |
|---|---|---|---|
| 0 | 5,099 | **767 (15.0%)** | 2,901 (56.9%) |
| 1 | 2,393 | **69 (2.9%)** | 2,008 (83.9%) |

**Slot 0 carries 5x the rate of 50/50 sites**, and only 2.4% of slot-0 sites are cleanly inverted.
So the excess is not mislabelled haplotypes -- it is sites called HAPLOID whose reads come from
BOTH parent alleles. For those the site should be ploidy 2 and there is no strand to name; the
`a|.` it emits is a haplotype claim about a locus that has two.

That is why every strand hypothesis failed: **the strand derivation is correct.** Cross-tabulating
the crossing bits against the assigned strand accounts for all 11,700 records exactly -- 5,429 +
5,134 under diploid parents, and 1,137 under haploid parents correctly inheriting
`parent_nested_strand` rather than an allele index.

### Where it comes from, and what is still missing

`copies = bit(trav_first) + bit(trav_second)` over `child.parent_crossing`. A bit that reads 0
spuriously undercounts copies, makes a diploid child haploid, and hands it the one traversal whose
bit survived. The mask is built once at descent from `pr.travs`; the settled pair is read through
`trav_arena`, and a re-score can add an allele to that space.

**Not yet isolated.** Testing whether a settled traversal lies beyond the mask's highest set bit
found 1,240 such pairs (12.2%), but they skew to strand 1 (1,053 against 187) -- the opposite of the
slot-0 excess. So that is not the mechanism either.

### The fix does not have to wait for the mechanism

The condition is self-detectable at run time and needs no truth: **a site called haploid whose
confident supporting reads split near 50/50 is not haploid.** On chr20 that is 767 of 5,099
slot-0 sites and 69 of 2,393 slot-1 sites. The safe action is to refuse the strand there --
`nested_strand = -1`, no `a|.` -- which turns a false haplotype claim into an honest absence.
That is strictly better for assembly than the present behaviour, where ~3,500 chr20 records carry a
label written identically to the ~95%-correct ones.

Fixing `copies` properly is better still, and the measurement above is the gate for it: the 50/50
population should go to zero.

**Status: root cause narrowed to the copies/ploidy determination for nested children; the specific
mask defect is NOT yet isolated. `src/` is at `b9a9c3653`, unmodified.**
