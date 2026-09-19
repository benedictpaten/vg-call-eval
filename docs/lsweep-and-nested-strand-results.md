# The L-sweep, and what verifies a nested haploid strand

Two questions, both run on chr20 with chr6 held out.

## 1. Can links that SPAN a junction see what adjacent links cannot?

The chain has no redundancy: stage 1 is `o[m+1] = o[m] ^ (d[m] < 0)` and stage 2 orients each block
against the previous one, so one wrong sign flips everything downstream. That is why 13 parity
changes make 24.72 Mb of mis-phased sequence on chr20.

A single-site bypass had already come back null, but it skips a few hundred bases and so stays
inside the ~2 kb mis-specified region, using the same mis-aligned reads. Reads span a median 21 kb,
so a link from a clean site ~10 kb before to a clean site ~10 kb after should be carried by reads
whose alleles at BOTH endpoints are fine. **Nothing in the algorithm computes it.**

Instrumented `phase_link(rel[m+1-L], rel[m+L])` for L in {2, 5, 10, 20, 50}.

**A correction that mattered.** The first pass compared each span against the parity implied by the
signs of `d`. That is the wrong baseline wherever stage 2 relinked -- there the orientation came
from the K x K vote, not from `d` -- so it silently excluded every relink event, which is the
population the whole idea targets. Re-run against the settled `o[]` dumped after stage 3:

| events with >= 1 disagreeing spanning link | chr20 | chr6 (hold-out) |
|---|---|---|
| **relink events** | **2 of 4** | **2 of 4** |
| cascade events | 1 of 9 | 4 of 15 |
| all events | 3 of 13 | 6 of 19 |
| junctions flagged | 769 / 61,012 (1.26%) | 1,078 / 162,518 (0.66%) |
| lift | 18x | 48x |

**The 50% recall on relinks replicates exactly.** And the split is mechanistically right: a spanning
link helps where stage 1 had no evidence, and does not where the reads are uniformly wrong.

At the cascade errors the null is thorough. Filtered to spans enclosing an ODD number of erroneous
junctions (an even number cancels and the test is vacuous -- these junctions cluster, so this
matters), chr20 gives **0 disagreements out of 6-8 informative tests at every L**, including L = 10
(median span 11.6 kb, 30 shared reads) and L = 20 (25.9 kb, 21 reads). Nor is anything else
anomalous there: shared reads, evidence per read, magnitude and span all sit on the background.

So at a cascade error the reads are wrong over spans far longer than the mis-specified region -- a
read crossing it is mis-assigned along its whole length. No link-based scheme reaches those.

**The actionable consequence**: stage 2 already looks `--phase-relink 10` sites each side, and these
spanning links reach 50. Widening `--phase-relink` is the existing-parameter version of the idea,
and it is measured in mis-phased bases in the sweep beside this document.

## 2. Does a nested haploid site's inherited strand agree with its own reads?

A nested haploid site is excluded from the phasing chain -- `phase_sites` skips `ploidy != 2` -- so
its strand comes only from `nested_strand_of` against the parent's settled pair. Parity with the
parent is guaranteed by construction. Correctness is never checked, and no gate can see it:
whatshap never assesses a half-missing record and the genotype does not move.

Checked against each site's own reads, held out by `read_strand_log_odds`' leave-one-out:

| | sites | opinionated placements | agree |
|---|---|---|---|
| chr20 | 10,883 | 245,291 | **74.37%** |
| chr6 (hold-out) | 8,768 | 216,160 | **84.92%** |

Against the ~95% the het self-check reaches where the strand IS checkable. And on chr20 the two
slots are wildly different:

| slot | chr20 | chr6 |
|---|---|---|
| 0 | **65.60%** (160,818) | 84.44% (82,961) |
| 1 | **91.08%** (84,473) | 85.23% (133,199) |

chr20's slot 1 is at the het ceiling while its slot 0 is 25 points worse; chr6 is balanced. 12.8% of
chr20's sites have a MAJORITY of their reads disagreeing.

**The likely mechanism.** `phase_haploid_slot` has three `return 0` fallbacks -- no phase entry,
ploidy != 1 or nested_strand < 0, and a PhaseCall about a different genotype. It returns 0 for
"strand 0" and for "cannot tell", and nothing downstream can distinguish them, so every unanswerable
site lands in slot 0 looking like a haplotype claim. A run distinguishing a real nested strand from
the fallback is pending.

The earlier count asymmetry has the same shape: chr20 emits 2,105 `x|.` against 1,012 `.|x` while
top-level phased hets are 49.9/50.1. From the cascade dump, 11,657 nested haploid children split
55.0/45.0 after the cascade, and the imbalance lives entirely in the flipped-parent group (59.1% on
strand 0 against 48.5% for unflipped parents) -- the cascade itself mirrors correctly, so the skew
is in the pre-cascade assignment.

## Update: the fallback is NOT the explanation

Re-run distinguishing a real nested strand from `phase_haploid_slot`'s fallback:

| class | sites | placements | agree |
|---|---|---|---|
| real strand 0 | 6,751 | 160,782 | **65.61%** |
| real strand 1 | 4,104 | 84,473 | **91.08%** |
| fallback (no strand) -> slot 0 | 26 | 36 | 0.00% |

Only 26 sites take the fallback, so it explains nothing; restricted to real strands the figure is
74.38%, unchanged. **The asymmetry is in genuine strand assignments.**

Two further facts, from the cascade dump (11,657 children):

- **The skew grows with nesting depth**: 52.3% on strand 0 at generation 1, 58.7% at 2, 66.6% at 3,
  65.4% at 4, 71.4% at 5. A cascade that were symmetric would not drift with depth.
- Children under a **ploidy-1 parent** sit at 64.7% on strand 0 against 53.6% under a ploidy-2
  parent. The ploidy-1 path inherits the parent's own strand, so a deep chain inherits through
  several levels.
- Every parent is in `phase_index`, so `frame_flipped` never silently defaults for a child whose
  parent was flipped. That hypothesis is out.

A candidate remains unconfirmed: `nested_strand_of`'s diploid branch tests
`parent_trav_first == carrying` before `parent_trav_second == carrying`, so a chain carried by BOTH
parent strands returns 0 rather than -1 -- and `g_nest_both`, the counter that exists for exactly
that case, reports 0 because the first test fires first. For a ploidy-1 child this should be
unreachable (a chain crossed by both parent copies ought to be genotyped at ploidy 2), so it does
not obviously explain the numbers, and it is recorded as the next thing to check rather than as a
diagnosis.

**None of this is visible to any gate in use**: whatshap never assesses a half-missing record and
the genotype does not move, so neither switch error nor F1 can see it. 3,117 chr20 records carry
these strands.

## Update 2: a real defect found, and a residual that is still open

### Refuted, with data

| hypothesis | verdict |
|---|---|
| `phase_haploid_slot`'s three `return 0` fallbacks pool unanswerable sites into slot 0 | **no** -- 26 of 10,883 sites |
| the parent's (or child's) `order_arbitrary` coin flip | **no** -- 0 for every child and every parent |
| a homozygous parent, where `nested_strand_of`'s first test wins | **no** -- those score *better* (82.4% vs 70.8%) |
| `frame_flipped` silently defaulting for a parent it cannot find | **no** -- the lookup hits 100% of the time |
| the 64-bit crossing mask overflowing on `tb` | **no** -- max traversal index is 1 and 16 |
| my own measure including reads from the other haplotype | **partly** -- restricting to reads that support the called allele lifts 74.4% -> 79.9%, but the slot gap survives |
| reference-vs-alt sequence attracting mismapped reads | **no** -- at the SAME carrying traversal (1), slot 0 scores 72.56% and slot 1 scores 95.72% |

### The defect

**1,442 children sit under a parent whose settled pair is homozygous, and hold strand 1.**
`nested_strand_of` cannot produce that: with `trav_first == trav_second` the first test fires and
returns 0. So those strands were derived against a parent pair that was heterozygous **at the time
the barrier derived them** and is homozygous in the record that ships.

2,358 children in total sit under a now-homozygous parent -- a parent that carries the chain on
BOTH copies, where no strand is meaningful at all -- and every one of them holds a real strand, not
the fallback.

`nested_strand` is derived once, in the barrier, against the parent's pair at that moment. The
read-phasing cascade carries a later **swap** of that pair. Nothing carries a later **change** of it.

### The residual, still unexplained

The defect above does not account for the main asymmetry. Restricted to heterozygous parents, where
the rule is well defined:

| | sites | placements | agree |
|---|---|---|---|
| slot 0, parent het | 5,697 | 107,816 | **70.77%** |
| slot 1, parent het | 2,596 | 59,496 | **92.12%** |

A 21-point gap with no mechanism yet. The per-site distribution is a broad spread, not a spike at
zero (slot 0: 2.4% fully inverted, 49.7% clean; slot 1: 0.5% and 81.4%), so this is not a sign error
on a subpopulation -- those sites' reads are genuinely mixed. Next thing to look at is how
`carrying` is computed for the ploidy-1 group in `resolve_generation`, since `relate_to_parent`
itself is symmetric and clean.

None of this is visible to F1 or switch error.
