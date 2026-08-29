# Rebuilding nested calling around inference, not VCF

The deferred-calling machinery has grown by iteration until several mechanisms decide the same
facts and none of them can be checked against the others. This is the plan to replace it with the
smallest design that produces the same answer, organised around what is being inferred rather than
around what will be printed.

Authorised by the repository owner as a rebuild: large-scale removal and replacement is in scope,
and the measure of success is code removed, not accuracy gained.

## The algorithm, as the owner states it

For each top-level chain at ploidy 1 or 2:

1. Genotype the top-level chain as today, using linkage, giving chosen haplotype traversals (one or
   two) for each snarl in the chain.
2. For each parent snarl and its haplotype traversals: (a) symbolic allele comparison to decide the
   genotypes to REPORT, which is an emission concern only; (b) compute the copy count of each nested
   chain, which may exceed 2 -- count and log those.
3. Recursively genotype each nested chain at its copy count, returning to step 1.

## Settled decisions

- **Greedy, permanently.** A parent settles before its children and is never revisited. Child
  evidence never informs a parent. Today is already greedy in this direction, so nothing is lost.
- **A nested chain's phase is expressed in its parent's frame.** Strand identity composes up the
  tree. The owner calls this the essence of the inference.
- **The output is each strand as a sequence of recombinations of the GBWT haplotypes**, including
  recombinations inside nested chains. The mosaic is the product; VCF is a projection of it.
- **Windowing and phasing behaviour does not change.** Reimplement for clarity if it comes out
  cleaner, gated on byte-identity. No parallelism inside a chain's decode: independent windows put
  the two decodings in unrelated frames, and reconciling them after the fact is a vote that gets
  noisier as the panel grows -- it would manufacture a switch at every seam, which is the artifact
  the seam pin exists to prevent. Cutting at a saturated gap is exact but `rho` reaches 1 only
  asymptotically, so there are no such gaps to cut at.
- **Copy count above 2 is counted, not represented.** The hidden state is an ordered PAIR of panel
  haplotypes; triples would be 42,875 states at m=35, a third `rho`, and a different GL layout.
- **Both conditionings are testable from one binary**: a child conditioned on its parent's hard
  called pair, versus on the parent's posterior over pairs with the copy count still hard-called.

## What the measurements already settled

**Linkage between nested chains is worth nothing.** Decoding each nested chain alone against
decoding all of a parent's children together, on three contigs:

| | chr20 | chr6 | chr17 |
|---|---|---|---|
| ALL / SNV / Indel F1 | 0 | +1.8e-6 | -1.3e-5 |
| SV F1 | -3.3e-4 | 0 | +8.4e-4 |

Every movement is one call, and the signs disagree across contigs. chr17 carries 273 multi-crossing
children where the other two carry none, so the case with the best chance of mattering does not.

**Therefore sibling chains have no transition between them, and their relative order and spacing are
unobservable.** That retires inter-chain ordering, orientation and distance entirely.

**Read retention across the recursion already exists**, in the read source rather than the caller:
descent is synchronous on the fetching thread and the cache is keyed by node-ID window, so a child's
reads hit the parent's resident window. Going back to the reads after settling was measured at
**+57.8%** (22.4M against 14.2M) because a scattered second pass re-fetches whole windows to reach a
handful of sites.

**A chain-unit read query would break that.** A top-level chain's node-ID span routinely exceeds the
4,096-node window, and a straddling query bypasses the cache and caches nothing. Scoring stays
snarl-unit inside a chain-unit loop.

**"Score only what the reads touch" buys 3-5%.** 96.3% of nested chains have discriminating reads
(94.2% with off-reference chains admitted). Worth doing as a cheap skip at descent -- it avoids the
`panel_alleles` query too -- but it cannot drive the architecture.

## Distances

- **Top-level chain**: reference positions, unchanged. The reference is a real coordinate there.
- **Within a nested chain**: the length of the settled traversal through the preceding snarl. Chain
  snarls are adjacent by construction, so this is intrinsic and needs no reference. Per-strand
  distances genuinely differ here, because the two haplotypes take different-length paths -- unlike
  the between-chain case, where 90% were exactly the reference difference.
- **Between nested chains**: none. There is no transition.
- **Parent to child**: the child's offset along the parent's settled traversal, per strand, from
  `Entry::frame_offset`. Linkage should decay across a large parent snarl, so containment is not the
  same as zero separation. To be tested against a no-transition arm and a fixed-distance control,
  which separates "decay matters" from "any transition matters".

**The switch that exists is not the switch step 7 names.** `VG_LINKAGE_FRAME_GAPS` defaults to **0**
-- reference-position gaps -- with 1 = one scalar frame offset shared by both strands and 2 = per
strand. So the parent-to-child frame distance is OFF today, and the plan's "no transition" and
"fixed nominal" arms do not exist yet. Separately, the plan asks for the within-chain distance to be
the length of the settled traversal through the preceding snarl, while the code differences
`frame_offset` along the PARENT's traversal. Close, not the same; step 7 has to choose which.

## The plan

Steps V1-V3 are validations the plan itself needs before it can be executed.

| | step | gate |
|---|---|---|
| V1 | **DONE, and it refuted the figure.** 250 sites was said to retain 55% pair correlation, needing ~1,237 for 0.05. Measured on chr20: mean -log P(no switch) = 0.0121 per step, so **4.8% retained at 250 sites and 247 needed for 0.05**. The shipped margin is almost exactly right. (Per strand rather than per pair it is 22% retained and ~494 needed; the pair is what the window decodes, so 250 is correct on the measure that matches the model.) **Step 4 is therefore deleted.** | measured; instrument arm byte-identical |
| V2 | **DONE, settled by reading.** Descent already visits EVERY child of every snarl: a child no called allele reaches is not skipped, it sets `retain_only` and descends anyway, precisely because linkage may move the parent onto it. So collect's enumeration is already broad enough for greedy recursion. The one narrowing is the reference gate (12,486 chains on chr20), which is NOT what step 8 removes and
is answered separately at 8b. | none needed |
| V3 | Restate step 5's gate. "Zero disagreements" is unachievable on current code because the `ploidy == 2 && parent_trav >= 0` population is a known defect. The gate becomes zero AFTER that fix, making it a prerequisite rather than an aside. | plan text only |
| 0 | Instruments: pin declines, read-evidence split, copy-count histogram | published, no logic change |
| 1 | **DONE, 217 lines.** Retire the stage-14 frame instrumentation, THEN delete `frame_reversed`, `frame_end`, `frame_total`, `n_reads`. The three frame fields are not dead as claimed: `frame_end`/`frame_total` are read in the cross-parent branch of that instrumentation, into a `frame_gap` the next line discards. `frame_offset` STAYS -- it is the parent-to-child distance. | byte-identical |
| 2 | **DONE.** Each nested chain decoded alone. Retires the POOLING; the deletions it enables (`unpositioned`, the anchor, `chain_index`) land at 6c, because a nested chain still enters a contig chain that the grouping then discards. | matches the measured per-chain arm exactly on chr20, chr6 and chr17 |
| 2b | **DONE.** The diploid group key is the chain's own boundary pair from the graph, not the alignment column: 6,481 groups partitioned identically across three contigs, 37 differing, all of them chains the alignment could not identify and isolated as singletons. This is what actually retired the alignment column as a grouping key. It also makes `align_rank` **provably** inert in the diploid comparator rather than merely unused -- the key is `(parent, chain)`, so every entry in a group belongs to one chain, the rank is constant across it, and the branch can never fire. Only the per-strand pass still reads it. | byte-identical; F1 unchanged or better on three contigs |
| M2 | **DONE. INERT on all three contigs, byte-for-byte.** Base and arm produce the same VCF and the same mosaic on chr20 (115,343 records), chr6 (296,777) and chr17 (134,888) -- not a matched F1, the same file, six times over. The alignment of the parent's two settled traversals orders nothing anything reads. Original statement of the step: The per-strand haploid bucket is keyed on `(phase set, strand)`, so unlike a diploid group it spans parents AND spans chains: it is the one surviving consumer that orders ACROSS chains, and M1 did not measure it. `VG_LINKAGE_NO_ALIGN_ORDER` already produces the arm and, given 2b, now affects nothing else -- so this is three contig runs and no code change at all. Inert retires `align_rank`, `set_align_rank`, `symbolic_columns`, `chain_backward` and the descent-time alignment together. Caveat: the switch drops the rank and the chain index in one move, so a result that MOVES has to be split before either can be read. | accuracy, chr20 / chr6 / chr17 |
| 3 | **FOLDED INTO 6c.** Splitting it out standalone requires keeping traversals alive longer -- more state, not less -- and 6c restructures the same call sites. Doing it twice would be churn for a performance gain this rebuild is not for. | -- |
| ~~4~~ | ~~Expand the margin~~ -- **deleted**, V1 shows 250 is correct | -- |
| 5 | **DONE.** `relate()` landed as a CHECK. 19,979 derivations on chr20, 52,800 over three contigs, **zero** disagreements and zero unanswerable. | byte-identical |
| 6a | **DONE.** The strand pass consumes the derivation instead of reading the stored field. The step-5 check passed and was still insufficient: it validated the derivation from the parent's ENTRY, the substitution takes it from the parent's PHASE CALL, and the `ploidy == 2` guard before consulting `trav_second` lived only in the first. The gate caught it on the first run. | byte-identical |
| 6b | **DONE.** `Entry::parent_trav` deleted, with `set_parent_trav`, the descent-time computation and the parameter from `record`/`respecify`. Three positional argument shifts in the tests, all caught by them and none visible to chr20. | byte-identical |
| M3 | **DONE. INERT on all three contigs, byte-for-byte.** The chain index, and the backward flip that reverses it, order nothing anything reads either. Together with M2 this retired `symbolic_columns`, the barrier's alignment block, `set_align_rank` and three `Entry` fields -- 387 lines. | byte-identity, three contigs |
| M4 | **DONE, and the answer is NO.** Pooling haploid nested chains across a (phase set, strand) is worth two structural variants on chr20 and two on chr6, none lost anywhere; small variants are a wash. The pooling stays, and 6c is re-scoped accordingly -- see the section above. | accuracy, three contigs |
| SC | **DONE.** `record()`'s ten trailing parameters became one `SiteContext` with named initialisers. They were all bool, int and size_t -- which convert to one another silently -- so a call site that dropped one still compiled with every later argument shifted. That happened three times on this branch and the compiler caught none of them. | byte-identical |
| H | **DONE, ~90 lines.** Housekeeping, before 6c. Eight empty conditional bodies and one wholly dead loop remain in `linkage_model.cpp` where 6d removed the counters but left their conditions, and their comments describe instruments that no longer exist. ~40 lines. Worth clearing before the restructure carries them through it. | byte-identical |
| 6c | **RE-SCOPED by M4, not started.** The generation loop becomes a recursion for the DIPLOID groups only; `respecify`, `Entry::generation` and `max_generation` go with it; the per-strand haploid pass stays as a final contig-wide phase because M4 says it earns its keep. `Entry::ploidy` stays -- it is the arity of the stored likelihoods, not a cache. | byte-identical |
| 7 | **DONE, and every arm came back a wash.** The parent-to-child distance is unmeasurable -- three arms spanning no-link to perfect-link, and on SV both alternatives change sign between chr20 and chr6. The hard-clamp conditioning is worth one false positive across two contigs. Neither is adopted: the default already spends nothing on the distance, and the conditioning arm's value is the machinery it shows to be unnecessary rather than the switch itself. See the two sections above. | accuracy, chr20 and chr6 |
| 8 | **DONE, and it is NOT inert -- the earlier row had that wrong.** `chain_reported_inline` does return false without `atomize_blocks`, and `--atomize-blocks` is ON by default under `--read-likelihood`: the rule holds back 391 chains on chr20, every one of which is now genotyped, recorded and phased. The line is still suppressed, at the render hand-off, beside the off-reference population that already had that shape. An emission rule was deciding what got inferred. `VG_CALL_INLINE_SKIPS_DESCENT` restores the old behaviour so both arms come from one binary. | scored against the old arm on three contigs |
| 8b | ~~Remove the reference gate~~ -- **declined, and it was already measured.** Admitting off-reference chains (`VG_CALL_NO_REF_NESTED`) adds 12,486 chains on chr20, never improved recall in any arm, and changes the record set, so it cannot be gated on record identity the way step 8 can. The gate stays on by default and the flag stays as an arm. | -- |
| 9 | **DONE.** The mosaic already reported each strand as a sequence of segments on panel haplotypes, so a recombination -- including one inside a nested chain -- was already a segment boundary. What it could not say is WHICH, and the nested ones are what nested calling is for. `PhaseCall::depth` and two columns, `nested_sites` and `max_depth`, appended at the end of the row with the version bumped to 4 so a positional consumer reads it unchanged. | VCF byte-identical; mosaic identical once the two columns are stripped |

Steps 1-3, H and 6 should be byte-identical. With step 4 deleted and 8b declined, **only M2 and
step 7 can move numbers** -- M2 by measurement rather than by change, and step 7 through the two
conditioning arms and the three parent-to-child distance arms. Everything else in the rebuild is
gated on producing the identical answer.

**M2, M3 and M4 go before 6c**, and M4 is not an optimisation of the order -- it is a precondition.
6c says "replace the barrier with the recursion", and a recursion settles one chain at a time. The
per-strand haploid pass does not: it pools every haploid nested chain on one (phase set, strand) into
a single contig-long run, across parents and across chains. Until that pooling is shown to be worth
nothing -- as the diploid pooling was -- there is no recursion to write, only a recursion with a
whole-contig pass bolted to the end of it. M2 and M3 cost nothing and delete code; M4 decides whether
6c is a restructure or a redesign.

## Three ordering keys, three inert answers, and the frame goes with them

The comparators carried three keys an entry could LACK, and each needed an all-or-nothing decision
per group to stay a strict weak ordering. Each was measured on chr20, chr6 and chr17, and each came
back byte-for-byte inert:

| switch | what it drops | result |
|---|---|---|
| `VG_LINKAGE_NO_ALIGN_ORDER` | the alignment of the parent's two settled traversals | identical, x3 |
| `VG_LINKAGE_NO_CHAIN_ORDER` | the snarl's index within its chain, and the backward flip | identical, x3 |
| `VG_LINKAGE_NO_FRAME` | the offset along the parent's settled traversal | identical, x3 |

With all three gone both comparators are TOTAL on every entry, so the all-or-nothing machinery has
nothing left to decide and goes too. And the frame's other use -- forming a distance -- was already
covered: its every arm had been measured worth nothing, including the default, which never entered
the block at all.

**~850 lines removed across the two deletions**, and `site_gap` loses two thirds of itself. What
survives in it is the case that matters: a pair where either site is unpositioned gets a uniform
transition, because differencing an anchor against a real coordinate is not a distance.

The pattern across every measurement in this rebuild is worth stating plainly. **Nothing about WHERE
a nested chain sits is observable** -- not its order among its siblings, not its order within itself,
not its distance from its parent, not its distance from the previous chain. The only thing the decode
needs of a chain is its IDENTITY, so that it can be told apart from its siblings, and that comes free
from the graph.

## M4 says no, and that is what bounds 6c

The haploid analogue of step 2: decode each haploid nested chain alone instead of pooling every chain
on one (phase set, strand) into a contig-long run.

| | chr20 | chr6 | chr17 |
|---|---|---|---|
| ALL F1 | -1.5e-5 | 0 | +8.3e-6 |
| SNV | -6.8e-6 | 0 | +5.8e-6 |
| JointIndel | -4.2e-5 | 0 | +1.5e-5 |
| **SV** | **-1.1e-3 (TP 434 -> 432)** | **-5.8e-4 (TP 939 -> 937)** | 0 |

Small variants are a wash and the sign changes with the contig, as it did for every arm so far. SV
does not: two structural variants lost on chr20, two on chr6, none gained anywhere. The sign never
reverses.

**So the pooling stays.** And the reason it differs from the diploid case is visible in the run: 352
of chr20's 649 strand groups are groups of ONE, and 393 of chr17's 728. A haploid nested chain is
routinely a single site with no internal linkage at all, so pooling is the only thing that gives it a
neighbour; a diploid nested chain is bigger and already has its own structure. The two populations
are not symmetric and the same answer does not apply to both.

**6c has to be re-scoped.** The plan says "replace the barrier with the recursion", and a recursion
settles one chain at a time. That is available for the diploid nested groups -- each is already
decoded alone -- and it is NOT available for the haploid ones, whose decode is a contig-wide sweep
that must run after every parent is phased. So the restructure is:

- the generation loop over a flat pending list becomes a recursion over the tree, for the diploid
  groups only;
- `respecify` goes, because a recursion knows the copy count before it genotypes and can `record()`
  once at the right arity;
- `Entry::generation`, `max_generation` and `PendingRecord::generation` go with the loop;
- the per-strand pass stays exactly as it is, as a final phase, and `Entry::ploidy` stays because it
  describes the arity of the stored likelihoods rather than caching a derivable fact.

That is a smaller and less elegant 6c than the one written above, and it is the one the measurements
allow. It is a restructure of two large functions with a byte-identity gate, and it should be started
from a clean tree with the owner present rather than at the end of an unattended run.

## Step 7, chr20: the parent-to-child distance is worth nothing either

Three arms over the step from a parent to its nested children, spanning the whole range the model
can express -- no transition at all, the reference-position difference the default uses, and
containment read as zero separation, which is the strongest link there is:

| chr20 | uniform (no link) | reference-position gap (default) | minimal (zero separation) |
|---|---|---|---|
| ALL F1 | 0.972403 | **0.972410** | 0.972405 |
| SNV | 0.985156 | **0.985190** | 0.985183 |
| JointIndel | **0.928397** | 0.928314 | 0.928314 |
| SV | 0.532701 | 0.533020 | **0.533675** |
| records | 115,159 | 115,343 | 115,375 |

Every movement is a handful of calls -- 5 of 91,493 true positives on ALL -- and the winner changes
with the class. This is the shape M1 had: **the parent-to-child transition is not measurable.** The
parent still conditions the child, through the context message; what is worth nothing is expressing
containment as a DISTANCE on top of that message. Which retires the last surviving distance in the
nested tree: between chains there is none, and now above them there is none either.

**chr6 confirms it, by disagreeing.**

| chr6 | uniform (no link) | reference-position gap (default) | minimal (zero separation) |
|---|---|---|---|
| ALL F1 | 0.977469 | 0.977475 | **0.977481** |
| SNV | 0.988020 | 0.988036 | **0.988041** |
| JointIndel | **0.939489** | 0.939460 | 0.939469 |
| SV | **0.586195** | 0.584746 | 0.584096 |
| records | 296,688 | 296,777 | 296,807 |

On SV both alternatives change sign between the contigs: uniform is -3.2e-4 on chr20 and +1.4e-3 on
chr6, minimal is +6.6e-4 on chr20 and -6.5e-4 on chr6. ALL stays inside +/-7e-6 everywhere. Two
contigs, two arms, and the winner is different every time.

**So the parent-to-child distance is unmeasurable, and the default already spends nothing on it.**
Mode 0 lets `site_gap` fall back to the reference difference, which costs no frame, no measurement
and no branch -- adopting either arm would ADD code to express a distinction the data does not
support. The default stands.

What this does retire is `VG_LINKAGE_FRAME_GAPS` and the diploid frame-gap block behind it, about 80
lines that are off by default and whose every arm is now known to buy nothing: mode 1 and mode 2
measure a distance along the settled traversal for the same step these three arms just showed is not
a step worth measuring.

**The hard-parent conditioning arm did not run.** It reported "0 context messages replaced" at every
generation, so the delta was never substituted and the arm measured nothing. Re-run with the decline
reasons counted apart -- no PhaseCall, wrong message width, wildcard strand -- and with a guard the
first version lacked: `WILDCARD` is `(size_t)-1`, so `first * m + second` wraps and can land inside
the array, writing a delta at a pair the panel never chose.

## Step 8, scored: emission was deciding what got inferred, and it cost a structural variant

Under `--atomize-blocks` -- ON by default with `--read-likelihood` -- a chain whose variation an
enclosing difference block has already spelled out gets no record of its own. That rule was applied
by SKIPPING the descent, so the chain was never genotyped, never recorded, and could not inform its
own parent's phase. It now descends, and only the line is suppressed, at the render hand-off beside
the off-reference population that already had exactly this shape.

391 chains on chr20, 137 on chr6, 311 on chr17:

| | chr20 | chr6 | chr17 |
|---|---|---|---|
| ALL / SNV / Indel / Deletion | identical, same TP/FP/FN | identical | identical |
| SV | identical (434/445/322) | identical (939/724/609) | **+8.4e-4: TP 466 -> 467, FN 344 -> 343** |

One structural variant gained, nothing lost on any contig. Small, but one-directional, and it is the
change the design argues for independently: what a block already printed is a fact about the OUTPUT,
and it has no business deciding what the model gets to see.

`VG_CALL_INLINE_SKIPS_DESCENT` keeps the old behaviour, which is what lets everything else in the
run be gated on byte-identity while this one step deliberately moves.

## The conditioning arms: a delta is as good as the posterior

A child chain is conditioned on a message from its parent, and that message is by default the
parent's POSTERIOR over ordered haplotype pairs -- the normalised alpha*beta harvested from the
parent's own decode. `VG_LINKAGE_HARD_PARENT` replaces it with a delta at the pair the parent
SETTLED on, which is what the greedy story says out loud: a decided parent has no residual
uncertainty to pass down.

Clamped 1,744 to 2,306 messages a generation on chr20 and 2,990 to 3,418 on chr6, with nothing
declined:

| | chr20 | chr6 |
|---|---|---|
| ALL / SNV / Indel / Insertion / Deletion F1 | identical to six decimals, same TP/FP/FN | identical, same TP/FP/FN |
| SV | -3.3e-4 (one extra FP, 445 -> 446) | identical, 939/724/609 |

**One false positive across two contigs.** The parent's residual uncertainty carries essentially
nothing to its children, and the owner's greedy story is empirically true here rather than merely
convenient.

**What that licenses is the largest deletion still on the table**, and it is not the arm itself --
it is the machinery that exists to produce the message the arm replaces. `posteriors_with_context`,
the sparse mask threaded through `segment_posteriors`, the alpha/beta harvest and its normalisation,
`wanted_parents`, and `parent_context` itself, which holds 16.4 MB on chr20 at generation 0 alone.
All of it computes a distribution that a one-hot built from the parent's `PhaseCall` reproduces to
six decimals.

It is not taken here: it moves one call, it is a redesign of the decode's interface rather than a
deletion behind a gate, and it should be a decision made deliberately. But it is measured, and it is
where the next real simplification is.

## A gate that could not fail

Every `cmp` in the measurement scripts was `$SCRATCHPAD/cmp`, a compiled one-off probe from earlier
work, because those scripts prepend the scratchpad to PATH to reach a pinned `vg`. It ignores its
arguments and exits 0, so **every byte-identity check reported PASS unconditionally** -- including
for two arms that differ by 26 kB. Three results were stated as byte-identical before this surfaced.

All of them were re-checked with `/usr/bin/cmp` and all of them hold: M2 is inert on three contigs,
and both housekeeping gates pass. But one claim was flatly wrong -- "this switch reaches nothing"
about a switch that changes 184 records -- and the tell had been in every log for an hour: the probe
printed its own output, which was filtered out as leftover instrumentation instead of being
explained.

## Where this run ended

**657 lines net removed from `src/`** (293 added, 950 deleted), every one of them gated:

| gate | scope |
|---|---|
| chr20 / chr6 / chr17 VCF byte-identical | the alignment deletion, the frame deletion, SiteContext, the housekeeping, the `ref_ploidies` fix, six switches |
| mosaic identical once step 9's two columns are stripped | the same, plus step 9 |
| `test/t/18_vg_call.t`, 317 assertions, PASS, run alone | all of it |
| full `vg test`, 12,547,761 assertions | all of it |

Step 8 is the one deliberate behaviour change, scored separately: one structural variant gained on
chr17, nothing lost on any contig, with `VG_CALL_INLINE_SKIPS_DESCENT` keeping the old arm so
everything else could be gated on identity.

**What is left, in the order it should be taken:**

1. **The conditioning simplification.** A delta at the parent's settled pair reproduces its posterior
   to six decimals, which makes `posteriors_with_context`, the sparse mask through
   `segment_posteriors`, the alpha/beta harvest and `parent_context` unnecessary. It moves one call
   and it is a redesign of the decode's interface, so it is a decision rather than a gated step --
   but it is the largest deletion still available.
2. **6c**, re-scoped by M4: a recursion for the diploid groups, `respecify` and the generation loop
   with it, and the per-strand haploid sweep left standing.
3. `ref_offsets` uses `operator[]` at 19 read sites from worker threads, the same defect just fixed
   for `ref_ploidies`.

## Expected size

Roughly 2,600-3,200 lines removed against 900-1,200 added, plus 31 functions, ~97 struct fields, 18
named concepts and ~103 counters. The three largest blocks counted directly:
`run_deferred_descent` 715 lines, the per-strand nested pass 816, the grouping and its vetoes the
rest.

## What this is not

It is not a performance or accuracy project. Everything measured today that survived was a deletion:
inter-chain linkage is worthless, and the off-reference population costs slightly more than it gives.
The value here is the code that stops existing, and that should not be read as anything else.
