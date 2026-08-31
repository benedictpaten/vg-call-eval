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
| V3 | **DONE.** Step 5's gate was restated: zero disagreements AFTER the `ploidy == 2 && parent_trav >= 0` fix, making it a prerequisite rather than an aside. It then read zero on 52,800 derivations. | plan text only |
| 0 | **DONE, then retired at 6d.** Pin declines, the read-evidence split and the copy-count histogram were published, answered their questions, and came out again with ~174 lines. | published, no logic change |
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
| 6c | **TWO INCREMENTS DONE, both byte-identical on three contigs; the rest NOT started.** Done: (1) at depth > 0 the contig runs are not built and every live site is grouped with its parent straight from the entries -- the veto that walked those runs had never fired on any contig, and its scope was wrong anyway, so a site that cannot be grouped is now decoded alone rather than chained to sites under unrelated parents; (2) the by-contig collection and its per-contig sort no longer run below the top level at all -- `deferred_nested` is gathered directly and sorted on the same key. **Not done: `respecify` and the generation loop.** `respecify` (87 lines) is retract-plus-record in principle, but it PRESERVES `explained_share` and five other fields that `record` would reset, and both take the collector's mutex -- so a naive unification deadlocks or silently changes 2,378 chains' `explained_share`. The loop is thin now that (1) and (2) have landed: converting it to a recursion changes decode ORDER, not content, and buys shape rather than anything measurable. | byte-identical, three contigs, plus TAP |
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

## The alpha/beta harvest is gone, and it cost one false positive

A nested chain is conditioned on a message from its parent. That message was the parent's POSTERIOR
over ordered haplotype pairs, harvested from its own decode by an alpha/beta pass with a sparse mask
and held at m*m doubles per parent -- 16.4 MB on chr20 at generation 0, 28.1 MB on chr6. It is now a
DELTA at the pair the parent settled on, built where it is consumed.

| | chr20 | chr6 | chr17 |
|---|---|---|---|
| every small-variant class | identical to 6 dp | identical | identical |
| SV | 0.533020 -> **0.532692** (FP 445 -> 446) | identical | identical |

**226 lines removed** -- `parent_context`, `wanted_parents` and their masks, both harvest blocks,
`posteriors_with_context`, and `alpha_out`/`beta_out`/`harvest_mask` from `window_posteriors` and
`segment_posteriors`.

One structural-variant false positive on one contig, nothing on the other two. It is the same
magnitude every arm here has moved, but unlike the others it never reverses sign -- there is no
contig where it gains.

**No switch, and that was decided rather than overlooked.** Every other output-moving change on this
branch keeps both arms in one binary -- `VG_CALL_INLINE_SKIPS_DESCENT`, `VG_LINKAGE_PER_CHAIN_STRAND`,
`VG_CALL_NO_REF_NESTED`. This one does not: the posterior path is gone, `git revert` of it conflicts
with the 6c increments that followed, and restoring it behind a flag was offered and declined. The
whole value of the change is the 226 lines, and a flag keeps them. Simpler wins; the cost is on the
record above, and the way back is this document plus commit `e7f7a27f9`.

**A use-after-free flattered the first attempt, and that is the part worth remembering.** The deltas
were held in a deque declared inside the loop that filled it, while the pointer vector indexing them
was read after that scope closed. It never crashed. On chr6 the freed memory was intact and the
answer correct; on chr20 it was not, and the arm appeared to cost NOTHING. The decision to delete 226
lines was very nearly made on that number. No gate could have caught it -- freed-but-intact memory
produces plausible output -- and what caught it was reading the scope while planning the next step.
When an arm looks better than expected, that is a reason to re-read it.

## 6c, first increment

At depth > 0 every live site is grouped with its parent and those groups replace the contig chains
wholesale, so the runs built at each generation were scaffolding for a fallback never once taken: 0
sites with no parent key, 0 whose parent has no live entry, 0 on a ploidy difference, on all three
contigs. They are not built there any more, and the group loop reads the entries directly.

The veto changes scope with them. It rejected a WHOLE contig run if any one site could not be
grouped; a site that cannot be grouped is now decoded alone instead -- which preserves what the veto
protected without chaining that site to sites under unrelated parents.

## Whole genome, 24 contigs: confirmed

Two fresh arms, `faedeb9e4` (this session's starting point) against the current tree, same inputs,
same scheduler. A fresh baseline rather than the `wgs-tt` arm on disk, which predates the session and
would have credited today with two days' work.

| autosomes, chr1-22 | base | final | |
|---|---|---|---|
| ALL F1 | 0.972941 | **0.972942** | +8.2e-7 |
| SNV | 0.984909 | 0.984909 | -1.5e-7 |
| JointIndel | 0.927538 | **0.927542** | +4.4e-6 |
| Insertion | 0.917687 | **0.917690** | +3.4e-6 |
| Deletion | 0.940780 | **0.940785** | +5.4e-6 |
| SV >=50 bp | **0.562112** | 0.561987 | -1.25e-4 |

All 24 contigs give the same picture: ALL 0.970507 -> 0.970508, SV 0.559939 -> 0.559817.

Small variants come out very slightly ahead -- 7 fewer false positives across the autosomes, 4.04 M
true positives unchanged. SV comes out slightly behind: 14,066 -> 14,064 true positives and 12,511 ->
12,517 false positives, which is **8 calls out of about 27,000**.

**And it answers the question the three-contig work left open.** chr20's single SV false positive is
NOT representative -- it is a coin flip. Per contig:

- **11 of 24 contigs are identical on every metric.**
- 8 contigs move on SV, and the signs disagree: chr19 +1.2e-3, chr17 +8.4e-4 and chr5 +3.7e-4 gain,
  while chr3 -1.8e-3, chr8 -7.3e-4, chr1 -5.8e-4, chr21 -5.2e-4, chr7 -4.3e-4 and chr20 -3.3e-4 lose.

That is the same shape every arm in this rebuild produced: a handful of calls, moving in both
directions, with the aggregate landing wherever the individual contigs happen to sum. The rebuild is
accuracy-neutral genome-wide, which is what it was supposed to be.

## What it cost, measured

Baseline `faedeb9e4` against the final tree, both pinned, run alone, and each verified byte-identical
to the arm it was scored as. chr20 twice per arm because peak RSS on this workload is noisy.

**Wall clock: no resolvable change.**

| | base | final | |
|---|---|---|---|
| chr20 | 176.3, 167.4 -> **171.9** | 171.3, 172.5 -> **171.9** | +0.0% |
| chr6 | 335.2 | 337.5 | +0.7% |

The base's own two chr20 runs differ by 8.9 s, so nothing here clears its own noise. Step 8 isolated
(`VG_CALL_INLINE_SKIPS_DESCENT=1`, 172.5 s) is indistinguishable as well.

**Peak RSS: unreadable.** chr20 base 4.13 and 3.49 GB across two runs of ONE binary; final 3.65 /
3.94 / 3.84. chr6 9.01 against 9.12. The spread within an arm exceeds any delta, exactly as the
deferred-descent work found when it tried to size retention this way.

**What is exact, because it counts objects instead of measuring them:**

| | chr20 | chr6 |
|---|---|---|
| linkage layer | 48.9 -> **44.0 MB** | 108.2 -> **97.0 MB** |
| per site | 233.9 -> **209.9 B** (-24.0) | 232.4 -> **208.4 B** (-24.0) |
| retained records | 413.3 -> 410.5 MB | 718.9 -> 711.7 MB |
| the linkage pass itself | 8.47 -> **8.02 s** (-5.3%) | 18.35 -> **17.80 s** (-3.0%) |

`Entry` lost exactly 24 bytes -- `align_rank`, `chain_index`, `chain_backward`, `frame_offset[2]`
and their padding -- and `PendingRecord` 16. The identical -24.0 B/site on two contigs is the
deletion appearing as arithmetic rather than as a measurement. Roughly 8 MB on chr20 and 18 MB on
chr6, against 3.8 and 9.1 GB peaks: **0.2%**.

**So this bought no speed and no meaningful memory, and the case for it was never that.** The
linkage pass did get 3-5% faster -- the per-parent O(|h1| x |h2|) alignment and the traversal walks
for the frame were real work -- but that pass is 5% of the run, so 0.45 s of 172. Read scoring
dominates and none of this touched it.

One number is worth keeping for its own sake. Step 8 descends into 391 extra chains on chr20 and
fetches **4,111 FEWER reads** (14,219,698 -> 14,215,587, cache hit 99% either way). The extra
genotyping is free because those windows are already resident, which is the claim the whole
collect-then-settle design rests on -- observed here rather than asserted.

## Accuracy, before and against after

Every small-variant class, all three contigs: identical to six decimals with identical TP/FP/FN.

| chr20 | F1 | precision | recall | TP | FP | FN |
|---|---|---|---|---|---|---|
| ALL | 0.972410 | 0.978673 | 0.966227 | 91,493 | 2,003 | 3,198 |
| SNV | 0.985190 | 0.995048 | 0.975526 | 73,181 | 355 | 1,836 |
| JointIndel | 0.928314 | 0.925869 | 0.930772 | 18,312 | 1,648 | 1,362 |
| Insertion | 0.917328 | 0.914699 | 0.919971 | 8,932 | 896 | 777 |
| Deletion | 0.938743 | 0.936206 | 0.941295 | 9,380 | 687 | 585 |
| SV >=50 bp | 0.533020 | 0.493743 | 0.579085 | 434 | 445 | 322 |

| chr6 | F1 | precision | recall | TP | FP | FN |
|---|---|---|---|---|---|---|
| ALL | 0.977475 | 0.982903 | 0.972107 | 272,156 | 4,742 | 7,809 |
| SNV | 0.988036 | 0.996188 | 0.980017 | 219,117 | 815 | 4,468 |
| JointIndel | 0.939460 | 0.938183 | 0.940741 | 53,039 | 3,927 | 3,341 |
| Insertion | 0.928627 | 0.926578 | 0.930686 | 26,129 | 2,207 | 1,946 |
| Deletion | 0.949913 | 0.949113 | 0.950715 | 26,910 | 1,543 | 1,395 |
| SV >=50 bp | 0.584746 | 0.564642 | 0.606335 | 939 | 724 | 609 |

The one moving cell in the whole matrix is chr17 SV: 0.532955 -> **0.533795**, TP 466 -> 467,
FN 344 -> 343, from step 8.

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

**What is left.** Nothing in the table. All three items below landed; the whole-genome confirmation
is running.

**`respecify` is gone**, 126 lines, byte-identical. It is retract-plus-record with the six fields it
silently preserved passed explicitly. Two things had to be got right for that to be equivalent:
`explained_share` is captured BEFORE the ploidy swap releases the original CallInfo (the alternate
does not copy it), and a chain that had no entry still gets the hard 1.0 the old fallback passed
rather than its real share -- the two old paths disagreed, the layer hands that value to
`apply_linkage_quality` where it discounts GQ, and it moved 5 chr20 records. Invisible to F1, which
does not read quality fields; caught by byte-identity.

**`ref_offsets` and `ref_ploidies` read through `ref_offset_of` / `ref_ploidy_of`.** Eighteen reads
used `operator[]`, which inserts on a miss, from OpenMP worker threads.

**The generation loop stays a loop, and now says why.** A recursion is not available: each group's
decode depends only on its parent's settled state, so level-order and depth-first would agree --
except that the per-strand haploid pass pools a whole contig's nested haploid sites on one strand,
which M4 measured worth two structural variants on chr20 and two on chr6. Depth-first would hand it
one subtree at a time. What the loop lost instead is the duplicated `has_entry` query and a counter
still named after `respecify`.

The conditioning simplification that used to head this list is **done** -- see the section above.

## Expected size

Roughly 2,600-3,200 lines removed against 900-1,200 added, plus 31 functions, ~97 struct fields, 18
named concepts and ~103 counters. The three largest blocks counted directly:
`run_deferred_descent` 715 lines, the per-strand nested pass 816, the grouping and its vetoes the
rest.

## What this is not

It is not a performance or accuracy project. Everything measured today that survived was a deletion:
inter-chain linkage is worthless, and the off-reference population costs slightly more than it gives.
The value here is the code that stops existing, and that should not be read as anything else.

## Design pass: one decode path for both nested populations

Authorised after the genome-wide measurement showed the two treatments are indistinguishable
(ALL +6e-8, SV -3.7e-4; 37 truth SVs lost against 30 gained, P = 0.46). The choice is therefore
free, and it is being made on simplicity.

### Target state

`resolve_generation` groups every live nested site of the current generation by (parent, chain),
**regardless of ploidy**. A group's ploidy is its sites'. Its context is a delta built from the
parent's PhaseCall -- over PAIRS for a diploid group, over SINGLE haplotypes for a haploid one -- and
the existing chain loop decodes both, because it already dispatches on `chain_ploidy`:

```
if (chain_ploidy == 1)      haploid_posteriors(sites)              <- gains the context
else if (context != null)   segment_posteriors(sites, ..., context)
else                        posteriors(sites)
```

One pointer serves both messages: `haploid_posteriors` accepts one of size `m`, `window_posteriors`
one of size `m*m`, and each already ignores a wrongly sized message. No dispatch on ploidy is needed
to build or pass it.

### What moves rather than merges

**Strand derivation.** A haploid group's strand is a property of the group -- one parent, one chain,
one carrying traversal -- and it is needed twice: to choose which of the parent's two haplotypes
conditions the chain, and to fill `PhaseCall::nested_strand`. It has to run before the grouping. It
is also the genuinely nested-specific logic and the part worth keeping intact: it carries the
`ploidy == 2` guard that chrX's haploid interior needs (0.94939 -> 0.93643 when it was missing) and
the nested-haploid-parent case that the identity match cannot find.

### What is deleted

The per-strand pass is 754 lines. Roughly 250 of them -- `by_strand`, the bucket sort, the site
build, the `haploid_posteriors` call, the settle loop, the `haploid_phasing` call and the PhaseCall
emission -- duplicate what the chain loop does. Those go, with `deferred_nested` and the hold-out
that fills it. The strand derivation, its counters and the `hap_contradicted` fixup stay.

### Increments, each with its own gate

| | step | gate |
|---|---|---|
| D1 | The chain loop's haploid branch takes the context. Inert: no haploid chain has one yet. | byte-identical |
| D2 | Strand derivation moves ahead of the grouping; the per-strand pass consumes the result instead of deriving it. | byte-identical |
| D3 | Behind a switch, route nested haploid sites into the (parent, chain) grouping with ploidy 1 and a single-haplotype delta, instead of into the strand buckets. | accuracy, three contigs then genome-wide -- expected neutral |
| D4 | Make it the default and delete the per-strand decode. | byte-identical against D3 with the switch on -- **done**, and it exposed three mosaic-only regressions the VCF gate could not see |
| D5 | The parent POSTERIOR, for both populations: harvest the per-strand marginal (`m` = 35 doubles per parent, against the `m*m` = 1,225 of the pair message deleted earlier) and use it in place of the delta. | accuracy -- **done and reverted**: byte-identical on chr20 and chrX, so the delta is sufficient |

D5 is the experiment the owner asked for and the one none of the five conditioning arms actually
ran: every one of them handed the child the parent's ARGMAX, which discards what the Li-Stephens
process computed about how certain that haplotype is. A constant-strength prior was swept (w = 0.5
to 1.0, 35x to infinity in prior odds) and was flat, which says calibration is not the issue -- but a
constant cannot emulate a per-site confidence, so it does not answer the question.

## D3 found a real bug on chrX, and it was the uniform rule that exposed it

Four contigs, switch off against on. The autosomes behave as the genome-wide conditioning arm
predicted -- a wash. chrX does not.

| chrX | off | on | |
|---|---|---|---|
| ALL F1 | 0.956487 | **0.958307** | +1.8e-3 |
| SNV | 0.968728 | 0.969484 | +7.6e-4 |
| JointIndel | 0.917183 | **0.922148** | +5.0e-3 |
| SV >=50 bp | 0.456872 | **0.476471** | **+2.0e-2** |
| small-variant FP | 4,409 | **4,095** | -314 |
| SV FP | 318 | **281** | -37 |

Fifty times any other movement in this rebuild, and the 1,120 records chrX "loses" are almost all
false positives. Autosomes for comparison: chr20 ALL -1.5e-5, chr6 +2.2e-5, chr17 +3.9e-6.

**The cause is structural.** chrX's 40,301 nested sites were bucketed by (phase set, strand) on a
contig that HAS no strand: `nested_strand_of` correctly returns -1 for a top-level haploid parent, so
they came out "on neither strand", unplaced, and whatever handled them then was emitting hundreds of
spurious calls. Under the uniform rule they are haploid children of a HAPLOID parent, the ploidies
match, and they join their parent's group and are decoded as part of a chain like anything else.

Two lessons, both about how this was gated rather than about the code:

**Three autosomes could not see it.** All 44,139 no-strand sites on chr20/chr6/chr17/chrX are chrX's,
and 95,339 more are chrY's; the 22 autosomes have zero. Every gate in this session ran on three
autosomes.

**The owner predicted it from the design, not the data.** "Surely chrX is just a special case of the
current code in which the top-level chain is (mostly) haploid?" -- which is exactly what it is, and
asking that question is what replaced a haploid special case with the ploidy-match rule that fixed
it. The rule I had written, "a haploid group does not contain its parent", would have excluded chrX's
parents for no reason and kept the bug.

## D4: deleting the pass was free, but three of its jobs were not the decode

The per-strand pass is unreachable once the grouping handles nested haploid chains, so D4 deleted it
along with `deferred_nested`, the two collectors that filled it, the dead ploidy-difference decline
counter, and the emission-masking arm behind `VG_LINKAGE_HAPLOID_PARENT` together with the two
`Site` fields nothing writes any more. Net **-731 lines**.

It was not free. The pass did three things besides decoding, and the grouping replaced none of them.
All three were caught by the unit tests and **none of them by the VCF**, which came out
byte-identical on chr20, chr6 and chrX with identical small-variant and SV F1 on all three.

**Haplotype slots are indexed by STRAND.** A haploid decode returns one haplotype and the render put
it in `hap_first` unconditionally; the mosaic reads `strand == 0 ? hap_first : hap_second` and calls
the other slot empty (`.`) rather than unexplained (`*`). So a chain on the parent's second strand
was reported carried on the strand it is not on and unexplained on the strand it is.

| chr20 mosaic | before | after |
|---|---|---|
| strand 0 named | 11,264 sites | 10,199 |
| strand 0 empty | 0 | 1,065 |
| strand 1 named | 8,958 | 10,023 |
| strand 1 **unexplained** | **1,066** | **1** |

20,222 named sites either way: a pure relocation, which is the signature to check for. 98 sites on
chrX, the same shape.

**The phasing was unconditioned while the posteriors were conditioned.** `haploid_posteriors` took
the entering message and `haploid_phasing` did not, so the genotype was settled against the parent
and the mosaic then named a haplotype chosen in complete ignorance of it. Both take it now, applied
to the Viterbi start of the chain's first window, falling back to uniform where the message and the
emission share no state. Seven chrX sites that had fallen to the wildcard now name a haplotype.

**A missing strand has two causes that want opposite renderings.** Under a haploid parent there is
no strand to choose because there is only one, and the haplotype is nameable -- that is all 44,139 of
chrX's, reproduced exactly by the new counter. Under a diploid parent whose settled pair does not
reach the chain there is no strand because the sample has no copy of the locus, and naming a
haplotype asserts a mosaic path through sequence the parent record does not carry.

**And a correction that belongs in the record.** The entering message is a prior in *form*, and both
this file and the header described it as one the emission could argue back against. It cannot: the
collector supplies a point mass, so a state the message excludes is unreachable. What the reads can
do is recombine away from it over the chain -- which means on a **one-site** chain the message
decides the answer outright. That is not a defect (79, 71 and 97 extra genotypes moved on chr20,
chr6 and chr17, and 2,599 on chrX, where it took SV F1 up 2.0e-2), but it is a stronger operation
than "prior" suggests and D5 should be read with it in mind.

**Gating note.** Three of these were invisible to VCF byte-identity and to F1 alike, because the
haplotype slots feed the mosaic only. The unit tests were the only gate that caught them -- and two
of the three fixtures had been recording nested sites at generation 0, a state the caller cannot
produce, which the deleted pass serviced and the grouping (gated to `generation > 0`) cannot. Any
future change to nested rendering wants the mosaic's per-strand accounting in the gate set, not just
the VCF.

## D5: the parent's posterior changes nothing, and the machinery is not kept

D5 was the last item: condition a child on its parent's POSTERIOR over panel haplotypes rather than
on the argmax of it. Built, measured, and **reverted** -- the arm is byte-identical to D4 on chr20
and chrX, VCF and mosaic, across 994 + 5,682 = 6,676 conditioned haploid groups. Not within noise:
not one record.

**Scoped to the haploid population, because the diploid half was already answered.** The diploid arm
ran the other way round earlier in this rebuild: the parent's exact posterior over ordered PAIRS was
the default and `VG_LINKAGE_HARD_PARENT` replaced it with the delta, which came out identical to six
decimals on every small-variant class on chr20 and chr6 at the cost of one SV false positive. A
strand marginal is a strictly lossier summary of that same distribution, so there was nothing there
to find. The haploid population is the one that had no parent message at all before D3.

**The first version of the arm was wrong, and the way it was wrong is the interesting part.** "The
parent's marginal on this strand" is not a strand-specific object. A top-level diploid site starts
from a uniform 1/(m*m), its emission is symmetric in the pair because `genotype_index` sorts, and its
two per-strand switch probabilities are equal, so the joint is invariant under swapping the strands
and the two marginals come out the same. Measured rather than assumed once the suspicion was there:
**673 of 887 diploid parents on chr20 have the two marginals equal to 1e-9, and 201 of 241 on chrX**
-- so it is most of them, not all, and the asymmetric quarter comes from steps where the parent's two
traversals have different lengths.

What IS strand-specific is the **allele** the parent settled on for the traversal the chain hangs
off. So the message became the posterior over which panel haplotype is carrying that allele: the
marginal restricted to haplotypes spelling it, renormalised. The argmax of exactly that distribution
is the delta, which makes the arm a clean softening rather than a different message.

**And that is why it is a null.** Once the parent's allele is conditioned on, the residual spread is
over haplotypes that all spell the same allele at the parent site -- and haplotypes agreeing there
tend to agree at the child, because that is what linkage is. The distribution's shape inside an
allele class carries nothing the child can use.

So the delta is not an approximation waiting to be improved. It is sufficient, now measured from
both directions:

| | message | result |
|---|---|---|
| diploid nested | exact pair posterior -> delta | identical to six decimals, one SV FP |
| haploid nested | delta -> allele-restricted posterior | **byte-identical** |

The plumbing is reverted rather than shipped behind `VG_LINKAGE_PARENT_POSTERIOR`: it is 203 lines,
it re-introduces a per-site retention that has no bound (which site is a parent is not known until
the next generation descends -- the deleted `wanted_parents` mask existed to answer exactly that),
and its value is the measurement, which is recorded here.

## Ploidy as a parameter: the premise, measured

"If done correctly the diploid/haploid property should be a simple parameter that simply gates
ploidy." Measured against the code rather than argued, and the answer is half yes.

| pair | diploid | haploid | identical lines | longest common run |
|---|---|---|---|---|
| `build_emission` / `haploid_emission` | 47 | 26 | 21% | 8 |
| `window_posteriors` / `window_haploid_posteriors` | 170 | 134 | 40% | 5 |
| `window_phasing` / `window_haploid_phasing` | 135 | 130 | 37% | 8 |
| `posteriors` / `haploid_posteriors` | 11 | 13 | **69%** | -- |
| `phasing` / `haploid_phasing` | 22 | 15 | 27% | -- |

**The cores are not duplicated.** They share their shape -- build emissions, apply the constraint,
forward, backward, reduce to genotype space -- and almost none of their code, because m against m*m
changes every index and the reduction is a different function. The "identical" lines are mostly
closing braces; no shared run exceeds 8 lines in 800. Merging them means a ploidy branch in the
innermost loop of the hot path, bought with the deletion of boilerplate. Not taken.

**The wrappers are duplicated, and they are where the damage was.** Both drift bugs this layer has
had were a parameter added to one ploidy and forgotten on the other: `haploid_posteriors` had no
`alpha_in` at all, and the haploid phasing had none while the haploid posteriors did -- which is how
the mosaic came to name haplotypes chosen in ignorance of the parent the genotype had just been
settled against. Neither is reachable through a single entry point. So `posteriors` and `phasing`
now take the ploidy and pick their own state space, `haploid_posteriors` and `haploid_phasing` are
gone, and so is the collector's adapter loop, since `phasing` returns `Phase` at both ploidies.

Line count is a wash (-3) and that is the design: the win is one door, not fewer lines.

Generalising to arbitrary ploidy was considered and rejected as speculative -- `-d` is rejected
outside {1, 2} at the CLI, at `graph_caller.cpp:720` and `:4883`, and asserted in
`snarl_caller.cpp:512`.

Inert: chr20 and chrX byte-identical, VCF and mosaic. `vg test` 12,547,771 assertions;
`18_vg_call.t` 318 alone.

## A gate that missed, again, and in a new way

D4 shipped failing `18_vg_call.t` 139. The test pattern-matches the wording of the
`[vg call] nested strands: ...` progress line; D4 replaced that line; and the gating had pinned the
binary and run TAP *before* the replacement went in, checking only VCF and mosaic byte-identity
afterwards -- which a stderr line does not touch.

The rule that would have caught it is not "check stderr too". It is that **a pinned binary
invalidates on any later source edit**, however inert the edit looks, because "inert" was only ever
established for the outputs that were compared. Fixed in `2c041db50`, with a second assertion added:
the old awk was vacuously true if the report never printed, so a report that stops being emitted now
fails loudly instead of passing every run.
