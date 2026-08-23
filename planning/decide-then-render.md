# Decide-then-render, haplotype-frame linkage, and parent conditioning

Successor to `nested-traversal-space.md`, which is implemented and validated on chr20 through
`4371c9b67`. That plan moved the linkage layer into the genotyper's own allele space; this one moves
the *ordering* of the pipeline to match, so that a genotype is decided before a record is rendered
and linkage runs along a haplotype rather than along the reference.

Designed by three independent passes under different framings (risk-first, value-first,
minimal-diff), each scored by separate feasibility, gate-quality and completeness reviewers, then
synthesised. The reviewers scored feasibility 5/10 across all three first drafts, and the corrections
that produced are in the section below — several planning premises turned out to be contradicted by
the code, including two of my own.


Working branch `read-likelihood-genotyping`, base `4371c9b67`. Precedent for format and tone: `planning/nested-traversal-space.md`.

## What is wrong now

Three things, in the order the target model puts them.

**Records are rendered before the genotype is decided.** `emit_variant` writes the line during the parallel read sweep from the pre-linkage genotype, and linkage afterwards patches the text it finds by `(contig, POS)`. On chr20 that costs, all [V]: 1,465 sites settle on a traversal the record has no ALT for, so the change is dropped — their false-positive rate is 61.7% against 11.0% where the layer *can* act, a 4.24x enrichment over GQ-matched records, and 455 of the contig's 2,008 false positives sit there. 1,362 of the 1,465 are top-level, whose render inputs are not retained, so "add the ALT at render time" reaches only 103 of them. 4,490 records reach a final genotype of 0/0 and there is no mechanism to retract a top-level line. `is_symbolically_reference` is evaluated at `src/graph_caller.cpp:2098` against the pre-linkage genotype and never revisited; the comparison itself is correct and recursive, only its timing is wrong. The recall mirror — a site whose pre-linkage genotype collapsed to reference, which linkage then settles on a non-reference genotype — produces no line and is not counted; nested chains have a re-render path (511 gained on chr20) and top-level sites have none, so for them it is architecturally forced [A].

**Linkage is a chain along a haplotype, and the code runs it in the reference frame.** Sites are sorted by reference POS (`src/linkage_model.cpp:2103-2109`) and the transition gap is a reference-POS difference (six sites: `:353, :489, :728, :853, :935, :1046`, all verified this session). Nested sites are held out of the diploid runs entirely (`src/linkage_model.cpp:1557-1571`, verified) and re-pooled by `(contig, strand)`, so a nested site links against same-depth sites under unrelated parents [A], and 12,516 chr20 children are skipped outright for having no reference path through them.

**The parent's genotype conditions nothing.** It selects the ploidy, the strand-group key, and whether the entry exists. The child's posterior is bit-for-bit identical for any two parent genotypes implying one copy on the same strand [A]. Target step 3 is the one part of the model the code does not implement at all.

Against that, four commits (`c6ed1ed3d..4371c9b67`) are pushed and chr20-validated: every genotyped snarl now reaches the linkage layer, ploidy and strand are derived once from which parent traversal carries the chain, and four incoherence classes that all came from treating "no VCF line was written" as "nothing to record" are closed. chr20 sits at ALL F1 0.97048 against a pre-refactor 0.96996, FP 2,008 against 2,093, every nested site on exactly one strand (6,716 of 6,716), no bare haploid GTs, no unphased records. None of the four is validated beyond chr20.

## Corrections to the backlog, verified this session

Six planning premises are contradicted by the code. Each changes a stage, so they are recorded before the stages rather than inside them.

**The retention instrumentation does not exist.** `git status` is clean apart from submodules; `grep -n "RETAIN\|retained_bytes\|would_retain\|TEMPORARY" src/graph_caller.cpp src/graph_caller.hpp` returns nothing. `src/graph_caller.cpp:30-53` is the descent-depth histogram (`g_descent_skipped_no_ref`, `g_descent_skipped_no_copy`), and `:2286-2330` is the `linkage_collector->record` call. The only retention meter in the tree is `LinkageCollector::bytes()`, reported at `src/graph_caller.cpp:699`. Stage 6 writes the sizing from scratch; it is not free.

**Retaining "the allele strings plus a bit" cannot render a record.** Every FORMAT field comes from `snarl_caller.update_vcf_info` (`src/graph_caller.cpp:2207` → `src/read_likelihood_caller.cpp:345-536`), which maps emitted alleles back to matrix columns by *structural comparison of `SnarlTraversal` objects* (`:390-399`), enumerates GL over the emitted allele count and looks each entry up in `info->genotype_lls` keyed by the sorted matrix-column multiset (`:444-452`), and derives QUAL by renormalising the all-reference posterior over that same map (`:508-520`). AD is `Number=R` over the emitted alleles. Rebuild the ALT list and every one of those has to be re-derived from the CallInfo. The minimum retention is the CallInfo plus the traversals plus the snarl — which is exactly what `PendingRecord` already holds (`src/graph_caller.hpp`, fields verified: `snarl, ref_path_name, ref_offset, travs, ref_trav_idx, genotype, ploidy, call_info, record_key, parent_record_key, parent_crossing, crossing_known, parent_trav, generation, dropped, emitted, buffer_thread, buffer_index`) and what the barrier already consumes at `src/graph_caller.cpp:3992`. There is no cheaper route; the design question is only whether that retention fits, not what to retain.

**`apply_linkage_change` is not only patch machinery.** `src/graph_caller.cpp:977-1030` recomputes GQ as the phred-scaled complement of the posterior, applies the `explained_share` discount, caps at GQI (the comment records that the discount alone does not make `GQ <= GQI` hold), rewrites GQN, and re-labels the `lowconf` FILTER that was set from the pre-linkage GQN. `test/t/18_vg_call.t:207-208` already asserts `GQ <= GQI`. A render pass fed the settled genotype must reproduce all of that or the deletion silently reverts post-linkage quality to per-site quality.

**`-1` in `haplotype_allele` is the escape state, not exclusion.** `src/linkage_model.cpp:139-155`: `ai < 0` gives `marginal[ai] * escape`, and both negative gives `overall * escape * escape`. Setting a non-carrying row to `-1` makes it a free wildcard at escape cost — a *weaker* constraint, not a restriction. The hard-constraint mechanism already in the tree is emission zeroing (`src/linkage_model.cpp:648`, `emissions[t][a * m + b] = 0.0`). Parent conditioning must use that.

**The `by_key` staleness is live, and the depth-2 strand story is not settled.** `src/linkage_model.cpp:1900-1908` is an 8-sweep loop whose own comment says it exists "so a nested site whose parent is itself nested finds its parent already placed". Inside it, `by_key[pc.record_key] = ParentPhase{... pc.trav_first, pc.trav_second ...}` (~`:1984`) executes *before* the block that assigns those fields (~`:1990-2006`), and `PhaseCall::trav_first`/`trav_second` default to `-1` (`src/linkage_model.hpp:422-423`). Separately, `linkage_phased` accumulates across generations (comment at `src/graph_caller.cpp:617-619`), so `by_key` is also seeded at `:1862-1866` with *correct* values from earlier generations. Consequence: a depth-2 site reading the seeded entry gets `parent.trav_first == parent.trav_second` (a nested parent has ploidy 1) and matches the first branch, so strand 0 always; a depth-2 site whose parent's in-loop insert has already overwritten that entry reads `-1`, matches nothing, and becomes *strandless*. Which one fires depends on the order of entries within `pending`, which derives from insertion order in the parallel sweep. So the measured 1,447/1,116 skew (ratio 1.30) may be a mixture of two mechanisms and may itself be run-dependent. Stage 3 measures the split before it fixes anything. `ParentPhase` is also a plain aggregate with no default member initialisers, so adding a field and updating only one of the two construction sites compiles clean and value-initialises the other to 0 — reproducing exactly the bug being fixed.

**The Poisson depth divisor is self-cancelling as written, so the "2x" is not supported.** `src/allele_likelihood.cpp:940-942` sets `depth_rate = local_read_rate(ranges) / effective_ploidy` where `effective_ploidy` is the site's ploidy, and `expected_reads` (`:53-61`) sums one term *per copy in the genotype*. So a diploid site expects `(local/2) × 2 = local × (len + L - 1)` and a nested haploid site expects `(local/1) × 1 = the same`. The copy count appears once in each place and cancels; today `DR` should read ≈1.0 at both classes, and changing the divisor to the region ploidy would *create* a 2x gap at nested sites (`DR` → ≈2.0), not close one. Also: `params.depth_ploidy` is assigned nowhere in `src/` (only the default `= 2` at `src/allele_likelihood.hpp:626` and two reads), so the header's claim that `-d`/`--ploidy-regex` feed it is already false. And `depth_weight` defaults to `0.1` (`src/subcommand/call_main.cpp:355`, plumbed at `:1503`), so the term is **on by default** and any change here moves default output. Stage 16 is therefore a measurement, and the code change is conditional on its sign.

Three smaller ones. The mosaic already declares a version (`src/graph_caller.cpp:1102`, `#mosaic-version 2`), asserted by `test/t/18_vg_call.t:467,480,487`. The VCF already uses `.` for "the strand that carries nothing" (`doc/read-likelihood-genotyping.md:553-557`, which also records why `*` was rejected there), so spelling the mosaic's empty strand `.` is consistent rather than novel. `LinkageModel` never reads `Site::ploidy`: model selection is `chain_ploidy = entries[indices.front()].ploidy` at `src/linkage_model.cpp:1600`, consumed at `:1644` and `:1729`, and the header at `:176` says "A whole chain shares one ploidy". Folding a ploidy-1 site into a diploid chain is a modelling change, not plumbing.

Finally, two circulating chr20 record counts conflict — 105,251 and 116,965 — and no gate should be written against either until one run settles it. Stage 1 reports it.

---

## Measured while planning, and two corrections to the brief

**Retention was priced before the plan was written, and the answer brackets rather than settles it.**
Temporary instrumentation summed the bytes a `PendingRecord` would hold for every snarl reaching
`emit_variant` on chr20, with nothing retained: **219,607 snarls, 296.6 MB, ~1,350 B each.** The
snarl count matches the 219,593 recorded linkage sites, which is the cross-check stage 7 asks for.
Two caveats, both load-bearing:

* The nested/top-level split from that run is **wrong** and must not be quoted. The test was
  `nested_context.active || current_generation > 0`, and `nested_context.active` is set only for
  `copies == 1` chains, so ploidy-2 nested chains were counted as top-level (192,207 reported against
  a true 165,408). Only the total is sound.
* It disagrees with the older 3.18 kB/chain anchor by 1.8x, because the two use different sizing
  models — 32 B per visit and 48 B per map node here, something heavier there. On chr20 the two give
  297 MB and ~620 MB. Projected to the largest contig at ~3.4x, 1.0 GB against 2.1 GB.

So the decision band in stage 7 (retain unconditionally / retain-and-release / abandon) is narrower
than the disagreement between the two estimates, which is exactly the case stage 7's own caveat says
estimates cannot settle. **Stage 7 stands, but its purpose is narrowed to arbitrating the two sizing
models against `LinkageCollector::bytes()`; the real answer is stage 8's measured peak RSS.**

**The denominator is packing density, not a per-contig ceiling.** `docs/wgs-performance.md:75-78`
gives the per-contig worst case as 6.1 GB (chr3), and line 4 gives the machine as 32 GB with several
contigs packed at once under a budget. Adding ~1–2 GB per contig therefore costs *concurrency*, and
so wall clock, rather than making any contig infeasible. The "24 GB budget" figure the plan text uses
below is not documented anywhere and should be read as the 32 GB machine under `--budget-gb`.

**Two backlog items in the brief were wrong, and the plan's corrections section is right to reject
them.** Both were passed to the plan as `[A]` and are now verified refuted:

* *Restricting the child's panel rows by setting them to `-1`* would **loosen** the model, not
  restrict it. `src/linkage_model.cpp:139-155` gives a `-1` row `marginal[ai] * escape` — a free
  wildcard at escape cost. The hard-constraint mechanism already in the tree is emission zeroing
  (`src/linkage_model.cpp:648`). Stage 19 must use that, and any description of parent conditioning
  as "restrict the panel rows" is misleading about what the code would do.
* *The Poisson depth term does not miscalibrate by 2x with handed-down ploidy.*
  `src/allele_likelihood.cpp:942` sets `depth_rate = local_read_rate / effective_ploidy`, and
  `expected_reads` (`:53-61`) sums one term **per allele in the genotype** before multiplying. The
  division and the summation cancel. Stage 20 keeps its measure-first shape, but the defect it was
  written to chase is not there as described.

---

# Phase I — instruments and small fixes (chr20-gated)

## 1. Byte-reproducible record order

**Goal.** Make two runs of one binary produce identical bytes. Six later gates are "output must not move", and that gate is not evaluable today: 72 chr20 record pairs sharing a position swap between runs [V].

**Changes.** `src/graph_caller.cpp:768-772`: the `std::sort` comparator over `all_variants` returns on `(contig, position)` only and `std::sort` is not stable, so records sharing a POS come out in per-thread-buffer concatenation order. `add_variant` files under `make_pair(make_pair(var.sequenceName, var.position), dest)` (`:597`). Extend the buffer key to `(contig, position, id)` and order on all three; `out_variant.id = print_snarl(snarl, false)` (`:2175`) is intrinsic to the site. Two consequential edits the key change forces: `linkage_changes.find(v.first)` and `linkage_phasings.find(v.first)` (`:815`, `:828`) are keyed `pair<string,size_t>` and must follow. Extract the comparator to a non-static, header-declared function so it is linkable from a unit test.

Order does not affect which record a patch lands on: `write_variants` already re-hashes the ID column per line and matches by `record_key` (`:790-800`), lazily and explicitly because "several records can share a `(contig, position)`, and every patch below must land on its own record".

**Scope limit.** `out_variant.id` is `print_snarl` only on the FlowCaller path. `VCFGenotyper` sets it from the input VCF (`:3120`), commonly `.`, so this stage does not make `vg call -v` or `vg deconstruct` reproducible. Say so in the commit rather than claiming the general property.

**Gate.** Baseline `4371c9b67`. (i) Two chr20 runs, same flags, `-t 8`: `cmp` on the VCF returns 0, and separately on the mosaic. (ii) `grep -v '^#' | sort | md5` identical to `4371c9b67`'s — order moved, content did not. (iii) Duplicate `(contig, POS, ID)` triples counted and reported as 0; without this the key looks total while it is not, and stages 8, 9 and 12 inherit an unproven premise. (iv) Report the chr20 emitted-record count, to settle 105,251 vs 116,965. (v) `vg test` (830) and `test/t/18_vg_call.t` (302) pass.

**Tests.** Unit test on the extracted comparator: it must be a strict weak ordering that separates two records with equal `(contig, position)` and different IDs, and sorting every permutation of a small vector must yield one order. This fails today for the right reason — the current comparator returns false in both directions, so output is a function of input permutation. `src/unittest/` has no `graph_caller.cpp` (verified: only `allele_likelihood.cpp`, `allele_likelihood_scoring.cpp`, `linkage_model.cpp`), so this stage creates it. **Do not add a fourth TAP determinism test:** `test/t/18_vg_call.t:370, 374, 402` already assert thread-count independence and already pass, which is direct evidence the in-tree fixtures do not produce ties. Record that in the commit; the chr20 `cmp` is the gate.

**Output moves:** yes (order only). **Reversibility:** one comparator plus two map key types; clean revert.

## 2. The mosaic stops asserting things it cannot support

**Goal.** Two defects that make the strand gates in stage 3 and the phasing gates in phase II unmeasurable. The mosaic spells "no sequence on this haplotype here" and "the panel cannot name a haplotype here" with one character, which is why #43's metric is unusable — raw wildcard segments rose 437 → 616 across the refactor while the count that means only the second thing fell 463 → 239. And a nested site's panel haplotype is copied from the parent and never checked against the child's own settled allele, so the mosaic can name a haplotype that demonstrably does not carry it. That second item is the one backlog bug no candidate plan scheduled.

**Changes.** (a) `src/graph_caller.cpp:1171-1173` writes `*\t*` under the single condition `hap == LinkageModel::WILDCARD`, so the two meanings are indistinguishable *at the writer*. They are set apart upstream — the deliberate empty strand of a nested haploid site carries `pc.nested_strand >= 0` (`src/linkage_model.cpp:1975-1982`), an unexplained strand does not — so the discriminator must be carried to the writer on `PhaseCall`. Emit `.` for the empty strand (matching the VCF's existing meaning, `doc/read-likelihood-genotyping.md:553-557`) and keep `*` for panel-unexplained. Bump `#mosaic-version` from 2 to 3 (`src/graph_caller.cpp:1102`) and update the header note at `:1114` and the column table at `doc/read-likelihood-genotyping.md:433`. Note the `gbwt_node` columns already use `.` for unresolvable (`:435`), so a consumer splitting on column index is safe and one grepping for `.` across columns is not.

(b) After step three moves a nested site's genotype (`src/linkage_model.cpp:2151-2185` updates traversals and alleles but not `hap_first`/`hap_second`), validate: if the haplotype named on a strand does not carry the settled compact allele in `hap_arena`, drop that strand to WILDCARD and count it. Propagate `PhaseCall::order_arbitrary` from parent to child, which today is set only in the diploid chain path.

**Gate.** Baseline stage 1. (i) The split reconciles against an independently produced number: the two new counts must sum to today's single wildcard count, and the panel-unexplained half must equal the separately reported figure (239 on chr20). Disagreement means one of the two code paths is wrong. Report both counts and the sum; do not gate on equality between *sites* and *records*, which are different populations — normalise to segments before comparing. (ii) The new counter "nested sites whose inherited haplotype does not carry the settled allele" is measured (unknown today, predicted > 0) and the panel-unexplained count then **rises by exactly that number**. Gate on the exact rise, not on the count falling: this stage deliberately claims less than before. (iii) VCF byte-identical to stage 1 — the mosaic is a separate file and no VCF field is touched. (iv) Harness consumers updated in the same change: `scripts/wgs/concat_mosaic.sh`, `scripts/wgs/mosaic_vcf_agree.py`, `scripts/wgs/nested_strand_check.py`, `scripts/tier2/mosaic_switches.py`, `scripts/tier2/phasing_benchmark.py`, each asserting the version.

**Tests.** Unit: a nested site whose parent's chosen haplotype carries allele A with the child settling on B — assert the emitted hap is WILDCARD and the counter increments. Fails today: the hap is asserted regardless. TAP on the nested fixture: both characters appear where both cases exist. Fails today because `.` never appears. `test/t/18_vg_call.t:480` (`version == 2`) and `:487` (exact header key set) must be updated in the same commit — they are the tests that prove the format changed.

**Output moves:** mosaic yes, VCF no. **Reversibility:** two commits; (a) is the format change with external consumers, keep it separate from (b).

## 3. Strand below depth 1: attribute first, then fix

**Goal.** Fix the depth ≥ 2 strand defect and the two mechanisms behind it. Roughly 166 chr20 records carry the wrong haplotype [V], and where the parent sits on strand 1 the child inherits `hap_first = WILDCARD` while the VCF asserts `a|.` — the two outputs contradict each other for the same record. The stage is split because, per the corrections above, the mechanism is not settled.

**Changes, in order.**

3a, instrumentation only, no behaviour. Report, per depth ≥ 2 site: whether its parent's `by_key` entry came from the accumulated seed (`src/linkage_model.cpp:1862-1866`) or from an in-loop insert; the resulting strand (0, 1, or none); and the sweep index at which it was placed. This partitions the population between "always strand 0" and "strandless", and it tells us whether the 1,447/1,116 skew is stable across two runs. If the strandless class is large, the headline claim "depth ≥ 2 strand is always reported as 0" is wrong and both later gates need re-derivation before the fix lands.

3b, the fixes. (i) Move the in-loop `by_key` insert (~`src/linkage_model.cpp:1984`) below the block that assigns `pc.trav_first`/`pc.trav_second` (~`:1990-2006`). (ii) Add `int8_t nested_strand` and `bool order_arbitrary` to `ParentPhase` (`src/linkage_model.cpp:1843-1859`) and populate them **at both construction sites** — `:1862-1866` and the moved in-loop insert. The struct is a plain aggregate; a partial update value-initialises to 0 and silently reproduces the strand-0 skew, so the compiler will not catch this and a reviewer must. (iii) In the strand derivation (`:1946-1952`), when the parent's own `ploidy == 1` and its `nested_strand >= 0`, take the parent's strand rather than running the identity match, and take the hap from whichever slot is not the wildcard.

3c, separately: `VCFOutputCaller::resolve_linkage()` (`src/graph_caller.cpp:601-604`) calls `resolve_linkage_generation(0, true)` and entries above generation 0 are skipped (`src/linkage_model.cpp:1512`), so any path reaching `write_variants` without `run_deferred_descent` drops every nested site from linkage, phasing and the mosaic. Loop generations there, re-reading `max_generation()` each pass and passing `last` only on the final one, as `run_deferred_descent` already does (`:3788`, `:4092`). **The fix is in `resolve_linkage()`, not `LinkageCollector::resolve()`** — that method (`src/linkage_model.hpp:478`) is called only from `src/unittest/linkage_model.cpp:362, 380, 399, 400, 421, 422`, so looping it changes test behaviour and leaves the shipped hole open.

**Gate.** Baseline stage 2, all chr20. (i) Records at depth ≥ 2 with `nested_strand == 1` rise from the 3a figure (predicted 0) to non-zero. (ii) A floor so that "place nothing" cannot pass: nested sites placed on exactly one strand must not fall below 6,716 of 6,716, and the placed/unplaced split at depth ≥ 2 must be reported against 3a's. (iii) Zero records where the VCF names an allele on a strand while the mosaic gives that record `*` on both strands — measurable only because stage 2 split the character, and using `scripts/wgs/nested_strand_check.py`, which walks per-strand wildcard intervals; **not** `mosaic_vcf_agree.py`, whose own docstring says it checks arity and phase-set correspondence. (iv) The half-called split: report it, but do not gate on the asserted 1.00 ± 0.10 band. On n ≈ 2,563 a binomial null gives σ ≈ 25, i.e. ±0.04, and the strand-0/strand-1 asymmetry can also come from the Viterbi's own ordering convention and `order_arbitrary`; a band derived from 3a's measurement, or a report, is honest and 1.10 is not. (v) Accuracy A/B: ALL F1 within 0.0005 of 0.97048, FP not above 2,060. Flipping which side of `a|.` an allele sits on does not change an unordered genotype — but it *does* change which `by_strand` group the site joins and therefore which haploid chain links it (`src/linkage_model.cpp:2011`), so a larger move means the regrouping is doing more than relabelling, and that is the finding. (vi) chr20 phase-block count reported either side (22 today).

**Tests.** Unit, extending the existing case at `src/unittest/linkage_model.cpp:523`: a diploid parent, a nested child placed on strand 1, a grandchild off that child; assert the grandchild's `nested_strand == 1` and its hap is not WILDCARD. Fails today two ways — the identity match cannot return 1 when `trav_first == trav_second`, and the strand-1 inheritance sets `hap_first = WILDCARD`. Unit: a parent and child resolving in one pass; assert the child takes the parent's strand. Fails today — the parent's in-loop entry carries `-1/-1`, no branch matches, and the site lands in `unplaced_no_strand`. Unit for 3c: `resolve_linkage_generation` over a two-generation entry set produces a PhaseCall for both generations; fails today at the `:1512` filter. The collector is free of graph types, so none of these needs a fixture graph.

**Output moves:** yes (3b, 3c may). **Reversibility:** four commits (3a instrumentation, then i/ii/iii together since they share a struct, then 3c). 3a is discarded after 3b lands.

## 4. Ploidy provenance

**Goal.** Three ways a site's copy number can be wrong with no visible symptom. A nested record called at ploidy 2 inside a haploid `--ploidy-bed` interval is a well-formed diploid GT; chr20 exercises no ploidy region at all, so only chrX can show it.

**Changes.** The mechanism is not the one the backlog names. `child_ploidy` already caps at the parent's ploidy (`src/graph_caller.cpp:4643`), so a child cannot exceed its parent through the normal path. The route to a diploid GT in a haploid interior is `call_snarl_internal(..., copies >= 1 ? copies : 2)` at `:4696`: a retained chain (`copies == 0`) is handed a hard-coded 2. (a) Replace that literal with the region's ploidy. (b) `src/graph_caller.cpp:4366-4367` prefers `ploidy_override` over `ploidy_at` unconditionally; the override is a copy count *within* the region's ploidy, so compose as `min(ploidy_override, ploidy_at(...))`. (c) `ploidy_at` (`src/graph_caller.cpp:568`, declared `src/graph_caller.hpp:166`) reads `interval_start` only, so a snarl straddling a BED boundary takes one edge's answer; take the minimum over the snarl's reference interval, on the ground that over-calling a haploid region invents a haplotype. This is a signature change across the four `ploidy_at` call sites (`:3365`, `:4366`, and two more), not two one-liners.

(d) The >64-traversal hole: `child_crossing_mask` returns unknown for a parent with more than 64 candidate traversals, and that parent plus its whole subtree keeps the pre-linkage ploidy, reported as `crossing_unknown` (`src/graph_caller.cpp:3900`, `:4144-4146`). **Measure first.** If chr20's count is small, the correct action is to drop those subtrees and count the drop, which is conservative and visible; building a wider mask for a rare population is the wrong trade. Stage 11 deletes the mask entirely, so this stage should not widen it.

**Gate.** Baseline stage 3. chrX with `scripts/wgs/chrX.par.bed`. (i) Records inside the haploid interior carrying a diploid GT: measure first; if the pre-change count is 0 the gate passes 0 → 0 and the change is not validated — say so rather than claiming it. If non-zero, it must reach 0. (ii) Snarls straddling a PAR boundary: report the count and assert each takes the haploid answer. (iii) Records inside the PAR carrying a haploid GT stay at 0 — this must not make the PAR haploid. (iv) chr20 byte-identical to stage 3: chr20 has no ploidy regions, so `ploidy_at` returns the fallback and `min` is the identity; any movement means (c) changed the fallback path. (v) chrX phase-block N50 reported, expected to change only at PAR boundaries — a ploidy change cuts a chain (`src/linkage_model.cpp:1553-1571`). (vi) `crossing_unknown` measured on chr20 and chrX and reported; not gated to 0 here.

**Tests.** Unit on `ploidy_at`: a snarl interval spanning a 1↔2 boundary returns 1, in both directions of travel — the value is that it pins both, since one direction passes by luck today. TAP with a `--ploidy-bed` marking a haploid region containing a nested fixture: no diploid GT inside it at any depth; fails today via the `:4696` literal. **`test/t/18_vg_call.t:235-236`** (`OUT_WINDOW_HAP == 0`, "`--ploidy-bed` leaves calls outside the window at the `-d` ploidy") is *contradicted* by a min-over-interval rule: a snarl starting outside the window but overlapping it becomes haploid. That assertion must be re-baselined in this commit, and the commit message must say which behaviour changed and why.

**Output moves:** yes on chrX, no on chr20. **Reversibility:** three commits; keep (c) last and re-run the chr20 no-op check after each. `min` is not trisomy support — ploidy is clamped to {1,2} in the collector (`src/linkage_model.cpp:1216`, `:1348`) and the model's state is an ordered pair; say so in the commit.

## 5. `-L` GL fold, made layout-aware

**Goal.** `merge_similar_alleles` folds GL with an i-major index against a caller that writes it colexicographically, so merged three-allele records get two of six GL entries transposed. Off by default, and when on, the failure is invisible.

**Changes.** `src/graph_caller.cpp:1827`: `gl_index = i*n - i(i-1)/2 + (j-i)` is i-major, and the comment nearby asserts `PoissonSupportSnarlCaller` is the only GL writer. It is not: `ReadLikelihoodSnarlCaller::update_vcf_info` writes GL in `AlleleReadLikelihoods::enumerate_genotypes` order (`src/read_likelihood_caller.cpp:444-452`), which `src/allele_likelihood.cpp:236-249` documents as colexicographic and notes differs from the Poisson caller. At n=3 the two disagree at indices 2 and 3 — `(1,1)` against `(0,2)`. But `src/snarl_caller.cpp:872-915` really is i-major and is still live, so **indexing unconditionally through `enumerate_genotypes` would corrupt merged Poisson records.** The fold must be layout-aware. `merge_similar_alleles(graph, site_traversals, site_genotype, sample_name, out_variant)` has no `snarl_caller` argument and `out_variant` does not name the writer, so the layout has to be passed in from the call site at `src/graph_caller.cpp:2226`, or the two writers unified first (a separate commit).

The second `-L` defect — `trav_to_allele` built at `:2078-2118` before `merge_similar_alleles` renumbers at `:2226`, so every merged record's patch is declined — is **not** fixed here. Phase II deletes patching, and the defect goes with it. Fixing it now is work thrown away.

**Gate.** Baseline stage 4. (i) Measure first: the count of merged records violating `test/t/18_vg_call.t:202`'s GT-indexes-max-GL invariant on a chr20 `-L` run. The fold is a max-marginal over collapsed classes, so a transposed read need not violate the argmax — if the pre-change count is 0 the chr20 gate is vacuous and the unit test is the only evidence. Say so. (ii) The same invariant on a `-L` run with `PoissonSupportSnarlCaller` must not regress — this is a gate condition, not a note in the risk section, because it is the way this change breaks something. (iii) Default (no `-L`) chr20 byte-identical to stage 4.

**Tests.** Unit on the fold: a three-allele colexicographic GL with distinguishable values and a merge of allele 2 into 1, asserting the folded vector, constructed so the transposition changes the *answer* and not merely a label. This requires extracting the fold — `merge_similar_alleles` is behind `protected:` in `src/graph_caller.hpp` — so the extraction is part of this stage, not a footnote. Fails today, exactly, because indices 2 and 3 are read transposed. Second unit case at the Poisson layout, asserting the fold is unchanged there.

**Output moves:** `-L` only. **Reversibility:** two commits (extraction, then the layout fix), both behind a default-off flag.

---

## 4. Ploidy provenance — MEASURED ZERO, NOT IMPLEMENTED

Every one of stage 4's four sub-items measures no observable effect, so none was implemented. The
mechanisms are real; their consequences are not reachable in the output.

chrX with `scripts/wgs/chrX.par.bed`, 114,207 records, on `a8b9ea448`:

| | |
|---|---|
| haploid interior records | 105,355 |
| of those with a genuinely diploid GT | **0** |
| PAR records | 8,852 |
| of those with a bare haploid GT | **0** |
| records whose REF straddles a PAR boundary | **0** |
| `crossing_unknown` (>64 parent traversals), chrX and chr20 | **0** |

The first measurement needed correcting before it meant anything: a naive "diploid GT" test counted
9,214, all of them `a|.` — nested haploid records in the deliberate half-called form, matched because
`.` sat inside the character class. Requiring both sides to be called gives 0.

Why zero. The hard-coded `2` at `src/graph_caller.cpp:4768` does hand a retained chain ploidy 2 in a
declared-haploid region, but the barrier then re-renders it at the ploidy the settled parent implies,
which for a haploid parent is at most 1. So the wrong ploidy is chosen and then corrected, and the
output never shows it. The same masking covers `ploidy_override` overriding `ploidy_at`.

**This is a Phase II prerequisite, not a closed item.** Phase II changes the barrier's role — the
record is built once, after the decision — so a ploidy chosen wrongly at descent is no longer
corrected downstream. These three edits (the literal, composing the override as a minimum, and
`ploidy_at` taking the minimum across the snarl's interval rather than reading its start) must land
*with* stage 9 or 10, where a gate can see them. Doing them now would be three unverifiable changes
to a numerical pipeline, and (c) is a signature change across four call sites at that.

`crossing_unknown` being 0 on both contigs also settles sub-item (d): there is nothing to widen the
mask for, and stage 11 deletes the mask anyway.

**What this run did validate**, incidentally and worth keeping: chrX is the only available contig with
mixed ploidy and more than one chain, and stages 1–3 are clean on it — all nested sites placed on
exactly one strand, no bare haploid GT anywhere in the PAR, no uncheckable masks. chr20 cannot test
any of that.

# 6. Whole-genome run 1

**Goal.** Validate the four pushed commits and stages 1–5 in one run. Nothing is validated past chr20 since `a27149728`, and the pushed work nearly doubled the linkage site count (117,210 → 219,246) and grew the arena 15.8 → 29.0 MB with chr20 wall +22% and the linkage pass 18.8 → 41.1 s. chrX with `--ploidy-bed`, chrY's haploid chains, and the acrocentric and centromeric contigs exist only at genome scale. This run also establishes the baseline every later comparison uses.

**Changes.** None. `scripts/wgs/schedule_wgs.py` to run, `scripts/wgs/bench_wgs.py` to score (aardvark `JointIndel` row, not `Indel`; plus truvari).

**Gate.** Against the stored `a27149728` figures in `docs/wgs-results.md`: autosomal ALL F1 not below 0.9699, SNV not below 0.9833, SV not below 0.5470 by more than 0.0020 — and additionally not below the level chr20 established for the four pushed commits, so a run that loses the measured +0.0005 fails rather than passing on the pre-refactor floor. Rates recomputed from summed TP/FP/FN, never averaged per contig. Per-contig peak RSS at matched thread count, warm page cache, ≥3 repeats, against the 0.7 GB noise floor and the 24 GB budget; summed CPU (not wall) within 10% of the stored 459.9 minutes. chr20's contribution must reproduce stage 5's chr20 figures **exactly** — after stage 1, identical code has no run-to-run noise, so "within harness noise" is the wrong tolerance.

Four structural invariants added to `bench_wgs.py` as hard assertions, each confirmed able to fail by running it against a pre-`4371c9b67` VCF where all four were violated (292, 75, 440, and 5,433-of-5,892 respectively): no nested record with a bare haploid GT; no record with no phase; no nested site on both or neither parent strand; no GT naming an allele index beyond its record's ALT list. The mosaic-accounting invariant is stated as **entries** versus **emitted records**, not "the mosaic accounts for exactly the emitted record count" — stages 9 and 13 make the latter false by design.

Also refit the scheduler's memory model (`peak GB ≈ 2.25 + 11.2e-6 × emitted_records`, `docs/wgs-performance.md:87`) and record the residual. It currently predicts 6.21 GB for chr1's 353,741 records against a measured 5.7 GB — a 0.5 GB residual, half the tolerance stage 6 wants to apply to it.

**Output moves:** no. **If it fails:** the failing subsystem is identified before phase II is written, which is why the run is here. Bisect on chr6 (2.7x chr20, ~10 min) rather than spending another genome.

---

# Phase II — decide-then-render

## 7. Price the retention

**Goal.** Answer the memory objection in the `LinkageCollector` header with a number, per contig, before writing the design around it.

**Changes.** New instrumentation (it does not exist — see corrections). For every snarl reaching `emit_variant`, accumulate the bytes a `PendingRecord` would hold, split by nested versus top-level and **keyed by contig**, reported under `--progress`. Report three subtotals, because they price three different decisions: top-level with both ploidies' `genotype_lls`; top-level with only the primary ploidy, since a top-level snarl's ploidy comes from the contig or the BED and the barrier never revises it, so the second answer is dead weight there and stage 8 will not retain it; and the nested figure, as a self-check.

**Denominator.** `vg call` runs one contig per process (`scripts/wgs/schedule_wgs.py`), and `docs/wgs-performance.md:74-78` states it: "peak memory scales with the contig rather than the genome ... per contig the worst case measured is 6.1 GB (chr3)". The decision quantity is therefore retention on the largest contig — chr1 at 353,741 records, measured 5.7 GB, roughly 3.4x chr20 — **not** a genome-summed figure no process ever holds.

**Gate.** One chr20 `--progress` run, output byte-identical to stage 5's. It must print the three subtotals and the top-level snarl count (expected ~165,408). Self-check: the nested figure must reproduce the already-measured 3.18 kB/chain and 87 MB for 27,404 chains within 10%, and the total must be consistent with `LinkageCollector::bytes()` at `src/graph_caller.cpp:699`. If it does not, the sizing model is wrong and no other figure from it may be trusted. Cross-check against the measured chain figure before trusting the estimate at all: 165,408 top-level snarls at 3.18 kB is ~0.53 GB on chr20 and ~1.8 GB projected to chr1, so the "1 kB per snarl / 165 MB" arithmetic in circulation is optimistic by roughly 3x.

**Decision rule, fixed before the number is seen, and per-contig.** Retain unconditionally if the projected chr1 delta keeps chr1 under 7.0 GB against a 24 GB budget with the refitted model from stage 6. Between 7.0 and 12 GB, retain but release each top-level chain's record as its generation resolves, and stage 8's gate becomes the released-early peak. Above 12 GB, retention is dead and phase II is rewritten around widening the ALT list at emission — which loses the recall mirror and the POS/REF/ALT renormalisation for top-level sites permanently, and that loss must be written into the design note rather than left implicit. The threshold is the maintainer's to set (open decision 1); what matters is that it is set first.

**Gate can fail:** if the estimate and `bytes()` disagree by more than 2x, neither can forecast chr1 and stage 8 must be run purely as a measurement — its byte-identical gate makes that safe.

**Tests.** None; this stage produces a number and a decision. It is not committed to the branch.

## 8. Retain the render inputs at top level, change nothing else

**Goal.** Land the whole memory cost as one commit whose gate is that output must not move, isolating the one risk a reviewer cannot check by reading.

**Changes.** The top-level emit is the first branch (`src/graph_caller.cpp:4486-4496`); the nested branch is the `else if (ploidy_override >= 0)` at `:4498`, and only it stages a `PendingRecord` (`:4536-4558`). Stage one for every snarl that reaches genotyping, in all three emit branches (`:4489`, `:4521`, `:4571`), and have `run_deferred_descent` ignore generation 0. Keep `set_want_alt_ploidy` off for top-level snarls. `pending_records` is sized per thread only inside `set_defer_nested_descent` (`:3760-3768`), so sizing must become unconditional, and the render queue must be a **second container**: `run_deferred_descent` moves `pending_records` out and clears it (`:3792-3797`), and `children_of` would otherwise bucket every top-level record under `parent_record_key == 0`.

`PendingRecord::travs` is held by value and the record is stored only after descent finishes with it (`:4702-4707`); the top-level branch has no descent loop reading `travs` afterwards, so the move is safe, but use the same staging discipline in both branches so they read alike. The use-after-move that cost 2,494 chr20 records (`906812957`) lived exactly here.

**Gate.** chr20 byte-identical to stage 7's (`cmp` = 0) — nothing reads the new records. Retained top-level count equals the top-level snarl count. Peak RSS at matched threads, warm cache, ≥3 repeats, against stage 7's, inside the band stage 7's rule accepted; report the arena figure too. If RSS exceeds `a27149728`'s chr20 peak of 3.91 GB, stop: stage 7's projection was wrong, and the route decision reopens here rather than three commits later.

**Tests.** TAP: `pending_record_count` under `--progress` is non-zero for a run with no nesting at all; fails today (0 for a graph with no child snarls). Confirm the count is not vacuous by deliberately skipping one per-thread queue in a scratch build and checking the assertion fires.

**Output moves:** no. **Reversibility:** one commit; the byte-identical gate means reverting cannot change output either.

## 9. Render after linkage resolves, from the same genotype

**Goal.** Move the emission point without moving the decision, so the largest code motion in the plan has a near-byte-identical gate and the behaviour change (stage 10) is a one-argument diff.

**Changes, and they are inseparable.**

(a) **`record()` must move out of `emit_variant`.** This is the point every candidate plan got wrong. `emit_variant` calls `linkage_collector->record(...)` at `src/graph_caller.cpp:2290`, guarded by `suppress_linkage_record`, which the barrier sets around its own emit (`:3991-3995`) precisely so the barrier calls `respecify` instead. If `record()` stays inside `emit_variant` and `emit_variant` moves to a post-sweep pass, the barrier resolves an empty collector and every nested revision, retraction and gain disappears. Move `record()` to the genotyping site. Its two arguments that came from the variant are both derivable there: `record_key` is `std::hash<string>{}(print_snarl(snarl, false))`, already computed that way at `:4545`, and `position` becomes `get_ref_position`.

(b) The position source therefore changes from the post-flatten `out_variant.position` to the pre-flatten reference position. Under patching, the post-flatten value was the `(contig, POS)` key the patch had to find; with no patch key, position reverts to what the model uses it for — ordering (`src/linkage_model.cpp:1540-1546`, `:2103-2109`) and the six transition gaps. `flatten_common_allele_ends` advances POS by the prefix every allele shares, typically 0 for a SNV and 1 for an indel.

(c) Suppress the sweep's emit for records that have a retained `PendingRecord`, and add a pass after `resolve_linkage` that walks them in generation order and calls `emit_variant` with `pr.genotype` — the same pre-linkage genotype. Parallelise it over the queue; the buffers are already per thread, and a serial render of six figures' worth of records is a new cost with no benefit.

(d) `PendingRecord` has no field for `nested_context.active`, which `record()` reads (`:2323`). Add one. Missing it is not a compile error — it writes a plausible small value, which is the failure mode the previous refactor hit at its stage 1.

**Scope.** `emit_variant` has six callers, verified: `:3379` (`VCFGenotyper`/`vg call -v`), `:3992` (barrier), `:4489`/`:4521`/`:4571` (`FlowCaller`), `:5155` (`NestedFlowCaller`, `-A`). `VCFGenotyper` and `NestedFlowCaller` have no `LinkageCollector`, and `NestedFlowCaller` passes a genotype-dependent `trav_to_flat_string` (`:5142-5152`), so its allele strings are not a function of the traversal alone. **Those two paths keep sweep emission.** State it in the commit.

**Gate.** Baseline stage 8. Not byte-identity — position now comes from a different source, and pretending otherwise produces a gate that will be waived. Total genotype differences below 50, **and every differing position must be one whose flatten prefix is non-zero** — dump the prefix per site and check, do not assume. ALL F1 equal to 0.97048 to four decimals; record count equal to stage 8's; patch declines equal to stage 8's, since the same patches are being applied by a different key. Above 50 differences, the motion is not faithful and the cause must be found before stage 10.

**Risk and the audit that addresses it.** A render subtly different from the emitter it replaces produces a plausible VCF. Grep every construction of `out_variant.samples[...]["GT"]` and every consumer of `Change::allele_i/allele_j` and `PhaseCall::allele_first/allele_second`, and check each goes through one application point. Second: `last_emitted` is thread-local and read immediately after the emit that filled it (`:4649-4665`) to build children's crossing masks. Moving the emit out of the sweep breaks that; until stage 11 deletes the masks, the sweep must compute the mask from `travs` directly, which `child_crossing_mask` already accepts. Third: `record()` moving to the genotyping site changes the collector's `entries` insertion order from sweep order to genotyping order. Genotype resolution is order-insensitive (chains sort on position then `record_key`), but the nested strand pass iterates `entries` in insertion order and reads `by_key` built during the same walk — the very interaction stage 3 fixed. This is the most likely way this stage fails, and stage 3 landing first is what makes it detectable.

**Tests.** TAP, permanent: no GT names an allele index beyond its record's ALT list — one `awk` line. Confirm it discriminates by planting a synthetic violation, or by pointing at the 48-record failure the previous refactor produced at this exact stage. TAP: `FORMAT/PS` is present on every record carrying a phased GT; confirm it fails by removing the phase application. Unit test on the render function against a hand-built site with a stored expected line, so a future change has something to fail against.

**Output moves:** yes, bounded. **Reversibility:** large but self-contained; keep stage 8 in place across any revert so the retention measurement is not lost.

## 10. Render from the settled genotype

**Goal.** The behaviour change the phase exists for. Six defects stop being representable rather than being fixed one at a time.

**Changes.** Hand the render pass the settled pair instead of `pr.genotype`. The collector stores it as compact indices `Entry::final_i`/`final_j`; add an accessor decoding them to traversals through `trav_arena`, since a compact index means nothing outside the collector. Everything downstream then follows because the ALT list at `src/graph_caller.cpp:2084-2118` is built by iterating the genotype it is given: no settled traversal is unrenderable; `is_symbolically_reference` (`:2098`) and `wants_line` (`:2267`) are evaluated on the settled pair, so a settled-reference site writes nothing and a settled-non-reference collapsed site writes a line for the first time; QUAL (`src/read_likelihood_caller.cpp:508-536`) is a pure function of the genotype passed in and is computed from the settled one; POS/REF/ALT and the arity of AD/GL/GQI are normalised against the settled allele set.

**And the part no candidate plan scheduled:** the render must reproduce `apply_linkage_change`'s post-linkage quality arithmetic — posterior-derived GQ, the `explained_share` discount, the cap at GQI, the GQN rewrite, and the `lowconf` FILTER re-labelling (`src/graph_caller.cpp:977-1030`). Without it, quality silently reverts to per-site quality on every changed record.

`LinkageCollector::Change` loses `called_i`/`called_j` (they existed to verify a patch landed on the right record) and `Entry::allele_offset`, `allele_arena`, `traversal_to_allele`, `vcf_allele_of`, `render_phase_pair` lose their only consumer — with them goes the whole `phase_fallback` population (5,015 chr20 records phased on the line's alleles rather than the model's).

**Gate.** Baseline stage 9, structural first because a wrong answer cannot fake these. (i) `unrenderable` 1,472 events at 1,465 positions → 0, and the counter deleted (`src/linkage_model.cpp:1707`, `:2167`). (ii) Records with a final GT of 0/0 or 0|0: 4,490 → 0 without `-a`. (iii) No record with GT 0/0 and QUAL > 0. (iv) Records gained because the settled genotype is non-reference where the called one was not: measured for the first time; report it, do not net it against the losses. (v) No GT past the ALT list. (vi) `GQ <= GQI` on every record — `test/t/18_vg_call.t:207-208` already asserts this and it is the check that the GQ arithmetic was carried over; confirm it fails by omitting the `explained_share` discount, which the code comment says leaves ~5% of records violating it. (vii) GL length is `na(na+1)/2` and GT is the argmax of GL (`test/t/18_vg_call.t:196`, `:202`) — both become hard constraints on the render once the ALT list is rebuilt.

Then accuracy, with the prediction stated in advance and the arithmetic shown: 336 of the 545 *judged* unrenderable sites are false positives and the layer lands at 11% where it can act, so acting on all of them should remove ≈275. **FP must fall from 2,008 to at most 1,800; below 1,900 is the abandon threshold.** ALL F1 at least 0.9710; FN not above 3,537. If FP does not fall, the 4.24x GQ-matched enrichment did not mean what it was read to mean, and that is the most important negative result this plan can produce — record it rather than tuning around it. Note the gate's own weakness: only 545 of the 1,465 fall in confident regions, so 63% of the population is unjudged and not at random. Diagnostic if it fails: FP rate at the 545 judged positions specifically. If that falls from 61.7% toward 11% while total FP rises, the cost is in the previously-unjudged population and the honest response is to report it and stop.

Do **not** gate on "QUAL == 0 iff GT is all-reference": `src/read_likelihood_caller.cpp:508-520` leaves QUAL at 0 whenever `have_ref` is false, on legitimate non-reference calls, independent of linkage.

**Tests.** Unit: for a site whose posterior argmax differs from the called pair, the settled-traversal accessor returns the argmax's traversals; fails today because no such accessor exists. TAP: a fixture where linkage moves a site onto a traversal the pre-linkage record had no ALT for; assert the emitted record carries that ALT. Fails on stage 9 by construction — that is the 1,465 population. TAP: a fixture whose pre-linkage genotype collapses to reference and whose settled genotype does not; assert a record exists. Fails on stage 9: no record is ever created. TAP: no record carries GT 0/0 without `-a`; confirm it fails against stage 9's output, where it is violated 4,490 times, and rebuild the fixture if it does not.

**Existing tests this breaks, all of which must be re-baselined in this commit:** `test/t/18_vg_call.t:120-121` (0/0 counts with `-a` versus `-v`), `:279` (non-0/0 count in a cyclic graph is 3), `:295` (non-0/0 counts match across GBWT enumeration), `:1259` (0/0 produces no non-ref), `:131-132` (GT comparison across two modes).

**Output moves:** yes, substantially. **Reversibility:** one argument. Revert restores stage 9 exactly, which is the whole reason stage 9 exists.

## 10 result: rendered from the settled genotype, and three ordering bugs the gate could not see

Implemented as three commits' worth of change: record at the genotyping site rather than at emission
(10a), normalise the locus name (10b), render from the settled pair (10c). Structural gate, against
stage 9:

| gate | stage 9 | stage 10 | target |
|---|---|---|---|
| (i) unrenderable settled genotypes | 1,472 | **496** | 0 |
| (ii) final GT of 0/0 or 0\|0 | 4,490 | **1,383** | 0 |
| (v) GT past the ALT list | 0 | 0 | 0 |
| (vi) GQ > GQI | 0 | 0 | 0 |
| records | 116,966 | 115,618 | — |

(i) and (ii) fall by 66% and 69% but **not to zero, and the gate as written called for zero.** The
residue is nested chains, which still emit inline from the pre-linkage genotype because the barrier
itself emits: only the top-level render pass moved. Driving both to zero needs the barrier to stop
emitting, which is stage 11's deletion, not stage 10's reordering. Recorded as a partial pass rather
than a pass.

Accuracy, against the stated prediction:

| class | stage 9 | stage 10 | delta |
|---|---|---|---|
| ALL | 0.97048 | **0.97231** | +0.00183 |
| SNV | 0.98436 | 0.98523 | +0.00087 |
| Insertion | 0.90843 | 0.91504 | +0.00661 |
| Deletion | 0.93266 | 0.93724 | +0.00458 |
| JointIndel | 0.91840 | 0.92390 | +0.00550 |
| SV | 0.51875 | 0.52099 | +0.00224 |

**The prediction was wrong about the mechanism, and that is the finding.** FP was predicted to fall
from 2,008 to at most 1,800, with 1,900 as the abandon threshold; it fell to 1,972 — a 36-record
improvement against a predicted 275. The gain is almost entirely recall: TP 91,154 → 91,454 and
FN 3,537 → 3,237, both moving by 300. So the 4.24x GQ-matched enrichment at unrenderable sites did
**not** mean those sites were false positives waiting to be removed; rendering from the settled
genotype instead *recovers true calls that the pre-linkage ALT list could not express*. 1,840 records
gained, 3,174 lost. On the gate's literal terms FP missed its threshold, and on the gate's own
"record it rather than tuning around it" instruction the result stands: every class improved, and the
mechanism is now known to be recall.

### Re-gated with phasing working, and the result confirms the invariant

The accuracy above was first measured while phasing was entirely absent (bug 2 below). Re-run with it
restored, chr20 comes out **identical to five decimal places, with identical TP/FP/FN**: ALL 0.97231,
91,454/1,972/3,237. Unphased records fell from all 116,952 to 144; record count unchanged at 115,618.

That equality is the evidence, not a formality. Aardvark's GT comparison is order-insensitive, so
identical counts across a run that gained phasing is a direct measurement of the property the phased
output is supposed to have — phasing re-orders a genotype without re-deciding it — at chr20 scale
rather than on a 70-record fixture. `nested_strand_check` agrees: `no_strand` 0,
`strand_not_recoverable` 0, strand correct on 180 of 286 decisive sites (62.9%) against a best trivial
control of 56.6%, unchanged from before the render moved.

### Three bugs, all of them ordering, none visible to the structural gate

Each was found only by a test or a fixture, and each would have passed every gate above.

**1. Locus versus path name.** `record_site` keys on `get_ref_position`'s answer, which is the base
path name `CHM13#0#chr20`; emission reduces it to `chr20`. Every patch lookup therefore missed and
**all 116,952 records came out unphased** — with no counter reading wrong, because the patch index was
consistent with itself. Fixed with `PathMetadata::parse_locus_name`.

**2. `emitted` is a snapshot, and the phase alleles with it.** `linkage_phased` holds `PhaseCall`
*copies* taken during resolution. `set_allele_map` writes `emitted` into the live `Entry`. With the
record now built after the decision, every copy says unemitted forever, so the patch-index build
(`if (!pc.emitted) continue`) skipped all 70 fixture records and reported `phasing: 0 sites phased` —
phasing absent, not mis-numbered. The same snapshot breaks the mosaic in the other direction: built
during resolution it would have admitted the ~100k sites that never become records. Both moved into a
new `finalise_linkage_outputs()` that runs after the render and reads the emitted set live, in one
pass (`emitted_records()`) rather than a scan per key — 219,600 entries against 219,600 phase calls is
not a lookup that can be answered individually.

**3. `--no-phased` is not "the same run without PS".** Where the linkage layer runs, turning phasing
off turns *nested calling* off with it (`src/subcommand/call_main.cpp:1810`, deliberate and documented:
a nested site's ploidy comes from its parent's phased genotype). So `test/t/18_vg_call.t`'s
permutation check compared a nested run against a non-nested one — valid only while the two emitted
the same records, which decide-then-render ends: 70 records against 63. The test was rewritten to put
`--no-nested` on both sides, where the invariant holds exactly (63 = 63, all joining, zero
non-permutations), and to **join on the snarl ID instead of pairing by line order**, because POS is not
an identity: one snarl legitimately moved from POS 10 to POS 9 between arms, and `paste` by line
number then reported 41 of 70 records broken. This is the same class as stage 10b's locus bug and the
third position-as-identity error in this phase.

This also confirms stage 11's open item (b) from the other direction: the configuration where the
collector is armed but nesting is off is live, and it is reached by `--no-phased`, not only by
`--no-nested`.

**Tests.** `test/t/18_vg_call.t` at 304 (was 303: the rewritten check gained a record-set assertion so
a join cannot silently drop rows). `vg test` 835 cases, 12,547,453 assertions.

## 10b. The residue was one defect, not two, and the fix is subtraction

Stage 10's gate items (i) and (ii) did not reach zero, and the reason turned out to be simpler than
the "partial pass" note above suggested. **Stage 10 converted only the top-level branch of
`call_snarl_internal` to staging.** The nested branch still called `emit_variant` inline during the
sweep (`src/graph_caller.cpp:4886`), and the barrier then called it a *second* time to replace that
line whenever the parent-implied ploidy differed.

So a nested chain's line was written twice, both times before its own generation resolved. An
independent audit established the ordering precisely: the selector `pr.generation != gen + 1` means a
generation-*g* chain is emitted during barrier iteration *g−1*, and `resolve_generation` skips
`entries[i].generation > generation`, so the site is *provably* unresolved when its line is written.

Once a line exists the only way to change it is a patch, and a patch has exactly two things it cannot
do. It cannot add an ALT, and it cannot withdraw a line. Those two impossibilities **are** gate items
(i) and (ii):

| | count | mechanism |
|---|---|---|
| (i) unrenderable | 496 | the settled traversal has no ALT on the already-written line, so the genotype is left as called |
| (ii) hom-ref | 1,383 | the settled pair maps to VCF allele 0 in both slots, so `apply_linkage_change` writes `0/0` and the line stays |

Both `++unrenderable` sites sit behind the same guard — `if (!e.emitted) continue;`
(`src/linkage_model.cpp:1771`, `:2292`) — so the counter can only fire for a site whose line already
exists. **It is therefore not a defect to be fixed but a count of how many records still take the old
path**, and it is structurally zero the moment nothing is written before the barrier. Item (ii)
follows the same way: the ALT list is built by iterating the genotype handed in
(`src/graph_caller.cpp:2262`), so an all-reference genotype leaves `alt` empty and
`wants_line = genotype_snarls || !out_variant.alt.empty()` is false. **Emission cannot write a `0/0`
line at all** — every one of the 1,383 was created by a patch.

**Two claims in the stage 10 note above are wrong and are corrected here.** The barrier's emit is not
"load-bearing twice because `respecify` needs the post-flatten position and allele map that only
emitting produces": `record_site` already passes a *pre-flatten* position and `no_allele_map`, with
`set_allele_map` supplying the map at render time, and its own comment records why that is sound —
"what the model uses position for is ordering and the transition gaps". Both are therefore deferrable
at the barrier exactly as they are at the sweep. And the residue is not two defects needing separate
work; it is one, and the fix removes code rather than adding it.

**The change.** The nested branch stages like the top-level branch; the barrier revises the staged
record in place (`pr.genotype`, `pr.ploidy`, `pr.call_info`) and respecifies the entry with the
pre-flatten key and no allele map; every surviving chain is handed to `render_records` after the
generation loop; one render pass writes every line, once, from the settled genotype. Deleted with it:
the barrier's `emit_variant` call, the `blank_buffered_line` replacement dance, the buffer handles on
a nested record, and the `landed` bookkeeping.

**One trap avoided, worth recording.** `crossing_known` is seeded from `parent_alleles.valid`, which
is gated on `emitted_this_call` — so stopping nested emission looked like it would make every
grandchild `crossing_unknown` and silently exempt it from revision. It does not:
`child_crossing_mask` sets `*known = true` unconditionally on entry and only clears it above 64
traversals, so the seed is overwritten and has no effect. The mask has been in *traversal* space on
both sides since stage 1; the `parent_alleles.valid` coupling is vestigial from when it was in
emitted-allele space, and the comment above it is stale.

### 10b result: both gate items zero, accuracy flat

chr20, against stage 10:

| gate | stage 9 | stage 10 | 10b |
|---|---|---|---|
| (i) unrenderable settled genotypes | 1,472 | 496 | **0** |
| (ii) final GT of 0/0 or 0\|0 | 4,490 | 1,383 | **0** |
| (v) GT past the ALT list | 0 | 0 | 0 |
| unphased records | 116,952 | 144 | **0** |
| records | 116,966 | 115,618 | 115,038 |

Both items reach zero, which was stage 10's gate as written, and every record is now phased.

| class | stage 9 | stage 10 | 10b | vs 10 |
|---|---|---|---|---|
| ALL | 0.97048 | 0.97231 | 0.97222 | −0.00009 |
| SNV | 0.98436 | 0.98523 | 0.98525 | +0.00002 |
| Insertion | 0.90843 | 0.91504 | 0.91457 | −0.00047 |
| Deletion | 0.93266 | 0.93724 | 0.93677 | −0.00047 |
| JointIndel | 0.91840 | 0.92390 | 0.92338 | −0.00052 |

Accuracy is flat: TP +16, FN −16, FP +35, ALL F1 down 9e-5. Stage 10's gain survives in full
(+0.00174 against stage 9). The scoring is deterministic on a fixed VCF, so the sign is real rather
than noise, but the magnitude is 0.009% and the change was made for coherence, not for accuracy.

**One thing to look at later, recorded rather than chased.** The 496 previously-unrenderable sites now
render on their *settled* genotype instead of keeping the called one, and that trade came out +16 TP
against +35 FP. At top level the same substitution was the source of stage 10's whole gain, so the
settled genotype being *worse* on average at exactly the sites where the patch could not express it is
a real asymmetry, concentrated in the indel classes (Insertion FP 878→887, Deletion 669→679). Worth a
targeted comparison of called-vs-settled at those sites; it bears on stage 19's conditioning, which is
about making the nested genotype better rather than merely renderable.

**Counter semantics corrected in passing.** `was_gained` was `!pr.emitted`, true for every record once
nothing is emitted, which reported "0 revised, 2,950 gained". Reading `!has_entry()` instead restores
the 2,514 revised / 518 gained / 411 retracted split exactly, which is also a check that the change
moved no chain between those populations.

## Open decision 3, resolved: unify the emission path rather than scope the deletion

Stage 11's item (b) asked what to do about the live configuration where the collector is armed but
nesting is off, since the patch path is the only mechanism there. It is reachable by `--no-phased`,
not only `--no-nested`: `src/subcommand/call_main.cpp:1810` turns nested calling off wherever linkage
runs and phasing does not.

**Decision: make retention and rendering conditional on the collector being armed, not on nesting.**
The alternative -- scoping the deletion so `apply_linkage_change` survives for that one configuration
-- keeps two emission paths alive, which is the exact thing stage 11 exists to prevent ("so the two
cannot drift"). With staging armed wherever the collector is, `run_deferred_descent` in a non-nested
run simply resolves generation 0 over an empty pending set and the render pass writes every record
from the settled genotype, which is the same rule rather than a second one.

This changes `--no-nested` output and therefore needs its own gate, run on chr20 before the deletion
lands: hom-ref records to zero and unrenderable to zero there too, on the same reasoning that took
them to zero under nesting. It is not a free deletion and is not treated as one.

After 10b the patch path is already dead in the default configuration -- every generation reports 0
genotypes changed, because a Change is only produced for a site with a line and no site has one at
resolution time. So `apply_linkage_change` is currently unreachable under `--read-likelihood --phased`
and reachable only in the configuration above. `apply_phasing` is still load-bearing everywhere: the
phase is applied by patching lines after the render, and moving PS and the phased separator into the
render is the second half of stage 11.

### 11a result: unifying the path improves the non-nested arm too

chr20, `--no-nested --phased`, decide-then-render against the old inline path:

| | old path | decide-then-render | delta |
|---|---|---|---|
| records | 105,251 | 103,904 | |
| hom-ref records | 3,131 | **0** | |
| unphased records | 416 | **0** | |
| unrenderable | 0 | 0 | |
| ALL F1 | 0.96468 | 0.96658 | **+0.00190** |
| SNV | 0.97809 | 0.97897 | +0.00088 |
| Insertion | 0.90580 | 0.91245 | +0.00665 |
| Deletion | 0.92814 | 0.93335 | +0.00521 |
| JointIndel | 0.91428 | 0.92009 | +0.00581 |

TP 89,897 → 90,210, FP 1,789 → 1,758, FN 4,794 → 4,481. Every class improves, and the shape matches
stage 10's gain at top level -- mostly recall, precision moving slightly the same way. So the
unification is not merely the precondition for the deletion; it is worth doing for this arm on its
own, and it makes the nested-versus-not comparison the prose rests on a fairer one, since both arms
now emit by the same mechanism and differ only in nesting.

The nested arm is untouched by this: where `nested_calling` is true the new condition is the old one.

## 11. Delete the patch machinery

**Goal.** Remove the path decide-then-render replaces, so the two cannot drift.

**Changes.** Delete `apply_linkage_change` (`src/graph_caller.cpp:852`+, ~210 lines, minus the GQ arithmetic stage 10 moved), `apply_phasing` (`:1270`+, moving its PS and phased-separator rendering into the render rather than dropping it), `blank_buffered_line` (`:3805-3813`), the tombstone branch in `write_variants` (`:781-790`), the `(contig, POS)` patch indices and their upsert loops, `change_declined`/`phase_declined` and their report (`:845-850`), `EmittedAlleles` (whose `num_alleles` is already write-only), the buffer handles on `PendingRecord`, and `child_crossing_mask` plus `NestedContext::parent_crossing`, `PendingRecord::parent_crossing`/`crossing_known`, `Entry::parent_crossing` (8 bytes a site) and the `crossing_unknown` counter — at decide time the parent's settled traversal is in hand and `crossings_of_child(travs[settled], child)` answers directly, with no 64-bit ceiling, so the >64-traversal hole closes here.

**Two things this stage must not do.** (a) `LinkageCollector::respecify` (`src/linkage_model.hpp:506-521`, `src/linkage_model.cpp:1319-1360`) is **not** patch machinery. It replaces an entry's likelihoods and ploidy — rebuilding the compact allele space and the GL vector at the parent-implied ploidy — and its header says "This is what makes the coherence guarantee structural rather than reported." Only its `(contig, position)` arguments belong to patching. Keep the ploidy/GL respecification under a name that says what it does; deleting it re-introduces the 440-record "both parent strands" class that `90d2bdce3` closed. Its four unit tests (`src/unittest/linkage_model.cpp:81-106` helper, `:432`, `:954`, `:979`) survive with it.

(b) **There is a live configuration that still needs the patch path.** `nested_calling` can be false (`src/subcommand/call_main.cpp:1735`, `:1744`, `:1835`) while the collector is armed by `linkage_weight > 0.0` alone (`:1848`), and `set_defer_nested_descent(true)` is reached only under `nested_calling` with a `FlowCaller` (`:1976-1980`). In that configuration nothing is retained, the render pass is empty, and `resolve_linkage()` inside `write_variants` plus `apply_linkage_change` is the only mechanism. Either make retention conditional on the collector being armed rather than on nesting — which contradicts stage 8's stated gating and must be decided there, not here — or scope this deletion to what only the deferred path used. **Resolve this before stage 9 is written** (open decision 3); it is load-bearing for the whole phase.

**Gate.** chr20 byte-identical to stage 10's. A deletion that changes output means something deleted was still live. `git diff --stat` net negative by at least 350 lines with fewer than 30 added. `vg test` (830) and `test/t/18_vg_call.t` (302) pass with no test deleted except the ones this stage names. Retracted-chain count identical to stage 10's: the mask and the direct crossing test must agree everywhere the mask was computable, so dump and compare per record rather than trusting totals — a disagreement is a defect in one of them. `crossing_unknown` reported non-zero on stage 10 and the mechanism gone here; report the count of newly-settled subtrees separately, because parents above 64 traversals now get their ploidy settled and their subtrees move, which means this stage is not purely mechanical after all.

**Tests.** No new tests for deleted code. Unit: derive a child's copies from a parent with 100 candidate traversals where the carrying one is index 80; fails on stage 10, where the mask returns 0 with `known == false`. And the four output-level invariants that replace backlog #66's nine barrier regression tests — written as invariants precisely because this stage deletes most of the paths those fixes lived in: no bare haploid GT, no unphased record, no nested site on both or neither parent strand, every nested site on exactly one strand. Confirm each fails against `a27149728`'s output, where all four are violated.

**Output moves:** no for the deletions, yes for the >64 subtrees. **Reversibility:** split into two commits — the pure deletions, then the mask replacement — so a regression can be attributed.

## 11b: the deletion, and a silent regression it uncovered

`apply_linkage_change` (216 lines), the `linkage_changes` index, the `Change` struct, the
`change_declined` counters and the `unrenderable` counter are gone. Net −132 lines across the two
files. Every generation reported **0 genotypes changed** on chr20 before the deletion, in both the
nested and the `--no-nested` arm, so the code was provably dead: a Change is only produced for a site
that already has a line, and no site has one at resolution time.

**The regression.** `apply_linkage_change` did two jobs, and only one of them stopped being
necessary. Alongside rewriting the GT it computed the post-linkage *quality*: GQ as the phred
complement of the HMM posterior, discounted by the explained-read share, capped at GQI, GQN blanked,
a stale `lowconf` cleared. When the genotype patch stopped being produced, that arithmetic stopped
running with it — and **nothing in the accuracy gate could see it**, because GQ is not used as a
filter here, so a run whose quality has silently reverted to the per-site value scores an identical
F1 and is differently calibrated. The counter told the story only in hindsight:

| | records carrying posterior-derived quality |
|---|---|
| stage 9 | 9,980 |
| stage 10 | 3,524 |
| 10b | **0** |

The cap at GQI is the load-bearing part — about +0.003 AUC and 1–2% fewer surviving false calls,
against +0.0001 to +0.0009 for the share discount alone — and it is what makes `GQ <= GQI` hold at
all. This is exactly the item stage 10 named as "the part no candidate plan scheduled", and it was
missed there and again in 10b.

Restored as `apply_linkage_quality`, driven by record key from a new
`LinkageCollector::moved_quality()` (record key → posterior, explained share) filled at the two places
the model moves a genotype. It touches GQ, GQN and FILTER only; the genotype is not its business any
more, because the line already carries the settled one.

**Deliberately not widened.** The posterior now exists for every settled site, not just moved ones, so
posterior-derived quality *could* apply to all 219,600. That is a different change with its own
measurement (AUC, and surviving false calls at matched recall) and it is not made here — the
restoration reproduces exactly the population the arithmetic used to cover.

**Counter renamed to match what it counts.** "N genotypes changed" was the number of patches
produced, which now reads 0 at every generation — a counter reporting that linkage changed nothing,
about the pass that decides every genotype in the output. It is "N genotypes moved by linkage" and is
incremented where the model's answer differs from the called one, whether or not anything is patched.

**Verified.** chr20 after the deletion: 115,038 records, unrenderable 0, hom-ref 0, GT past ALT 0,
unphased 0, ALL F1 0.97222 with identical TP/FP/FN — the VCF differs from the previous commit's only
in GQ and GQN, which is the repair. The quality pass is live: 15,068 genotypes moved, GQ changed on
5,910 records, GQN on 9,504, `GQ > GQI` 0, and the cap visibly acting in both directions (2 → 6 capped
at GQI 6; 43 → 19). `vg test` 835 cases / 12,547,462 assertions; TAP 304.

**Still to delete, and why not yet.** `apply_phasing` is the last patch pass. It is load-bearing:
the phase is applied by rewriting a rendered line's GT separator and appending PS. Moving it into the
render is the remaining half of stage 11 and needs the phase lookup built between the barrier and the
render rather than after it — which is possible, because the `emitted` filter that forced
`finalise_linkage_outputs` to run late exists for the mosaic, not for phasing: the render knows
whether it is writing a line. `render_phase_pair`, `vcf_allele_of_traversal` and the whole
`phase_fallback` population go with it, since the render knows the traversal→allele map it has just
built.

### 11c result: phasing moved into the render, output byte-identical

`apply_phasing` (161 lines), the phasing patch in `write_variants` (30), the `linkage_phasings` index,
`phase_declined_allele` and `vcf_allele_of_traversal` are gone. **chr20 output is byte-identical to the
patch-based run** — the right gate for a refactor whose whole claim is that it changes nothing but
where the work happens.

Phasing was a patch for exactly one reason: a PhaseCall names a *traversal* pair, the VCF needs allele
numbers, and the map between them did not exist until the record was built. Rendered instead, that map
is `trav_to_allele`, complete and correct, a few lines above in the same function. No fallback is
needed and none is taken: **zero phase refusals of any kind** on chr20, against 144 unphased records at
stage 10 and a 192,045-site `render_phase_pair` fallback population.

Two mistakes on the way, neither visible on the 70-record fixture, both worth recording:

- `trav_to_allele` is a `std::map` keyed *by* traversal, not a vector indexed by it. Bounds-checking a
  traversal index against `map::size()` — the number of alleles the record carries — refused 101,947
  phases and left 7,189 records unphased. `operator[]` would have been worse than the wrong bound: on a
  miss it inserts a default 0, writing allele 0 for a traversal the record does not carry into the very
  map `set_allele_map` is then handed.
- Adding PS before `update_vcf_info` put it second in FORMAT where the patch had appended it last.
  Same fields, same values, every line different — which costs the byte comparison for no gain.

`build_render_phases()` fills the lookup between the barrier and the render, deliberately with **no
`emitted` filter**, unlike the mosaic. That filter is what forced the bookkeeping to run after the
render in the first place; a render-time lookup does not need it, because a site with no line never
looks itself up.

**Left standing on purpose.** `render_phase_pair` and `PhaseCall::allele_first/allele_second` are now
written and never read, and the misleading `phase_fallback` report is removed, but the conversion
itself stays: `allele_*` carries the *compact* pair before that block overwrites it with the VCF one,
so removing it needs a rename rather than a deletion, and doing both at once inside the phasing path
was not worth the risk in the same commit that moved it.

---

# Whole-genome result

24 contigs under `schedule_wgs.py` into `work/wgs-dtr`, both prior arms left intact. 4,970,004
records. Scored through the same `bench_wgs.py` path as the arms it is compared against.

**Autosomes**, summed counts (averaging per-contig F1s would weight chr21 like chr1):

| class | inline | post-linkage descent | decide-then-render | delta vs inline |
|---|---|---|---|---|
| ALL | 0.97034 | 0.97032 | **0.97288** | +0.00254 |
| SNV | 0.98372 | 0.98371 | **0.98492** | +0.00119 |
| Insertion | 0.91046 | 0.91042 | **0.91811** | +0.00765 |
| Deletion | 0.93348 | 0.93350 | **0.94105** | +0.00757 |
| JointIndel | 0.91948 | 0.91946 | **0.92720** | +0.00772 |

TP 4,025,074 → 4,039,888 (+14,814), FP 94,189 → 88,172 (−6,017), FN 151,914 → 137,100 (−14,814).
**Both** precision and recall improve, which chr20 alone did not show: there the gain was recall-only
with FP nearly flat. Genome-wide FP falls 6.4% and FN 9.8%.

Every one of the 22 autosomes improves. The largest gains are chr16 (+0.00543), chr15 (+0.00497),
chr2 (+0.00396) and chr19 (+0.00386) — the contigs with the most segmental duplication, which is
where collapsed parents and nested children are most common.

**Invariants hold on all 24 contigs**: zero records with a hom-ref genotype, zero unphased records,
zero genotypes past the ALT list. Peak RSS tops out at 6.05 GB on chr1.

## Runtime and read I/O

The read I/O goal is met exactly. Eliminating the five-sweep penalty was the entire reason for the
single-sweep design, and it is gone to within 0.1%:

| | inline | post-linkage descent | decide-then-render |
|---|---|---|---|
| reads fetched | 609.1 M | 903.3 M (+48.8%) | **609.7 M (1.001x)** |
| CPU user | 6.26 h | 11.45 h | **9.96 h (+59%)** |
| wall clock | 3.20 h | 4.89 h | 6.03 h |
| linkage pass, serial | 1,064 s | 950 s | 1,762 s |

Wall clock is not comparable across separately scheduled runs -- the scheduler packs contigs
concurrently and the contention differs -- so **CPU user is the number to read: +59% against inline**,
and against the deferred arm it is a 13% saving with the read penalty removed as well.

The linkage pass accounts for only ~0.2 h of the +3.7 h. The rest is the sweep and the render doing
more work per site, and the layer holding more sites.

**The per-contig spread points at the same cause as the accuracy regression.** Slowest ratios against
inline: chrX **4.66x**, chr1 2.36x, chr4 2.26x, chr3 2.29x, chr2 2.02x, chr7 1.40x. chrX is the worst
by a wide margin and is exactly the contig with 7.4x site inflation.

So the excess line-less sites in the layer cost **both** accuracy (on the haploid contigs) and
runtime (everywhere, worst where they are densest). That reframes the fix: a site is recorded so its
nested children can inherit a strand, but **a site with no children does not need to be in the layer
at all**, and on chrX almost none of the 700,000 extra sites can have any. Filtering on "has a nested
child that needs this strand" would cut the density inflation, restore the transition model to the
density it was tuned at, and remove most of the added CPU -- one change against three symptoms, and
it removes the cost rather than working around it.

The information is available at the right moment: `record_site` runs during the sweep, and descent
runs immediately after in the same call, so whether any child was staged is known before the record
has to be final -- the same "stage now, complete after descent" discipline `travs` already uses.

## The chrX regression: two haploid bugs, both fixed

**Resolved.** Two bugs in the haploid path, one of them introduced by this phase's own depth->=2
strand fix. Neither had anything to do with linkage, symbolic collapsing, or site density -- all three
of which were proposed here first and are wrong. What found them was elimination: the deficit is flat
across `--linkage-weight` 0.01 / 1 / 2 (TP 82,435 / 82,697 / 82,689), so linkage cannot be
responsible; and only 39 of 4,247 lost records had a new record nested inside them, so collapsing was
not relocating the difference into a child. That left the rendering and the ploidy paths.

**Bug 1: a strand was assigned without requiring a diploid parent** (`src/linkage_model.cpp:2026`).
`strand = 1` was conditioned on `parent.ploidy == 2`; `strand = 0`, three lines above it, was not. A
haploid top-level site has ploidy 1 and no strand of its own, so the nested-parent branch misses and
the identity match fires unconditionally, handing its children strand 0. The renderer turns a strand
into `a|.`, which asserts the locus is diploid and the other strand carries nothing -- false on a
haploid contig, where the correct rendering is a bare `a`.

| chrX non-PAR | haploid-shaped | `a|.` |
|---|---|---|
| inline | 103,064 | 2,333 |
| before the fix | 93,444 | 8,056 |
| after | 102,669 | **0** |

The inline arm's 2,333 were wrong too, just fewer. chrY carried 10,994 of them and now carries none.

**Bug 2: an unreached child was descended at a hard-coded ploidy 2** (`src/graph_caller.cpp:4879`).
`child_ploidy` returns `min(copies, ploidy)`, so the fallback exceeded the bound the same call had
just applied, and on `-d 1` it named a ploidy the contig does not have. Output-neutral in fact --
`set_want_alt_ploidy` keeps the other ploidy's answer and the barrier respecifies anyway, so the
initial choice only decides which is primary -- and kept as a correctness fix rather than a measured
one.

**Result.** chrX 0.94939 -> **0.95666**, so it goes from the one regressing contig to +0.00727 against
inline. With both fixes in, **all 23 scoreable contigs improve**:

| autosomes + chrX | inline | decide-then-render | delta |
|---|---|---|---|
| ALL | 0.96989 | **0.97253** | +0.00264 |
| SNV | 0.98329 | 0.98460 | +0.00130 |
| Insertion | 0.91023 | 0.91788 | +0.00765 |
| Deletion | 0.93335 | 0.94091 | +0.00757 |
| JointIndel | 0.91909 | 0.92687 | +0.00778 |

TP +15,316, FP −6,844, FN −15,316. chr20 is byte-identical across both fixes, so the autosomes are
provably untouched by them.

**A test gap this exposed.** The existing depth->=2 strand test could never have caught bug 1: its
parent is a het *diploid* site, so the strand it hands down is real. A haploid-parent case is now
asserted and confirmed to discriminate -- it fails against the restored bug and passes with the fix.

**And a claim retracted.** The `pow(rho, weight)` density dependence documented below is real
arithmetic and remains worth fixing, but it is NOT the cause of anything measured here: at weight 1,
where the transition composes correctly, chrX scored 0.93726 against weight 2's 0.93643. It is a
latent mis-parameterisation, not this regression.

## chrX regresses, and the mechanism is identifiable

**chrX: 0.94939 → 0.93643, −0.01296**, FN +3,143 against FP −1,067. It is the only contig that gets
worse, and it is worth more than the note it would get as an outlier.

Including chrX the whole-genome figure is still strongly positive — ALL 0.96989 → 0.97212 (+0.00223),
TP +11,671, FP −7,084, FN −11,671 — so this is a real cost inside a larger gain, not a wash.

What was lost is specific: of 7,566 chrX records that disappeared, **3,461 carried a bare haploid
`1`** and 1,858 a bare `0` (the latter correctly, being reference). 4,021 were gained. The loss
concentrates in non-PAR: 5,892 lost against 2,076 gained, versus 1,674/1,945 in the PAR.

The mechanism is the **site-to-record ratio at ploidy 1**:

| | sites in the layer | records | ratio |
|---|---|---|---|
| chr20 | 219,600 | 115,038 | 1.91 |
| chr1 | 662,724 | 369,226 | 1.79 |
| **chrX** | 816,291 | 110,804 | **7.37** |
| **chrY** | 354,607 | 64,554 | **5.49** |

The ratio itself is expected and is not the defect: a haploid site needs only *one* reference allele
to collapse and write no line, where a diploid site needs two, so P(no line) is p rather than p², and
haploid regions therefore produce far more line-less sites. What is not expected is the consequence.
Those sites are recorded so that nested children can inherit a strand — but they also **enter the
linkage chain**, and on a haploid contig that chain is a single-strand HMM in which they now outnumber
real records 7:1. Their emissions and their transition gaps steer the Viterbi path for the records
that do get emitted. chrX's linkage moved 15,116 genotypes against the prior arm's 6,170.

**Correction: it does not predate this phase.** The claim above was that recording collapsed sites
arrived with the older stage 2, and the deferred arm disproves it — `wgs-defer` has **114,383** sites
on chrX against inline's 114,412, and scores 0.94931 against 0.94939. The 7x inflation arrived in
*this* plan's work, and the regression is ours.

**Correction: the first diagnostic was misconfigured and its result should be ignored.** A chrX probe
at `--linkage-weight 1` came out at 0.89823, which looked like a refutation of the density
explanation. It was run with `-d 2` where `call_wgs.sh` uses `-d 1` plus the PAR bed, so it changed
base ploidy as well as the weight and measured neither. Both whole-genome arms did use the same
configuration — `chrX: haploid with diploid PAR (--ploidy-bed)` appears in both schedule logs — so the
comparison itself is sound; only the probe was wrong.

**What is established.** The lost calls are not relocated. Of 4,247 chrX records with a non-reference
genotype that disappeared, only **39** have a new record nested inside them; **4,208 were lost with
nothing taking their place**. So this is not symbolic collapsing moving a difference down into a
child — it is the settled genotype changing from alt to reference. A decision change, not a
representation change, which puts the linkage model at the centre.

**And the obvious explanation is wrong.** A collapsed site is not an uninformative site: its panel row
is over *traversals*, and those traversals genuinely differ from each other — they differ inside the
child chains. Such a site says which panel haplotype the sample follows, which is exactly the evidence
linkage wants. Adding informative sites to a chain should not degrade it.

## The transition model's switch rate depends on site density

`switch_probability` computes `rho = rho_min + (1-rho_min)(1-exp(-gap/scale))` and returns
`pow(rho, weight)`, with `weight` defaulting to 2. **Raising a per-step probability to a power breaks
composition**: subdividing a span cannot preserve its total switch probability.

Measured on chrX's mean inter-site gap, which fell from 1,348 bp to 189 bp:

| | one step of 1,348 bp | 7.1 steps of 189 bp |
|---|---|---|
| weight 1 | switch 0.12701 | 0.13235 (0.96x — telescopes, as it should) |
| **weight 2** | switch 0.01613 | **0.00277 (5.83x stickier)** |

At weight 1 the distance term telescopes exactly and only the `rho_min` floor accumulates. At weight 2
it does not telescope at all, so **the effective switch rate is a function of how many sites are in
the chain rather than of genomic distance.** Any change to what enters the layer silently re-tunes the
linkage strength: `weight = 2` was fitted at ~114k sites on chrX and is now applied at 816k.

This is a defect independent of chrX, and it means the autosomes were re-tuned too — 1.9x denser, so
mildly stickier, which evidently helped there.

**Whether it explains the chrX loss is still open** and the test is a probe at weight 1 with the
correct `-d 1 --ploidy-bed`, where the model is density-independent. The candidate mechanism is
rigidity: with the per-step switch probability 5.8x smaller, the Viterbi path is pinned by whichever
haplotype pair was favoured earlier and cannot switch when a later site disagrees, so real haploid alt
calls get pulled to the pinned haplotype's reference allele. That is consistent with chrX's linkage
moving 15,116 genotypes against 6,170, but consistency is not evidence and the probe decides it.

The principled fix, if confirmed, is to make `weight` act on the survival probability so that it
composes: `1 - rho = pow((1-rho_min) * exp(-gap/scale), 1/weight)`, for which N steps multiply to
exactly the one-step answer over the same span. It changes what the number 2 means, so it needs a
re-fit — which the deferred list already anticipates.

**The question to answer, and not by guessing:** should a site that will produce no VCF line
participate in the linkage chain at all, or only in the strand inheritance its children need? Those
are separable — a site can be in the collector for its children without being a node in its parent
chain's HMM. The measurement is a chrX arm with line-less sites excluded from the chain but retained
for strand, against this one. It bears directly on stage 15, which rebuilds the chain per haplotype,
and on stage 19's conditioning; doing it there rather than as a patch here is the right sequencing,
because stage 15 changes what a chain *is*.

Not deferred out of convenience: chrX is ~5% of the genome and the loss is 3,143 true calls, against
14,814 recovered elsewhere.

---

# Phase II complete

Stages 1–12 are landed and gated. The invariant the phase was built for holds by construction rather
than by inspection: **nothing is written before it is decided.** There is one emission path, one
genotype behind every record, and no patch machinery at all.

| chr20 | stage 9 | now |
|---|---|---|
| unrenderable settled genotypes | 1,472 | **0** |
| records with GT 0/0 or 0\|0 | 4,490 | **0** |
| unphased records | — | **0** |
| GT past the ALT list | 0 | 0 |
| records | 116,966 | 115,038 |
| ALL F1 | 0.97048 | **0.97222** |
| SNV | 0.98436 | 0.98525 |
| Insertion | 0.90843 | 0.91457 |
| Deletion | 0.93266 | 0.93677 |
| JointIndel | 0.91840 | 0.92338 |

`--no-nested`, which the unification brought onto the same path, moved 0.96468 → 0.96658 with hom-ref
records 3,131 → 0 and unphased 416 → 0.

## 12. Documentation for phase II

**Changes.** `doc/read-likelihood-genotyping.md:564-565` still documents the three deleted `nested_*` FILTERs and how to interpret them (verified stale today, before any of this work). Rewrite that section around decide-then-render: the genotype is settled before the record exists, so there is nothing to flag. Also stale and not in any candidate plan's scope: `:210-211` ("linkage re-decides genotypes afterwards" — the architecture stage 10 inverts), `:370` (the merge header asserting DP/QUAL/GQ/FILTER are computed over the pre-merge allele set), `:389-435` and `:433` (mosaic format and column table, changed by stage 2), `:646` and `:779` (GQN blanking and `lowconf` clearing, whose owner moved in stage 10). In the eval repo: `docs/nested-calling-design.md` (eleven `nested_*` references including two results tables), `docs/coverage.md:233` (`apply_linkage_change`), `planning/nested-traversal-space.md:275`. Mark superseded rather than deleting, per that repo's convention.

**Gate.** No `nested_diploid`/`nested_haploid`/`nested_unreachable`/`respecify`/`apply_linkage_change` string survives outside a section marked as history — grep. Every source line cited in the rewritten sections resolves in the tree at this head, checked by a script in the eval repo so it can be re-run; confirm it fails by running it against `4371c9b67`, where `:564-565` cites FILTERs the source no longer emits. The four eval documents quoting pre-refactor whole-genome figures are **not** regenerated here — they need run 2, and updating them from chr20 would be worse than leaving them stale.

---

## 7-8 result: retention priced exactly, route decided

**407.25 MB on chr20**, walked from the retained objects rather than estimated: 222,623 records,
5,097,064 traversal visits, 991,557 genotype likelihoods, ~1,830 B a record, 14x the collector's own
29 MB arena. Landed byte-identical as `d151dc3d7`.

Stage 7's estimate-arbitration was abandoned in favour of measuring, for a reason worth keeping:
**peak RSS cannot resolve a delta this size on one contig.** Six runs of a single binary spread
3.39-4.42 GB, wider than the retention being priced and wider than the documented 0.7 GB noise floor.
The two prior estimates were also wrong in opposite directions -- 297 MB (27% low) and ~620 MB (52%
high) -- so a third estimate would have added a guess, not an answer.

**Route: retain unconditionally.** Projected by record count, chr1 (3.02x chr20) is ~1.23 GB, taking
its measured 5.7 GB peak to ~6.9 GB and chr3's 6.1 GB worst case to ~7.3 GB. Against the documented
32 GB machine packing contigs under a budget, that is about four concurrent contigs instead of five:
a throughput cost near 20%, not a feasibility limit. And it is transient -- stage 10 renders and can
release, so the peak spans only the sweep-plus-barrier window. The plan's 7.0 GB threshold against a
"24 GB budget" is not used; that budget figure is undocumented, and the real constraint is packing
density under 32 GB.

Two corrections to the stage-8 text itself:

* The retained population is **not** "top-level snarls" and the gate's "retained top-level count
  equals the top-level snarl count" is unsatisfiable. It is snarls reaching an emit branch with no
  ploidy override: 165,408 top-level plus 26,799 children that `call_top_level_snarls` reaches by
  recursing on failure. Those belong in the render container -- the barrier does not revise them
  either.
* Staging in "all three emit branches" is right, but the `parent_child_trav_sets` branch is
  unreachable on the default path: `nested` is set only by `-A`, `--top-down` or `--bottom-up`, never
  by `--nested`.

**And the move point has exactly one right answer, of three candidates.** Moving `travs` into the
record at emit time costs 12,302 chr20 records, because descent runs after every emit branch and
reads `travs` to see which children the called alleles reach -- the same failure as `906812957`, five
times larger. Moving it after symbolic descent, where the nested branch completes its own staging,
breaks four `-A` and `--top-down` tests, because the `-A` recursion builds each child's
`ChildTraversalSets` from `travs[allele_idx]` further down. Only after that recursion is correct. The
nested branch survives completing early only because `pending_this` is never set on the `-A` path --
correct by accident of configuration, not by construction. Anything Phase II moves near here needs
that ordering checked rather than inherited.

# Phase III — haplotype-frame linkage

Read the code comment at `src/linkage_model.cpp:1548-1556` before starting this phase. It records that letting nested sites into the diploid runs took chr20 from 22 phase blocks to 9,460 and N50 from 248 Mb to 1.08 Mb, "and the switch rate only looked flat because short blocks make switch error cheap". Unifying depths is exactly that move. Every gate below carries a phase-block condition for that reason, and phase III is the part of this plan most likely to be reverted.

## 13. Per-strand transition, called with one value

**Goal.** Land the arithmetic generalisation with both strands given the same value, so it is reviewable against a byte-identical gate before any distance changes.

**Changes.** `src/linkage_model.cpp:161-183` (`transition_apply`) takes one scalar `rho` and forms `stay*stay*in[a][b] + stay*jump*(row[a]+col[b]) + jump*jump*total`, applied symmetrically to the ordered pair; `viterbi_step` (`:234-300`) derives one `S`/`J` pair for both slots. Generalise to `(rho_a, rho_b)`: `stay_a*stay_b*in[a][b] + stay_a*jump_b*row[a] + jump_a*stay_b*col[b] + jump_a*jump_b*total`, same O(m²). `viterbi_step` is the harder half because the leave-one-out maxima (`Top2`, `:180-210`) are already per-axis, so the change is confined to how the four candidate terms are scored. Route all six `switch_probability` call sites (`:353`, `:489`, `:728`, `:853`, `:935`, `:1046`) through the pair form, passing the same value twice; the last three are the haploid path and are equally affected.

**Gate, stated once so it cannot be read two ways.** With `rho_a == rho_b` the new expressions reduce algebraically to the old ones, but the existing grouping `stay * jump * (row[a] + col[b])` **cannot** be preserved once the two coefficients differ — it must split into two differently-coefficiented terms, which is a guaranteed re-association. So byte-identity is not achievable and claiming it produces a gate that gets waived. The gate is: at most 5 differing genotypes on chr20, each individually characterised as a near-tie whose two candidates differ by less than 1e-9 in log space; more than 5, or any difference not of that form, means the algebra and not the rounding. Plus: `vg test` passes with all 28 existing `[linkage_model]` cases unmodified — they are the real gate here.

**Tests.** Unit: `transition_apply` with `(rho, rho)` against a test-local reference implementation of the old single-rho form on a random state vector, asserting agreement to 1e-12 rather than bit-identity. Confirm the test discriminates by perturbing one term's factorisation. Unit: with `rho_a = 0` and `rho_b = 1`, the result equals the closed form (first strand unchanged, second uniform); fails today because the function cannot express two rhos, so the test does not compile.

**Output moves:** marginally. **Reversibility:** two functions, mechanical.

## 14. Record haplotype-frame coordinates, consume nothing, and decide whether the frame matters

**Goal.** Supply the distance the frame needs, and give the phase an off-ramp before its riskiest stage.

**Changes.** At descent (`src/graph_caller.cpp:4693-4740`) the parent's `travs` and the child's snarl are both in hand, so the bp offset of the child's start along each carrying parent traversal and the child's own span are computable there — per traversal, because a diploid parent's two traversals give two offsets. `Entry` holds no lengths (verified: position, contig, arena offsets, `num_alleles`, `called_i/j`, ploidy, generation, `retracted`, `final_i/j`, `emitted`, `nested`, `parent_record_key`, `parent_crossing`, `parent_trav`, `explained_share`, `record_key`, `start_node`, `end_node`), so this is a new parameter on `record()` — already a 17-argument function called from the parallel hot path — plus three fields. For a top-level site the offset is the reference position and the span is the snarl's reference length, so the frame is uniform at every depth. Compute the offsets in `graph_caller.cpp` and hand them in as plain integers: `LinkageCollector` is deliberately free of graph types, and losing that property makes the unit tests unwritable without a fixture graph.

Then instrument, inert: for every adjacent site pair on chr20, the reference-POS gap, the haplotype-frame distance, the ratio distribution, how many adjacent pairs reorder under the haplotype frame, how many reorder under the consensus rule (children visited, then bp length, then start node id) when a diploid parent's two traversals disagree about child order, and the distribution of `|offset_first − offset_second|`. Also: how many of the 12,516 no-reference-path children become orderable.

**Gate.** chr20 byte-identical to stage 13's — nothing reads the new fields. Invariants, each printed and each able to fail: every nested entry's offset lies inside its parent's span; a parent's children's spans do not exceed the parent's span on the carrying traversal. Do **not** gate "zero entries with no computable offset": the retained population (`copies == 0`, `g_descent_skipped_no_copy`, ~296 on chr20) is descended into and has `carrying_trav = -1` because that is set only under `copies == 1` (`src/graph_caller.cpp:4664-4677`), so there is no called parent traversal to measure along. Count them instead. Arena growth ~12 bytes a site, ~2.6 MB on chr20's 219,246 sites, reported by `bytes()`, netted against stage 11's 8-byte-a-site saving. CPU (not wall) within +3% of stage 13: summing node lengths along a traversal is graph work inside the parallel sweep, and if it is expensive it shows as a regression with no output change.

**The off-ramp, stated before the measurement.** If fewer than 1% of adjacent pairs reorder and 99% of gap ratios sit inside 1.05, the frame change buys nothing measurable and stage 15's distance half is dropped — the value of phase III is then entirely in the pooling fixes and stage 16. Only stating this first makes it takeable.

**Tests.** Unit: a parent spanning 1,000 bp with two children at offsets 100 and 600 on one traversal and 100 and 400 on the other; assert both offsets round-trip and the disagreement is reported. Fails today: `record()` has no such parameters. Unit: a child whose offset exceeds its parent's span is rejected and counted, not stored; confirm by passing one deliberately.

**Output moves:** no. **Reversibility:** additive.

## 13-14 result: the arithmetic lands byte-identical, and the frame matters for distance only

**Stage 13** generalised `transition_apply` and `viterbi_step` to one switch probability per strand,
both given the same value. `viterbi_step` needed less than the plan expected: its four candidates
already *are* the four stay/jump combinations, so each takes its own axis's coefficient and the
leave-one-out top-2 maxima do not change at all.

Better than its own gate: the gate allowed 5 differing genotypes on chr20, each to be characterised as
a sub-1e-9 near-tie, on the reasoning that the re-association is unavoidable. It is unavoidable, and
chr20 still comes out **byte-identical** -- it falls below printed precision everywhere. The gate was
right to allow for it and right that a waived byte-identity claim would have been worthless.

Which test does the work is not the obvious one. The reduction test -- pair form with `(rho, rho)`
against a local copy of the old arithmetic -- **cannot** catch a swap between the two stay factors,
because passing one value makes them equal. Perturbing `stay_a * jump_b * row[a]` to
`stay_b * jump_b * row[a]` is caught only by the asymmetric test (`rho_a = 0`, `rho_b = 1` against the
closed form), which is also the test that could not be written before this change.

**Stage 14** measured the frame instead of assuming it. chr20:

| | |
|---|---|
| adjacent sibling pairs | 25,098 |
| reorder in the haplotype frame | **152 (0.606%)** |
| gaps within 1.05 of the reference gap | **91.05%** |
| children whose two parent traversals disagree on offset | 11,035 (mean 1,699 bp, max 55,751) |
| children with no offset on a called traversal | 2,714 |

The off-ramp required reordering under 1% **and** 99% of gaps inside 1.05. Only the first holds, so it
is not taken.

**The split reshapes stage 15 more usefully than either verdict would have.** Order barely changes;
distance changes materially for about 9% of adjacent pairs. So:

- **Keep** the per-haplotype distances. That is what stage 13 now makes expressible.
- **Drop the consensus ordering rule.** The plan reserved design effort for sorting by the longest
  allele's order with the shorter breaking ties, against the case where a diploid parent's two
  traversals disagree. At 152 pairs in 25,098 that rule cannot pay for its complexity: a simple
  deterministic order (reference-frame offset, then start node id) is indistinguishable in effect and
  cannot introduce an ordering bug of its own. Stage 15's risk budget goes to the distances.

Nothing is stored yet -- no `Entry` field, no `record()` parameter. Adding a 20th argument to a
function in the parallel hot path before knowing whether the answer justified it would have been the
wrong order, and the measurement did change the design.

## 15. One chain per haplotype, at every depth, with per-strand distances

**Goal.** Dissolve nesting as a separate linkage stage. This is the agreed design and the largest change in the plan.

**Changes.** Four things, and the plan is honest that the third is not plumbing.

(a) Fix the pooling defects first and gate them alone, because they do not depend on the frame: `by_strand` is keyed `(contig, strand)` (`src/linkage_model.cpp:1884`, filled `~:2011`), so a nested site links against same-depth sites under unrelated parents, and on a multi-chain contig two unrelated chains' strand 0 become one haploid chain. Key on `(phase_set, strand)`, which is on the parent's `PhaseCall` and is exactly the unit within which phase is comparable. And `if (idxs.size() < 2) continue;` (`~:2101-2103`) drops small groups entirely — not cosmetic, because `freq_prior` defaults to 5.0 (`src/subcommand/call_main.cpp:353`) and acts on a chain of one, so skipping changes the genotype. The diploid path already fixed this same defect for singleton chains and its comment records 258 chr20 sites going missing from the mosaic.

(b) Replace the reference-POS order with traversal order under the consensus rule, and the single gap with the per-strand distances from stage 14 through stage 13's pair form.

(c) **Fold nested sites into the parent's chain — and this requires a modelling change, not a sort-key change.** Nested entries are held out of `chainable` before runs are built (`src/linkage_model.cpp:1557-1571`, verified), so reordering alone is a no-op: there is nothing in the chains to reorder. Once folded, a chain mixes ploidy 1 and ploidy 2 sites, and `chain_ploidy = entries[indices.front()].ploidy` (`:1600`) selects the model for the whole chain at `:1644` and `:1729`. `Site::ploidy` exists but **`LinkageModel` never reads it** (verified), and `posteriors()` indexes `genotype_index(ai,bi)` against buffers sized n(n+1)/2 while a ploidy-1 `Site`'s `genotype_ln_likelihood` has only `num_alleles` entries — an out-of-bounds read. So this needs a real single-copy emission under the two-haplotype latent state, specified before it is written. Budget it as its own commit and its own review.

(d) Reconsider the chain-cutting rule (`:1553-1571`) only for the *nested* case: a nested ploidy-1 site is no longer a ploidy change inside a chain once there is one chain per haplotype. A *regional* ploidy change must still cut — across chrX's PAR boundary there is no haplotype correspondence to carry.

**Gate.** Baseline stage 14, structural before any F1, because a wrong consensus order or a wrong per-strand distance produces a well-formed VCF with slightly worse genotypes. (i) Phase blocks: chr20 autosomal count at 22 and N50 at 248 Mb. Say explicitly how these remain the same measurement after the hold-out is deleted — if "block" changes meaning, the gate compares two different quantities and must be re-derived from stage 14's run instead. (ii) Nested sites placed on exactly one strand: 6,716 of 6,716, re-measured on stage 14's output rather than assumed from `4371c9b67`, because stages 3, 4 and 10 all move that population. (iii) Determinism: two runs byte-identical (stage 1's instrument), and the consensus order independent of which parent traversal is enumerated first — assert by permuting the traversal list in a unit test. If the tiebreak is not total the sort is unstable and stage 1's gain is lost. (iv) Group counts before and after for (a); sites linked across a phase-set boundary → 0; nested sites in groups of one, and the genotype changes at those sites, counted — a non-zero count proves (a) was not cosmetic. (v) **Switch error** on chr20 against `docs/tier2-phasing.md`'s recorded 2.30% (34-haplotype) and 3.43% (4-haplotype), via `scripts/tier2/mosaic_switches.py`. This is the one number that can fall while every other gate here passes, and no candidate plan gated it. (vi) Accuracy on chr20 **and chr6**: ALL F1 not below stage 14's; JointIndel up by at least 0.0005 against 0.91840, since that is where the emission is flat and the panel has something to add. If JointIndel does not move, the frame change bought nothing measurable and (b) should be reverted rather than kept on principle. (vii) chrX for (a)'s two-unrelated-chains case, which chr20 cannot exercise. (viii) CPU and per-contig peak RSS against the 24 GB budget — one chain per haplotype per contig at every depth makes every site at every depth resident at once.

**Tests.** Unit: two nested sites on strand 0 with different parent phase sets do not influence each other's posterior; fails today (one chain, mutual influence). Unit: a single-site group receives a posterior differing from its raw likelihood when `freq_prior > 0`; fails today (skipped, so posterior equals likelihood exactly). Unit: a diploid parent whose two traversals visit children in opposite orders produces one deterministic sequence, unchanged under input permutation; fails today (no consensus rule). Unit: two strands with different deletion content across the same pair of sites receive different switch probabilities; fails today because the gap is a reference-POS difference identical for both strands, so the test asserts an inequality that cannot hold. Unit: a diploid chain with a nested ploidy-1 site in the middle is linked to both neighbours with phase continuous across it; fails today (the site is in `by_strand`, not the chain).

**Existing tests this invalidates:** `src/unittest/linkage_model.cpp:403` ("The collector sorts by reference position, not arrival order") and `:597` ("The phasing comes back in reference order even with nested sites in it") both encode the premise this stage changes. Re-baseline both here and say what replaces them.

**Output moves:** yes. **Reversibility:** four commits in the order (a), (b), (c), (d); tag stage 14 as the fallback. (c) is the one that is not cheap to revert.

## 15(b) needs a specification too, and here it is

The plan singled out (c) as "not plumbing" and needing a specification before it is written. (b) needs
one as well, and the reason only becomes visible once stage 13 exists.

**The problem.** The distance a haplotype has travelled between two adjacent sites depends on the
indel content it carries *between* them — which is exactly what the Viterbi is deciding. A transition
probability that depended on the source state would be a legitimate HMM, but it would need one `rho`
per source state, O(m) values a step, and that destroys the factorisation into row and column sums
that keeps `transition_apply` and `viterbi_step` at O(m^2). Stage 13's pair form takes two scalars for
precisely this reason: two is affordable, m is not.

**So (b) is only implementable as an approximation, and the choice has to be stated rather than
absorbed into the code.** The distances are derived from the *called* alleles and fixed before the
Viterbi runs:

    gap_a = ref_gap + (length of the called first allele at the source site - its reference length)
    gap_b = ref_gap + (length of the called second allele at the source site - its reference length)

Taken from the source site, since that is the sequence being traversed on the way out of it. This is
what a read-based phaser does — the distance is treated as a property of the assembly, not of the path
being scored — and it is exact wherever linkage does not move the genotype, which is 219,600 minus
15,068 sites on chr20, about 93%. Where linkage does move a site, the distance used is the one implied
by the per-site call rather than by the settled one.

**What this buys, bounded in advance.** Stage 14 measured 9% of adjacent gaps differing by more than
5% between frames. `scale` is 10 kb, so at typical gaps a 5-10% change in distance moves `rho` by
5-10%, and `weight = 2` squares that to 10-20% on a minority of steps. That is a perturbation to the
transition, not a new source of information, so a large accuracy gain would be surprising and should
be treated as suspicious rather than welcome.

**The kill criterion is the plan's own and is kept:** JointIndel must rise by at least 0.0005 on
chr20, or (b) is reverted rather than kept on principle. Stating the approximation first is what makes
that criterion meaningful — otherwise a null result is ambiguous between "the frame does not matter"
and "the approximation threw the signal away".

**Cost.** Two `int32` on `Entry` (~8 bytes a site, ~1.8 MB on chr20), two more arguments on `record()`,
and two fields on `Site`. The plan already flagged `record()` as a 17-argument function in the parallel
hot path; this takes it to 19.

## 15(a) and 15(b) results: one no-op kept, one revert

**15(a) landed with no measured effect, and half of it is unevidenced.** Keying `by_strand` on
(phase set, strand) is a real fix -- a strand only means something inside a phase set -- but chr20
**cannot exercise it**: PS is per chain, so on a single-chain contig the old and new keys are the same
key. chrX, which has a phase-set boundary at the PAR, shows 25 sites the old key would have pooled
across it. Both contigs are byte-identical and chrX's F1 is unchanged to five decimals. The singleton
half fires on **neither** contig (0 groups of one), so nothing measured supports it at all; it rests on
the argument that `freq_prior` acts on a chain of one. Kept as correctness, not as a measured gain, and
both counters print so the next contig that does exercise them says so.

**15(b) is REVERTED on its own criterion.** chr20 against 15(a):

| class | 15(a) | 15(b) | delta |
|---|---|---|---|
| ALL | 0.97231 | 0.97220 | −0.00011 |
| **JointIndel** | 0.92390 | 0.92329 | **−0.00061** |
| Insertion | 0.91504 | 0.91442 | −0.00062 |
| Deletion | 0.93724 | 0.93678 | −0.00046 |
| SNV | 0.98523 | 0.98526 | +0.00003 |

The criterion was JointIndel up by at least 0.0005 or revert rather than keep on principle. It fell.

**The negative result is interpretable, and that is what writing the approximation down first
bought.** Stage 14 had already bounded the upside: 9% of adjacent gaps differ by more than 5% between
frames, which at `scale` = 10 kb is a 10-20% perturbation to `rho` on a minority of steps. A
perturbation to the transition is not new information. So the finding is narrow and clean: **the
distances really are different, and feeding those differences into this transition model does not
improve genotypes.**

### What this does to the rest of stage 15

(b) was the distance half. With it out, and with (a) measuring nothing on either contig, the case for
(c) -- folding nested sites into the parent's chain -- is materially weaker than when the phase was
planned, because (c) was to be the thing that made the per-strand distances *reach* the nested sites.

(c) is also the expensive part: `LinkageModel` never reads `Site::ploidy`, and `posteriors()` indexes
`genotype_index(ai,bi)` against buffers sized n(n+1)/2 while a ploidy-1 `Site`'s
`genotype_ln_likelihood` has only `num_alleles` entries -- an out-of-bounds read. So it needs a real
single-copy emission under the two-haplotype latent state, specified and reviewed, and the plan says as
much.

**Recommendation, on the evidence rather than on the plan's momentum:** do not build (c) yet. Two of
phase III's three landed sub-stages measured nothing and one measured worse. That is not an argument
that nesting-in-the-chain is wrong; it is an argument that the *reason* given for it -- carrying
per-strand distances to nested sites -- has been measured and does not pay. If (c) is built it should
be justified by something else, and the honest candidate is the pooling argument: a nested site
currently links only against same-strand nested siblings, never against the top-level sites around it,
and that is a real modelling gap independent of distance. That is a different claim, needs its own
prediction, and should not inherit (b)'s.

## 15': Order from the traversal alignment, distance from the traversals, no reference below the anchor

**Supersedes (b), (c), (d) and the first draft of this stage.** That draft had a units error and is
recorded below as a rejected alternative, because the way it fails is instructive.

**Why.** A *covering reference* is coming: the current reference plus non-overlapping contigs covering
the graph's nested nodes, so snarls on paths the reference never visits become callable. A nested site
then has **no reference position**, so anything that orders or spaces nested sites by
`Entry::position` is undefined rather than approximate.

### Two corrections that shape the design

**Loops: subsequent visits are MASKED, and the unit stays the chain.** A parent traversal may enter the
same nested chain more than once. Each visit sits at a different point along the haplotype and would
have its own distance, which makes the natural unit a (chain, visit) -- but that also makes a chain
produce more than one record, which is stage 17's copy-number question. **Decided: mask visits after
the first.** The unit stays the chain, stage 17 is not pulled forward, and the copy-number
representation stays capped.

This is already the behaviour on both sides, and the point of the decision is to make it deliberate
rather than incidental:

- Ploidy already masks. `child_ploidy` does `if (crossings > 1) { crossings = 1; }`
  (`src/graph_caller.cpp`), so `copies` counts *traversals that cross* rather than crossings summed
  over traversals -- which is the correct ploidy semantics, and worth stating because the two are easy
  to conflate. A chain crossed twice by ONE traversal is one haplotype carrying two copies, not two
  haplotypes carrying one each, and it must not be genotyped at ploidy 2.
- Distance already masks. `traversal_offset_span` returns the first crossing and stops
  (`src/graph_caller.cpp:3637-3645`). Under this decision that is correct rather than a limitation, and
  should be documented and asserted as such instead of left as an implementation detail.

**The masked population is tiny, measured on the whole-genome run:** the multi-crossing warning fires
**0 times on chr20** and **242 times on chrX**. So masking costs essentially nothing on the autosomes
and little on the haploid contigs, which is what makes deferring stage 17 cheap. Replace the
per-occurrence warning with a counter so the size of the deferred question is reported every run rather
than inferred by grepping a log -- it is currently gated on `--progress` and printed once per
occurrence, which is the wrong shape for a population that needs sizing.

**The two chosen traversals are alignable, and the alignment gives the order.** The parent's two
settled traversals are two paths through one snarl, sharing its boundaries and whatever nodes they
have in common. Aligning them yields a merged sequence of chain *visits*: a visit both traversals make
is one column; a visit only one makes slots in between its neighbouring shared anchors. That is a
single order both haplotypes agree on, derived from the graph rather than from a heuristic tie-break,
and with visits after the first masked there is exactly one column per chain. It replaces the "sort by
the longest allele, shorter breaks ties" rule, which stage 14 had already measured as barely load-bearing at the
sibling level (0.606% of adjacent sibling pairs reorder).

### ORDER and DISTANCE are computed separately, and that is the whole fix

The rejected draft made one number do both jobs, and that is precisely what broke it.

**Order: a lexicographic snarl-tree key.** `(top-level anchor, offset within parent,
recursively down the tree)`, compared componentwise. This preserves subtree containment *by
construction* -- no arithmetic claim is being made, so none can be violated. Within one parent,
sibling order comes from the traversal alignment above.

**Distance: pairwise, between consecutive sites only.** For two adjacent sites C1 (under parent P1)
and C2 (under P2), the distance along one haplotype is

    remaining span of P1's settled traversal after C1
      + anchor-frame gap P1 -> P2
      + offset of C2 within P2's settled traversal

Computed for the pair, never as a difference of two absolute coordinates. This cannot reorder anything,
because it is not the sort key.

**Why the rejected draft was wrong, and why it would not have been caught.** It defined
`frame(child) = frame(parent) + offset`, summing a *reference* position with a *haplotype-walk* length.
The offset is bounded by the traversal's length, not by the parent's reference span, so an insertion
makes it arbitrarily larger: P1 at POS 1,000 carrying a 40 kb insertion with a child at offset 30,000
gets frame 31,000, while P2 at POS 5,000 with a child at offset 50 gets 5,050 -- inverting the true
order. Every cross-parent step is mis-scaled by the parent's net indel content. And **the sort hides
it**: the group is sorted by the same key then used for spacing (`src/linkage_model.cpp:2262-2273`), so
gaps come out non-negative by construction and no statistic reveals the inversion. The only symptom is
the F1 number -- which is exactly how 15(b) failed.

### Required fixes, from an adversarial review of the draft against the source

1. **Measure the cross-parent population before building.** Stage 14's 0.606% / 91.05% are
   **sibling-only** -- computed within one parent's `children_of` (`src/graph_caller.cpp:4888-4951`,
   counters at `:3657`) -- while the comparison this stage changes happens across a whole
   `(phase_set, strand)` group, one contig on one haplotype, where cross-parent adjacency is the common
   case. The neutrality prediction currently rests on a statistic about a different population. Measure
   inversion count and gap-ratio agreement over adjacent pairs *within a group*, plus each parent's net
   settled-versus-reference excess.

2. **One `int32`, keyed by the traversal it was measured along -- not two keyed by strand.** A haploid
   nested parent has only `trav_first`, so only slot 0 is fillable, yet its child's strand is
   `parent.nested_strand`, which can be 1: the child would read an unwritten slot, sort to the group
   head, and hand its neighbour a spurious multi-megabase gap. That is the shape of the 448-site
   regression already recorded at `src/linkage_model.cpp:2044-2053`. Tag the frame by its traversal and
   resolve trav->strand at placement exactly as `parent_trav` already does (`:2145-2152`); for a
   `copies == 1` child `parent_trav` *is* the tag, so one value suffices.

3. **`copies == 2` nested sites are not in the nested path at all.** `nested_context.active =
   (copies == 1)` (`src/graph_caller.cpp:5003`), so a two-copy child has `Entry::nested == false`,
   lands in `chainable`, and is sorted and spaced by `entries[].position`. The gate "zero nested sites
   ordered by a reference position" would compile and pass while these sites keep reading one -- and
   under a covering reference they would sort on a meaningless value. Split the overloaded `nested`
   flag (it currently does three jobs: cut the chain, hold out of the diploid run, route to the
   deferred strand pass) with an explicit "descended" bit, then either give these sites a chain or drop
   them from scope explicitly.

4. **The 11,035-child figure does not justify the pair form.** It is het-parent-only, from *called*
   traversals, over all `children_of` at descent (`src/graph_caller.cpp:4925`). The barrier's
   `copies == 2` population is deferred nested chains from *settled* traversals, already counted as
   `carried_on_both` (`src/linkage_model.cpp:2091`). A homozygous parent tests one crossing bit twice
   (`src/graph_caller.cpp:4030-4033`), so its two offsets are identical and the pair form is vacuous
   there. Re-derive the het/hom and descended split of the actual population, or drop the claim.

5. **Three write sites, and a parent-side index the barrier does not have.** `set_parent_trav`'s return
   is discarded and is a silent no-op when no entry exists (`src/graph_caller.cpp:4047`,
   `src/linkage_model.cpp:1513-1522`); `respecify` (`:4158`) and the `record` fallback (`:4185-4188`)
   are the other two. A setter alone drops the frame for gained chains (~2,950, the `!has_entry`
   population), which then enter the layer with a frame of 0. And `traversal_offset_span` needs the
   *parent's* `travs`: generation-1 parents are top-level records in `render_records`
   (`src/graph_caller.cpp:5094`) which the barrier never indexes by `record_key`, so the largest slice
   of nested sites has nothing to walk from. Build a `record_key -> (container, thread, index)` map over
   `pending` + `render_records` at the top of the barrier, and pass the frame through `record()` as well
   as the checked setter.

6. **Return which boundary was entered, and carry a direction bit.** `traversal_offset_span` opens on
   `node == start || node == end` and discards which (`:3637-3645`). Depth 1 survives, but a grandchild
   under a reversed crossing is mirrored -- offsets measured from the wrong end -- and sibling order
   inside that subtree reverses. This is masked today only because v1 descends solely where the
   reference also goes, so every child has a reference path and is flipped onto it; **stage 16 removes
   exactly that gate**, so the construction would be well defined only on the population it is being
   retired from. Orientation *is* recoverable: `PendingRecord::snarl` and `::travs` are co-oriented by
   construction (`:4499`, `:4785`). Also note `*span` is the child's extent along the *parent*, not the
   child's own length -- do not use it as the latter.

7. **Retract "the fallback population is empty by construction."** It is false: three `continue` paths
   in the barrier leave a child with no settled parent traversal. Count them and give them a defined
   behaviour.

8. **Types.** A frame feeding `Site::position` (a `size_t`) through an unsigned gap must not be able to
   arrive negative or unset -- that wraps rather than failing. Keep the frame signed to its own last
   consumer and check at the boundary.

### Gate -- non-regression, and now with the right denominator

Stage 14's numbers cannot support the neutrality prediction (fix 1), so the gate is: the **new**
cross-parent measurement first, then chr20 and chrX ALL F1 and JointIndel not below stage 15(a), with a
null result a PASS. The payoff is the capability, not accuracy: 15(b) already measured the accuracy
upside of a better frame as negative on reference-anchored data. Determinism byte-identical across two
runs and independent of traversal enumeration order. Cost reported per contig, since the traversal walk
moves from the parallel sweep (where stage 14 measured it free) to the serial barrier.

## 16. Children with no reference path

**Goal.** 12,516 chr20 children are skipped for having no reference path, so REF and POS are undefined for them. In the haplotype frame they are orderable and linkable; only *rendering* needs a reference POS.

**Changes.** `src/graph_caller.cpp:4638` increments `g_descent_skipped_no_ref` and `continue`s. Descend instead, genotype, record with `emitted = false`, and render nothing. Stage 14's frame coordinate is defined for them — it is an offset along the parent's traversal, not a reference position — which is why this cannot precede stage 14. Their children inherit a strand from them as from any other nested parent. These sites have never been through `emit_variant`, which indexes `called_traversals[ref_trav_idx]` unchecked; the barrier guards this at `:3945-3958` and the render path needs the same guard, where for these sites the guard is the normal case rather than the exception.

**Gate.** Baseline stage 15. (i) No record is emitted for a no-reference-path site. Do **not** gate "the record count must not change by a single line": after stage 10, `wants_line` is a function of the settled genotype, and adding 12,516 sites to the chains changes the settled genotype of emitted neighbours. Gate instead on a bounded and reported neighbour-flip count, with the flips themselves attributable to a changed settled genotype. (ii) `g_descent_skipped_no_ref` 12,516 → 0, and the counter renamed rather than left describing a population that no longer exists. (iii) Linkage sites rise by ≈12,516 and `bytes()` by the corresponding arena. (iv) **Read I/O at +0%** — these chains are inside the parent's window and their reads are already resident, and the maintainer's constraint is exactly this. Measure fetched-read count, not wall time. (v) CPU within a stated bound for genotyping 12,516 extra chains. (vi) ALL F1 not below stage 15's; descendants of these children that *do* have a reference path gain a strand, count reported and greater than 0. A fall means 12,516 sites' worth of noise on one contig, which is a finding.

**Tests.** Unit: an entry with `emitted = false` and no reference position enters the chain, is phased, and produces no `Change`. The collector already supports unemitted entries, so confirm the test fails by asserting the *count* of such entries after a descent, which is 0 on stage 15. TAP: a chain reachable only from a non-reference allele; assert its children carry a strand. Fails on stage 15, where the chain is never visited.

**Output moves:** no new records, but neighbours may flip. **Reversibility:** one `continue` restored.

## 16 attempted and BLOCKED, with the diagnosis so far

Reverted, not landed. What is established, and what is not.

**Removing the skip works.** Replacing the `continue` at the no-reference-path test does descend into
the population: child calls go 30,416 -> 42,932, exactly +12,516, and the descent-depth histogram grows
at every level. So the chains are reached and genotyping is attempted.

**Nothing comes out.** Linkage sites stay at 219,600, retained chains stay at 30,416, records stay at
115,038, and read I/O is flat (14,231,576 against 14,251,354, -0.14%). So all 12,516 return early
during setup, before they are recorded or staged. That also means the read-I/O gate is *satisfied* --
these chains are inside the parent's resident window, exactly as the plan predicted -- but on a
population that produces nothing, so the number does not yet mean what it will mean.

**Two things were tried on the setup path.** The first was widening
`common_names.empty() && parent_child_trav_sets != nullptr` to admit a symbolically descended child.
That does not fire: `common_names` asks whether the reference path NAME is present in the snarl, while
the skip test asks whether the parent's reference TRAVERSAL crosses the child, and these children
generally have the name. The second was falling back to the parent's interval when
`get_ref_interval` returns -1 for a descended child, which is the right shape -- the
`use_parent_interval` branch already leaves `ref_trav` empty and takes the first traversal as a
pseudo-reference, and the read-likelihood genotyper already guards `ref_trav_idx >= 0` in both places
it uses it. Neither changed the counts.

**Where the diagnosis stopped.** A `no_reference_path` flag was threaded through `NestedContext` and
`PendingRecord` to mark these sites as never-render, and instrumentation shows it is **false at the
child's own entry to `call_snarl_internal`** -- 0 of 12,516 -- even though it is set 12,516 times in the
parent's descent loop, immediately before the call, on the same thread, after `saved` is taken. That
contradiction is unresolved and is the next thing to chase. One of the debug counters was also placed
in the wrong function (its anchor text occurs twice in the file), so part of that instrumentation was
measuring nothing, and it should be redone with the placement checked before anything is concluded.

**Why it is reverted rather than left.** With the skip removed and nothing produced, the change is
12,516 wasted child calls per contig for no output. Reverting also keeps the read-I/O and CPU baselines
honest for whoever resumes.

**What resuming needs.** Establish where the 12,516 return -- instrument every `return false` between
the ref-interval block and the genotyping call, with each counter's placement verified against the
enclosing function. The likely candidates are the traversal finder returning nothing without a
reference to anchor on, and the `assert` on `ref_trav`'s endpoints, which is inside
`if (!use_parent_interval)` and so should be skipped but is worth confirming. Only once a chain is
reachable and genotyped does the rest of stage 16 -- linking it, phasing it, and passing a strand to
its children while emitting nothing -- become testable.

## 17. The copy-number cap
 
**Goal.** A chain one traversal crosses twice — tandem duplication, cycle — is counted as one copy and its second copy is never scored.

**Changes and the representation question that gates them.** `child_ploidy` caps at 1 with a `--progress` warning (`src/graph_caller.cpp:3744-3752`, verified), and ploidy is clamped to {1,2} at `src/linkage_model.cpp:1216` and `:1348`. Lifting the cap collides with four representations no candidate plan named: `Entry` stores exactly two settled alleles (`final_i`, `final_j`); `Site::genotype_ln_likelihood` is defined only for ploidy 1 (allele-indexed) and ploidy 2 (triangular, `n_gt = k*(k+1)/2`); `PhaseCall` has exactly two strands; and `LinkageModel`'s state is an ordered pair of *panel haplotypes*, so two copies on one haplotype has no representation at all — the same objection that rules out three copies. Further, `nested_context.active = (copies == 1)` (`src/graph_caller.cpp:4677`) and `carrying_trav` is set only under `copies == 1` (`:4664-4669`), so at `copies == 2` the child stops being a nested site and loses its strand and parent linkage. And two copies of one snarl print the same ID, so `(contig, POS, ID)` stops being total — stage 1 regresses undetectably — and two entries collide on one `record_key`, which is the sole identity for `has_entry`, `retract`, `set_parent_trav` and the `by_key` map, all first-match lookups.

**So: measure before designing.** Stage 14 reports how many chr20 crossings are currently capped. If it is small — say under 100 sites — the correct outcome is to keep the cap, convert the warning into a counted and reported class, and record the measurement. That is an honest end for this stage and it costs nothing to discover. If it is substantial, the smaller build is linkage-and-mosaic only: score and phase both copies, emit one record, and defer the two-records-at-one-POS representation to the same later change as `--nested-pseudo-ref`. Which of the two gets built is open decision 6.

**Gate.** Conditional on the measurement. If built: capped count → 0; sites rise by exactly the previously-capped count; a copy index added to both the sort key and the `record_key`, with stage 1's determinism gate re-run; total copies above 2 refused and counted rather than clamped; accuracy reported at the affected positions specifically, since a few hundred sites are invisible in a contig F1. Two-sided: a large accuracy move on a rare population is a bug, not a result. Also check identifiability before the accuracy number — both copies of a tandem duplication share sequence, so if their likelihoods are identical by construction, the honest outcome is copy number 2 without phasing the copies apart, said plainly.

**Tests.** Unit: a traversal visiting a child twice reports 2; fails today (reports 1, sets `capped`). Unit: total copies of 3 is refused, not clamped; fails today (silently clamped). `test/t/18_vg_call.t:279` is a cyclic-graph fixture and is affected by both this stage and stage 10.

**Output moves:** yes, if built. **Reversibility:** the cap is one line.

---

# Phase IV — parent conditioning and the depth term

## 18. Measure whether conditioning has anything to remove

**Goal.** Before choosing a factor, find out whether restricting the child's panel to haplotypes carrying the parent's settled traversal removes any rows.

**Changes.** Read-only, over data already resident: per nested site, how many panel haplotypes carry the child chain (non-negative entries in the parent's and child's `hap_arena` rows) against how many of those also take the parent's settled traversal, which is `Entry::final_i`/`final_j` decoded through `trav_arena`. One thing the plan must not pretend is free: there is no `record_key` → entry index, and `has_entry`, `retract` and `set_parent_trav` are linear scans over `entries` under a mutex. Finding each nested site's parent entry needs that map — trivial to add, but it is work.

**Gate.** chr20 byte-identical to stage 17's. Decision rule fixed before the number is seen: if the median panel-row count falls by less than 10% under the restriction, the conditioning is a no-op on this panel and stage 19 is not built — record the measurement and stop. Above 25%, build it. In between, report and let the maintainer decide. Also report the count of nested sites where the restriction leaves fewer than 2 rows, and the distribution of conditioned carrier counts, not just the mean: `freq_prior` is an exponent over multiplicity (default 5.0), and the header records that at 4-haplotype panel sizes the prior is inert because multiplicity barely varies. Conditioning can push a 34-haplotype site into that regime. Run the same measurement on a 4-haplotype chr20 graph in the same pass — it is cheap and it prevents building a rich-panel-only feature without knowing it.

**Tests.** None; this stage produces a number.

## 19. Condition the child's panel on the parent's settled traversal

**Goal.** Make the parent's genotype enter the child's model. Of the three candidate factors: rescoring the child's reads under the parent's sequence is excluded by the reads constraint; pinning the child's latent haplotype to the parent's strand is delivered by stage 15 putting the child in the parent's chain on the parent's strand, so it is not separate work; that leaves restricting the state space, which is also the only candidate that makes `freq_prior` conditional.

**Changes, and the mechanism is not the obvious one.** Do **not** set `haplotype_allele[h] = -1` for non-carriers: `-1` is the escape state (`src/linkage_model.cpp:145-155`, verified — `marginal[ai] * escape`, and `overall * escape * escape` when both are negative), so that makes a non-carrier a free wildcard rather than removing it. It is a *weaker* constraint than the status quo. Use emission zeroing, on the model of the constraint path already in the tree at `src/linkage_model.cpp:648` (`emissions[t][a*m+b] = 0.0`). Then the multiplicity that `freq_prior` exponentiates is computed over the surviving rows, which is the point.

**Gate.** Baseline stage 18. (i) The falsification: two parent genotypes both implying one copy on the same strand but with disjoint carrier sets must give the child *different* posteriors. Today they are bit-identical, so the test fails by construction before and is the cleanest demonstration that step 3 was actually done. **But it is not sufficient alone** — it is satisfied by any conditioning, including one that decodes the wrong parent traversal, since `final_i`/`final_j` are compact indices that must go through `trav_arena`. So add: a unit test pinning the surviving rows to the specific traversal the parent settled on, and a counter for child sites whose surviving rows are inconsistent with the parent's settled traversal, which must read 0. (ii) Nested sites whose conditioned panel is emptied fall back to the unconditioned panel plus the wildcard, **not** to an all-wildcard panel — an all-wildcard panel destroys the frequency prior and presents as a precision gain while wrecking recall, which is the failure mode the `escape` header warns about. If the emptied count exceeds 10% of nested sites, the restriction is too aggressive and the weaker fallback must be measured before the stage is judged. (iii) Report the fraction of posterior mass on the wildcard at nested sites before and after: `escape` (default 1e-2) was tuned against unconditioned panels, and zeroing rows effectively changes it. (iv) Accuracy at **nested positions specifically**, on chr20 and chr6: FP falls, FN does not rise, ALL F1 does not fall; the overall figure is dominated by top-level sites this cannot touch. State a floor from stage 18's measured reduction rather than asserting one now. A fall means conditioning on the parent over-restricts the child, which is a real possible answer for the one stage in this plan with no structural defect behind it — nothing is broken here, the model is merely less informed than the target says. The accuracy gate should be allowed to kill it.

**Tests.** The falsification and the traversal-identity test above. Second: where every panel haplotype takes the parent's settled traversal, the conditioned result is bit-identical to the unconditioned one. Third: a site whose conditioned panel is empty produces the same posterior as the unconditioned site; confirm by removing the fallback and observing the posterior collapse onto the wildcard.

**Output moves:** yes. **Reversibility:** one masking step at one call site.

## 20. The Poisson depth divisor — measurement first

**Goal.** Settle whether the depth term's per-haplotype rate is wrong at nested single-copy sites. Per the corrections above, the arithmetic does not support the "2x" claim as stated, and `depth_weight` defaults to 0.1, so anything done here moves default output.

**Changes, stage 20a: measurement only, no code change.** `DR` is set unconditionally (`src/read_likelihood_caller.cpp:273`) and emitted as a FORMAT field (`:377-381`) whether or not the term is armed, so the baseline is free. Dump the `DR` distribution on chr20, split by site class: top-level diploid, nested ploidy-1, chrY, non-PAR chrX. The prediction from reading the code is that nested ploidy-1 and top-level diploid both centre near 1.0, because `depth_rate = local_read_rate / effective_ploidy` (`src/allele_likelihood.cpp:940-942`) and `expected_reads` sums one term per copy (`:53-61`), so the copy count cancels. If that is what the data shows, **the diagnosis is wrong and no code change is made** — record the measurement and close the item. If nested ploidy-1 sites centre near 0.5, the copy count is being applied twice somewhere the reading missed and 20b proceeds.

**Changes, stage 20b, conditional.** If it proceeds: the divisor becomes the region's ploidy rather than the site's copy number. This is not one argument. `AlleleLikelihoodCalculator::compute` (`src/allele_likelihood.hpp:640-642`) is a pure virtual reached only through the `SnarlCaller::genotype` virtual (`src/read_likelihood_caller.cpp:83`), which has three implementations, plus `src/unittest/allele_likelihood_scoring.cpp`. And `params.depth_ploidy` is assigned nowhere in `src/`, so it cannot serve as the channel without being wired first. `scale_depth_rate` (`src/read_likelihood_caller.cpp:321`, `:335`) exists solely to rescale the rate for the alternate ploidy; under the new convention it is deleted, not re-conventioned — and note that the both-ploidies scoring at `:310-336` itself survives, since the read-likelihood matrix does not depend on ploidy and scoring each chain once at both ploidies while its reads are resident is exact.

**Gate for 20b.** Default flags, not a tier-2 arm — the term is on by default. Median `DR` at nested ploidy-1 sites moves at least halfway to 1.0 from 20a's measured value. chrY byte-identical: there region ploidy equals site ploidy, so nothing may change, and this catches botched plumbing. Called-traversal length distribution at nested haploid records reported either side. chr20 ALL F1 not below stage 19's, with JointIndel and Deletion moving in the predicted direction by more than the 0.0002 noise level; if recall falls, the sign is opposite to the mechanism claimed and the change is reverted, not tuned. Count and report every site where the region ploidy is not knowable at `compute` time and a fallback is taken — a silent fallback to the site ploidy reproduces the bug for exactly the sites that matter. `--depth-quality` (off by default) consumes `DR` too, so A/B both.

**Tests.** Unit in `src/unittest/allele_likelihood.cpp`: at region ploidy 2 and site ploidy 1, `depth_rate == window_rate / 2`; fails today (`window_rate / 1`). Second case at region ploidy 1, site ploidy 1, asserting no change — there to pin chrY, not to catch the bug. Third: `expected_reads` for a one-allele genotype is half that for a two-allele genotype of the same allele at the same rate, pinning that the ploidy appears exactly once.

**Output moves:** 20a no, 20b yes on default output. **Reversibility:** 20a nothing; 20b one signature change plus a deletion, with the deletion as its own commit.

---

# 21. Whole-genome run 2, and the results pages

**Goal.** One genome-scale A/B for phases II–IV against run 1, and the regenerated results pages the eval repo's stale documents need.

**Changes.** No caller code. Regenerate `docs/wgs-results.md`, `docs/wgs-performance.md`, `docs/sv-residual-errors.md` and the whole-genome figures in `docs/nested-calling-design.md` from run artefacts. Also unscheduled by every candidate plan and stale by this point: `planning/vg-call-linkage-hmm.md` (the design document for the reference-frame linkage model stage 15 replaces — the single most stale artefact this work produces), `docs/tier2-depth-term.md` (stage 20), `docs/tier2-phasing.md` and `docs/tier2-parameters.md` (stages 15 and 19 change the state space `escape` and `freq_prior` act on), `docs/sv-fp-anatomy.md` (stage 10 removes ~275 FPs). Rewrite `doc/read-likelihood-genotyping.md`'s linkage sections for the post-stage-15 architecture — there is no separate nested linkage stage to describe any more — including `:221` (Li-Stephens over reference distance) and `:756` (a chain is a maximal run of one ploidy).

**Gate.** Against run 1, at matched thread count, warm cache, ≥3 repeats for RSS, CPU not wall, rates recomputed from summed counts. Autosomal ALL F1 above run 1's; SNV not below; SV not below by more than 0.0020. Per-contig peak RSS not above run 1's by more than stage 7's measured retention plus the 0.7 GB noise floor — stage 7's number is what makes this a gate rather than a hope — and the memory model refitted if its worst residual exceeds 1 GB. Every structural invariant from run 1 still zero, genome wide, with the mosaic invariant in its entries-versus-emitted form. chr20's contribution reproduces stages 19/20's chr20 figures exactly. The citation-resolution script from stage 12 run over the regenerated documents, so a stale line number cannot ship. Every new regression test demonstrated to fail against a named earlier commit, with that commit recorded beside the test; a test whose failing baseline cannot be named does not go in.

**Output moves:** no.

---

# Whole-genome runs: where, and how many

**Two, at stages 6 and 21.** There are exactly two questions worth a genome. Does what is already pushed hold at scale — four commits, chr20-only since `a27149728`, with the linkage site count nearly doubled and chrX/chrY/acrocentric contigs unexercised. And does the finished work hold at scale. Anything between those answers a question chr20, chr6 and chrX answer more cheaply.

Run 1 sits **between** the small fixes and phase II, not after phase II, for a specific reason: phase II is the largest change in the plan, and discovering a scale defect in already-pushed code while a much larger change sits on top of it makes the two indistinguishable. It sits after stages 1–5 rather than before because stages 2, 3, 4 and 5 all move output and a run before them would be re-run.

**A third run is conditional, not optional-by-default.** It is required if run 2 fails on a contig chr20 cannot represent — chrX's mixed ploidy through stages 4 and 15, chr1's scale through stage 8's retention — *and* the fix is confined to a single stage, so the re-run is a re-validation rather than a new question. It is not justified if chr20, chr6 and chrX all agree in direction and magnitude with a mechanism that explains the sign for every stage from 13 on. Phase III deliberately does not get its own genome run: it is the phase most likely to be revised, and chr20 + chr6 + chrX + chrY covers the multi-chain and pure-haploid cases chr20 cannot.

---

# Open decisions

**1. Retention threshold, and in what unit.** Stage 7 measures per-contig retention and projects chr1. The proposal is: retain unconditionally under 7.0 GB projected chr1 peak against the 24 GB budget, retain-and-release between 7.0 and 12 GB, abandon above 12. Retain-and-release is a different stage 8 (a lifetime that ends at generation resolution rather than at write time), and abandoning means widening the ALT list at emission, which permanently loses the recall mirror and the POS/REF/ALT renormalisation for top-level sites. The tradeoff is memory headroom on the largest contig against those two capabilities. The measured anchor is 3.18 kB per nested chain, giving ~0.53 GB on chr20's 165,408 top-level snarls and ~1.8 GB projected to chr1 before the alt-ploidy saving — but the projection depends on chr1's snarl count scaling with its record count, which is unverified. Only the maintainer can set the band, and it must be set before stage 7 runs.

**2. Which factor carries parent conditioning.** The plan builds state-space restriction (stage 19) and argues it is the only candidate that makes `freq_prior` conditional, with pinning delivered free by stage 15 and read-rescoring excluded by the reads constraint. The alternative reading is that pinning is the whole of step 3 and stage 19 should not be built. Stage 18 is designed to be neutral between them — it measures how many panel rows the restriction removes, which is informative either way — so the decision can wait for that number. The tradeoff: state-space restriction is the stronger intervention and the one that can over-restrict; pinning is weaker, safer, and possibly already done.

**3. Does retention follow nesting or the collector?** `nested_calling` can be false while the collector is armed by `linkage_weight > 0.0` alone, and stage 8 as written gates retention on `defer_nested_descent`. In that configuration nothing is retained, the render pass is empty, and `apply_linkage_change` is the only mechanism — so stage 11 cannot delete it. Option A: gate retention on the collector being armed, so decide-then-render covers every linkage-active configuration and stage 11 deletes cleanly, at the cost of retaining on runs that never nest. Option B: keep retention gated on nesting and scope stage 11 to delete only what the deferred path used, leaving two decision paths in the tree. This must be answered before stage 9 is written, not at stage 11.

**4. What QUAL should be once the record is built after the decision.** Today it is a declared function of the pre-linkage genotype and is never patched, which is why 4,490 hom-ref records carry non-zero QUAL. Computing it once forces a choice between the posterior-derived quality and the per-site likelihood ratio, and stage 10 cannot inherit an answer because there is no consistent answer to inherit. Related: today's GQ is capped at GQI on measured grounds (+0.003 AUC, 1–2% fewer surviving false calls). Stage 10 keeps the cap so the change is one thing at a time; whether the cap is still the right operation once the quality is computed rather than patched is a separate measurement.

**5. Should a settled-reference site still emit a line?** Stage 10 stops emitting ~4,490 chr20 records. That is right on the target model's own terms — whether a genotype implies a VCF record is answered by whether the symbolic alleles differ from the reference — but it changes record counts on every contig, and every downstream comparison in the eval repo is calibrated against them. `--genotype-snarls` exists to force lines for sites that do not want them. Is that a sufficient answer, or is a separate flag wanted? Confirm before stage 10 is written.

**6. Stage 17: build it, build the smaller version, or defer?** The three options are: full two-copy representation, which needs a copy index in the sort key and the `record_key` and a decision about two records at one POS; linkage-and-mosaic only, which scores and phases both copies and emits one record, strictly better than today and defers the representation question; or keep the cap and convert the warning to a counted class. Stage 14's capped-crossing count decides whether the question is worth answering at all. Related and separate: if the two copies of a tandem duplication are not distinguishable by the read matrix, calling the locus at copy number 2 without phasing the copies apart is the honest outcome — but that changes what the output means and needs sanctioning before it is built.

**7. The mosaic format change (stage 2) breaks any consumer outside this harness.** Five in-repo consumers are known and updatable. The version bump from 2 to 3 is the migration path, and someone has to own the possibility of an external consumer that reads `.` as an unknown token.

**8. Is stage 15's phase-block figure a hard gate or a reported cost?** Letting nested sites into the runs once took chr20 from 22 blocks to 9,460 and N50 from 248 Mb to 1.08 Mb. If the haplotype-frame design costs a small factor rather than two orders of magnitude, whether that is acceptable against the accuracy it buys is a judgement this plan cannot make.

---

# Deliberately deferred

**Per-generation re-reading of the reads.** Excluded by agreement; the arm that cost +48.8% read I/O does not come back. Every stage preserves scoring each chain once at both ploidies while its reads are resident (`src/read_likelihood_caller.cpp:310-336`), including stage 10, where the ploidy choice becomes a single selection at render time, and stage 16, which is gated at +0% read I/O for this reason. This is also why read-rescoring is ruled out as a conditioning factor in stage 19 rather than evaluated.

**`--nested-pseudo-ref`.** Stage 16 admits the 12,516 no-reference-path children into linkage, phasing and the mosaic and emits nothing for them. Giving them a REF and POS is a representation problem entangled with stage 17's two-copies-at-one-position question; both belong in one later change about representation.

**Widening the ALT list at emission.** The fallback if stage 7's measurement rejects retention. Stage 7 measures its cost so the choice rests on two numbers, but it is not planned as work, and if it becomes necessary phase II is rewritten and the lost capabilities go into the design note rather than being left as an implicit gap.

**Backlog #51, the 2.3x offsetting indel-pair enrichment.** Untouched. Nothing here bears on the mechanism, and folding it in would add an accuracy question to stages that already have one.

**Backlog #43 proper, phase-block fragmentation.** Stage 2 fixes its metric — that is the point of doing stage 2 and stopping. The investigation needs the corrected metric to exist first.

**Backlog #68, whether `--top-down`/`-A` double-descends at scale.** Cheap to check with the descent-depth histogram already in the tree, but on a path this plan does not otherwise touch, and `--top-down` is documented as measuring worse than the default on every axis including recall.

**Per-haplotype mosaic output tracing a path per strand.** Unblocked by the completed work and nearly free after stage 15 — one chain per haplotype per contig *is* a path per strand. It is a new output format with its own consumers and belongs after this work, not inside it.

**Re-fitting `weight`, `freq_prior`, `scale` and `escape`.** Every one is a default with a recorded measurement behind it, and stages 15 and 19 change the state space they act on: a conditioned panel has less multiplicity for a `freq_prior` exponent of 5, and haplotype-frame distances change what `scale`'s 10 kb means. Re-fitting is a legitimate follow-up; a stage that changes both the model and its parameters has no interpretable gate. If stage 15's accuracy gate fails, a `scale` sweep is the first response — a harness parameter sweep, not a code change, and not pre-committed to.

**Total copies above 2.** Stage 17 at most lifts the cap from 1 to 2 and refuses 3, because the model's state is an ordered pair of panel haplotypes. Making a third copy representable is a rewrite of the state space.

**Backlog #66's nine barrier regression tests, as tests of the barrier.** Stage 11 deletes most of the paths those fixes were made in, so tests written against them are work against machinery that is going. They are replaced by four output-level invariants that hold across the deletion and are confirmed to fail against `a27149728`.

---

# Where this plan is guessing

Marked [A] where a stage depends on a claim that is agent-reported or inferred rather than measured, with what to check and when.

**Before stage 3.** [A] That the 1,447/1,116 half-called skew is caused by the depth-2 strand mechanism at all. Two mechanisms are live (see corrections) and which fires depends on entry order in the parallel sweep, so the skew may be a mixture and may be run-dependent. Stage 3a exists to attribute it; if the strandless class dominates, gates (i), (ii) and (iv) all need re-derivation before 3b lands. [A] That ~166 chr20 records carry the wrong haplotype — this follows from the ratio, so it inherits the same uncertainty.

**Before stage 4.** [A] That any records inside a haploid `--ploidy-bed` interior currently carry a diploid GT. The mechanism is identified (`copies >= 1 ? copies : 2`), but the population is unmeasured; the gate has a measure-first abort for this reason. [A] That `crossing_unknown` is large enough to be worth acting on rather than dropping.

**Before stage 5.** [A] That any merged `-L` record currently violates the GT-indexes-max-GL invariant. The index transposition is verified; whether it changes the fold's *answer* on real data is not, because the fold is a max-marginal over collapsed classes.

**Before stage 7.** [A] That chr1's top-level snarl count scales with its record count, which is how chr20's per-snarl figure gets projected to chr1. If snarl density differs, the projection is wrong in an unknown direction. [A] The sizing model itself — a `PendingRecord` shell size and a per-map-node cost are estimates, fine for a 10x decision and not for a 1.5x one. If the routes land within 1.5x, stage 8's measured peak RSS is the answer and stage 7's estimate is not.

**Before stage 9.** [A] That the flatten prefix is small at every site class, which is what bounds the position change to "a handful of near-ties". A site class where the shared prefix is large would move more, and the gate is written to catch that (every difference must sit at a non-zero prefix) rather than to assume it. [A] That no other thread-local read inside `emit_variant` and its callees changes when the emit moves. The audit is specified; a stale read that happens to equal the correct value on chr20 would pass the gate and fail elsewhere.

**Before stage 10.** [A] That the layer's preferences at the 1,465 unrenderable positions are good — the 11% figure is measured on sites where the layer *could* act, which are not the same sites, and 63% of the 1,465 are unjudged and not at random. The FP floor and the judged-subset diagnostic exist because this could fail while the change is still right. [A] That the recall-mirror population is non-trivial at top level; the nested analogue is 511 and the top-level figure has never been measured, so no magnitude is predicted.

**Before stage 11.** [A] That the direct `crossings_of_child` test agrees with the mask everywhere the mask was computable. The per-record comparison in the gate is what checks it; disagreement means `travs` was moved out from under one of them, the class of bug `906812957` fixed.

**Before stage 14.** [A] That the haplotype frame differs from the reference frame enough to matter. This is the whole premise of stage 15's distance half and it has never been measured; the off-ramp is stated first so it can be taken. [A] That computing traversal offsets at descent is cheap enough for the parallel sweep — gated at +3% CPU, unmeasured.

**Before stage 15.** [A] That folding nested sites into the parent's chain does not reproduce the 22 → 9,460 block collapse the code comment records. The mechanism differs (one chain per haplotype rather than a ploidy change mid-chain), but the comment is a measured warning about exactly this move. [A] That `Site::ploidy` becoming load-bearing is a bounded change — the out-of-bounds indexing is verified, the size of the emission work to fix it is not. [A] That the consensus order (children visited, bp length, start node id) is the right *biological* order. It is total and deterministic, which is what the forward pass needs; a wrong order costs transition distances rather than correctness, which is the one comfort, and the plan proceeds on it without evidence.

**Before stage 16.** [A] That these chains' reads are already resident so admitting them adds no read I/O. The gate measures fetched-read count at +0% rather than assuming it.

**Before stage 17.** [A] That chr20 has any capped crossings at all. Stage 14 reports it and the stage may honestly end in "not on this data".

**Before stage 19.** [A] That restricting the panel improves anything. This is the one stage with no structural defect behind it: nothing is broken, the model is merely less informed than the target says. Stage 18's decision rule can kill it before it is built, and stage 19's accuracy gate can kill it after.

**Before stage 20.** [A] That `observed_reads` at a nested single-copy site counts reads from one haplotype rather than from the window's full pile. This is the pivot the whole item turns on, and the code reading says the copy count already cancels between divisor and sum — which means today's `DR` should read ≈1.0 at both site classes and the "2x miscalibration" is unsupported as stated. Stage 20a measures it, and the honest outcome may be that no code change is made.

**Throughout.** [A] The chr20 emitted-record count: 105,251 and 116,965 are both in circulation. Stage 1 reports it and no later gate should quote either figure until then. [A] Every standing baseline used as an absolute constant — 22 phase blocks, N50 248 Mb, 6,716 of 6,716 on one strand, 2,767 mosaic wildcards, 239 panel-unexplained, 511 nested gained — was measured at `4371c9b67`, and stages 3, 4, 10, 15 and 16 all move the populations they count. Each gate that uses one must re-measure it on the immediately preceding stage's output rather than citing `4371c9b67`.

---

# Left on the table

Written at the point the work was wrapped up, so that what is unfinished is unfinished on the record
rather than by omission. Ordered by what a successor would most want to know.

## 1. Stage 16, blocked mid-diagnosis (highest value, and the covering-reference enabler)

Fully described above. In one line: removing the no-reference-path skip does descend into all 12,516
chr20 chains, and every one of them returns early during setup, so nothing is recorded, staged or
emitted and read I/O is flat. The unresolved contradiction is that the `no_reference_path` flag reads
false at the child's own entry to `call_snarl_internal` although it is set immediately before the call
on the same thread. Reverted, with the next instrumentation step named. **This is the blocker for the
covering reference**, since off-reference nested chains cannot be called until it is fixed.

## 2. The reference dependency is 88% removed, not removed

Nested chains are ordered by a snarl-tree tuple and can be spaced along the parent's settled traversal.
What still reads a reference position:

- **Cross-parent steps, 12% of adjacent pairs in a group.** Forming them in the frame needs the
  distance from the earlier parent's END to the later parent's START, and only parent start positions
  are stored. One field on `Entry` (the parent's reference span) closes it.
- **The top-level anchor**, by design -- under a covering reference this becomes the covering contig's
  coordinate by the same mechanism.
- **`copies == 2` nested sites**, which are not in the nested path at all (`nested_context.active =
  (copies == 1)`), so they sit in `chainable` sorted by reference position. The gate "zero nested sites
  ordered by a reference position" would pass while these violate it. The overloaded `nested` flag
  needs splitting first -- it currently does three jobs.
- **Frames not written**: 3,389 chains the settled traversal does not cross and 1,202 with no layer
  entry, of ~30,015. An earlier claim in this document that the fallback population was "empty by
  construction" was wrong and is retracted above.

## 3. A traversal distance is a WORSE predictor than a reference distance, and nobody knows why

The most surprising result here, measured twice by unrelated derivations: spacing nested chain steps
along a traversal instead of along the reference costs ~0.0005 of JointIndel on chr20 (15(b), from
per-site called alleles: −0.00061; 15', from the parent's settled traversal: −0.00052). Same magnitude,
same direction. So it is not the labelling and not the derivation.

**No mechanism has been established for this.** It is the single most interesting open question in the
phase, because the traversal distance is by construction closer to the sequence the sample actually
carries. Candidate explanations, none tested: the reference distance is implicitly a better proxy for
recombination rate than physical distance on the sample's own haplotype; `weight`/`scale` are fitted to
reference distances and a different distance scale needs a refit; or ~7% of steps changing distance is
simply noise that happens to land negative twice.

## 4. `switch_probability` is density-dependent, which is a latent mis-parameterisation

`pow(rho, weight)` raises a per-step PROBABILITY to a power, so it does not compose: the same 1,348 bp
span is **5.83x stickier** divided into 7 steps than taken in one. The model's effective switch rate
therefore tracks site density, not genomic distance, and `weight = 2` was fitted at ~114k chrX sites
and is now applied at 816k.

This is real arithmetic and independently verified, and it is **not** the cause of anything measured
here -- chrX at weight 1 scored 0.93726 against weight 2's 0.93643. The composing form is
`1 - rho = ((1 - rho_min) * exp(-gap/scale))^(1/weight)`, for which N steps multiply to exactly the
one-step answer. Fixing it changes what the number 2 means and needs a refit, which is why it was not
done here.

## 5. The uniform jump

The transition redraws a strand's panel haplotype **uniformly over all m**, so there is no population
structure and no local haplotype frequency in the transition at all (frequency enters only through
`freq_prior` in the emission). For a recombination model that is the crudest available choice. It is
also what keeps the row/column factorisation valid and the step at O(m^2): any non-uniform jump breaks
the collapse to `row[a]`/`col[b]`/`total`. Untested, and noted as a structural observation rather than
a finding.

## 6. Phase IV, not started

Stages 18-20: whether parent conditioning has anything to remove, conditioning the child's panel on the
parent's settled traversal, and the Poisson depth divisor. All three are measure-first stages and all
are independent of stage 16, so they are the cheapest available next work. Nothing here blocks them.

## 7. Stage 17, deliberately deferred

Visits after the first are masked by decision: one copy for ploidy, the first crossing for distance.
Measured cost of deferring: **0 occurrences on chr20, 242 on chrX**. Representing a second copy is the
copy-number question and stays closed.

## 8. Smaller items

- **`transition_apply`'s pair form has no live consumer.** Both strands are passed the same value; its
  only intended consumer is `copies == 2` nested sites, which cannot reach it (item 2).
- **`render_phase_pair` and `PhaseCall::allele_first/allele_second` are written and never read.**
  Removing them needs a rename, because `allele_*` carries the compact pair before that block
  overwrites it with the VCF one.
- **`schedule_wgs.py`'s memory predictions do not model retention** (+50% per contig). Safe only
  because they were already conservative; should be refitted.
- **A one-record discrepancy** between the descended-no-reference counter (12,516) and the frame
  instrumentation's count of the same condition (12,517 measured from called traversals). Both atomic,
  so not a lost increment. Unexplained.

## What this phase got wrong, for calibration

Recorded because the pattern is more useful than any single error. Five claims made confidently here
and then measured false: that the barrier's emit was load-bearing for `respecify`; that the residual
gate items were two defects rather than one; that `unrenderable` was pure dead weight (it was
accidentally protecting haploid sites); that the site-density stickiness explained the chrX regression;
and that the fallback population was empty by construction. Plus three arithmetic errors inside one
instrumentation pass, caught only by reasoning through the null case afterwards -- which is now checked
at runtime.

The thing that worked was elimination against measurement, not hypothesis: the chrX regression was
found by ruling linkage out across three weights and ruling relocation out at 39 of 4,247, which left
only the rendering path.
