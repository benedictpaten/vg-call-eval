# Nested genotyping, as implemented

What `vg call` actually does today, read out of the code rather than out of the design notes: how a
traversal becomes a symbolic allele, how a genotype becomes one or more VCF records, and how the
caller recurses into nested chains and settles their ploidies.

This is a **reference description, not a plan**. Every design note it descends from
-- [nested-calling-design.md](nested-calling-design.md),
[planning/decide-then-render.md](../planning/decide-then-render.md),
[planning/symbolic-diff-decomposition.md](../planning/symbolic-diff-decomposition.md) --
carries stages that were withdrawn, superseded or re-measured, and reading them for current
behaviour means reading history. This page describes the code as it stands.

**Pinned to vg `b136a179a` (branch `read-likelihood-genotyping`, 2026-08-25).** Line numbers are
cited so each claim can be checked, and so the citations rot visibly rather than silently when the
code moves.

**Configuration.** Symbolic collapsing and nested descent are on by default under
`--read-likelihood`, and decline under `--no-nested` and `--no-phased` (a nested site's ploidy comes
from its parent's *phased* genotype, so phasing off means nesting off — `call_main.cpp:1945`). Block
emission is on by default and off under `--no-atomize-blocks`; it refuses to combine with
`-a/--genotype-snarls`, `--legacy`, `--bottom-up` or `--top-down`. The three-phase structure below is
armed wherever the **linkage layer** is, not only where nesting is (`call_main.cpp:2085`), so a
non-nested run resolves generation 0 over an empty pending set and still renders every record from
the settled genotype. One rule rather than two.

---

## Part A — symbolic alleles, and what a genotype gets reported as

### 1. Resolve the site, or decline

`resolve_site` (`symbolic_allele.cpp:38`) requires `into_which_snarl(start)` to return a snarl whose
boundaries match this one, accepting the **forward or the reversed** pairing: `flip_snarl` reverses
any snarl whose reference path runs backwards and the caller then works on that reversed copy, whose
start node is the original end node. Requiring start-to-start turned symbolic collapsing silently off
for 7.4% of chr20's sites.

Unresolvable means projection degenerates to a bare node list with no symbols. Every consumer gates
on `symbolic_site_resolvable` separately rather than assuming, because the failure is silent: with no
symbols, every child reads as "not matched" and a whole subtree gets dropped instead of delegated.
On chr20 this population is 9,279 sites, reported by the run itself.

### 2. Project each traversal

`symbolic_allele` (`symbolic_allele.cpp:95`) walks the visits and replaces each excursion through a
child chain with one `SymbolicStep`. Three things carry the design:

- The symbol is the **chain's** boundary pair (`chain_bounds_of`), not the child snarl's, so two
  traversals entering and leaving by the same boundaries carry the same symbol *however* they cross
  it. That is the entire point.
- Only a genuine child of *this* site may be symbolised. A site that is itself a member of a longer
  chain sees the enclosing chain's bounds; collapsing on those swallows the site's own interior,
  makes every allele compare equal, and erases the variant.
- A chain entered but never left is emitted as a plain node step, never as a symbol that eats the
  tail — the one error mode that would make unrelated alleles compare equal.

Optionally it returns visit ranges partitioning the traversal contiguously, one range per step, which
is what lets a step range be turned back into sequence. A chain symbol's range is half-open at the
exit, because the exit boundary node belongs to whatever step comes next.

### 3. Two consumers, two questions

**"Is this allele the reference *here*?"** — `is_symbolically_reference` (`graph_caller.cpp:2084`),
used in `emit_variant`'s allele loop (`graph_caller.cpp:2710`). A called traversal whose symbolic form
equals the reference's collapses onto **allele 0** instead of becoming a long ALT; its differences
belong to the child chains' own records. This is the SNV-recall mechanism: 55,222 of 142,707 autosomal
SNV false negatives were sitting inside a large allele vg itself had emitted.

**"Where exactly does it differ?"** — `symbolic_diff` (`symbolic_allele.cpp:215`), edit distance with
**substitution at cost 1**. That is a disambiguation, not a different metric: under insert/delete-only,
`[a,b]` against `[b,b]` has two minimal alignments of equal cost and nothing chooses between them.
Remaining ties break deterministically (diagonal, then deletion, then insertion), because this
function decides how many records a snarl emits, and an unstable tie-break would make the output
depend on nothing the caller controls. O(|ref| × |alt|); a pair too large degrades to one block
spanning both alleles, which is the pre-existing whole-allele behaviour, and sets `out_degraded` so
the population is counted rather than assumed empty.

### 4. Build the site record

`emit_variant` (`graph_caller.cpp:2643`) iterates the genotype it was handed:
`STAR_ALLELE_MARKER` → `*`, `MISSING_ALLELE_MARKER` → `.`, symbolically-reference → allele 0,
otherwise dedup by allele string. It keeps `trav_to_allele`, because GT, GL and the linkage layer are
written in VCF allele numbering and only this function knows the mapping from traversals to it.

### 5. Decompose into difference blocks

`emit_block_records` (`graph_caller.cpp:2216`), placed **after** the site record is finished on
purpose: every field a block does not redefine is inherited from a record that has already been
through `update_vcf_info`, flattening and merging.

- Project the reference and each called haplotype; align; then **cluster all haplotypes' blocks by
  overlap of reference step range, with touching counting as overlap**. That merges a deletion on one
  haplotype abutting an insertion on the other — they cannot be separate records, because the two
  alleles would have to disagree about the same reference span.
- Per cluster: cut the reference span; cut each haplotype's aligned span via `alt_before_ref`; borrow
  one anchor base to the left if any allele would be empty (one base, matching what
  `flatten_common_allele_ends` leaves, so the two agree on how an indel is spelled); dedup by string;
  and carry the site's phase across. **Three GT shapes carry a phase set and two of them are
  haploid**: a diploid pair `a|b` transfers its orientation by slot (mapped through
  `trav_to_allele`, since the site's GT is in site-allele space and the block's slots are in
  genotyper order); a nested chain called at ploidy 1 transfers its **strand**, `a|.` or `.|a`,
  which is the only place the VCF can say which haplotype the allele sits on; and a genuinely
  haploid locus keeps `PS` as a block label with no orientation to inherit. Only a slash-separated
  GT is unphased, and only there is `PS` dropped. Testing for the diploid pair first is what made
  every haploid block come out as a bare `a` with `PS` erased — 3,452 lines genome-wide, from 1,517
  snarls, losing the strand and the phase set the unsplit record would have carried.
- `AD` and `GL` are looked up through `site_of_block`, so **every block of a snarl reports the same
  evidence**. That is honest only about arity; `INFO/SB` marks the replicated set so a consumer can
  avoid double-counting it.
- **Every refusal returns −1, meaning "the site record stands."** A case this does not understand
  degrades to the behaviour it was going to have anyway, not to a wrong record.

---

## Part B — recursing on nested chains

Three phases. **Reads are touched in the first only.** The predecessor design descended once per
generation and paid +48.8% reads fetched (903M against 607M genome-wide) for it; the observation that
removes that cost is that the expensive object — the per-read per-allele likelihood matrix — does not
depend on ploidy, and ploidy is the only thing an ancestor's genotype determines about a chain.

### Phase 1 — collect, descending inline while the reads are resident

`call_snarl_internal` (`graph_caller.cpp:5217`); the descent loop at `graph_caller.cpp:5853`. For each
non-trivial child of the current snarl:

1. **Reference gate.** The reference traversal must cross the child, else REF and POS for its record
   are undefined. (`--nested-pseudo-ref` is where that would be handled.)
2. **Exactly-once gate** — `chain_reported_inline` (`graph_caller.cpp:1989`). Under block emission, if
   every called haplotype's crossing of the chain falls inside a difference block, that block's ALT
   already spells the route through it, and the chain's own record would report the same variation
   twice. The subtlety that matters: test each haplotype's **own** projection, not the reference's
   chain step. A haplotype that *deletes* the chain puts the reference step inside a block while
   crossing zero times; testing the reference step fired 60× too often and deleted 399 records
   linkage was still entitled to move.
3. **Ploidy from the parent, not the contig** — `child_ploidy` (`graph_caller.cpp:4517`) counts how
   many *called* parent traversals cross the child, requiring start then end **in order**
   (`crossings_of_child`, `graph_caller.cpp:4392`; testing the boundaries independently counts a
   traversal that touches both on unrelated excursions), capped at the parent's ploidy. A chain
   crossed by one parent allele and deleted by the other is haploid *there*, whatever the contig
   says. A single allele crossing a chain twice is capped rather than modelled, and said so in the log.
4. **Zero copies is not a skip.** The chain is descended with `retain_only`: genotyped now, emitted
   never. The parent is not settled, and linkage may still move it onto an allele that reaches this
   chain (296 such on chr20). Going back to the reads at the barrier to find out is precisely the
   cost this design exists to remove.
5. **Genotype at both ploidies.** `set_want_alt_ploidy(true)` keeps the *other* ploidy's entire answer
   in `alt_ploidy_info`, so re-ploidying at the barrier is a **revision** rather than a re-call.
6. **Stage a `PendingRecord` and write nothing.** The `travs` move happens only after the descent loop
   *and* the `-A` recursion have finished reading them (`graph_caller.cpp:5998`) — one block earlier
   broke four `-A`/`--top-down` tests; at emit time it cost 12,302 chr20 records.

Context travels in a saved-and-restored thread-local `NestedContext`, which is safe because descent is
synchronous on the calling thread. `parent_crossing` is a mask over the parent's own **candidate
traversals**, not its VCF alleles — a nested parent's allele numbering does not exist yet, because its
record has not been built — with `crossing_known` beside it so an uncomputable mask (more than 64
alleles) is never read as "nothing crosses".

Retention is measured, not estimated (`graph_caller.cpp:5113`): the run walks the objects and reports
MB, traversal visits and genotype likelihoods. Peak RSS cannot resolve a delta this size — six runs of
one binary on chr20 spread 3.39 to 4.42 GB, wider than the retention itself.

### Phase 2 — settle at a barrier, one generation at a time

`run_deferred_descent` (`graph_caller.cpp:4733`). No reads. Per-thread queues are merged; `children_of`
indexes by parent key; `record_by_key` indexes **both** containers, because generation-1 parents are
top-level records living in `render_records`, which is the largest slice of nested sites by far.

For `gen = 0..max`: resolve that generation, then for each pending chain at `gen+1` whose parent has
just settled —

- Recompute copies from the settled pair, tested in **traversal space** against `parent_crossing`.
  Testing the compact allele index instead retracted 3,615 chains against a true 190 on chr20; the two
  agree only when every allele at the parent is panel-carried.
- `parent_trav` and the copy count come from the same two booleans, so the count and the identity
  cannot disagree (−1 neither settled traversal carries it, −2 both do).
- **copies == 0** → `drop_subtree`, iteratively over an explicit stack. The settled parent does not
  carry the chain, so the sample has no copy of it, and nothing inside a sequence the sample lacks
  exists either.
- **copies == ploidy, and already in the layer** → nothing to revise. The `has_entry` half is not
  redundant: a chain no called parent allele reached is staged but deliberately *not* recorded
  (2,713 of them on chr20), so for those a matching ploidy proves nothing.
- **otherwise** → swap in `alt_ploidy_info`/`alt_ploidy_best`, `respecify` the layer entry (or
  `record` it, for the retained population that was never filed), then recompute the children's
  crossing masks against the parent's now-current traversals.
- `set_frame` measures where the chain sits along each settled parent traversal, indexed by traversal
  **order** rather than by strand: a haploid parent has `trav_first == trav_second`, so both slots are
  written and neither can be read unset.

The generation bound is re-read after each pass, because gaining a chain can create a layer entry at a
generation the collector has never held, and a snapshotted bound would leave it — and everything below
it — outside every resolve pass: emitted but never settled, never phased, absent from the mosaic.

### Phase 3 — render

`render_retained_records` (`graph_caller.cpp:4681`). Phases are built first, since every generation has
settled by now. Then, parallel over queues, each record asks the collector for `settled_traversals` and
hands **that** pair to `emit_variant`.

This is the load-bearing property of the whole redesign. The ALT list, the symbolic-reference test that
decides whether a line exists at all, QUAL, and the arity of `AD`/`GL`/`GQI` are *all* built by
iterating the genotype passed in — so handing in the settled genotype makes every one of them agree
with the call by construction, instead of being patched towards it afterwards. There is no patch pass,
and the two things a patch structurally cannot do stop being possible rather than being flagged:

| what a patch could not do | was | is |
|---|---|---|
| add an ALT the line lacks | 496 `unrenderable` events | impossible by construction |
| withdraw a line that settled on the reference | 1,383 records left at GT 0/0 | impossible by construction |
| express a ploidy the parent contradicts | three `nested_*` FILTERs | deleted |

Thread-locals, not records, are the hazard here: `emit_variant` reads `nested_context` and
`current_generation`, which in a batch pass hold whatever the last snarl this thread genotyped left
behind. A stale `nested_context.active` would file a top-level site into the nested strand population
silently, so the context is reset per record rather than trusted.

---

## What the structure buys

- **Coherence flags: 7,458 → 0** genome-wide, and by construction rather than by filtering.
- **Reads fetched** back to the single-sweep figure, against +48.8% for the descend-per-generation arm.
- **One place a line is written**, and one genotype it is written from, for top-level and nested
  records alike.
- Accuracy is a wash against the arm this replaced — precision up, recall down, in every class — so
  the change stands on the guarantee and the I/O, not on F1.

## Where it is still loose

Two known asymmetries, neither of them accidental:

1. **`AD` and `GL` are replicated across the blocks of a snarl** rather than being apportioned per
   record. `INFO/SB` makes the replication recoverable, which is a disclosure, not a fix. This is the
   one field set that is knowingly not per-record.
2. **`crossing_unknown` chains are left exactly as the sweep left them.** A parent with more than 64
   candidate traversals, or one that emitted nothing on that invocation, yields a mask the barrier
   refuses to interpret — correctly, since reading an unknown mask as "no allele crosses" silently
   exempted these chains from revision. They are counted and reported, not resolved.

Beyond those, `--nested-pseudo-ref` (a record for a chain the reference does not cross) remains
unbuilt: on chr20, 4,428 ALTs carry a chain the reference does not cross, covering 57.2% of the bases
involved — variation with no record of its own, reachable only inside the allele that contains it.
