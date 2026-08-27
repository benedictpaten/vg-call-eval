# Ordering nested chains by the traversal alignment, not by the reference

Nested chains in a diploid snarl are currently ordered, and spaced, by reference position. This is
the plan to order them by the alignment of the parent's two settled traversals instead, and to space
them per strand along those traversals.

## The motivation, measured

The reference cannot order a chain it does not visit. On chr20:

```
9,022 chains a called allele reaches but the reference does not
by parent:  1=474  2=148  3=56  4=47  5=30  6=26  7+=171
478 parents carry two or more, which the reference cannot order
```

952 parents carry at least one such chain and **478 carry two or more** -- for those, reference
position supplies no relative order at all, not merely a worse one. 171 parents carry seven or more.
Block emission independently reports 4,543 ALTs carrying a chain the reference does not cross,
covering **57.6% of the bases** in those ALTs: variation reachable only inside the allele that
contains it.

**These chains are not descended into today.** The descent loop skips 12,486 children for having no
reference path, because REF and POS for their records would be undefined. So the ordering problem
does not arise yet -- it arises exactly when `--nested-pseudo-ref` is built, and those 478 parents
are what it will immediately need an order for.

That makes this a **prerequisite for pseudo-ref calling** rather than an independent improvement,
which is a firmer motivation than accuracy and should be how the work is judged.

## What the alignment gives that the reference cannot

Two chains the reference misses are not unorderable -- they are unorderable *by the reference*.
Aligning the parent's two settled traversals orders them whenever both haplotypes cross both, which
is the common case. So the alignment does not merely replace reference ordering; it covers a
population reference ordering cannot reach.

## Do not average the two distances

`transition_apply(in, m, rho_a, rho_b, out)` already takes a switch probability **per strand**, and
says why:

> The two haplotypes of a diploid sample recombine independently, so the distance each has travelled
> since the previous site is its own -- which is what the haplotype-frame work needs and what a
> single scalar cannot express.

The call site passes `rho, rho` with "the per-strand distances that will make them differ come with
the haplotype frame." Averaging the two traversal distances would collapse them back into the one
scalar the function was written to escape. Pass both.

## Two prior measurements that constrain the test

**Ordering is predicted inert on today's data.** The frame instrumentation reports **0 of 5,540
same-parent adjacent pairs reorder** under the traversal key on chr20, and the shared-chain
out-of-order population is 25 occurrences on 17 sites (0.023%). Where the reference *does* visit
both chains, it already agrees with the traversals.

**Traversal spacing has been tried once and lost.** Spacing steps along a traversal instead of the
reference cost ~0.0005 of JointIndel on chr20 (0.92390 -> 0.92338), reproduced by two independent
derivations. That was a *single* traversal-derived gap -- which is what averaging would reproduce.
The untested variant is per-strand.

## Plan: three stages, each with its own gate

**Stage A -- per-strand plumbing, no behaviour change.**
`Site::gap_to_previous` becomes two values; `site_gap` returns a pair; `window_posteriors` and
`window_phasing` compute `rho_a` and `rho_b` from them and pass both to `transition_apply`. Where
only one gap is known, use it for both, which is today's behaviour.
*Gate: byte-identical on chr20.* Both gaps are equal until stage B, so anything else is a bug.

**Stage B -- fill the two gaps from the frames.**
`Entry` already carries `frame_offset[2]`, `frame_end[2]` and `frame_total[2]` per settled traversal
slot, written by `set_frame` at the barrier for 26,241 chains on chr20. For adjacent siblings under
one parent, `gap[slot] = frame_offset[slot]_next - frame_offset[slot]_prev`.
*Gate: accuracy* -- aardvark and truvari on chr20 and chr6, against a same-binary baseline via a
runtime switch. This is where the -0.0005 lived, and where the per-strand form gets its first test.

**Stage C -- alignment-based ordering.**
Project the parent's two settled traversals with `symbolic_allele`, align them with `symbolic_diff`
-- the machinery `tally_haplotype_diff` already runs -- and match each child chain to its symbol by
boundary nodes (`Entry` carries `start_node`/`end_node`; `chain_bounds_of` gives a symbol's bounds).
Order children by position in the alignment.
*Gate: count reorderings against reference order.* Predicted zero on today's data, in which case the
change is byte-identical and its value is capability rather than accuracy.

Staged this way, B is the only step with a known risk and A and C should both be inert. Running them
together would leave an accuracy result uninterpretable, which is the trap the pin experiment fell
into by bundling a model change with a speedup.

## Details to settle before coding

- **The 17 out-of-order sites.** Where the alignment cannot match a chain the two haplotypes share,
  fall back to reference order for that group and count it. Do not guess an order.
- **Per-strand gaps are well defined in the groups.** The per-parent groups contain only children
  whose ploidy equals the parent's, so both parent traversals cross every child and both strands
  have a distance. Ploidy-1 children go through the per-strand haploid chains, where there is one
  strand and the existing frame fallback already applies.
- **Chains with no reference position at all** have no `Entry::position` to sort by, so stage C is
  what makes them orderable. Until pseudo-ref calling exists there are none in any chain, so stage C
  cannot be tested on them -- only on the population where both orders exist and agree.

## What this does not buy

No accuracy today. Every chain currently in a linkage chain has a reference position, and reference
distance has twice measured as the better predictor. The payoff is the 478 parents whose children the
reference cannot order, and those are unreachable until `--nested-pseudo-ref` descends into them.

---

## RESULT: built in seven commits; ordering is free, distances are free except on SV

vg `245058e97..c711e97fa` on `read-likelihood-genotyping`. Four of the seven are defect fixes, three
of those found by an adversarial review panel rather than by the author.

| | commit | gate | outcome |
|---|---|---|---|
| A | `245058e97` | byte-identity | per-strand `rho_a`/`rho_b` plumbing — **passed** |
| B | `88a47ab03` | accuracy | frame spacing behind `VG_LINKAGE_FRAME_GAPS`, **off by default** |
| C | `af95aba9a` | reorder count | **0 of 18,097**, byte-identical |
| B′ | `a6fb73e45` | accuracy | start-to-start convention — supersedes B's numbers |
| D | `d51972cd4` | byte-identity | rank carried into `record()` — 897 chains recovered |
| E | `1832d794d` | byte-identity | comparator was undefined behaviour — 204 cycles removed |
| F | `c711e97fa` | byte-identity | frame replay (1,255 → 0 lost) + `site_gap` unknown distance |

Final gates: chr20 **and** chr6 byte-identical, VCF and mosaic; `test/t/18_vg_call.t` **315/315**;
652 `linkage_model` and 86 `symbolic_allele`/`symbolic_diff` assertions.

### The answer the plan was written to get

**Ordering costs nothing.** 0 of 18,097 same-parent adjacent pairs reorder under the alignment key,
on chr20 and again on chr6, so the output is byte-identical. Where both orders exist they agree
exactly — which the plan predicted from the earlier 0-of-5,540 instrumentation and which now holds
with the alignment key rather than the frame key.

**Distances cost nothing on small variants and something on SV.**

|  | chr20 m1 | chr20 m2 | chr6 m1 | chr6 m2 |
|---|---|---|---|---|
| ALL / SNV / JointIndel / Ins / Del | **0** | **0** | **0** | **0** |
| SV F1 | −3.1e-3 | **−2.2e-3** | −5.5e-4 | −5.5e-4 |

Every small-variant class moves by exactly zero, to five decimals, across 1,586–2,318 differing VCF
lines. SV loses on both contigs, as +3 to +4 false positives.

**And the reason the effect is small at all is measured, not argued: 90.2% of the traversal distances
on chr20, and 88.0% on chr6, are EXACTLY the reference difference.** Two adjacent siblings share a
boundary node, so wherever the parent's settled traversal follows the reference through them the two
measures coincide to the base. Only the residual tenth can move anything, and what it moves is SVs.

### Two measurement errors this plan's own formula would have avoided

The plan says `gap[slot] = frame_offset_next - frame_offset_prev` — start to start. `Entry`'s doc
comment says `offset(next) - frame_end(prev)` — end to start. The implementation followed the comment,
and it was wrong twice:

- **It changed the convention as well as the metric.** `site_gap`'s reference fallback differences two
  `position` values, and a snarl's position is where it begins. So the arm tested two things at once,
  and the first reported result — small variants losing 5e-6 of ALL and 2e-5 of JointIndel — was the
  convention, not the metric.
- **It made most steps unmeasurable.** A span is inclusive of its closing boundary node and adjacent
  siblings share that node, so end-to-start is negative for every adjacent pair: 10,103 of chr20's
  18,235 same-parent steps read as "no frame" and were exempted.

That second error also hid the per-strand population. End-to-start found **85** steps where the two
strands' distances differ; start-to-start finds **1,650** on chr20 and **1,899** on chr6 — nineteen
times more, because the refused adjacent-sibling pairs are exactly the ones where the two haplotypes'
spans differ. With them included, **the per-strand form beats the single scalar** (chr20 −2.2e-3
against −3.1e-3, one true positive recovered), which is the first evidence the split earns anything.
The earlier claim that the per-strand form had almost nothing to act on was an artefact.

A separate direction bug — frame offsets run in traversal order while sites are sorted in reference
order — is real but inert: 92 steps on chr20, 130 on chr6, byte-identical on both. Kept because it
looked like the explanation for the negative steps and was not.

### Three defects the review found, and one it correctly rejected

- **`set_align_rank` and `set_frame` wrote to entries that did not exist yet.** The barrier computes
  these facts ~200 lines before a revise block that `record()`s a chain the sweep never filed, and
  both setters are fail-quiet. 897 chains lost their rank; 1,255 slot-writes lost their frame, and
  were *misattributed* to "the settled traversal does not cross". The rank now travels as an argument
  to `record` (as `parent_trav` already did); the frames are staged and replayed.
- **Both sibling comparators were undefined behaviour.** "Compare on this key only if both operands
  have it, else fall through" is intransitive: 204 cycles over small triples with the alignment rank,
  and **81 in the frame-and-position form that predates it**, so this was latent before this work. A
  sentinel is transitive but displaces the unplaced chain and moved 241 chr20 lines; the key is now
  chosen once per group, which is transitive, moves nobody, and degrades to the previous behaviour.
  Only 10 groups on chr20 and 17 on chr6 fall back entire.
- **`site_gap` called 1 bp "no distance".** At the shipped parameters 1 bp gives rho = 1e-12 and
  stay = 1.0 — two sites about which nothing is known asserted to be *perfectly linked*. `SIZE_MAX`
  gives rho = 1.0 and a uniform transition. Guarded on a new `Site::unpositioned` rather than on
  `position == 0`, since zero is a legitimate coordinate.
- **Rejected:** that the affected chains were the reference-invisible population. They are not — those
  never reach the layer. They are chains no *called* parent allele reached during the sweep.

### What remains, for the pseudo-reference work

**Reference-invisible chains are not genotyped at all.** The descent gate `continue`s before
`call_snarl_internal`, so nothing is computed for chr20's 12,486 of them. A blast-radius map
established that removing the gate alone is a no-op — `call_snarl_internal` returns false at the
`common_names.empty()` check — and that relaxing that check *without* also relaxing
`use_parent_interval` hard-aborts: `get_ref_interval` is handed the parent's path, both step maps come
back empty, and `assert(start_steps.size() > 0 && end_steps.size() > 0)` is live, with no `-DNDEBUG`
anywhere in the Makefile. Both must be relaxed together.

The load-bearing fact for that work: **the genotyper needs no reference.** `ref_path_name` and
`ref_range` are parameters `read_likelihood_caller.cpp` never reads again — reads are reached through
merged node-ID ranges. A no-reference site is fully genotypable.

Two gaps must be closed first, and the larger is not the one this plan named:

- **Intra-chain snarl order: 3,650 pairs on chr20.** `align_rank` identifies a *chain* —
  `chain_bounds_of` returns the enclosing chain's boundary pair — while the entries are *snarls*, so
  two siblings in one chain share a rank and the alignment orders neither. The frame offset is
  per-snarl and can close this without the reference, giving a two-level key: the chain by the
  alignment, the snarl by its offset along the settled traversal.
- **The 154 ambiguous chains** (shared but unmatched, or crossed twice). With the fallback removed,
  the all-or-nothing rule sends any group containing one of them back to reference order entire —
  which is what removing the reference forbids.

And the position must be **supplied as absent, never borrowed from the parent**:
`have_reference = (prev.position > 0 && e.position > 0)` gates the frame-derived gap, so a borrowed
non-zero position makes that test true, the frame is never consulted, and the model differences two
identical values.
