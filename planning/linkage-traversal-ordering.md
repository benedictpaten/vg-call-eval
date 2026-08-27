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
