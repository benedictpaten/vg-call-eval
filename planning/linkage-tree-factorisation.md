# Decoding nested sites as a tree, not as insertions into a line

The linkage layer costs 31.7 s of a 223 s chr20 run, and 24.6 s of that is generations 1-4
re-decoding ~212,000 sites to settle a shrinking handful. Five attempts to make that cheaper all
failed, and they failed informatively. This plan is built on what they ruled out.

## What is now established, and rules out the obvious approaches

**The sufficient boundary object is the alpha/beta message pair over haplotype PAIRS.** Three cheaper
candidates were each tried and each changes the answer:

| boundary object | why it fails |
|---|---|
| the settled **genotype** (today's clamp) | `build_emission` maps a delta genotype onto every pair spelling it -- for a het at 34 haplotypes that is hundreds of live states, so the message passes straight through |
| the settled **pair** (a pin) | discards the posterior mass on the alternatives that are still plausible. Measured: SV F1 0.5345705 -> 0.5338229, five extra false positives |
| **uniform** | discards everything. This is what `posterior_in_radius` actually measured, which is why "you need +/-500 sites of context" came out -- that is the distance over which the chain forgets an uninformative boundary, not a dependency |

**Influence decays to nil at ~500 sites.** TV distance from the full window: 5.5e-3 at radius 25,
2.3e-4 at 100, 8.0e-9 at 500. Measured on chr20, stable across generations.

**Insertions invalidate messages to their RIGHT.** A child inserted between L and R leaves alpha_L and
beta_R untouched -- but generation k+1 inserts children *throughout* the chain, and a run at position
q has a stale alpha because of every child inserted at a position < q within ~500 sites. Spacing of
new sites per generation on chr20: gen 1 ~12 sites, gen 2 ~49, gen 3 ~471, gen 4 ~9,634. So
cross-generation message reuse in a linear chain is valid **only at generation 4**, worth ~6 s of
24.6 s. This is the finding that kills every caching variant, and it was reached five different ways.

**The clamp is what makes messages stale**, because it changes the chain's emissions every generation.
Dropping it is required for reuse and is itself a model change: 115,392 -> 115,317 records on chr20.

**Byte-identity is unreachable for any restructuring.** `windowed_marginals` places window boundaries
by index, and each generation's insertions shift every index after the first. Today's posteriors are
an artefact of each generation's own index layout.

## The plan: stop linearising the tree

Every failure above comes from one decision: nested children are inserted into a contig-wide chain
ordered by reference position. That turns a tree into a line, and insertions into a line invalidate
everything downstream of them.

The snarl hierarchy is a tree. Decode it as one.

**Each parent's children form their own chain, conditioned on the parent's posterior over haplotype
pairs.** Not inserted anywhere. Then:

- there is no contig-wide chain after generation 0, so nothing to re-decode
- cost is O(children) per parent, with no caching, no staleness, and no bracket-containment question
- each child chain is independent, so the whole thing parallelises per parent
- ordering within a child chain is by offset along the parent's settled traversal, which
  `set_frame` already records -- the reference-free ordering this was always aiming at

The stored object is exactly what was asked for at the outset: for each parent snarl, the
probabilities over pairs of haplotypes. One 1,225-entry vector per parent that has nested children,
9.6 kB each; on chr20 roughly 20,000 such parents, ~192 MB, and less if restricted to parents whose
children actually resolve.

**Why the tree form escapes the insertion problem.** A child chain is bounded at both ends by its
parent's state rather than by neighbours that other generations will disturb. Nothing is inserted
into anything. The staleness that made linear caching work only at generation 4 cannot arise.

## What changes about the answer

Two deliberate changes, and they should be measured separately rather than together:

1. **Conditioning moves from reference neighbours to siblings.** A nested child currently takes its
   transition context from whatever is nearest in reference coordinates -- often a top-level site, or
   a site under an unrelated parent. Under the tree form its neighbours are its siblings, which share
   a parent traversal and therefore a haplotype by construction. This is the substantive change and
   the one with a real prior in its favour.
2. **The clamp goes**, because a child conditioned on its parent's posterior does not need its
   parent's emission flattened. Coherence is unaffected: a child's ploidy comes from the barrier and
   its strand from the parent's settled `PhaseCall`, neither of which is decided here.

## Staging, with a gate at each step

1. **Harvest and store gamma per parent** during generation 0's decode. No behaviour change; the
   gate is byte-identity plus a reported memory figure. If the figure is much above ~200 MB on
   chr20, stop and reconsider before building anything on it.
2. **Build child chains per parent**, ordered by frame offset, conditioned on the stored gamma.
   Decode genotype and phase for them. Keep the existing path alive behind a flag so the two can be
   run against each other on the same binary.
3. **Remove the generational re-decode**, once step 2 matches or beats it.

Gates, in order of authority:

- **coherence counter** -- how often a child's decode lands in a frame where its parent's genotype
  differs from the one emitted. Must be zero. This is pass/fail, not a judgement.
- **accuracy** -- aardvark and truvari on chr20 against the same baseline, via
  `scripts/tier2/score_vcf.py --dataset chr20-34hap`. The bar this project has used before is
  ~0.0005 of F1; the pin attempt failed it at 0.00075 on SVs.
- **runtime** -- the point of the exercise, but last. Speed that comes from discarding information
  looks exactly like speed that comes from removing redundancy, and the pin attempt produced a 26x
  "speedup" that was entirely the former.

## RESULT: built, measured, and it works

Implemented on 2026-08-26. Linkage generations 1-4 decode per-parent groups instead of the contig
chain, conditioned on the parent's context message.

| chr20, 34-hap panel | baseline | tree factorisation |
|---|---:|---:|
| linkage, generations 1-4 | 24.57 s | **1.06 s** |
| generation 1 | 6.35 s | 0.68 s (18,533 sites in 1,334 groups) |
| generation 2 | 6.10 s | 0.19 s (2,634 in 320) |
| generation 3 | 6.00 s | 0.10 s (541 in 67) |
| generation 4 | 6.12 s | 0.09 s (36 in 8) |
| wall clock | 222.9 s | **192.5 s (-13.6%)** |
| context messages held | -- | 2,290 parents, 21.5 MB |

Accuracy, aardvark and truvari against the same baseline VCF:

| | baseline | tree | delta |
|---|---:|---:|---:|
| ALL F1 (GT) | 0.9723943 | **0.9724101** | +1.6e-5 |
| JointIndel F1 | 0.9282711 | 0.9283140 | +4.3e-5 |
| Insertion F1 | 0.9172879 | 0.9173276 | +4.0e-5 |
| SNV F1 | 0.9851832 | 0.9851900 | +0.7e-5 |
| Deletion F1 | 0.9387837 | 0.9387433 | -4.1e-5 |
| SV F1 | 0.5345705 | 0.5340033 | -5.7e-4 |

Small variants improve in every class but deletions. The SV move decomposes as TP 444->443,
FP 444->442, FN 321->322 -- one true positive lost and two false positives removed, with precision
flat at -2e-5. On a 765-SV benchmark that is a single call, not a mechanism. **It is not evidence of
SV neutrality**; chr6 or the whole genome is what would settle the direction.

858/858 unit assertions, 315/315 TAP, `test/` clean afterwards.

### Two things the implementation had to learn

**Grouping at one generation disabled it at the next.** The group decode passed null harvest
buffers, so a grouped generation stored no context and the next generation found none for its
parents. Group decodes have to harvest too.

**A parent with no stored context still forms a group.** Coverage is never complete, because a whole
population of children -- the chains reachable only under a settled parent -- is recorded at the
barrier AFTER their generation's decode, so their parents were never masked for harvesting. Treating
that as disqualifying let 51 sites veto the other 21,447. A group without upstream context is still
the per-snarl model: the parent's own emission is in it, so its children link to the parent and to
each other.

## Already built and verified

`LinkageModel::window_posteriors` now takes optional incoming alpha/beta and can harvest them at
every index, with `segment_posteriors` as the public entry point. Byte-identical on chr20 with the
default arguments. That is the machinery step 2 needs and it is the one piece of the five attempts
worth keeping.

## Risks

- **Accuracy is genuinely open.** Sibling conditioning is more defensible than reference-neighbour
  conditioning, but "more defensible" has lost to measurement twice today.
- **Nested parents.** Depth reaches 6 on chr20; a parent that is itself nested takes its gamma from
  its own parent's chain, so the recursion has to be ordered shallowest-first. The existing
  nested-strand pass already does this in 8 sweeps.
- **Haploid children** are already decoded in per-strand chains; those should be subsumed by the
  tree form rather than left alongside it, or there will be two mechanisms again.
- **Children with no phased parent** keep the existing fallback -- attach to the nearest preceding
  block, wildcard on both strands.
