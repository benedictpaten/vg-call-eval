# One allele space: the parent's genotype, not the parent's VCF record

The linkage layer currently learns about a site only if that site wrote a VCF line. That is a
presentation decision standing in for a genome fact, and it costs the thing nested calling exists to
produce: a haplotype path through a nested snarl.

## What is wrong, measured

`record()` is called from inside `emit_variant`, gated on `added` — a line was buffered. A parent
whose called alleles all collapse to the reference symbolically emits no line, so it never enters the
layer, so it has no phase, so its children have no strand to inherit.

chr20, instrumented over all four outcomes for a nested site:

| outcome | sites |
|---|---|
| placed on one strand | 478 |
| **no phased parent** | **289** |
| diploid under the settled parent | 3 |
| unreachable under the settled parent | 0 |
| crossing mask uncheckable | 0 |

The VCF carries exactly 292 strandless haploid records, which is 289 + 3. So 99% of them are there
because the parent was never phased, and none because the strand was undecidable. Genome-wide the
autosomes carry 3,060 such records.

**The naive fix does not work.** Calling `record()` for an unemitted parent under the emitted allele
numbering enters it as `0/0`: symbolic collapsing maps every one of its called traversals to allele 0,
because that is what collapsing means. A homozygous-reference site has no strand distinction, so the
child still inherits nothing. The parent's two haplotypes are distinguished *only* inside the child
chains — the collapsed numbering has thrown away precisely the information the child needs.

## The design

One allele space, and it is the genotyper's, not the VCF's: **a site-local compact set of distinct
traversals**.

* `record()` takes the called traversal pair, the panel as traversal indices, and the genotype
  likelihoods keyed by traversal indices — all three of which the caller already has in that form.
  The remap into VCF alleles disappears at the source, along with the two-numbering hazard the code
  warns about twice ("Passing the wrong one wrote past the end of the genotype vector and corrupted
  the heap").
* **Every genotyped snarl records**, emitted or not. An `emitted` flag on `Entry` says whether there
  is a line to patch. A collapsed parent becomes a genuine heterozygous site in traversal space,
  phaseable against the panel, and its children inherit a real strand.
* Sites that emitted a line additionally store a traversal→VCF-allele map, one `int8` per compact
  allele, used only at patch time.

### Why this is affordable

The obvious objection is memory: traversal counts run to ~35 on a 34-haplotype panel, and a
triangular GL vector over 35 alleles is 630 entries against today's 3–6.

It does not arise, because `build_emission` only ever indexes the GL vector at `(ai, bj)` where both
come from `allele_at(site, hap)` — a panel-carried allele — and the constraint needs only the called
pair. So the *reachable* allele set is exactly

> {the called pair} ∪ {traversals some panel haplotype carries}

which is the same order of size as today's emitted-allele set. The compact space is built from that
union; traversals that are neither called nor panel-carried cannot be spelled by any panel pair and
are dropped. Expect the arena to grow slightly, where collapsing previously merged two distinct
traversals that spell the same sequence, and not by orders of magnitude. Today's figure is 14.1 MB
for chr20's 117k sites, so this is checkable rather than assumed.

### Ploidy and strand become one computation

Today the child's ploidy comes from `child_ploidy` over the parent's *pre-linkage* genotype at descent
time, its strand from `parent_slot`, and its coherence from a third comparison of the crossing mask
against the settled pair. Three derivations of one fact, which is why they can disagree.

With the parent's settled traversals in hand there is one:

```
crossed = { s in parent's settled strands : parent's traversal on s crosses this child }
ploidy  = |crossed|
strand  = the single element of crossed, when |crossed| == 1
```

Ploidy and strand agree by construction, and `parent_slot` — whose own comment concedes it "agrees
with allele_first only by luck" — is not needed.

### What this deletes

Because the disagreement becomes unrepresentable rather than merely rare:

* `NestedIncoherence`, `apply_nested_filter`, and the three `##FILTER=<ID=nested_*>` header lines
* the `final_diploid` / `final_absent` / `mask_unknown` counters and their report
* `Entry::parent_slot`, and `parent_crossing`'s 64-allele ceiling as a correctness concern
* the panel remap added at the barrier for review finding 2, since nothing needs remapping
* `EmittedAlleles`' role as the carrier of an allele numbering for descent

## Staging, with a gate on each

1. **Collector in traversal space, emitted sites only.** `Entry`, `record`, the arenas,
   `build_emission`, `resolve_generation`, `Change`/`PhaseCall` rendering through the new map.
   *Gate:* chr20 A/B against `a27149728` — genotypes and phasing near-identical, arena size
   within a small factor. A large accuracy move here is a bug, not a result.
2. **Record unemitted sites.** Children inherit strands from parents that wrote no line.
   *Gate:* strandless haploid records on chr20 fall from 292 toward 0; mosaic wildcard segments fall.
3. **One derivation for ploidy and strand**, then delete the incoherence machinery above.
   *Gate:* no coherence FILTER can be emitted because the code no longer contains one; every nested
   haploid record carries a strand.
4. **Tests and documentation.** Unit tests on the compact-space construction and on strand derivation
   from a settled parent; the end-to-end fixture that has blocked task #56; regenerate the results
   pages.
5. **Then, separately:** per-haplotype mosaic output tracing a path per strand through nested snarls.
   This design is what makes it possible — the parent's per-strand traversal is exactly what a path
   needs and exactly what is not recorded today.

## Risks

* **A settled pair with no VCF allele.** In traversal space the Viterbi can reach a traversal the
  emitted record has no ALT for — impossible today, because the space *was* the emitted alleles. The
  patch cannot add an ALT, so such a change is skipped and counted rather than applied. This is new
  behaviour and needs the count reported, not hidden.
* **Emitted-site genotypes may move**, because the model now sees panel-carried traversals that
  collapsing previously merged. That is the point, but it makes stage 1 a real A/B rather than a
  refactor.
* **Arena growth** is bounded by argument, not yet by measurement. Stage 1's gate settles it.
* Reverting is cheap per stage: each is a separate commit with its own gate, and `a27149728` is the
  validated fallback.
