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

## Stage 1 result

Implemented and gated on chr20 against `a27149728`. Two structural results arrived before any
accuracy number and stand on their own:

| chr20 | baseline | stage 1 |
|---|---|---|
| crossing masks the sweep could not compute | 381 | **0** |
| coherence disagreements | 3 | **0** |
| arena | 14.11 MB | 15.77 MB (+11.7%) |
| records | 116,958 | 117,048 |
| genotypes changed | 8,855 | 7,742 |
| bare haploid (strandless) GTs | 292 | 298 |
| records left unphased | 0 | 6 |
| wall / peak RSS | 168 s / 3.91 GB | 169 s / 3.18 GB |

Accuracy holds, which is what the gate asked. Recall rises in every class and precision falls by
slightly more, so the F1s move a little and mostly down:

| chr20, aardvark GT | baseline | stage 1 | delta |
|---|---|---|---|
| ALL | 0.96996 | 0.96982 | −0.00014 |
| SNV | 0.98424 | 0.98428 | +0.00004 |
| JointIndel | 0.91642 | 0.91571 | −0.00072 |
| Insertion | 0.90783 | 0.90640 | −0.00143 |
| Deletion | 0.92945 | 0.92956 | +0.00012 |
| SV (truvari ≥50 bp) | 0.51768 | 0.51301 | −0.00467 |

ALL goes 91,138/2,093/3,553 to 91,148/2,130/3,543 TP/FP/FN: ten more true calls bought with
thirty-seven more false ones. That is the "genotypes may move" risk landing as a small, *directional*
cost rather than as noise, and it is the one thing here to re-check genome-wide rather than accept
from one contig. Two residual classes are also worth naming rather than burying: 550 settled
genotypes name a traversal with no ALT and keep their called genotype, and 426 records are phased on
the line's alleles instead of the model's. Both are the same underlying fact -- the model can prefer
an allele the emitter did not write -- and stage 2 is where that stops being invisible.

The mask population is gone because a mask over *candidate traversals* does not need the parent to
have emitted anything — which also closed task #67 without it being worked on. Every coherence
bucket reads zero with no special-casing, which is the precondition stage 3 needs.

**The gate found a real defect, and the way it surfaced is worth recording.** Small-variant scoring
returned all zeros: aardvark aborted in region generation on 48 chr20 records whose GT named an
allele the record had no ALT for (`.|2` on a record with one ALT). The cause was the haploid nested
regenotyping path, which built its `Change` — and, through `nested_regenotyped`, its `PhaseCall` —
straight out of compact indices with no render step. It was the one site of four that stage 1 missed,
and nothing typed it: a compact index is a plausible small integer.

Three things came out of it, all kept:

* the render at that site, so all four sites now go through `vcf_allele_of`;
* a guard in `apply_linkage_change`/`apply_phasing` that declines any patch naming an allele past the
  ALT list, with declines counted and reported by reason rather than silently dropped;
* `render_phase_pair`, because declining a *phase* patch costs the record its strand and its phase
  set. When the settled pair has no ALT the phase now names the pair the line actually carries. That
  invariant — a phased GT is always a permutation of the line's own — is what `apply_phasing` has
  always checked, and stage 1 had started violating it on 1,627 chr20 records.

The regression test is one awk line and belongs in the suite permanently: **no GT may name an allele
index beyond its record's ALT list.**

**Stage 1 also fixed a mis-call the test suite had pinned.** `18_vg_call.t`'s nested haploid fixture
asserted the parent came out homozygous for a chain-spanning deletion, with a comment conceding the
call was wrong and only pinned because it was what the code did. In traversal space the parent comes
out heterozygous, which is right — the reads are from one deleted and one crossing haplotype — so the
chain is reachable on one strand, the nested site is called at ploidy 1, and it names its strand. The
assertions were rewritten to the correct behaviour. This is the "emitted-site genotypes may move"
risk below arriving as intended rather than as damage.

## Stage 2 result

Implemented: `record()` no longer runs inside the `add_variant` branch, so every genotyped snarl
enters the layer and `emitted` says whether there is a line to patch. `PhaseCall` carries the same
flag, so an unemitted parent is phased -- its children read their strand off it -- without entering
anything that counts or patches records. That last part is load-bearing: the mosaic's site counts are
index arithmetic over the vector it is handed, so a collapsed site inside a run would inflate it and
break the invariant that the mosaic accounts for exactly the emitted records. It gets a filtered
vector.

The gate asked for strandless haploid records to fall from 292 toward zero. They do:

| chr20 | baseline | stage 1 | stage 2 |
|---|---|---|---|
| bare haploid (strandless) GTs | 292 | 298 | **18** |
| nested haploid records *with* a strand | 1,650 | 2,155 | **2,509** |
| nested sites with no phased parent | 289 | 288 | **19** |
| records with a strand the panel does not explain | 463 | 511 | **341** |
| collapsed sites phased with no line | — | — | 101,864 |
| linkage sites / arena | 117,148 / 14.11 MB | 117,210 / 15.77 MB | 219,246 / 29.04 MB |
| wall / peak RSS | 168 s / 3.91 GB | 169 s / 3.18 GB | 205 s / 3.05 GB |

And accuracy improves in every class, over stage 1 *and* over baseline -- which also retires the
precision worry stage 1 raised, since false positives now land below where they started:

| chr20, aardvark GT | baseline | stage 1 | stage 2 | vs baseline |
|---|---|---|---|---|
| ALL | 0.96996 | 0.96982 | **0.97041** | +0.00046 |
| SNV | 0.98424 | 0.98428 | **0.98434** | +0.00010 |
| JointIndel | 0.91642 | 0.91571 | **0.91816** | +0.00174 |
| Insertion | 0.90783 | 0.90640 | **0.90829** | +0.00046 |
| Deletion | 0.92945 | 0.92956 | **0.93241** | +0.00296 |
| SV (truvari >=50 bp) | 0.51768 | 0.51301 | **0.51903** | +0.00135 |

ALL false positives: 2,093 baseline, 2,130 stage 1, **2,014** stage 2.

**The mosaic wildcard half of the gate needs its metric corrected rather than reported as missed.**
Raw wildcard segments went 437 -> 588 -> 616 and wildcard *sites* 2,402 -> 2,964 -> 2,868, so the
count rose against baseline. The cause is not worse phasing: 859 more records gained a strand than at
baseline, and a nested haploid record written `1|.` puts a `*` on its empty strand, because the
mosaic spells "no sequence on this haplotype here" and "the panel does not name a haplotype here"
with the same character. The count that means only the second thing fell below baseline, 463 -> 341.
Task #43 should not use the raw wildcard count as its metric; the two cases need distinguishing in
the mosaic format first.

**Two costs, both real.** The arena is 29.04 MB against 15.77 MB and the linkage pass 42.8 s against
18.8 s, because the site count nearly doubled -- 219,246 against 117,210, the difference being every
snarl that collapses to the reference. Wall clock is +22%; peak RSS did not move. Whether recording
*every* collapsed snarl is necessary, or only those with children, is the obvious lever if this
becomes a problem at whole-genome scale.

**A coherence class came back, and it is stage 3's to remove.** "Diploid under the settled parent"
was 0 in stage 1 and is 440 here. Nothing regressed: recording collapsed parents means far more
children now *have* a checkable parent (3,969 nested sites in the first generation against 765), so a
disagreement that was previously invisible is now counted. Deriving ploidy and strand from one
computation is what makes it unrepresentable.

## Risks

* **A settled pair with no VCF allele.** In traversal space the Viterbi can reach a traversal the
  emitted record has no ALT for — impossible today, because the space *was* the emitted alleles. The
  patch cannot add an ALT, so such a change is skipped and counted rather than applied. This is new
  behaviour and needs the count reported, not hidden.
* **Every write into a VCF allele field is a place the refactor can be incomplete.** Stage 1 renders
  compact → VCF at each of the four `Change`/`PhaseCall` construction sites, and missing one is not a
  compile error and not a test failure: it writes a plausible small integer into a GT. Stage 1 did
  miss one — the haploid nested regenotyping path, which wrote compact indices into both the genotype
  patch and, via `nested_regenotyped`, the phase patch. Audit by grepping for assignments to
  `allele_i`/`allele_j`/`allele_first`/`allele_second` and checking each one passes through
  `vcf_allele_of`. The invariant that catches it downstream is cheap and belongs in the test suite:
  **no GT may name an allele index beyond its record's ALT list.**
* **Emitted-site genotypes may move**, because the model now sees panel-carried traversals that
  collapsing previously merged. That is the point, but it makes stage 1 a real A/B rather than a
  refactor.
* **Arena growth** is bounded by argument, not yet by measurement. Stage 1's gate settles it.
* Reverting is cheap per stage: each is a separate commit with its own gate, and `a27149728` is the
  validated fallback.
