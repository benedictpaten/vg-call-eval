# One allele space: the parent's genotype, not the parent's VCF record

> **Superseded by `planning/decide-then-render.md`.** Stage 1 of this plan landed; the rest was
> replaced by the decide-then-render phase, which removes the patch machinery this plan was written
> around rather than extending it. `apply_linkage_change`, `apply_phasing`, `Change`, the `nested_*`
> FILTERs and the `unrenderable` counter no longer exist.

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
from one contig. Two residual classes are also worth naming rather than burying: settled genotypes
that name a traversal with no ALT and keep their called genotype, and records phased on the line's
alleles instead of the model's. **The figures first recorded here -- 550 and 426, and later 509 and
1,781 -- were all wrong**, because the report was inside a block guarded on the generation having
nested sites while the counters are also incremented from the diploid chain sweep, which runs
regardless. The true chr20 figures are 1,472 and 5,015, a 2.9x and 2.8x undercount. See the section
on what those sites turned out to be. Both are the same underlying fact -- the model can prefer
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

## Stage 3 result

The barrier now records *which traversal* of the parent's settled pair carries a chain, and both
ploidy and strand read off that one value: the chain is carried by this traversal, so it has one copy
and sits on that traversal's strand. The traversal rather than its index is what makes it hold --
`record` sorts the called pair and the Viterbi then orients it, so an index recorded at descent means
nothing by the time the child is placed, which `parent_slot`'s own comment conceded ("agrees with
allele_first only by luck").

Deleted, because the disagreement is now unrepresentable rather than merely rare: `NestedIncoherence`
and its three kinds, `apply_nested_filter`, the three `##FILTER=<ID=nested_*>` header lines, the
separate mirror check for ploidy-2 children, `final_diploid`/`final_absent`/`mask_unknown`, and
`Entry::parent_slot`. **Net −215 lines** (184 added, 399 removed).

Gate, both halves:

| chr20 | stage 2 | stage 3 |
|---|---|---|
| `nested_*` FILTER headers in the output | 3 | **0, and none emittable** |
| nested sites whose parent settled off their traversal | — | **0** |
| nested sites carried on *both* parent strands | 440 | 440 |
| nested sites with no phased parent | 19 | 19 |
| bare haploid (strandless) GTs | 18 | 18 |
| records whose strand the panel cannot explain | 341 | **239** |
| mosaic wildcard sites | 2,868 | **2,767** |

Accuracy is unmoved from stage 2 to four decimal places on every class except SV (0.51903 →
0.51838, four calls' worth), which is what deleting a FILTER-only mechanism should do: it never
touched a genotype.

**Two things the single derivation settled that the three could not.** The "parent settled off the
traversal this chain hangs from" bucket reads **zero** in every generation -- so the 440 that used to
be reported as a coherence disagreement were never inconsistent at all. They are chains the parent
carries on *both* its settled traversals: a genuinely diploid locus where the record names one allele
because it was genotyped at ploidy 1 and the barrier had no ploidy-2 answer kept for it. The old code
called that `nested_diploid` and filtered it; it is better described than flagged, and the mosaic now
names both haplotypes there instead of marking the strand unexplained -- which is where 101 of the
recovered wildcard sites come from.

The deleted mirror check was testing `allele_first`, a VCF allele, against a mask over candidate
traversals. It was the last of this refactor's two-numbering bugs and had been reporting on the wrong
axis since it was written.

**The residual 18 strandless records split cleanly**, which is the point of naming the classes: 11
have no phased parent at all, and 7 are the both-strands case whose GT cannot say "on both" at ploidy
1. Neither is a coherence failure. Closing the first needs the parent to be reachable at all; closing
the second needs a ploidy-2 answer retained for every chain, which is a genotyping decision rather
than a phasing one.

## Every incoherence class is zero, and they were all one bug

Stage 3 left four residual classes and this closes all of them. They turned out to be one mistake
wearing four hats: **the existence of a buffered VCF line was standing in for whether the linkage
layer needed updating.** Harmless while only line-bearing sites entered the layer; wrong from the
moment stage 2 let collapsed sites in.

| chr20 | stage 3 | fixed |
|---|---|---|
| carried on both parent strands | 440 | **0** |
| carried on neither | 0 | **0** |
| no phased parent | 19 | **0** |
| bare haploid (strandless) GTs | 18 | **0** |
| records with no phase at all | 75 | **0** |
| nested sites placed on exactly one strand | 5,433 of 5,892 | **6,716 of 6,716** |
| records | 117,097 | 116,965 |

So a mosaic path exists for each strand through every nested snarl, which is what the exercise was
for. And accuracy rises in every class against both stage 3 and the pre-refactor baseline:

| chr20, aardvark GT | baseline | stage 3 | fixed |
|---|---|---|---|
| ALL | 0.96996 | 0.97040 | **0.97048** |
| SNV | 0.98424 | 0.98433 | **0.98436** |
| JointIndel | 0.91642 | 0.91816 | **0.91840** |
| Insertion | 0.90783 | 0.90829 | **0.90843** |
| Deletion | 0.92945 | 0.93241 | **0.93266** |
| SV (truvari >=50 bp) | 0.51768 | 0.51838 | **0.51875** |

ALL goes to 91,154/2,008/3,537 TP/FP/FN against a baseline 91,138/2,093/3,553 -- better on all three
axes at once, so the 132 records cascading retraction drops were false positives.

### The four causes

1. **440 "carried on both parent strands".** The barrier *did* reach these, *did* compute
   `copies == 2`, *did* find the retained ploidy-2 answer and *did* re-render with it -- then threw
   the result away, because that answer collapses to the reference, so no line was buffered and
   `last_emitted.buffer_thread < 0 || !wrote` returned before `respecify`. **The theory recorded in
   the stage-3 section above -- that no ploidy-2 answer had been retained -- was wrong**, and was
   disproved by reproducing the class on a fixture where adding `--genotype-snarls`, whose only
   relevant effect is to force a line to be wanted, revised both chains correctly off
   `alt_ploidy_info`.
2. **19 "no phased parent".** Retraction never cascaded: dropping a chain left its descendants
   emitted and pointing at an entry that no longer existed. This also explains the class's shape --
   zero at generation 1, because a generation-1 child's parent is top level and never a pending
   record, so it could never be the thing retracted.
3. **75 records with no phase.** Two independent bugs in `respecify`: it never updated
   `Entry::emitted`, and it never updated `Entry::position` -- while re-emitting at a different
   ploidy changes the emitted allele set, and `flatten_common_allele_ends` advances POS by the prefix
   every allele shares. The patch indices are keyed on (contig, POS), so those patches were not
   declined, they were never looked up.
4. **A hole stage 2 opened**: the `copies == 0` retraction was conditional on a line existing, so a
   line-less chain the settled parent does not carry stayed in the layer at a contradicted ploidy.

### What this says about the design

The user's premise -- "the parent's traversals are always fully decided before the child is
genotyped" -- is **not** what the code does, and every one of these bugs lived in the gap. Children
are genotyped *and emitted* during the read sweep from the parent's pre-linkage genotype; the barrier
only retro-fits, re-deriving `copies` and repairing what disagrees. Selecting between two precomputed
answers is sound (the ploidy-2 answer really is always there), but the *retro-fit* has to be
complete, and four separate places assumed a line was the unit of repair.

The structural alternative -- retain every nested chain, emit nothing until the ploidy is settled,
and `record()` once at the settled ploidy with no `respecify` at all -- would make all four
impossible by construction rather than fixed one at a time, and would delete `respecify`,
`blank_buffered_line`, the tombstone branch in `write_variants`, and the buffer handles on both
`PendingRecord` and `EmittedAlleles`. It costs no extra read passes, since the `CallInfo` is already
retained. That is the version worth building if this area is touched again.

Counters now split by whether a line exists: a figure mixing records with entries that were never
records cannot be read as a defect count, which is how 440 sites looked like a coherence
disagreement when 433 of them had no line at all.

## What the unrenderable genotypes are, measured

The class where the model settles on a traversal the record has no ALT for was carried for three
stages as a curiosity with a reported figure of a few hundred. Both halves of that were wrong.

**The count was undercounted 2.9x by a reporting bug.** `++unrenderable` is incremented from the
diploid chain sweep as well as the nested one, but the report sat inside
`if (phasing_out != nullptr && !deferred_nested.empty())`, so every pass without nested sites
incremented and never printed. chr20's true figure is **1,472 events at 1,465 distinct positions**,
against a reported 507. The phase-fallback counter was undercounted the same way, 5,015 against 1,781.

**And they are not a curiosity.** Every one is a genotype change the layer wanted and could not
apply -- the `best == before` test comes first, so a site only reaches the render attempt if the
model disagreed with the call. Joining the dumped positions against aardvark's per-record verdict:

| group | judged | FP rate |
|---|---|---|
| the model wanted to move it and **could not** | 545 | **61.7%** |
| the model wanted to move it and **did** | 1,003 | 11.0% |
| every judged record | 90,908 | 1.8% |

The unmatched comparison is confounded, since the layer only touches low-confidence sites, so the
control is GQ-matched: **336 false positives observed against 79.3 expected at matched confidence,
a 4.24x enrichment**, consistent across every GQ bin (63.5% against 20.7% at GQ<1, 64.0% against
14.3% at GQ<20). **455 of chr20's 2,008 false positives -- 22.7% -- sit at these 1,465 positions.**
Mostly indels (1,144), then SVs (263), few SNVs (122).

So the earlier speculation that the frequency prior was overruling the reads, and that dropping the
change was therefore a feature, is **not supported**: where the layer can act it lands at 11% FP,
and where it cannot the surviving call is wrong 62% of the time.

### Which fix, and why the obvious one does not reach it

**1,362 of the 1,465 are top-level diploid sites, not nested chains.** Adding the ALT at render time
needs the retained `CallInfo`, which exists only for nested chains, so it covers 103 of them. The
option that reaches the bulk is to widen the ALT list *at emission*: emit an allele for every
panel-carried traversal, not only the called ones, so every choice the layer can make is renderable
by construction. The data is already there at that point -- the genotyper scores every candidate
traversal, which is exactly why `AD` does not sum to `DP`. The costs are a larger VCF and the return
of records carrying an ALT no genotype names, which is the `INFO/NGT2` problem the single-sweep
design removed.

Two limits on the measurement, both real. Only 545 of the 1,465 sites fall inside the benchmark
confident regions, so 63% are unjudged and not at random. And "FP" says the emitted call is wrong; it
does not prove the model's preferred allele is right, only that the status quo is bad. The 11% figure
where the layer can act is the best available evidence that its preferences are good.

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
