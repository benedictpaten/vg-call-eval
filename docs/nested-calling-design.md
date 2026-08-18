# Nested calling by symbolic alleles, with ploidy propagation

Design, implementation and results. Stages 0-3 and 5 are built and measured behind `vg call
--nested`; Stage 4 is deliberately last and not started. Every claim is cited so it can be
re-checked.

**Genome-wide result, autosomes, HG002 against T2T-Q100.** Both arms run with the same binary and
scored through the same `bench_wgs.py` path. The default arm reproduces the published numbers to the
last digit, which is what makes the comparison clean:

| autosomes | default | `--nested` | delta |
|---|---|---|---|
| SNV F1 | 0.9752 | **0.9833** | +0.0081 |
| SNV recall | 0.9567 | **0.9742** | +0.0175 |
| SNV precision | 0.9945 | 0.9926 | -0.0019 |
| Indel F1 | 0.9147 | **0.9189** | +0.0042 |
| ALL F1 | 0.9626 | **0.9699** | +0.0073 |
| SV F1 | 0.5134 | **0.5478** | +0.0344 |

**SNV false negatives fall from 146,786 to 87,373 -- 59,413 recovered**, against the 55,222 swallowed
SNVs predicted from the offline analysis. SNV recall of 0.9742 clears PanGenie's 0.9659 and lands
close to the 0.9737 projected from the swallowed count; the graph's own ceiling is 0.9828. SVs improve
on both axes at once, 1,071 more true calls and 263 fewer false ones.

**Phasing, autosomes, measured with whatshap against the same phased truth:**

| | pairs | switch % | hamming | hamming % | blocks | block N50 |
|---|---|---|---|---|---|---|
| default | 2,442,552 | 2.41% | 1,185,383 | 48.53% | 22 | 248.384 Mb |
| `--nested`, chains fragmented | 2,502,583 | 2.35% | 780,795 | 31.20% | **9,460** | 1.079 Mb |
| `--nested` + per-strand nested chains | 2,510,608 | 2.40% | 1,212,244 | 48.28% | **22** | **248.386 Mb** |

**The switch percentages in this table were previously quoted a hundred times too low** -- 0.0240%
against 0.0241% and so on. `phasing_benchmark.py` recomputes `all_switchflip_rate` as switches divided
by assessed pairs, a fraction, and printed it under a column headed `switch %`; this table copied the
figure and attached a percent sign. The script now converts, and prints the hamming rate beside it for
the reason below. The comparison between arms is unaffected -- all three rows were mislabelled the
same way -- and [tier2-phasing.md](tier2-phasing.md) was never wrong: it reports 2.30% for chr20 at 34
haplotypes, which is what 2.40% here should always have read as.

Block structure matches the default -- 22 blocks, N50 within 1.4 kb -- while phasing 68,056 more
variants than it does. Switch error is 2.40% against 2.41%, and that comparison only means anything
now the blocks are the same length: the fragmented run's apparently better 2.35% was the trivial win
short blocks always give, which the benchmark script warns about in its own docstring.

**What the two columns together say about long-range phase.** At 2.4% per adjacent het pair, the
relative orientation re-randomises every forty sites or so, and every switch flips everything
downstream of it -- so blockwise hamming sits at 48% whatever the block length. A 248 Mb block is a
statement about which sites share one `PS`, not about phase that is trustworthy across a chromosome.
`tier2-phasing.md` argues hamming is uninformative over long blocks, and is right; the corollary worth
stating alongside it is that the phase is then only locally meaningful, so nothing should read the
mosaic as chromosome-scale truth. At 2.4% the caller is where a 34-haplotype panel puts it -- that doc
compares it against the 0.5-2% statistical phasers reach from thousands -- and the lever is the panel.

Hamming still moves the wrong way between the last two rows: 1,212,244 against the default's
1,185,383, 2.3% worse. Restoring long blocks forces one global orientation per chromosome where
9,460 short blocks could each be locally correct and globally meaningless -- which is why the
fragmented run's 780,795 flatters it. On 68k more phased variants the per-variant gap is close to
flat, but it is not an improvement and should not be reported as one.

**Runtime and memory are unchanged**, chr20, same graph and reads: default 136.5 s / 3.64 GB against
`--nested` 135.9 s / 3.44 GB. Genome-wide both arms take the same ~55 minutes.

## Why

`vg call` emits one record per top-level snarl, and a `SnarlTraversal` runs from snarl start to snarl
end through every interior node. Nested variation is therefore already baked into each top-level
allele string, and two traversals differing only *inside* a nested chain become two long, nearly
identical ALTs. Measured on HG002 against T2T-Q100, autosomes:

| | |
|---|---|
| SNVs vg misses that sit inside a large allele vg itself emitted | **55,222** of 142,707 (38.7%) |
| the same, among variants PanGenie called | 71.2%, against a **0.6%** rate among variants vg calls correctly |
| indels likewise | 10,735 of 78,274 (19.6%) |
| vg SV "false positives" that are same-length substitutions | 2,219, against PanGenie's 29 |
| of those, differing at <=10 bases | **90.6%** — including a 4,710 bp allele differing at 3 |

One record can swallow many variants: across chr1+chr6+chr20, 349 large records account for 8,168
missed small variants, a mean of 23.4 each and a maximum of 435. The median swallowing record is
4,420 bp of REF; the largest is 210,215 bp.

The failure is not that these variants are unreported. It is **all-or-nothing at snarl scope**:
choosing one traversal commits to every nested variant inside it at once, so a slightly wrong long
allele loses all of them together. The 0.6% control rate proves the scorer credits nested variants
when the enclosing allele is right — aardvark compares by local haplotype — so the 71.2% are cases
where the allele is wrong.

Two further facts bound the work. **39.9% of vg's SNV false negatives are not in the panel at all**
and are irreducible. Of the rest, recovering the swallowed set would move SNV recall from 0.9570 to
about 0.9737, past PanGenie's 0.9659; the graph's own ceiling is 0.9828.

Full workings: `docs/small-swallowed.md`, `docs/fn-representable.md`, `docs/sv-nocall.md` and
`docs/sv-delta.md` in the companion evaluation repository.

## The idea

A **leaf snarl** contains no nested snarls; a **non-leaf snarl** contains nested chains.

For a non-leaf snarl, rewrite each traversal as a **symbolic allele**: the ordered sequence of steps
where every maximal excursion through a nested chain is replaced by a single symbol naming that
chain, rather than the concrete path taken through it. Allele identity at this snarl is then decided
by comparing symbolic sequences.

Two consequences follow:

- **A traversal symbolically equal to the reference traversal is the reference allele here.** Its
  difference lies entirely inside nested chains, and belongs to those chains' own records. The
  4,710 bp allele differing at 3 bases stops being a top-level substitution; the three SNVs become
  three nested calls.
- **A traversal that skips a chain, or crosses different chains, is symbolically distinct** and stays
  a genuine top-level call. A real deletion is still reported as a deletion.

The second half is **ploidy propagation**. Once the snarl is genotyped, each child chain's ploidy is
the number of called parent alleles that cross it: both alleles cross it, ploidy 2; one crosses and
the other deletes it, ploidy 1, and the chain is called haploid; neither crosses it, ploidy 0, no
descent, star allele. Recursion continues down the hierarchy with that ploidy.

**The read-likelihood model does not change.** What changes is how traversals are collapsed into
alleles for emission, and which snarls are visited at what ploidy.

## What already exists

This is not a new representation. The decomposition provides it:

| piece | where |
|---|---|
| a `Visit` carrying *either* a `node_id` *or* a `Snarl` — "each step is given as either a node or a child Snarl" | `deps/libvgio/deps/vg.proto:283` |
| `SnarlManager::chains_of(parent)` — the child chains to symbolise | `src/snarls.hpp:513` |
| `NetGraph`, in which "each chain and unary child snarl is treated as an ordinary node" — the symbolic view, already built | `src/snarls.hpp:245`, `net_graph_of` at `:521` |
| `into_which_snarl(visit)` — maps a visit to the child snarl it enters | `src/snarls.cpp:923` |
| a loop already projecting a traversal's visits onto child snarls | `src/graph_caller.cpp:3071` |

The deprecated `NestedFlowCaller` used Snarl-carrying Visits already — it had the representation and
lacked the ploidy discipline. Its Visit form is also why `--bottom-up -T` aborts in the GAF emitters.

### The defect this also fixes

Recursion today is a side-effect of emission. `FlowCaller::call_snarl_internal` ends with
`ret_val = trav_genotype.size() == ploidy && added` (`graph_caller.cpp:3749`), and `emit_variant`
returns `true` even when it emitted nothing, because the `genotype_snarls || !alt.empty()` gate at
`:1893` is simply skipped for a hom-ref call and control falls through to `return true` at `:1976`.
So **a snarl genotyped hom-ref writes no record, reports success, and its children are never
queued.** `RecurseOnFail` only fires where the caller could not decide at all — in practice, sites
with no reads. Under this design, recursion becomes an explicit decision driven by propagated
ploidy, not a by-product of whether a line was written.

## Decisions taken

**Arg-max before marginalisation.** Several concrete traversals can share one symbolic allele. Under
arg-max, the model runs exactly as it does now over concrete traversals, and collapsing to symbolic
happens strictly afterwards, for emission and propagation. Every quantity stays defined on a single
concrete traversal pair: the mixture weights `w_h`, the expected count `λ_G`, and hence
`DR = N_eff/λ_G`.

Marginalising over a symbolic class breaks that at the **depth term**. Class members differ inside
nested chains, which may contain indels, so they differ in length and therefore in `λ`. A logsumexp
over genotypes carrying different `λ` is a sound mixture but leaves `λ_G` undefined for the reported
call, so `DR`, the depth term and any gate on them have no single value. Marginalisation is
therefore a later stage with its own answer for `λ` — candidates are the arg-max member's `λ`, a
length-weighted average over the class, or a `λ` recomputed from the symbolic allele's
reference-projected length.

Note this is only lossy where it matters least. Symbolically equal traversals cross the same chains,
so ploidy propagation is identical whichever member arg-max picks; and where the call is symbolically
hom-ref nothing is emitted, so the parent's GQ is unused. The case marginalisation would improve is
the het one — a real deletion competing against a *class* of symbolically-reference traversals whose
probability is currently split across its members.

**Propagated ploidy is already safe for the depth model.** `local_read_rate` is reads/bp over the
read source's fetch window and is deliberately not divided by ploidy; the division happens at point
of use. A nested site called at ploidy 1 inside a heterozygous parent therefore gets the right `λ`
without further work.

**Ordering nested sites for phasing rests on the reference, which constrains the off-reference
work.** A nested record takes its position from `get_ref_interval` on the reference path, exactly as
a top-level one does, and chains are built by sorting on `(position, record_key)` and cutting at
ploidy transitions. So a parent and its children interleave by reference coordinate, and two disjoint
nested chains are ordered against each other because both are anchored on the same reference.

That is well-defined *only* because v1 descends into chains the reference also crosses. The
restriction was adopted below for a narrower reason -- REF and POS are ill-defined for the record --
but it is also what makes an ordering exist at all, and that is the stronger reason. Under
`--nested-pseudo-ref` there is no reference interval: the existing fallback hands back the *parent's*
interval, so every off-reference nested site under one parent collapses onto a single position and
two disjoint chains there are mutually unordered, tie-broken only by a snarl-derived key. Stable is
not meaningful -- the transition model uses inter-site distance, and distances taken from a shared
parent position are fictional.

So off-reference nested calling needs an ordering answer as well as a REF/POS answer. Either order
sites by offset along the parent allele's traversal, mapped back onto the parent's reference span, or
keep off-reference nested sites out of the chains and emit them unphased.

A related property, pre-existing rather than introduced here but made more common by nesting:
inter-site distance is measured along the reference while the called haplotype's own path through a
snarl may be much longer or shorter. Not measured.

**Off-reference nested content: opt-in, off by default.** A chain crossed only by a non-reference
parent allele has no reference path through it, so REF and POS for a nested record are ill-defined.
v1 descends only into chains the reference traversal also crosses — which is where the swallowed
variants are, since they are on-reference by construction. `--nested-pseudo-ref` enables descent
elsewhere against a pseudo-reference. Off-reference calling proper is intended later work.

**Repeated traversal of a chain by one allele** (cycles, tandem duplication) gives per-haplotype copy
number above 1. v1 caps at the {1, 2} the rest of the system assumes and logs the occurrence.

**Nested sites in the linkage layer and the mosaic: measured, not assumed.** Including them is more
correct for phasing but changes chain construction and the mosaic's site accounting, which the
harness asserts must equal the record count. Stage 5 runs it both ways and decides on the numbers.

**Shared-flank trimming: not needed -- vg already does it.** This was proposed here as
`--trim-shared-flanks`, on the hypothesis that a genuine 300 bp deletion in a snarl whose boundary
anchors add shared flank is emitted at 400 bp and so fails truvari's size-similarity test, and that
trimming would recover part of the 561 vg-only false negatives rejected on similarity
(`docs/sv-unmatched.md`).

**The hypothesis was wrong and the option is redundant.** `flatten_common_allele_ends` already runs
unconditionally on every record, in both directions: it trims the shared suffix, then the shared
prefix, leaves one anchor base and adjusts POS. Measured on all 836 records with REF >= 50 bp in the
`--nested` chr20 output, **100% have zero trimmable shared prefix and zero trimmable suffix
remaining**. Trimming is already maximal, so a 300 bp deletion is already emitted at 300 bp.

That also explains why long records survive at all. Trimming can only remove sequence shared at the
*ends*; the records that remain have differences near *both* ends, so the shared sequence is in the
**middle**, where flank trimming cannot reach. The test fixture is exactly this shape -- a 22 bp
allele differing at base 1 and base 22 with twenty identical bases between them.

Reaching interior shared sequence needs atomisation, which was built and measured already:
`--atomize-substitutions` gave SV F1 0.4998 -> 0.5015 on chr20-4hap and was judged not worth it.
Symbolic collapsing has since taken a much larger bite out of the same population, which weakens
atomisation's marginal case rather than strengthening it.

## Stages

Each stage has a gate. A stage that misses its gate stops the sequence rather than being carried
forward on the assumption a later one will rescue it.

### Stage 0 — Offline validation, no C++

Establish that the swallowing records are non-leaf snarls with child chains. If they are leaf snarls
holding one large complex bubble, symbolic collapsing buys nothing there and this whole design is
aimed at the wrong population.

Method, all from existing artefacts: take the snarl hierarchy from `vg snarls` on chr20, and the `AT`
INFO field of each emitted record, which already gives each allele's traversal as a node path.
Project each `AT` onto the child chains and compare symbolic sequences.

**Gate**: the records predicted to collapse to reference must account for the bulk of chr20's share
of the 2,219 same-length substitution false positives and of the 55,222 swallowed SNVs. Report the
fraction of swallowing records that are non-leaf.

### Stage 1 — Symbolic allele encoding

`symbolic_allele(trav, snarl, snarl_manager)` returning a symbol sequence, with equality and hashing.
Symbol granularity is the **child chain**, not the child snarl, and symbols carry orientation.

Unit tests on synthetic hierarchies: chain present, chain absent (deletion), chain crossed twice,
reversed orientation, trivial chains, and a nested chain inside a nested chain. Also an early probe
that panel enumeration behaves at a nested snarl under a propagated ploidy, since that is a
precondition for Stage 3 and is cheaper to find out now.

### Stage 2 — Reference-equivalence collapsing at emission

Map called traversals to symbolic alleles at emission; a traversal symbolically equal to the
reference collapses to allele 0, and only symbolically distinct ALTs are emitted.

**Gate**: chr20's same-length substitution false positives largely disappear; small-variant F1 does
not fall.

**Measured: the first half holds and the second does not, and the reason matters.** On chr20, SV
false positives fall 367 -> 282 and SV F1 rises 0.4944 -> 0.5106; small-variant F1 falls 0.9646 ->
0.9596, losing 912 variants, 585 of them (64%) inside the 660 records the collapse dropped.

Those long alleles were not pure noise. aardvark compares by local haplotype, so it was crediting
small variants carried *inside* them; collapsing the record removes that correct nested content along
with the wrong top-level allele. **Stage 2 is therefore not separately shippable**, as an earlier
draft of this plan claimed. It is demolition, and Stage 3 is the rebuild: the two are halves of one
change and `--nested` must not ship with only the first.

### Stage 3 — Ploidy propagation and an explicit recursion contract

Replace emission-driven recursion with explicit descent: per child chain, ploidy is the number of
called parent alleles crossing it; ploidy 0 gives a star allele and no descent. Fix
`emit_variant`'s `return true` on a no-emission path, so "nothing was written" can no longer read as
"this site is resolved".

New flag `--nested` enables the mode. It does not replace `-A` or `--top-down`, both of which stay
until Stage 6 decides.

**Gate**: nested SNVs appear as their own records; the swallowed-SNV count falls toward zero; record
count and peak memory stay within the scheduler's `2.25 + 11.2e-6 * records` budget, refitting the
coefficient if the record count rises materially.

### Stage 4 — Marginalised symbolic likelihood — *deferred to last, not yet started*

Only if Stages 2–3 land, which they now have. Deliberately sequenced after Stages 5 and 6: it is the
one stage that changes the model rather than the output, it needs an answer for `λ` before it can be
written at all, and everything measured so far is arg-max, so it should be A/B'd against a settled
baseline rather than a moving one. Requires a defined `λ` for a symbolic class, per the decision above. A/B
against arg-max on the four tier-2 arms.

**Gate**: het-site accuracy improves without `DR` or the depth term becoming unreportable.

### Stage 5 — Linkage, phasing and mosaic for nested sites

**Partly measured, and it found an open defect.** Nested sites already flow through the linkage
layer, because linkage is on by default under `--read-likelihood`; the question was never whether
they participate but whether they should, and what breaks.

What works: ploidy propagation reaches the output. On chr20, 2,135 records are single-allele
genotypes -- nested sites called at ploidy 1 because only one parent allele crosses them -- against
114,831 diploid and 81 unphased. None carry a missing or star allele.

What breaks: **the mosaic's site-total invariant.** The mosaic must account for every emitted record,
and the harness asserts it. On the default it holds exactly (105,251 sites, 105,251 records); under
`--nested` it does not (116,789 against 117,047, a gap of 258). So 258 nested records reach the VCF
without reaching a linkage chain, and the mosaic no longer describes the whole call set.

**Found and fixed.** `LinkageCollector::resolve` skipped any chain of fewer than two sites. That is
right for genotype resolution -- a lone site has nothing to link to -- but it also skipped the
*phasing* output, so such a site never reached `phasing_out` and the mosaic never counted it. The
guard was invisible while chains were maximal runs of one ploidy along a contig, which made a
singleton chain vanishingly rare; propagated ploidy creates them in quantity, because an isolated
ploidy-1 site between diploid neighbours is a chain of one. The guard now skips only empty chains,
and both arms account for every record: default 105,251 = 105,251, `--nested` 117,047 = 117,047.

The default path is unchanged, verified field by field on GT/DP/AD/GQ/PS across all 105,251 chr20
records. (A first comparison appeared to show 104,773 differences; that was `bcftools view -H`
reformatting floats when round-tripping the bgzipped file against a raw read of the fresh one --
a comparison artefact, not a behaviour change.)

Still to measure: whether including nested sites *helps* phasing accuracy, which needs a whatshap
comparison against the default arm.

#### Original plan

Run both ways — nested sites in the linkage chains and out of them — and decide on measured phasing
accuracy. Whichever wins, the mosaic's site-total invariant must still hold.

### Stage 6 — Shared-flank trimming, then full evaluation

`--trim-shared-flanks`, measured *after* symbolic collapsing so it is credited only with what it
adds. Then the four tier-2 arms, the whole genome, and the PanGenie comparison refreshed.

**Gates**: SNV recall past PanGenie's 0.9659, toward the 0.9737 the swallowed count implies; SV F1
up; small-variant F1 not down; runtime and memory acceptable for the 24-contig laptop run.

### Stage 7 — Coherence between a nested call and its parent's final genotype — *partly done, and one confident change was withdrawn*

Nested descent decides which children to visit, at what ploidy, and on which parent strand, from the
parent's **pre-linkage** genotype. Linkage then rewrites parents and nothing revisited those decisions.
Investigating that turned up four things. Two are fixed, one is reported rather than fixed, and one was
implemented, measured, found to make the output worse, and reverted.

#### The mechanism this rests on: a crossing mask

Descent now carries a **crossing mask** into `LinkageCollector::Entry` -- one bit per parent VCF allele,
set where that allele crosses this child chain -- so `resolve` can ask questions about the parent's
*final* genotype that it previously could not ask at all. Eight bytes a site, placed beside
`parent_record_key` rather than beside the slot it was meant to supersede, since after a `uint8_t` it
would cost sixteen. chr20's retained linkage state goes 12.31 MB to 13.20 MB.

#### 1. Fixed: 255 of 2,135 nested sites have a ploidy their parent's final genotype contradicts

The population used to be an upper bound -- nested sites whose parent's genotype moved at all, 345 on
chr20. With the mask it is exact. Under the parent genotype linkage actually settled on:

- **88 sites are diploid**: both parent haplotypes cross the child, so the locus has two alleles there
  and the record names one.
- **167 are unreachable**: neither haplotype crosses it, so the sample carries no copy of the chain
  under its own parent record and the call has no haplotype to sit on.

Neither is dropped. Both are flagged, `FILTER=nested_diploid` and `FILTER=nested_unreachable`, with
header descriptions saying what the flag means. The read evidence at the child is real, and it is
evidence *against* the parent's new genotype as much as the parent is evidence against the child, so
deleting the record would be asserting that the parent won that argument. Neither aardvark nor truvari
filters on FILTER as this harness invokes them, so scoring is unaffected and the flags are information
rather than a silent exclusion.

**This check survives the problem that sank the strand derivation below**, and for a specific reason: it
asks about the *unordered* pair. Whether both called alleles cross, neither crosses, or exactly one
does never requires knowing which allele is on which strand.

Their strand handling is honest rather than falling back to a slot that means nothing: a diploid one
names the parent's haplotype on both strands, an unreachable one names neither, and both stay out of the
per-strand chains of step three, which link within one haplotype and have no place for a site on two or
on none.

**The gap that remains** is re-genotyping the 88 at ploidy 2, and it cannot be done from what is
retained: `Entry` holds the *haploid* likelihood vector, indexed by allele, because that is what the site
was genotyped at, and a diploid genotype needs the triangular vector that was never computed. It needs
either a speculative second genotyping of every candidate child at both ploidies -- affordable, since the
candidates are 1.8% of loci and the read fetch would be shared -- or a genuine re-call from the write
path. Not started.

#### 2. Fixed: the mosaic put nested segments in the wrong place

`write_mosaic` reads the phasing as one ordered sweep per contig and closes a segment only where the
haplotype changes. Nested sites break that by construction: placing one needs its parent already phased,
so they are appended after every chain, leaving the vector in two runs.

Out of order, a nested site 451 kb into chr20 shared a run with one 65 Mb along. chr20's mosaic carried
five segments spanning tens of megabases, one claiming 284 sites between `ref_start` 451,374 and
`ref_end` 65,512,343, with `start_node`/`end_node` anchors to match. **The site totals still added up**,
which is the invariant the harness checks, so the file looked complete while being wrong about where
several hundred sites were.

Fixed by sorting the phasing into `(contig, position, record_key)` at the end of `resolve` -- one
guarantee for every consumer rather than each having to know, and unit-tested. chr20 goes from 4,563
segments with five out-of-order rows to 5,111 with none, over the same 117,047 sites and the same 2,510
wildcard site-slots.

Sorting then exposed a latent hazard next door: the strand count was read off `phasing[i].ploidy`, the
*first* site on the contig. A diploid contig can now open with a nested haploid site, which would have
dropped that contig's entire second strand -- silently, since a one-strand mosaic is exactly what a
haploid contig is supposed to look like. It takes the maximum over the run now. chr20 opens at position
24 with a diploid record, which is the only reason this was never visible.

#### 3. Withdrawn: deriving the nested strand from the parent's phased allele pair

`parent_slot` is an index into the parent's *called traversal* order, recorded at descent. By the time
the parent is phased, `record` has sorted its allele pair and the Viterbi has oriented that pair against
the panel, so slot 0 and `allele_first` coincide only by chance. Using one to index the other looked like
a category error rather than a tunable heuristic, and the crossing mask makes the correct derivation
available: ask which of the parent's *phased* alleles crosses the child.

It was implemented, and it moves **636 of 1,851** placed nested sites on chr20 -- 34%, close to the coin
flip a wrong index predicts. Then it was measured against the phased truth, and **it is worse**. Paired
on the 142 sites decisive for both runs: the recorded slot is right 106 times (74.6%), the derived slot
64 (45.1%). Guarding the derivation against parents whose allele pair has no determined order (below)
recovers most of the gap but does not close it: 94 (66.2%). So it was reverted, and the mask is kept only
for the ploidy check of item 1, which does not depend on the ordering.

Reverted rather than shipped-and-flagged because a change that moves a third of a population needs to be
right, and "my reasoning says the index is wrong" lost to a direct measurement. **Why it lost is not
resolved.** The leading explanation is item 4 -- for many parents nothing determines the order of the
pair, so indexing it is no better founded than indexing traversal order -- but that accounts for the
45.1% to 66.2% step and not the remaining eight points. A second candidate not yet tested is whether the
VCF's allele order and the mosaic's haplotype order stay aligned through the fallback branches that build
a `PhaseCall`, since the measurement reads the strand from one and the frame from the other.

#### 4. Reported: heterozygous sites whose allele order nothing determines

Building a `PhaseCall` for a diploid site takes the allele each phased panel haplotype carries. Where
neither haplotype spells either called allele, the last-resort branch writes the pair in sorted order --
and that pair is then emitted as the phased `GT`, in the block, indistinguishable from a site the panel
actually oriented. `PhaseCall::order_arbitrary` now marks those, and the run reports how many there are.

They are the reason the derivation in item 3 has less ground under it than it appears, and they are worth
a look in their own right: a phased genotype that asserts an orientation nothing chose is a guess dressed
as a call, which is the thing this caller avoids everywhere else.

#### What could not be measured, and why it matters more than any of the above

Every phasing number published for `--nested` is computed on the **diploid records alone**: whatshap
refuses a VCF of mixed ploidy, so the 2,135 haploid nested records -- precisely the ones this stage is
about -- are excluded from every switch and hamming figure. The strand is not in the VCF either. A
haploid `GT` of `1` with a `PS` says which block the site is in and nothing about which strand of it, so
there is no field for a phasing tool to read.

`scripts/wgs/nested_strand_check.py` closes enough of that gap to have produced the verdict in item 3,
by recovering the strand from the mosaic -- a nested site is a wildcard on the strand it is not on, and a
wildcard breaks a segment -- and checking it against the phased truth in a local frame. Its limits are
why item 3 says "not resolved" rather than "the old code was right":

- Only **157 of 2,135** nested calls are decisive. 1,334 have no exact `POS`/`REF`/`ALT` match in the
  truth, 95 have a homozygous truth genotype that cannot tell a right strand from a wrong one, and 549
  have no single recoverable strand because wildcard segments merge across neighbouring sites.
- The frame has to be local. At 2.4% switch error per adjacent het pair a block has no single
  orientation -- blockwise hamming is 49% -- and a first version of the check used a block-wide majority
  and measured exactly that noise.
- The recovery itself checks out: it finds 572 sites differing between the two arms where the caller
  reports 636 strands moved, the rest being sites it cannot compare.

The way to settle this properly is to put the strand in the VCF -- `a|.` and `.|a` for a haploid record
inside a diploid phase set -- so whatshap can score it directly. That is a separate change with its own
scoring risk and is the obvious next thing to do.

**On chr20, what the shipped state reports:** 88 diploid and 167 unreachable of 2,135 nested sites, none
unchecked, so 255 records flagged; 180 heterozygous sites carrying an allele order the panel does not
determine; 621 sites with a strand it does not explain, which now includes the 167 that legitimately
have none. Runtime 166.47 s against the pre-change binary's 166.61 s on the same machine under the same
ambient load, and peak RSS 3.28 GB against 3.30 GB -- so free. Both absolute figures sit well above the
135.9 s recorded earlier in this document, which was measured on a quieter machine; the same-day
pre-change run is the only comparison worth making.

**Accuracy against the pre-change binary, chr20, shipped state against pre:**

| | ALL | SNV | Indel | SV |
|---|---|---|---|---|
| before | 0.9697 | 0.9841 | 0.9156 | **0.5177** |
| after | 0.9697 | 0.9841 | 0.9157 | **0.5164** |

Seven small-variant false negatives recovered (3,535 to 3,528) against four more structural false
positives (406 to 410). The SV cost is traceable: the 255 flagged sites are held out of step three's
per-strand chains, since a site on both strands or on neither belongs to no single-haplotype chain, so
four of them keep a per-site genotype linkage would otherwise have corrected. Not a wash in the caller's
favour and reported as it came out.

Structure is untouched: records 117,047 either way, one phase block with 116,983 sites carrying `PS`,
114,907 diploid and 2,140 haploid in both. With the withdrawn derivation still in, the diploid whatshap
comparison was bit-identical to the pre-change run -- 58,807 pairs, 1,627 switches, hamming 29,008 --
which bounds the collateral of everything that stayed.

## Testing

Unit tests accompany each stage as described. Beyond those:

- **A regression test that fails on today's build**: in `test/t/18_vg_call.t`, a synthetic graph with
  a SNP nested inside a snarl that also carries a deletion, asserting the SNP gets its own record and
  the parent emits no compensating substitution. This is the whole design in one case, and it should
  be written before Stage 2 so it is red first.
- The existing 272 checks in `18_vg_call.t` must pass unchanged.
- Harness assertions in the evaluation repository on the swallowed-SNV count and the same-length
  substitution false-positive count, so a recovered win cannot silently regress.

## Risks

- **Record-count growth against the memory model.** The scheduler packs contigs under
  `2.25 + 11.2e-6 * records` GB. Nested emission raises the record count; the coefficient is fitted
  and will need refitting, and the per-contig worst case (6.1 GB on chr3 today) may move.
- **REF and coordinates for off-reference nested records**, deferred by keeping `--nested-pseudo-ref`
  off by default, but it is the thing that will make off-reference calling hard later.
- **Interaction with the linkage layer**, which currently assumes one site per position with a single
  ploidy; nested sites break both assumptions and Stage 5 exists for it.
- **Trimming being credited with symbolic alleles' work**, addressed by ordering.

## Working practice: detecting when a subcommand has finished

Three times in the investigation that produced this plan, a watcher polled for a condition that could
never become true, and the work sat finished and unnoticed — once for about an hour. The failures
were: `pgrep -f "schedule_wgs.py"` from a shell whose own command line contained that string, so it
matched itself; waiting on a `BENCH_DONE` marker the program does not emit; and `pgrep -cx`, where
`-c` is not supported by BSD `pgrep`, so the loop exited immediately and reported success.

Rules for this work:

1. **Prefer running the long job in the background directly** and letting the harness notify on exit.
   A separate watcher process is the exception, for external state the harness cannot see.
2. **Wait on a PID, never on a pattern.** Capture the pid at launch and poll `kill -0 "$PID"`. A
   pattern match can match the waiter.
3. **Never wait on a completion marker without first confirming the program emits it**, by reading
   the source or a prior run's output.
4. **Check the artefact, not the process**: a finished run is one whose output file exists and whose
   exit status was captured, not one whose process has vanished — those differ when a job is killed.
5. **Verify flags on this platform.** BSD and GNU differ in `pgrep`, `head -n -N` and `sed -i`, all
   three of which have already cost time in this project.
