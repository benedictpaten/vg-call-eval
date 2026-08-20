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

#### 3. Fixed: the strand reaches the VCF, and it made the file more accurate

A nested haploid record used to be a bare `GT` of `1` plus a `PS`, which says which block the site is in
and nothing about which strand of it. The strand existed in the caller and in the mosaic and nowhere a
consumer would look -- and whatshap refuses a mixed-ploidy file, so the records nested calling creates
were excluded from every phasing figure ever published for it.

It is now written as `a|.` or `.|a`: the position in the pair is the strand, and `.` is the haplotype
that carries nothing, because the parent's other allele deletes the chain there. `.` rather than `*`,
which would say "absent because something deleted it" and say it more precisely -- but `*` is an ALT
allele, and adding one changes the arity of `AD`, `GL` and `GQI`, all written long before the strand is
known.

On chr20, 1,590 of 2,135 nested records carry a strand (285 `a|.`, 1,305 `.|a`). The other 604 correctly
do not: 255 are the flagged sites of item 1 and the rest have no reachable phased parent.

**It also removed 34 false positives**, which was not the point of the change and is the better argument
for it -- a bare haploid genotype on a diploid contig was being read as something it is not:

| | ALL F1 | FP | SNV F1 | Indel F1 | Indel precision | SV F1 |
|---|---|---|---|---|---|---|
| bare `1` | 0.9697 | 2,166 | 0.9841 | 0.9157 | 0.9090 | 0.5164 |
| `a\|.` | **0.9699** | **2,132** | **0.9842** | **0.9163** | **0.9100** | 0.5164 |

#### 4. Not resolved, and the measurements that failed to resolve it

Which strand a nested site belongs to is *determined* rather than estimated -- the parent's other allele
deletes the chain -- so it looked like the one thing here that could be got exactly right.
`parent_slot`, the index recorded at descent, indexes the parent's **called traversal** order; by the
time the parent is phased, `record` has sorted the pair and the Viterbi has oriented it against the
panel, so slot 0 and `allele_first` coincide only by chance. Deriving the strand from the phased pair
instead is the obvious correction, and it moves 636 of 1,851 sites.

**Three measurements, and the first two were wrong.** Recorded here in order because each correction
matters more than the conclusion:

1. **Site-level check against the phased truth, strand recovered from the mosaic.** Paired on 142 sites:
   recorded slot 74.6%, derived 45.1%. Read as "the derivation is much worse", and the derivation was
   reverted on it.
2. **Same check with the strand read from the VCF** once item 3 made that possible -- no recovery
   losses. Recorded slot 73.9% (116/157), derived 47.1% (74/157). Consistent with (1), which is what
   made it trustworthy.
3. **The constant-strand control on the same 157 sites, which should have been run first.** "Put every
   nested site on strand 1" scores **74.5%** -- better than the recorded slot's 73.9%. The recorded slot
   is 1 for about 82% of nested sites, so its apparent accuracy was the base rate of that skew and not
   information. The derived slot is roughly balanced, which is the whole reason it scored 47%: against a
   one-sided truth subset, a balanced guess loses to a constant. **Neither convention beats a constant,
   so this check cannot compare them.**

**The instrument that can** is relative phase, which a constant cannot game, and item 3 is what made it
available. whatshap over the call set including nested sites:

| | pairs | switches | hamming |
|---|---|---|---|
| recorded slot | 58,938 | 1,655 | 29,045 |
| derived from the phased pair | 58,934 | 1,661 | 29,087 |

Six switches apart on 58,900 pairs. **The two conventions are indistinguishable**, so there is no
evidence to change the behaviour and the revert stands -- on the absence of a difference rather than on
the difference the site-level check appeared to show.

**What the same numbers do say is that nested strands are badly assigned either way.** Against the
diploid-only comparison of 58,807 pairs and 1,627 switches, the nested sites add 131 assessed pairs and
28 switches: a **21% switch rate against the 2.77% baseline**. That is a real, bias-free signal, and it
points at the mechanism rather than the indexing. The candidates, in order of suspicion: the parent's
own phase is only 2.8%-switch accurate and the child inherits it; the 180 heterozygous sites of item 5
whose allele pair has no determined order at all; and the 82% slot-1 skew, which nothing in the design
predicts and which no one has explained.

#### 5. Reported: heterozygous sites whose allele order nothing determines

Building a `PhaseCall` for a diploid site takes the allele each phased panel haplotype carries. Where
neither haplotype spells either called allele, the last-resort branch writes the pair in sorted order --
and that pair is then emitted as the phased `GT`, in the block, indistinguishable from a site the panel
actually oriented. `PhaseCall::order_arbitrary` now marks those, and the run reports how many there are.

They are the reason the derivation in item 3 has less ground under it than it appears, and they are worth
a look in their own right: a phased genotype that asserts an orientation nothing chose is a guess dressed
as a call, which is the thing this caller avoids everywhere else.

#### The measurement, and what it still cannot see

`scripts/wgs/nested_strand_check.py` is the site-level check of item 4, kept together with the control
that invalidated it. It reads the strand from the VCF where item 3 put it and falls back to recovering it
from the mosaic -- a nested site is a wildcard on the strand it is not on, and a wildcard breaks a
segment. Its limits, since they bound anything built on it:

- Only **157 of 2,135** nested calls are decisive. 1,344 have no exact `POS`/`REF`/`ALT` match in the
  truth, 97 have a homozygous truth genotype that cannot tell a right strand from a wrong one, and the
  rest have no strand to test.
- The frame has to be local. At 2.4% switch error per adjacent het pair a block has no single
  orientation -- blockwise hamming is 49% -- and the first version used a block-wide majority and
  measured exactly that noise.
- **Its verdicts are not usable for comparing conventions at all**, per the control in item 4. It is
  retained because that is worth knowing, not because the percentages mean anything on their own.

The strand recovery from the mosaic did check out, for what it is worth: it found 572 sites differing
between the two arms where the caller reported 636 strands moved, the rest being sites it cannot compare.

**On chr20, what the shipped state reports:** 88 diploid and 167 unreachable of 2,135 nested sites, none
unchecked, so 255 records flagged; 180 heterozygous sites carrying an allele order the panel does not
determine; 621 sites with a strand it does not explain, which now includes the 167 that legitimately
have none. Runtime 166.47 s against the pre-change binary's 166.61 s on the same machine under the same
ambient load, and peak RSS 3.28 GB against 3.30 GB -- so free. Both absolute figures sit well above the
135.9 s recorded earlier in this document, which was measured on a quieter machine; the same-day
pre-change run is the only comparison worth making.

**Accuracy against the pre-change binary, chr20, cumulative over items 1 to 5:**

| | ALL | SNV | Indel | SV |
|---|---|---|---|---|
| before | 0.9697 | 0.9841 | 0.9156 | **0.5177** |
| ploidy check + mosaic order | 0.9697 | 0.9841 | 0.9157 | 0.5164 |
| **+ strand in the VCF** | **0.9699** | **0.9842** | **0.9163** | 0.5164 |

The middle row cost four structural false positives (406 to 410) for seven small-variant false negatives
recovered, and that cost is traceable: the 255 flagged sites are held out of step three's per-strand
chains, since a site on both strands or on neither belongs to no single-haplotype chain, so four keep a
per-site genotype linkage would otherwise have corrected. The strand emission then took 34 false
positives back off, which more than pays for it on small variants and leaves SV where it was.

Structure is untouched: records 117,047 either way, one phase block with 116,983 sites carrying `PS`,
114,907 diploid and 2,140 haploid in both. With the withdrawn derivation still in, the diploid whatshap
comparison was bit-identical to the pre-change run -- 58,807 pairs, 1,627 switches, hamming 29,008 --
which bounds the collateral of everything that stayed.

### Stage 8 — Genotype nested chains *after* linkage has settled their parents — *Stage 0 measured, gate missed on its letter*

Stage 7 made the ploidy incoherence exact and flagged it. The incoherence itself is still there, and it
is structural rather than incidental: descent decides a child's ploidy from the parent's **pre-linkage**
genotype, and linkage then rewrites parents. The proposed fix is to reorder — score the top-level snarls,
let linkage settle them, then descend. Greedy, and coherent by construction.

Two things had to be measured before touching the driver, because the design's cost and its payoff both
turn on numbers nobody had.

#### Descent depth: 6, not 2

Instrumented on the current default, chr20, 27,400 symbolic child calls:

| depth | child calls | share | at or below |
|---|---|---|---|
| 1 | 21,839 | 79.70% | 79.70% |
| 2 | 4,464 | 16.29% | 96.00% |
| 3 | 989 | 3.61% | **99.61%** |
| 4 | 98 | 0.36% | 99.96% |
| 5 | 8 | 0.03% | 99.99% |
| 6 | 2 | 0.01% | 100% |

**The gate was depth ≤ 3, and on its letter it is missed.** On its substance it is met: the gate was a
proxy for how many barriers a level-synchronised design pays for, 99.61% of descents sit at depth ≤ 3,
and the deep tail is 108 calls on a chromosome.

What the tail costs is one `resolve` pass per level. On chr20 that pass is 4.24–4.49 s of a 205 s run, so
six full-chain passes would be about 27 s — **+11% of runtime for the last 0.4% of descents**. Paying it
in full is not the answer and neither is capping the depth, which would leave 1,097 descents (4.0%)
incoherent. The answer is that a level-*k* pass does not need the whole chain: linkage decays over 10 kb,
`window_posteriors` already does windowed inference, and levels 2–6 add 4,464, 989, 98, 8 and 2 sites
respectively. Restricting each pass to the windows around its own new sites makes every level after the
first nearly free. That is now part of Stage 1 rather than a later optimisation, and it is the depth
measurement that made it so.

#### The children nothing has ever counted

Descent skips a child when no called parent allele crosses it (`copies <= 0`) and never revisits the
decision, so a parent that linkage moves onto an allele which *does* cross the chain leaves a call
nobody makes. It cannot be flagged — there is no record to flag — and it appears in no total. chr20:

| child descents skipped | | |
|---|---|---|
| no reference path through the chain | 12,359 | not a coherence problem; this is the `--nested-pseudo-ref` population |
| no called allele crosses it | 2,620 | |
| — parent linkage could not move | 0 | |
| — crossing mask unreadable | 0 | |
| — still uncrossed by the final genotype | 2,324 | correctly skipped, before and after |
| — **crossed by the final parent genotype** | **296** | 208 gaining one copy, 88 gaining two |

**All 296 hang off a parent whose genotype linkage actually rewrote.** That matters because a set mask
bit is not by itself evidence of a change: two traversals can flatten to the same VCF allele, and the
reference traversal crosses every chain descent considers, so bit 0 is set at essentially all of them.
Cross-checking against the parents linkage moved was meant to separate the real class from that
artefact, and it found the artefact to be empty here.

So the reordering reaches **more than twice** the population Stage 7 flagged: 255 records at a
contradicted ploidy, plus 296 decisions never made at all. 45.2% of descents emit a record (12,383 of
27,400), so the 296 are worth roughly 134 new records on chr20. Scaling by the flagged population, where
chr20's 255 is 1/29.2 of the genome's 7,458: about **8,700 gained descents and 3,900 new records
genome-wide**, against 7,458 flagged — 0.23% of the 5,041,066 records either way.

The largest population here is none of those. **12,359 children are skipped because the reference does
not cross them**, 4.7× everything else on this table, and reordering does nothing for them. That is the
off-reference decision recorded under [Decisions taken](#decisions-taken), now with a number against it
for the first time.

#### One thing the reordering gets for free

Stage 7's remaining gap was re-genotyping the 88 `nested_diploid` sites, and it is blocked by what is
retained: `Entry` holds the *haploid* likelihood vector because that is the ploidy the site was called
at, and a diploid genotype needs the triangular vector that was never computed. Deciding ploidy before
genotyping removes the problem rather than solving it — the vector computed is the right one.

#### Stage 1: `--nested-after-linkage`, implemented and measured on chr20

Default off. About 570 lines across `graph_caller.{cpp,hpp}`, `linkage_model.{cpp,hpp}` and
`call_main.cpp`, nearly all of it moving existing code: nothing in the HMM, the per-site genotyper, or
emission.

**How it is put together.** `resolve()` becomes `resolve_generation(k, last, ...)`. Sites of a later
generation do not exist yet -- the caller has not descended into them. Sites of an *earlier* generation
stay in the chains but are **clamped**: their emission becomes a delta at the genotype they settled at
and their phase is pinned to the pair already emitted for them, so they still carry transition context
for this generation -- a generation alone is far too sparse for a 10 kb decay -- while being unable to
move. `build_emission` already maps a non-finite likelihood to zero mass, so the clamp needed nothing
from the model; the phase pin generalises the mechanism the window seams already used.

Descent queues a `PendingDescent` instead of recursing, but **only where linkage can still move the
parent** -- a snarl with no linkage entry has a final genotype already, so its children are visited
inline exactly as before. On chr20 that keeps 11,166 of the descents out of the barrier. The queued
entry deliberately retains no traversals: the crossing mask against the parent's settled allele pair
gives both the copy number and the strand, which is all the child call takes.

One thing the barrier forces: the child's **strand** now comes from the settled allele pair rather than
from `parent_slot`, the traversal-order index, because the traversals are gone by then. Those two
conventions were measured against the phased truth in Stage 7 and are indistinguishable -- 1,655
switches against 1,661 -- so this is a forced choice that costs nothing known.

**The loop, chr20:**

| pass | sites settled | genotypes changed | seconds | then descended |
|---|---|---|---|---|
| generation 0 | 109,966 | 6,490 | 4.00 | 14,086 queued, 11,866 called, 2,220 not carried |
| generation 1 | 5,550 | 1,799 | 3.39 | 2,600 queued, 2,319 called, 281 not carried |
| generation 2 | 1,203 | 483 | 3.31 | 607 queued, 561 called, 46 not carried |
| generation 3 | 218 | 148 | 3.27 | 60 queued, 58 called, 2 not carried |

**The resolve cost is flat in the number of sites it settles** -- 218 sites cost the same 3.3 s as
5,550 -- which is the full-chain rebuild the Stage 0 depth measurement predicted, and the direct
argument for the windowed level-*k* pass. Total linkage time 17.2 s against 4.2 s inline.

**Coherence, which is the whole claim:**

| chr20 | inline | deferred |
|---|---|---|
| `nested_diploid` | 88 | **0** |
| `nested_unreachable` | 167 | **0** |
| children hanging off a parent linkage moved afterwards | 838 | **0** |

Zero with 4,782 opportunities, not zero because the check could not fire. No child was lost to a
missing parent either -- that counter reports nothing, which is to say zero of 14,086.

**What changed in the output:**

| | |
|---|---|
| records | 117,047 -> 116,960 |
| records lost | 198, of which 165 were `nested_unreachable` |
| records gained | 111 |
| **the genotype itself differs** | **747 (0.64%)**, including 94 ploidy flips 1<->2 |
| only the phase orientation differs | 21,862 (18.7%) |
| phase sets differing | 20 |

The genotype change is where the design aimed: 0.64%, against 255 flagged plus 296 never-made
decisions predicted in Stage 0. Scored against the benchmark, chr20 moves by a hair in the right
direction on every class -- ALL F1 (GT) 0.96986 -> 0.96990, SNV 0.98412 -> 0.98419, SV >=50 bp 0.51309
-> 0.51434 -- and all of it is precision: small-variant FP 2,134 -> 2,104 and SV FP 410 -> 407, against
TP 91,156 -> 91,135. Fourth-decimal changes, reported because the alternative is to leave the
impression that the reordering is free of them. The 165 vanished `nested_unreachable` records are that population
disappearing rather than being flagged, and the 111 gained records are consistent with 296 gained
descents at the 45% rate at which a descent emits a record at all.

**The phase orientation churn is large and it costs nothing.** Generation 0's Viterbi runs over a
chain missing the 7,081 deferred sites, so its path re-routes: 18.7% of records come out with their
strands the other way round. Only 20 phase sets differ, so those are re-routings *inside* a block --
exactly what a switch rate measures, and the reason this had to be measured rather than argued about.
Against the phased Q100 truth, chr20:

| | inline | deferred |
|---|---|---|
| assessed pairs | 58,937 | 58,947 |
| switches | 1,655 | **1,635** |
| switch rate | 2.8081% | **2.7737%** |
| blocks | 1 | 1 |
| block N50 | 66.208 Mb | 66.209 Mb |
| phased variants | 69,893 | 70,087 |

Twenty fewer switches, and the block structure is unchanged. On 1,655 switches the counting error alone
is about 41, so **the two are indistinguishable and the point estimate is slightly the better one**. The
churn is a re-labelling of which strand is which, not a loss of phase.

A second coherence gain fell out of preparing that comparison. A nested haploid record written as a
bare `GT=1` names no strand, so nothing can place it and whatshap cannot read it: **545 such records
inline, 279 deferred**. Deciding the parent's genotype first halves the population that has no strand to
inherit.

Measuring it needed the call set made strictly diploid first, which is worth recording because two
separate things blocked it: the whole-genome truth VCF carries haploid chrX and chrY, so it has to be
subset to the contig under test, and `--half-missing ref` rewrites `a|.` but leaves a bare `GT=1`, which
still trips whatshap's uniform-ploidy check. The bare-haploid records are dropped rather than guessed --
545 and 279 of about 117,000 -- and the half-missing ones get 0 on the empty strand, which is the
existing approximation and not a new one.

**Cost: runtime 219.1 s against 205.2 s, +6.8%. Memory is not yet measured properly and the first
reading should not be used.** The deferred run peaked at 3.86 GB against the shipped run's 2.94 GB, but
the control that would attribute that -- the same binary with the flag off -- was run across a machine
sleep, so both its wall time (1,540 s) and its 3.35 GB peak are unusable. What that spoiled reading does
suggest is that some of the increase is *not* deferral at all, since the inline path is byte-identical
in output and has no queues; `Entry` growing 8 bytes a site for the generation and settled genotype
accounts for under a megabyte, so it is not that either. This needs a clean back-to-back pair before any
number is quoted, and it matters: the whole-genome scheduler packs contigs on a fitted
`2.25 + 11.2e-6 * records` GB model, so a real 30% miss would exhaust memory on chr1.

**The bit-identity test the plan promised does not exist.** At `--linkage-weight 0` the linkage pass
returns before producing any phasing, so deferred descent would have no settled parent to read and
would drop every deferred child; the two orders cannot be compared that way. Nor are the generation-0
records guaranteed identical with linkage on, since a deferred generation-0 chain is 6% sparser than
the inline one. What is exact, and was verified, is that **inline mode is byte-identical to the shipped
run** -- 117,047 records with identical content and an identical mosaic -- which is the regression that
protects the shared write path. Re-verified with the final binary after every edit: 0 differing records
and 0 differing mosaic lines.

#### Stage 9: post-linkage descent costs half again as much read I/O, and why

Stage 8 shipped the coherence guarantee and the whole-genome run confirmed it holds. It also
confirmed the price, and the price is disqualifying:

| genome-wide, 24 contigs | inline | post-linkage descent |
|---|---|---|
| coherence flags | 7,458 | **0** |
| records | 5,041,066 | 5,037,529 |
| ALL F1 (autosomes) | 0.97034 | 0.97032 |
| SNV F1 | 0.98373 | 0.98371 |
| Indel F1 | 0.91952 | 0.91946 |
| SV >=50 bp F1 | 0.54854 | **0.54901** |
| **reads fetched** | **607,088,639** | **903,322,552 (+48.8%)** |
| gbz-base subprocess spawns | 58,229 | 78,823 (+35.4%) |

Every class shows the same signature -- precision up, recall down -- and lands in the fifth decimal
except SV, which gains in the fourth. **The accuracy case is a wash, so the change rests entirely on
the coherence guarantee, and at +48.8% read I/O it does not rest.** On a four-haplotype panel the
sign even flips: chr20-4hap ALL F1 0.950239 -> 0.950177, the same rich-panel asymmetry nested calling
itself showed when it shipped.

The cost is locality, not work. Deferred descent makes *fewer* child calls than inline (14,804
against 27,400 on chr20) but makes them in separate contig-wide passes -- a median of 1 + 4 sweeps,
up to 1 + 6 -- and every window each pass touches re-spawns `gbz-base` against a 22 GB database. The
cache hit rate stays at 99% throughout, so nothing thrashes; there are simply five sweeps where there
was one. Isolated that costs 7% (chr20 205.2 s -> 219.5 s); under the three-way packing the scheduler
uses it costs about 2x (chr20 205 s -> 416 s), because three processes each do five sweeps against one
database.

**The fix is to stop re-reading rather than to read less**, and what makes it possible is that the
expensive object -- the per-read per-allele likelihood matrix -- does not depend on ploidy, while
ploidy is the *only* thing an ancestor's genotype determines about a chain. Children are genotyped
over their own full traversal set and never constrained to the parent's allele, which is the
`--top-down` lesson, so scoring every chain once at both ploidies answers every question the barrier
can ask, at every depth. Two measurements sized that redesign before any of it was built.

**How many chains an unconditional sweep must visit**, counted from `vg snarls` on chr20 plus the
reference path's 2,031,992 nodes. The reconstruction validates exactly against the caller: 165,408
top-level snarls, the number the log reports.

| reference-crossed chains, chr20 | |
|---|---|
| top level | 165,408 |
| nested, depths 1-6 | 40,922 / 12,986 / 2,857 / 609 / 25 / 2 |
| **nested total** | **57,401** |

Descent considers 30,020 of those today -- 27,400 descended plus 2,620 skipped for want of a crossing
allele -- so 27,381 sit in subtrees it prunes, and an unconditional sweep visits 1.91x the nested
chains. That sounds disqualifying and is not: chain *visits* are the unit of work, and those go
195,428 -> 222,809, **+14.0%**, or at most +29 s at chr20's measured 1.049 ms per visit. Against a
five-sweep penalty of 14.5 s isolated and 211 s scheduled.

The count also produced a sharper rule than "everything the reference crosses". A chain can only be
reached by a genotype the barrier can settle on, and the barrier chooses among the parent's *candidate
traversals* -- a bounded set, typically two to four, already in hand. So descend into a child iff any
parent candidate crosses it: exact, and it prunes far harder than reference-crossing, which admits
every chain the reference passes through whether or not the sample's alleles could go there.

**What retaining a chain to the barrier costs**, measured by sizing the actual objects rather than
pricing them -- `LinkageCollector::bytes()` exists because that same argument was once arithmetic
there too. chr20, every nested chain descent makes: 27,404 chains, 65,409 traversals, 746,269 visits,
211,167 genotype likelihoods, **87.0 MB**, or 3.18 kB a chain. At the ceiling of all 57,401 chains
that is 182 MB on chr20 and roughly 550 MB on chr1 -- 5% and 8.6% of their peaks. Affordable, so the
fallback of retaining only for children of movable parents is not needed.

Two things that fell out of the measuring rather than the design. The per-chain figure is 3.18 kB and
not the 2.52 kB the first run reported, because the first reporting placement printed at the end of
the calling sweep, before the deferred passes had run, and so counted 11,166 chains of 27,404 -- an
understatement of 2.3x in the number the whole design rests on. And retaining the `CallInfo` itself
turns out to remove the obstacle that got the earlier re-genotyping attempt abandoned: `emit_variant`
can then be called unchanged at render time, where reconstructing `AD`, `DP`, `GL`, `GQ`, `GQI`,
`GQN`, `GP` and `DR` from a compact form would have needed the parallel rendering path task #45 was
dropped over.

#### Stage 10: one read sweep, and the coherence guarantee kept

Stage 9 established that post-linkage descent works and costs half again as much read I/O, and sized
the way out. This is that redesign, and it holds both properties at once.

**What makes it possible.** The expensive object is the per-read per-allele likelihood matrix, and it
does not depend on ploidy -- while ploidy is the *only* thing an ancestor's genotype determines about
a chain, because children are genotyped over their own full traversal set and never constrained to
the parent's allele. So scoring every chain once, at both ploidies, answers every question the barrier
can ask, at every depth, from one visit to the reads.

Three phases. Descent runs inline in the single sweep the reads are resident for, and each nested
chain keeps its genotyping rather than emitting immediately. Between generations the barrier asks what
the parent's *settled* genotype implies: `respecify` moves a chain to the other ploidy before its own
generation resolves, `retract` drops one the parent turns out not to carry, and a chain the settled
parent reaches for the first time is rendered from what was kept. On chr20 that is 459 chains revised,
302 reachable only under the settled parent, 190 retracted.

That 302 is worth its own line. The five-sweep arm found **296** of the same population by an unrelated
route -- post-hoc crossing masks against already-emitted records -- so two mechanisms with nothing in
common agree to within six.

**Whole genome, 24 contigs, against both prior arms:**

| | inline | five-sweep | single sweep |
|---|---|---|---|
| coherence flags | 7,458 | 0 | **0** |
| reads fetched | 607,088,639 | 903,322,552 | **609,856,118 (+0.5%)** |
| records | 5,041,066 | 5,037,529 | 5,037,820 |
| ALL F1, autosomes | 0.97034 | 0.97032 | 0.97031 |
| SNV F1 | 0.98373 | 0.98371 | 0.98370 |
| Indel F1 | 0.91952 | 0.91946 | 0.91941 |
| SV >=50 bp F1 | 0.54854 | 0.54901 | 0.54861 |
| peak RSS vs inline | 1.00 | -- | median 1.02, range 0.74-1.42 |

Accuracy is the wash it has been throughout -- fourth to fifth decimal, precision up and recall down,
SV the one class ahead. The single-sweep arm sits within 0.00005 of the five-sweep arm everywhere,
which is the result to expect and therefore the useful check: identical likelihoods and identical
settled genotypes, differing only in emission path, so a real gap would have meant a rendering bug.
chr20 also runs *faster* than the inline arm, 173 s against 205, because the retracted records are
work the inline arm did and then kept.

**Two measurements sized this before it was built.** An unconditional sweep must visit 57,401
reference-crossed nested chains on chr20 against the 30,020 descent considers -- 1.91x the chains but
only +14.0% of chain visits, which is the unit of work. And retention costs 72.4 MB, of which 51.9 MB
is traversals, measured with protobuf's own `SpaceUsedLong` rather than priced by hand.

**What the retention is not.** The first estimate of it was arithmetic -- 48 bytes a visit against
component counts -- and came to 87 MB while the byte total the design was justified on was wrong in
both directions at once: too low per object, and counting only nested chains. The measured figure is
smaller than the guess, and a later attempt to blame protobuf overhead for a 1.4 GB memory rise was
wrong by twentyfold. The rise was a dropped `ploidy == 1` condition making every top-level site carry
an alternate CallInfo it could never use; the retention counter could not see it because it counts
nested chains only, so 72 MB of retention sat beside 1.4 GB of resident memory with no apparent
contradiction. A counter covering part of a population is worse than none.

**And what the design does not fix.** It was argued that choosing the allele list with the genotype
would eliminate unused ALTs. Measured across all three arms they run 4,484 / 4,480 / 4,301 of about
120,600 -- ~3.7% either way. Those are multi-allelic records whose GT names a subset of the ALTs, not a
nested-calling artefact. The `NGT2="."` population that claim was generalised from is the narrow
ploidy-switch slice, about 2% of the total. A ploidy revision can no longer strand an ALT; an ALT can
still go unused.

`INFO/NGT2` is retired with this. It reported what a nested haploid site would call at ploidy 2 for a
caller that could not act on it -- the header said as much -- and the barrier now acts on it. The data
behind it stays and is load-bearing.

**`NestedIncoherence` and the three nested FILTERs stay**, though they now never fire. No nested FILTER
appears anywhere in the whole genome, which is exactly their value: the detector for a coherence
violation is armed and reports none. Retiring them would remove the only thing that would say if the
guarantee ever broke.

#### Instrumentation

143 lines across `graph_caller.{cpp,hpp}` and `linkage_model.{cpp,hpp}`, reported under `--progress` and
otherwise inert. Verified behaviour-neutral: chr20 re-called with it produces 117,047 records whose
content is identical to the shipped run's, and every pre-existing counter reproduces exactly (2,135
nested sites, 88 diploid, 167 unreachable, 10,248 diploid children, 8,799 genotypes changed, 255 records
flagged). The only diff is the relative order of two records sharing a position, which the output sort
has always left to thread arrival.

Kept rather than reverted: the depth histogram says how many barriers a Stage 1 build actually ran, and
the skipped-child counter is the check on whether the 296 became calls.

## Shipped as the default

`vg call` turns nested calling on wherever `--read-likelihood` runs, and phasing on wherever the
linkage layer runs. `--no-nested` and `--no-phased` restore the old behaviour. Both follow the pattern
`--linkage-weight` already used: the default declines where its preconditions are absent, and only an
explicit request errors.

**Whole genome, autosomes, against T2T-Q100, one build:**

| | old default | new default |
|---|---|---|
| SNV F1 | 0.9752 | **0.9833** |
| SNV recall | 0.9567 | **0.9740** |
| ALL F1 | 0.9626 | **0.9699** |
| Indel F1 | 0.9147 | **0.9191** |
| SV F1 | 0.5134 | **0.5467** |

5,041,066 records over 24 contigs in 59.3 minutes, at no runtime or memory cost.

**Tier-2, `readlik` arm, GT comparison, ALL F1 -- and the caveat the genome run cannot show:**

| dataset | before | after |
|---|---|---|
| chr20, 34-hap | 0.9645 | **0.9698** |
| chr6, 34-hap | 0.9689 | **0.9749** |
| chr20, 4-hap | 0.9507 | 0.9502 |
| chr6, 4-hap | 0.9601 | 0.9598 |

**The gain is a rich-panel effect.** On 34 haplotypes recall rises about 1.3 points for about 0.3 of
precision -- chr20 SNV F1 0.9780 to 0.9841. On 4 haplotypes recall barely moves (0.9329 to 0.9335)
while precision slips (0.9691 to 0.9675), so F1 lands marginally *down*. That is what the mechanism
predicts: nested calling recovers variants buried inside long collapsing ALTs, and a four-haplotype
panel enumerates few of them, so there is little to recover while the extra-records cost still applies.
The right trade for HPRC-scale graphs, and not a free one everywhere.

**Did it reach the thing it was built for?** The Stage 0 gate was the 2,219 same-length
substitution false positives and the population of SVs with no record anywhere near them. Measured
against `--no-nested` on the same graph, reads and scoring:

| | `--no-nested` | current default |
|---|---|---|
| same-length substitution SV false positives | 2,219 | **294** (PanGenie: 29) |
| truth SVs missed with no vg record within 100 bp | 1,230 (46.7%) | **258 (13.7%)** |
| truth SVs recovered / lost, per variant | -- | **+1,174 / -119** |

The gate is met. The residual gap against PanGenie is a different population, and the inverse
pathology -- one event written as several sub-threshold records -- exists at 156 loci, up from 98.
Both in [sv-residual-errors.md](sv-residual-errors.md).

**Genome-wide coherence, the first time these counters have run at scale:**

| | |
|---|---|
| sites phased | 5,041,066 |
| nested haploid sites | 54,551 |
| flagged `nested_diploid` | 2,458 |
| flagged `nested_unreachable` | 5,000 |
| children called at ploidy 2 | 237,814 |
| of those, parent moved by linkage | 10,767 |
| flagged `nested_haploid` (2 -> 1) | **0** |
| het sites with an undetermined allele order | 5,931 |

0.15% of records carry a coherence flag. The 2 -> 1 class is empty with 10,767 opportunities, so the
incoherence really is one-sided: 13.7% of haploid nested children against 0% of the diploid ones, a
population four times larger.

**Phasing, now that the nested records are readable.** Including them: 2,512,675 pairs, 61,169
switches, 2.43%, blockwise hamming 48.29%, 22 blocks, N50 248.386 Mb. Against the diploid-only
2,510,608 pairs and 60,319 switches, the nested records contribute **2,067 pairs and 850 switches -- a
41% switch rate against a 2.43% baseline**. Their genotypes are good and the haplotype they are
assigned to is close to a coin flip.

That, and the 5,931 heterozygous sites whose allele order nothing determines, are the two defects this
work is shipping in the default. Both are written up under
[Open, and now shipping in the default](#open-and-now-shipping-in-the-default) with what has already
been ruled out and where to look next.

## Open, and now shipping in the default

Two defects belong to this work rather than to some future project, and both are properties of the
**default** output as of vg 956864c18. Neither is a hypothesis: both are measured, and for the first
the obvious explanations have already been ruled out.

### 1. Nested haploid records are assigned to a haplotype barely better than chance

Genome-wide, including the nested records in the whatshap comparison for the first time:

| autosomes | pairs | switches | switch % |
|---|---|---|---|
| diploid records only | 2,510,608 | 60,319 | 2.40% |
| including nested | 2,512,675 | 61,169 | 2.43% |
| **the nested contribution** | **2,067** | **850** | **41%** |

So a nested site's *genotype* is good -- that is what the whole recall gain rests on -- and the
haplotype it is placed on is close to a coin flip. On chr20 alone the same figure is 21%; the genome
is worse.

**What has been ruled out.** The strand is derived from `parent_slot`, an index into the parent's
called-traversal order. Deriving it instead from which of the parent's *phased* alleles crosses the
child -- the obvious correction, since that is what a phase means -- was implemented and measured: the
two conventions are indistinguishable, 1,655 switches against 1,661 on ~58,900 chr20 pairs. So this is
not an indexing error, and the fix is not a better index. See item 4 of Stage 7 for the two earlier
measurements that appeared to say otherwise and why both were wrong.

**Where to look next**, in order of suspicion:

1. **The child inherits the parent's phase, and the parent's phase is itself 2.4%-switch accurate.**
   Test: condition the nested switch rate on whether the parent's own adjacent pairs switch. If the
   nested rate collapses to the parent's when the parent is locally stable, the child is faithfully
   inheriting a bad frame and the fix belongs upstream in the diploid phasing, not here.
2. **The 5,931 sites of item 2 below.** A nested site hanging off a parent whose allele pair has no
   determined order has no frame to inherit at all. Measure those separately.
3. **The 82% `parent_slot == 1` skew.** Nothing in the design predicts it and nobody has explained it.
   Find out what makes the crossing traversal almost always the second in `trav_genotype` order --
   most likely `ref_trav_idx` sorting late, but that is a guess and should be checked rather than
   assumed.

**Instrument**: `scripts/wgs/nested_strand_check.py`, with the warning in its own docstring that its
site-level percentages cannot compare two conventions -- a constant beats both. Relative phase through
whatshap is the measure that works, and it works only because the strand now reaches the VCF as `a|.`
or `.|a`.

### 2. Heterozygous sites emit a phased genotype whose orientation nothing chose

When neither phased panel haplotype spells either called allele, `LinkageCollector::resolve` falls
through to writing the allele pair in **sorted order**, and that pair is then emitted as a phased `GT`
inside the block, indistinguishable from one the panel actually oriented. `PhaseCall::order_arbitrary`
marks them and the run reports the count: **5,931 genome-wide**, 180 on chr20.

Small against 5,041,066 phased sites, and still wrong in kind rather than in degree: a phased genotype
that asserts an orientation nothing chose is a guess dressed as a call, which this caller avoids
everywhere else. It is also the most likely single contributor to item 1, since a nested child of such
a parent inherits a frame that was never decided.

**Options, none of them free.** Leave the pair unordered and emit the site unphased, which costs those
records their `PS` and fragments nothing but is honest. Or expose the flag in the record -- a `FORMAT`
field or a `FILTER` -- so a consumer can tell the two apart, which keeps the phase set intact. The
second is cheaper and strictly more informative; the first is what the rest of the caller would do.

## Retiring INFO/NGT2, and what the docs rest on

`INFO/NGT2` reported, for a nested site called at ploidy 1, the genotype its own reads would give at
ploidy 2. It existed to expose the incoherence this design removes: under single-sweep calling a
nested chain is genotyped at the ploidy its parent's settled genotype implies, so the alternative
ploidy is no longer a discrepancy to report but an input the caller consumes and discards. It was on
73,262 records genome-wide, all of them saying the same thing. Removed; `alt_ploidy_best` stays,
because the barrier reads it.

The whole-genome arm the results pages quote (`work/wgs-single`) was called *before* that removal, so
its header still declares the tag. Removing an INFO field cannot change a genotype, but that is the
kind of claim worth checking rather than asserting, so chr20 was re-called with the final binary and
compared: **116,945 records both ways, and the two bodies are identical as multisets once the tag is
stripped.** The only difference is the order of records that share a position -- the emission
buffer's sort does not tie-break beyond the position, and the barrier now inserts nested records in
a different sequence. So the arm on disk faithfully represents the code being pushed.

Every FILTER count on those pages is now `PASS`: **all 5,037,820 records genome-wide carry no
FILTER at all**, coherence or otherwise. The three nested coherence FILTERs remain in the header as
a live invariant check, and firing zero times is what they are for.

The results pages were refreshed against this arm, and one class of staleness was fixed at the
source rather than in the text. `bench_wgs.py` rewrites `docs/wgs-results.md` wholesale, so its
prose lives in the generator -- but the *numbers* inside that prose were typed in, and had gone
stale silently: the tables moved with each rescore while the sentences above them kept quoting the
previous arm, and the diff looked clean because only the tables changed. Those figures are now
computed (`scripts/wgs/bench_wgs.py`), and the quality-gate and false-positive-population tables
that had been hand-measured are now generated too (`scripts/wgs/sv_quality_gates.py`, whose ungated
row is asserted against truvari's own F1 so a mis-accounting cannot pass silently).

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
