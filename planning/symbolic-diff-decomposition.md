# Symbolic diff decomposition: one snarl, several variants

Today a snarl emits at most one record. The symbolic form feeds a single equality test
(`is_symbolically_reference`, graph_caller.cpp:1880), yielding one bit per haplotype -- same route
as the reference (allele 0) or not (one ALT spanning the whole snarl, concrete sequence throughout).
This plan replaces that bit with an alignment, and emits one record per difference block.

## The algorithm

**Step 1 -- symbolic projection** (implemented, `symbolic_allele.cpp`). Over the alphabet

    Sigma_S = {oriented nodes of S interior to no child chain} U {[[C_j, +/-]]}

project each traversal T to sigma(T) by replacing each excursion through child chain C_j with one
symbol carrying the chain's own boundary nodes and the crossing direction, never the route inside.

**Step 2 -- alignment.** Find a monotone partial matching M of sigma(R) against sigma(H_k) with
r_i == h_j for every (i,j) in M, minimising edit cost. The complement of M decomposes into maximal
gaps, each of which is one REPLACE block b = ([i,i'), [j,j')), either range possibly empty.

**Step 3 -- haploid emission.** One record per REPLACE block: REF from R's steps [i,i'), ALT from
H_k's steps [j,j'), chain symbols inside a block expanded to that haplotype's concrete route, POS
from the reference offset of step i. MATCH blocks emit nothing.

**Step 4 -- diploid join.** Cluster the two block lists by transitive overlap of reference step
range. Per cluster: one haplotype only -> 0/1; both, same range, same ALT -> 1/1; both, same range,
different ALT -> 1/2; overlapping unequal ranges -> one record over the union, each ALT re-expressed
across the union.

**The exactly-once rule.** Steps 2-4 are not self-contained: they change what descent must do one
level down. A chain symbol can appear inside a REPLACE block, where it is expanded to concrete
sequence. The rule:

    ploidy(C_j) = #{ k : sigma(H_k) places C_j as a MATCH step }

routed through the existing `child_ploidy` path (graph_caller.cpp:5078) rather than new plumbing.

**Scope of that rule, corrected.** It was first justified as preventing the chain from being
reported twice -- once expanded in a block ALT and once as its own record. That justification is
wrong for the dominant case. Symbolic descent passes `parent_child_trav_sets == nullptr`
(graph_caller.cpp:4967-4968, and the comment states it: "Children are genotyped independently, with
no parent traversal sets"), so in the child's own `call_snarl_internal`:

    if (common_names.empty()) {
        if (parent_child_trav_sets == nullptr || parent_ref_path_name.empty()) {
            return false;                          // graph_caller.cpp:4515-4522
        }
    }

`common_names` is the intersection of path names visiting the chain's two boundary nodes, so a chain
the reference does not cross has it empty and the child returns false before genotyping. It emits
nothing. The `use_parent_interval` branch at :4551 is dead on this path for the same reason. This is
stage 16's blocker, and the reason a covering reference is wanted.

**And there is a more direct gate than that one, with the population already counted.** Descent never
even reaches such a child: graph_caller.cpp:5674-5680 asks `child_ploidy` over the reference
traversal alone and `continue`s when it is zero, incrementing `g_descent_skipped_no_ref`. Its own
comment says so -- "v1 descends only where the reference also goes. A chain crossed only by a
non-reference allele has no reference path through it, so REF and POS for its record are undefined".
On chr20-34hap that counter reads **12,516**, against 30,416 child calls actually made. So roughly
a third of the child chains the called alleles reach are never given a record, and their sequence
reaches the VCF only inside their parent's ALT.

That is the off-reference burial population, measured, and it did not need a new counter.

Enumerating where a chain symbol lands in a REPLACE block:

| case | chain crossed by | today | duplicate? |
|---|---|---|---|
| (a) | H only | no reference path -> `return false` | no |
| (b) | R only | called at the ploidy its crossings give; block reports the deletion, chain reports its interior | no, complementary |
| (c) | both, same chain, unmatched by the alignment | ploidy 2, chain emits, and both block REF and ALT contain routes through it | **yes** |
| (d) | R and H cross different chains in one block | reduces to (a) for one, (b) for the other | no |

So the rule is a **loop-case guard** for (c), not a general anti-duplication measure. It still gives
the right answer everywhere -- (b) resolves because the rule is per-haplotype, so a chain crossed
identically by H1 and deleted by H2 is MATCH on H1 and gets ploidy 1 -- but stage 6's gate must be
written against (c) alone, whose size is expected to track the `g_child_multi_crossing` population
(0 on chr20, 242 on chrX).

**And the consequence that runs the other way, which is the larger one.** Today an off-reference
chain's internal sequence reaches the VCF only inside the parent's whole-snarl ALT; it has no record
of its own and cannot get one without a reference to call against. Decomposition does not make it
callable, but it does put that sequence in a much smaller allele -- the same "buried in a large
allele" mechanism symbolic collapsing addressed (55,222 of 142,707 autosomal SNV false negatives),
reached from the other side. This gain is **independent** of the split-into-many-records effect
measured at 1,000 records: it applies to every REPLACE block containing an off-reference chain,
whether or not the snarl splits. The population is unmeasured and belongs in stage 2. If it is large
it is a better argument for this change than the block count is.

Caveat on scope: the above is the shipped `--read-likelihood` symbolic-descent path.
`NestedFlowCaller` (`-A`) has its own recursion near :5350 which does thread traversal sets, so the
`return false` may not apply there; untraced, and immaterial because the flag refuses that path.

Separate caveat, found in review: MATCH-ness of a chain symbol does NOT imply the two traversals
resume on the same node. `symbolic_allele` picks the exit as the earliest later occurrence of
*either* boundary and sets `backward` from the *entry* node (symbolic_allele.cpp:82-104), so two
traversals can carry an identical `SymbolicStep` and leave at different boundaries. The ploidy rule
survives this; the tempting corollary that the step after a matched chain symbol is itself matched
does not. Do not assert it.

## What review changed about the picture

### F1 cannot adjudicate this change

Aardvark's precision denominator is the query **record count**. Verified on the shipped arm
(`work/tier2-chr20-hap32/results/aardvark-readlik/summary.tsv`):

    recall    = truth_tp/truth_total = 91470/94691 = 0.9659840956  == metric_recall
    precision = query_tp/query_total = 91765/93772 = 0.9785970226  == metric_precision

exact to ten decimal places. `truth_total` is 94,691 across all four 34-hap arms; `query_total`
swings 93,772 -> 96,594. So recall is fair and precision is not, because decomposition moves the
denominator mechanically. The differences this harness resolves are 0.002-0.01 F1
(`docs/tier2-chr20-results.md:21` records this arm moving 0.9699 -> 0.9722, i.e. 0.0023). A 1-2%
denominator move is the same order as the entire signal.

**Consequence:** GT F1 is reported with `query_total` printed beside it, and never gated on.
The fair gate is truth-side counters at fixed denominators, plus BASEPAIR query bases falling.
BASEPAIR has `truth_tp == query_tp` exactly in every arm checked (378,392/378,392 on `readlik`),
a single symmetric base count rather than two per-record tallies. Its diagnostic baseline is the
thesis: the caller claims 427,914 bases where truth has 390,636, and on insertions 135,480 against
95,512 true. That excess is the whole-snarl-ALT signature the change exists to remove.

Truvari alone measures representation, not accuracy: `score_vcf.py:135` passes `--sizemin 50
--sizefilt 50` with no `--refdist/--pctseq/--pctsize`, and `docs/tier2-sv-errors.md:715-724`
measures `truvari refine` at +0.081 to +0.096 F1 on the *same calls*. Gate the refined number only.

### The blocker is `update_vcf_info`, not the emitter

Two facts, both verified:

- Both support callers open with `assert(traversals.size() == variant.alleles.size())`
  (snarl_caller.cpp:302 Ratio, :790 Poisson). A hard abort, so the block path must be gated to
  `ReadLikelihoodCallInfo` or every non-`--read-likelihood` run dies on the first split record.
- Inside `ReadLikelihoodSnarlCaller::update_vcf_info`, `site_to_scored` is built by whole-traversal
  structural equality (`traversals_equal`, read_likelihood_caller.cpp:391-398). A block sub-range
  matches nothing, so `site_to_scored[s] == -1` for every block allele: **AD comes out all zeros**
  (:409-421), `all_mapped` goes false and **GL is dropped** (:427-473), QUAL falls to 0.

This is the silent path and the most likely way the change ships broken. It is also the largest new
work, the least testable, and the thing that decides whether the split can ship at all.

### A pre-existing bug this plan uncovered: flipped snarls have no symbolic projection

`flip_snarl` (graph_caller.cpp:210-214) is applied whenever the reference path runs backwards
through a snarl (:4561-4565, and again at :5350-5354), and `snarl` is held by value so the flipped
copy is what `emit_variant` and `symbolic_allele` receive (graph_caller.hpp:1045 says so).

`symbolic_allele` then computes `site_ptr = into_which_snarl(site.start().node_id(),
site.start().backward())`. For a flipped snarl that lookup *succeeds* -- the boundary index maps
`(end.node_id, !end.backward)` to the canonical snarl (snarls.cpp:926-934) -- but the guard at
symbolic_allele.cpp:41-45 then compares node ids: canonical start (= original start) against
`site.start()` (= original *end*). They differ, so `site_ptr = nullptr`, `is_child` is false at
every visit (:74-76), and sigma(T) degenerates to the plain node list with no chain symbols at all.

So symbolic collapsing -- and with it the whole nested-calling benefit -- is **silently inactive for
every snarl whose reference path runs backwards.** Consequences:

- The measured benefit of symbolic collapsing was measured over an unknown subset.
- For these snarls block decomposition would shred at node granularity, so the 372-block worst case
  concentrates here.
- Every child chain is a run of matched *plain nodes*, so the ploidy rule above reads ploidy 0 for
  every child of a flipped parent, turning double-reporting into **non**-reporting.

The population is unmeasured. 5,181 of 115,038 chr20 records (4.5%) have snarl-ID boundary nodes
running counter to node-id order, which is where flips live -- suggestive, not a measurement, since
canonical snarl orientation is not determined by node id order. Stage 2 measures it properly.

## Three decisions before any emitter code

### D1 -- alignment cost model and tie-break

The spec says "minimum number of edits, allowing insertions, deletions and matches", i.e. the
LCS/diff model with no substitution. **That is underdetermined, and the tie-break decides how many
variants come out.** Worked example, sigma(R) = [a,b] against sigma(H) = [b,b]:

- match r_1 with h_0: delete a, MATCH, insert b -> **two** REPLACE blocks separated by a match
- match r_1 with h_1: one REPLACE block a->b -> **one** record

Both alignments have LCS length 1 and cost 2. The spec cannot choose between them.

**Recommendation: substitution at cost 1.** Then a->b costs 1 against del+ins at 2, so the
one-block reading wins strictly and the degeneracy disappears. This is a disambiguator, not a
correction to the cost model -- it encodes "prefer fewer, larger blocks", which is the same
preference the aggregation rule already expresses. The tie-break must still be pinned (recommend
lexicographically-least under MATCH < SUB < DEL < INS) because ties remain, and an unstable
tie-break makes output irreproducible under `-t`.

This contract is an input to five other areas -- the block key, the child's POS and suppression,
AD and GQI, the offline pre-flight, and output determinism. One written contract, one unit test on
a deliberately tied pair, before anything downstream.

### D2 -- per-block AD and GL

The blocker above. Options, in increasing cost:

| option | what it means | why it may not be enough |
|---|---|---|
| omit AD/GL when N>1 | smallest change | three TAP gates treat *absence* as failure: 18_vg_call.t:196 `GL_BAD`, :202 `GT_MISMATCH`, :206 `GQ_ABOVE_GQI` all print (= count bad) on a missing field |
| replicate the site AD/GL onto every block | arity-correct, cheap | N records each claim the same evidence, so any consumer summing or averaging double-counts; F1 cannot see it |
| max-marginal fold over `genotype_lls` per block | the only principled answer | ~200 lines of new numerics whose only gate is a counter; redefines GQN and QUAL; depends on the D1 traceback |

**Recommendation: replicate for the experiment, behind the flag, with an INFO field naming the
snarl and block index so the replicated set is recoverable -- and do not ship it.** The justification
is narrow and must be stated as such: F1's blindness to quality fields makes replication safe for
*measuring whether decomposition helps*, and unsafe for production precisely because that blindness
means nothing would ever catch it. If the measurement says decomposition helps, the fold becomes
required work, not optional.

### D3 -- record identity

**Recommendation: keep the ID unsuffixed and add a block field to `BufferedRecordKey`.** Suffixing
(`>1>6.2`) needs a normaliser at four sites, and two of them abort rather than degrade:
`name_to_snarl.at(name)` (:2797) and `top_level_ref_info.at(...)` (:2877) both throw on an
unrecognised name inside `#pragma omp parallel for` (:2836), so a throw terminates. Worse,
`is_top_level` (:2740-2752) tests `chrom_of_name.count(print_snarl(*snarl))` -- an *unsuffixed*
parent name -- so every nested record would read as top-level.

The sort key needs the field regardless: `buffered_record_key_less` (:626-635) tie-breaks on
`a.id < b.id` and stops, and two blocks *can* share a POS (a deletion on one haplotype abutting an
insertion on the other), so without a block field the comparator is non-antisymmetric and
`std::sort` becomes input-order dependent -- the exact defect graph_caller.hpp:145-155 records as
72 differing record pairs between two chr20 runs.

## Stages

Each stage has a gate. A stage that misses its gate stops the sequence.

### Stage 0 -- the free calibration, no build

**Does aardvark credit a decomposed record as TP, or charge it as FP?** Everything downstream
depends on the answer and it is already computed:

    work/tier2-chr20/results/sweep-chr20-4hap-atomize0.vcf.gz   104,165 records,  0 duplicate IDs
    work/tier2-chr20/results/sweep-chr20-4hap-atomize8.vcf.gz   104,234 records, 39 duplicate IDs

Same binary, one minute apart (Aug 10 15:52 and 15:53), differing *only* by record decomposition:
104,234 - 104,165 = 69, and 39 IDs carry those 69 extra records. atomize0 was never scored -- there
is no `aardvark-*atomize0*` directory. Two `score_vcf.py` runs settle it.

**RESULT: PROCEED.** Scored both under fresh labels `s0-atomize0` / `s0-atomize8`, same session,
same scorer invocation. Decomposition alone -- 69 records across 39 IDs, nothing else different:

| | atomize0 | atomize8 | delta |
|---|---|---|---|
| GT ALL F1 | 0.948779 | 0.948778 | **-0.000001** |
| GT ALL recall | 0.932232 | 0.932274 | +0.000042 |
| GT ALL precision | 0.965923 | 0.965875 | -0.000048 |
| GT ALL query_total | 91,528 | 91,547 | +19 |
| GT ALL truth_tp | 88,274 | 88,278 | +4 |
| BASEPAIR ALL F1 | 0.921069 | 0.921069 | **+0.000000** |
| GT JointIndel F1 | 0.853129 | 0.853247 | +0.000118 |
| GT Insertion / Deletion F1 | -- | -- | **exactly 0** |

Three things this settles, and one it sharpens:

1. **Aardvark credits a decomposed record.** It is not structurally hostile to decomposition. GT F1
   moves by one part in a million and BASEPAIR F1 not at all, so representation change alone does
   not manufacture or destroy score.
2. **BASEPAIR is confirmed as the neutral ruler.** Flat to six decimals, with `truth_tp == query_tp`
   as always.
3. **Insertion and Deletion are exactly unchanged**, confirming `--atomize-substitutions` only ever
   touched same-length records -- so this calibration covers the substitution population only, and
   the length-changing blocks the new algorithm reaches are uncalibrated.

**The sharpened finding, which replaces "F1 is unfair" with something more useful.** The marginal
decomposed record is credited at **14/19 = 73.7%**, against the arm's 96.6% average precision. So
the denominator does move, but not for free -- added records must earn their TP, and they earn it at
a below-average rate. The mechanical exposure is therefore
`N x (0.966 - 0.737) / query_total`: for N = 2,000 added records that is about **-0.005 GT F1**,
squarely inside the 0.002-0.01 band this harness resolves.

**Revised gate for stage 7:** GT F1 is usable, not disqualified, but it must be reported with
`query_total` and the marginal credit rate beside it, and BASEPAIR remains primary. A GT F1 that
falls by less than the mechanical exposure is a *pass*, not a regression -- and that is a judgement
the numbers must be shown for, not folded into a single column.

This also disposes of the prior art as an effect-size prior. `docs/tier2-sv-errors.md:558-577`
reports `--atomize-substitutions` at SV F1 +0.0017 to +0.008 and rejects it, but its baseline was
overwritten in place: the chr20-4hap record count moved 104,165 -> 109,476 between the atomize8
scoring (Aug 10) and today's baseline (Aug 19). **5,311 records of revision drift against 69 records
of treatment -- a 77x confound.** The numbers in that table cannot be compared to anything current.
Process rule, absolute: re-run the baseline arm from the same binary in the same session under a
fresh label.

Note also why that prior art does not settle the question even if it were clean: it was restricted
to same-length biallelic records specifically to avoid needing an aligner, so it could touch no
length-changing block; it worked in base space with no notion of a chain, so it could delegate
nothing to a nested record; and it never touched the diploid join, where the 1/2 -> 1/1 collapse is
1,334 chr20 records -- larger than the 1,000-record multi-block population.

### Stage 1 -- three standalone zero-diff fixes

Each lands alone, each is independently useful, none changes output except the bug it fixes.

1. **`flatten_common_allele_ends` underflow** (:2364-2376). `max_flatten_len` is `size_t`;
   `min_allele_len == 0` makes `max_flatten_len == min_allele_len`, so `--max_flatten_len`
   underflows to `SIZE_MAX` and the backward pass reads `""[0-1-0]`. Out-of-bounds, not a wrong
   answer. Unreachable today because every allele carries both boundary nodes; reachable the first
   time a pure-insertion or pure-deletion block is rendered. `if (variant.alt.size() == 0) return;`
   does not save it -- a pure-deletion block has one (empty) ALT. Two lines.
2. **AD arity in the missing-allele fixup** (:2186-2193). The `*` fixup pushes to `alt`, `alleles`
   and `info["AT"]` but runs *after* `update_vcf_info` (:2142) wrote Number=R AD, so the record has
   two alleles and one AD entry. No test covers it. The AD-arity gate is the main tripwire for
   everything downstream and cannot be green until this is fixed.
3. **`symbolic_allele` step -> visit-range out-parameter**, with a unit test. The ranges partition
   `[0, visit_size)` contiguously (`i` advances by 1 or to `exit`, one symbol per advance,
   symbolic_allele.cpp:88-118) and a chain's range is `[i, exit)` with the exit boundary belonging
   to the *next* step. This is what makes per-block sequence extraction possible at all.

**Gate:** byte-identical output on the chr20 arm except for records exhibiting bug 2.

**RESULT: PASS.** chr20-34hap `readlik`, 115,038 records, **byte-identical** to
`work/tier2-chr20-hap32/results/readlik.vcf.gz` (`cmp` on the record bodies). Bug 2 changes nothing
on this data because the arm has **zero** records carrying `*` in any ALT -- the missing-allele path
is simply not exercised here, which is why the arity defect survived. GL needed no matching fixup:
it is already omitted whenever the genotype carries a marker, so the two were never inconsistent.

### Stage 2 -- counters only, no output change

Four counters, because four numbers currently rest on a Python proxy rather than the implementation:

- blocks-per-record histogram (reproduces the 1,000 / 372 figures from inside the caller)
- the flipped-snarl `site_ptr == nullptr` population -- **the number that says whether the existing
  nested-calling benefit was measured over 95% or 60% of snarls**
- chain-multiplicity disagreement between sigma(R) and sigma(H_k)
- blocks falling under 50 bp on a >=50 bp parent, which sizes the truvari size-filter exposure
- case (c): chain symbol present in both sigma(R) and sigma(H_k) but unmatched by the alignment --
  the only genuine double-reporting population, and the gate for stage 6
- called ALTs containing a chain the reference does not cross, and the share of allele length that
  chain accounts for. The *count* of such chains is already known -- `g_descent_skipped_no_ref` is
  12,516 on chr20-34hap -- so what this adds is the base share, i.e. how much of an allele's length
  is chain the reference never visits. That is what says whether shrinking the allele around it
  matters.

**Gate:** blocks-per-record from the caller agrees with the offline proxy to within the flipped-snarl
population. If it does not, the proxy was measuring something else and the motivating numbers are
wrong.

**RESULT: the gate FAILS, and it was right to exist.** chr20-34hap, flag off, output byte-identical:

    atomize: 115996 called ALTs projected, 115777 difference blocks,
             638 ALTs a diff would split (>=2 blocks)
    blocks/ALT: 0=1204 1=114154 2=488 3=96 4=20 5=16 6=2 7=2 8=1 9=2 10=3 11=1 12=2 14=2 15+=3
    9279 sites where projection is inert because the snarl does not resolve (flip_snarl)
    47 ALTs carrying a chain the reference also crosses but the alignment did not match
    4428 ALTs carry a chain the reference does not cross, 1456000 of 2545504 bases (57.2%)

**The offline proxy overstated the split population by 57%: 638, not 1,000.** The proxy diffed
INFO/AT node paths, which have no chain symbols, so a difference inside a child chain counted as its
own block. The caller's symbolic projection collapses exactly those, which is the whole point of
symbolic alleles. Every figure in this plan derived from the proxy has to be read at the lower
number, and the motivating "0.87% of records split" is really **0.55%**.

**The flipped-snarl population is 9,279**, against 115,996 ALTs projected -- about 7.4% of sites
have symbolic collapsing silently off. Real, moderate, and previously unknown.

**Case (c) is 47 ALTs.** The double-reporting population the exactly-once rule exists for is
0.04% of called ALTs. The correction to the rule's scope was right, and stage 6 is very nearly a
no-op guard rather than a fix.

**The off-reference finding is the strongest number here.** 4,428 ALTs carry a chain the reference
never visits, and those chains are **57.2% of those alleles' bases**. More than half the sequence in
such an allele has no record of its own and cannot get one without a reference to call against. But
see the arithmetic below before treating it as a win.

**One counter of mine was mis-specified and is fixed rather than reported.** It compared the whole
traversal's length against the smallest block's, and a single block spans everything except the two
boundary steps -- so it was measuring boundary-node length and called 74,885 of 115,996 ALTs
"exposed" to truvari's 50 bp filter. Only a genuine split can move a variant across that filter, so
the counter now requires `blocks >= 2` and the exposure cannot exceed 638.

### The arithmetic this forces, stated before the stage 7 run

638 ALTs split, into 1,623 blocks where they currently produce 638 records: about **+985 records**.
Stage 0 measured the marginal decomposed record being credited at 73.7% against a 96.6% average, so
the mechanical GT F1 cost of adding them is about

    985 x (0.966 - 0.737) / 91528 ~= -0.0025

**That is the same order as the entire signal this harness resolves**, and it is a cost the change
must overcome before any gain shows. Two consequences:

- A single-block ALT spans steps [1, m-1), i.e. everything but the two boundary nodes, so
  decomposition barely changes it -- `flatten_common_allele_ends` already trims those ends. The
  99.45% of ALTs that do not split get essentially nothing from this change. In particular the
  4,428 off-reference ALTs are **not** helped unless they also split, because their allele does not
  get smaller.
- So the case rests on 638 ALTs, and it has to beat -0.0025 on GT F1. BASEPAIR, which has no record
  denominator, is where a real improvement would show first and is the metric to believe.

### Stage 3 -- the alignment as a pure function

`vector<Block> symbolic_diff(const SymbolicAllele&, const SymbolicAllele&)` with the D1 contract,
and unit tests on: identical inputs (no blocks), pure insertion, pure deletion, the deliberately
tied pair from D1, a repeated chain symbol (loops), and a block whose ALT string equals its REF
(two distinct nodes spelling the same bases -- MATCH is step identity,
symbolic_allele.hpp:56-58, so this is reachable and today's string dedup at :1938-1982 is what
catches it).

**Gate:** unit tests green; the function is deterministic under repeated invocation and independent
of thread count.

**RESULT: PASS.** Full unit suite green at **12,547,745 assertions in 853 test cases**.
`symbolic_diff` landed in `src/symbolic_allele.{hpp,cpp}` with the substitution
cost model and the documented tie-break (prefer the diagonal, then deletion, then insertion, walked
backwards). 41 assertions in 10 cases under `[symbolic_diff]`, plus 77 in 12 under
`[symbolic_allele]` including the new visit-range partition tests. The load-bearing case is
`plain({1,2})` against `plain({2,2})`: it must come out as **one** block, and would come out as two
under an insert/delete-only model. Tested in both mirror orientations so the answer cannot come from
a left/right asymmetry.

Two things landed alongside, both needed by the emitter rather than by the diff:
- `symbolic_allele`'s visit-range out-parameter, whose ranges partition `[0, visit_size)`
  contiguously -- verified by a helper the tests run on every fixture, since a dropped visit would
  silently shorten an allele by a node and read as a real variant.
- `out_alt_before_ref`, giving the alt step index on arrival at each reference step. This is what
  lets two haplotypes express their alleles over one shared reference span, which the diploid join
  needs and which no amount of per-haplotype block data provides.

### Stage 4 -- per-block quality, per D2

**Gate:** AD arity equals allele count on every record; GL present on every record; the three TAP
quality gates green. A counter for records whose AD is all-zero must read zero -- that is the
`traversals_equal` failure mode and it is otherwise silent.

### Stage 5 -- the emitter, behind an off-by-default flag

The flag refuses, with an explicit error rather than silently degrading: `-a/--genotype-snarls`,
`--legacy`, `--bottom-up`, and a non-null `trav_to_string`.

`-a` matters more than it looks. `nested_calling` (hence `symbolic_manager`) is ON by default under
`--read-likelihood` and is NOT cleared by `-a` (call_main.cpp:1733-1745), so gating on
`symbolic_manager != nullptr` alone would atomize `vg call -a --read-likelihood` and break the
sample-independent record set that `wants_line` (:2215) contracts for. The harness depends on it
directly: `scripts/wgs/stage0_symbolic.py:103` reads `%ID` from a `-a -A` run as a snarl inventory
and `scripts/wgs/nocall_anatomy.py:59` asks a `-a` VCF whether a snarl exists at a position.

Empty block list falls through to the legacy whole-snarl path, which preserves the all-MATCH case
and the `*` missing-allele record. Note the residual: a snarl carrying `MISSING_ALLELE_MARKER` on
one haplotype and a genuine REPLACE on the other has a non-empty block list, so the fallback does
not fire and the `alt.empty()` fixup is unreachable on the block path. Handle it explicitly.

**Gate:** with the flag off, output is byte-identical to the previous binary. With it on, zero
records have AD arity mismatch, zero have all-zero AD, and `(CHROM, POS, ID, block)` is unique.

**Flag-off gate: PASS** -- byte-identical, 115,038 records (stage 1 result above).
**Existing tests: PASS** -- full unit suite 12,547,745 assertions / 853 cases, and
`test/t/18_vg_call.t` **304/304**. Note the TAP suite must be run as `cd test && prove t/18_vg_call.t`;
running it from the repository root gives `is: command not found` and a bogus FAIL, because
`bash-tap-bootstrap` is sourced by a relative path.

**Refusals.** `-a/--genotype-snarls` refuses with the intended message. The `--legacy`,
`--bottom-up` and `--top-down` refusals were originally placed after caller construction, where
`-k`'s own "pack file is required" check fires first and the refusal is never reached -- so they
were untestable without a pack and, worse, would have made a user wait for a 22 GB graph load to be
told the combination is invalid. Moved to option-validation time alongside the other
option-compatibility checks. The state-dependent refusals (nested calling declined, caller does not
emit VCF) stay late, because they cannot be answered earlier.

Note what byte-identity can and cannot cover: `add_allele_path_to_info` (:1815-1856) writes the
*whole* traversal into INFO/AT. Because step 0 and the last step are always MATCH (the boundary
visits cannot be symbolised), a "whole-snarl" REPLACE block spans steps [1, n-1) and its AT is a
strict sub-walk of today's. So "records that did not split are byte-identical" is achievable only
if AT stays whole-traversal on the block path; if AT becomes the block's range, that gate is
unachievable and the harness's AT consumers silently see truncated routes with no arity change.
**Decide AT before stage 5, not during it.**

### Stage 6 -- descent delegation, per the exactly-once rule

Implements `ploidy(C_j)` as the count of MATCH placements, guarded for the flipped-snarl case from
stage 2 (where every child is a run of matched plain nodes, so the naive rule reads ploidy 0 for
every child and drops the whole subtree).

**Gate:** total record count across the hierarchy does not fall; the counter for "chain suppressed
because it appears in no MATCH block" matches the case-(c) population measured in stage 2. Expect a
small number, possibly zero on chr20 -- a large one means the alignment is failing to match
identical symbols and D1's tie-break is wrong.

### Stage 7 -- validation

Arm `readlik`, dataset `chr20-34hap` (`work/tier2-chr20-hap32`), baseline
`work/tier2-chr20-hap32/results/readlik.vcf.gz` -- verified 115,038 records and current (Aug 23).
The 4-hap chr20 arm is 109,476 records; do not confuse them. Run via
`run_arms.py --readlik-extra "<flag>"` (`scripts/tier2/run_arms.py:195`, which exists for exactly
this), then `score_vcf.py --label <fresh-label>`. Not `param_sweep.py`, whose `--param` is
restricted to `{mismap-max, mismap-min}` (:177).

**Gate, in priority order:**

**RESULT (first pass; superseded by the final run below, which adds phase inheritance and the
stage 6 gate).** Both arms scored from the same binary in the same session, fresh labels
`s7-base` / `s7-atomize`, chr20-34hap:

| | baseline | atomized | delta |
|---|---|---|---|
| records | 115,038 | 115,686 | +648 |
| GT ALL recall | 0.965984 | 0.966143 | **+0.000158** |
| GT ALL precision | 0.978597 | 0.978516 | -0.000081 |
| GT ALL F1 | 0.972250 | 0.972290 | +0.000040 |
| GT ALL query_total | 93,772 | 93,929 | +157 |
| GT JointIndel F1 | 0.927672 | 0.928032 | **+0.000360** |
| BASEPAIR ALL F1 | 0.924542 | 0.924802 | **+0.000260** |
| BASEPAIR ALL query bases | 427,914 | 427,744 | **-170** |
| BASEPAIR ALL truth_tp | 378,392 | 378,441 | **+49** |

Every truth-side counter holds or rises; recall rises in every class except Deletion, which is flat.
BASEPAIR query bases fall while truth_tp rises -- the predicted direction, and the thesis of the
change stated as directly as this harness can state it.

**The predicted mechanical penalty did not materialise at the predicted size.** 648 records were
added but `query_total` moved only +157, because most block records fall outside the confident
region or below a size filter. So the -0.0025 exposure computed from stage 0 was an overestimate by
a factor of four, and GT F1 came out marginally positive instead.

**Invariants:** AD arity mismatch **0**, GL absent **0**, `(CHROM,POS,ID)` collisions **0**, SB
present on exactly the 1,137 block records. The output validates under `bcftools view`.

**The coordinate arithmetic is verified against the reference FASTA, not just internally.**
`bcftools norm --check-ref w -f chr20.fa` reports **mismatch_removed = 0** on both arms:

    baseline   total/split/joined/realigned/mismatch_removed: 115038/0/0/3626/0
    atomized   total/split/joined/realigned/mismatch_removed: 115686/0/0/3863/0

So every block record's REF string is the actual reference sequence at the POS the block was placed
at. That is the check that matters for stage 5's hardest part -- per-step offsets, the anchor-base
rule, and the shared-span extension in the diploid join are all wrong in ways that would show here
and nowhere else. The +237 realigned records are bcftools left-aligning the additional indels the
split produces, which is expected rather than a finding. All-zero AD went 97 -> 140, and the gate as written
("zero") was wrong because the baseline already had 97: the split re-expresses the same population,
72 of the 140 now sitting on block records while ordinary ones fell 97 -> 68. Not new breakage.

**Two defects found by looking at the output rather than the counters**, both fixed and re-run:

1. Block records carried `PS` beside an **unphased** `GT` where the site record was phased. Both
   self-contradictory and a scoring confound, since the baseline's records are phased and these
   would be compared on unequal terms. Block records now inherit the site's phase by mapping each
   genotype slot to the site allele its haplotype carries, and drop `PS` when they genuinely cannot.
2. My truvari size-exposure counter was measuring boundary-node length (see stage 2).

**Interpretation, stated plainly.** The direction is consistently positive on the metrics stage 0
established as fair, and every safety invariant holds. But the magnitude -- +0.00026 BASEPAIR F1,
+0.00016 GT recall -- is roughly **an order of magnitude below the 0.002-0.01 band this harness
resolves**. On chr20-34hap this change is real, safe, and below the noise floor of its own
evaluation. That follows from the population: 489 sites out of 115,038, or 0.43%.

| must hold | why |
|---|---|
| aardvark recall (fixed truth-side denominator) holds or rises | the only fair per-record metric |
| BASEPAIR `query_total` falls toward `truth_total`; `truth_tp` holds | the excess-claimed-bases thesis, stated directly |
| truvari **refined** SV F1 holds or rises | unrefined truvari measures representation, worth +0.08 F1 on identical calls |
| AD-arity, all-zero-AD, `(CHROM,POS,ID,block)`-collision counters all zero | the silent failures F1 cannot see |
| GT F1 reported with `query_total` beside it | visible, never gated |

## Left on the table, deliberately

- **Per-block GL by max-marginal fold.** Required to ship, not required to measure. See D2.
- **The flipped-snarl fix itself.** Stage 2 measures it; fixing it is a separate change with its own
  gate, and it should not ride along with decomposition because it moves the baseline.
- **`allele_core_length` re-stratification.** `merge_similar_alleles` gates on
  `allele_core_length(alleles) < allele_merge_min_len` (:1512-1514), and a block's core length is
  the block's span, so `-L` stops firing where it fires today. Every size-stratified figure in the
  harness re-bins for the same reason. Unmeasured; MAT-entry count is the only observable and the
  report does not collect it.
- **The mosaic/record-count invariant.** 18_vg_call.t:838-841 checks
  `sum(mosaic site counts) == record count`. One Entry per snarl but N records, so it becomes false
  in general -- and it will stay GREEN, because nest_hap.vcf has exactly 2 records and neither can
  split. Nothing forces the restatement.
- **`set_allele_map` last-writer-wins across blocks.** Its allele half writes `allele_arena`, read
  only into `PhaseCall::allele_first/allele_second`, which linkage_model.cpp:2708-2713 records as
  "written and never read" (grep confirms no reader). So every proposed treatment is behaviourally
  identical in the VCF and no test can distinguish them. Whoever revives `allele_*` inherits an
  arbitrary block's numbering.
