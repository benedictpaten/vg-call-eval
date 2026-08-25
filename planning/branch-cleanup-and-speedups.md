# What is left on the read-likelihood branch: a review for redundancy, simplification and speed

A second pass over the branch (`fork/read-likelihood-genotyping`, 23 source files, +16,606
lines against `origin/master`) after the two fixes in `docs/performance.md` landed. Everything
below is measured or read out of the code, not inferred from shape.

Method: a whole-run `sample` profile of chr20 at 20 ms across every thread, self-time attributed
by walking the call tree; a clone detector over branch-owned lines only, with every candidate
checked against the merge-base so upstream duplication is not reported as ours; a
declared-but-never-called sweep over the branch headers; and the 24 whole-genome logs for
counters that are constant.

## Where the time goes now

chr20, `-t 5`, after the key index and the clamped-chain guard. 197.2 s wall, 319 thread-seconds
sampled:

| | share of thread time |
|---|---|
| blocked on a `gbz-base` child | 33.5% |
| parked in OpenMP (a phase on one thread) | 30.9% |
| allocator | 8.3% |
| `panel_alleles` GBWT walk | 5.8% |
| GAF parse and protobuf construct/destroy | 5.3% |
| linkage decode | 2.8% |
| **read scoring -- the likelihood model itself** | **2.1%** |
| window deliver | 0.4% |
| symbolic projection and atomize | 0.1% |
| other | 10.7% |

Two things stand out. **64% of the thread budget is blocked**, half of it on a child process and
half on phases that run on one thread. And of the third that is working, the model the whole
branch exists to compute is 2.1%, while getting reads into memory to feed it -- GAF parsing,
protobuf `Alignment` construction and destruction, and the allocator traffic both cause, which
the uncategorised bucket is mostly more of -- is around 12%. **The caller spends five times
longer building read objects than scoring them.**

## Speedups

### S1. Prefetch the next read window -- the largest single item

33.5% of the thread budget is a worker thread asleep in `waitpid`. Sites are visited in node-ID
order, so the window after this one is known before this one is consumed: spawn its query into a
second temp file and `waitpid` only when the current window runs out.

The headroom is measured, not assumed. Total CPU across `vg` and every child runs at a median
357% of a possible 1000% through the calling phase, and `-t 8` -- the crude version of the same
overlap -- takes chr20 to 173.6 s from 197.2 s, the calling phase falling 81.0 s to 58.8 s.

Bounded by two things worth saying out loud: the children are I/O-bound on the 21 GB read
database, so the ceiling is the disk's tolerance for concurrent reads; and `schedule_wgs.py`
already packs contigs onto the machine, so a whole-genome run gains less than a single-contig
one. This is worth the most in the A/B loop.

- **Size:** up to ~40 s of chr20 in isolation; less under the scheduler.
- **Risk:** moderate. Per-thread state, a second temp file per thread, and a spawn that must be
  reaped on every exit path including exceptions.
- **Verify:** byte-identical VCF and mosaic; `GAF-Base: N subprocess queries` unchanged.

### S2. Snarl decomposition, 46 s on one thread

23% of the run, four of five threads parked, and it is upstream code
(`IntegratedSnarlFinder::find_snarls_parallel`). `vg call -r snarls.pb` skips it, but the caller
builds its finder with `extra_node_weight` biasing the decomposition towards the reference path's
endpoints and `vg snarls` does not, so a cached file is not obviously the same decomposition.

- **First step is a measurement, not a change:** dump the caller's own snarls once, call chr20
  with `-r`, and check the VCF is byte-identical. If it is, this is 46 s per invocation for a
  one-off cost -- which matters most in the A/B loop, where the same contig is called repeatedly.
- **Risk:** low, because the gate is byte-identity.

### S3. Cut linkage chains at pinned sites

Generations 1-4 still decode the whole contig because each has at least one non-nested site, and
one such site keeps its whole 192,045-site chain alive. A per-site pin zeroes every state but
one, so it severs the path: a free site's answer depends only on the sites between its bracketing
pins. Decoding only those sub-chains would take generations 1-4 to near zero.

The care needed is in `phasing()`, whose windows are chained by a sequential seam pin, so
"skip a window" is not obviously safe there the way it is in `posteriors()`, whose windows are
independent by construction. A pin that `window_phasing` declines to apply -- the case where the
pinned pair cannot spell the constrained genotype -- does not sever, so the rule has to key on
pins that survive rather than on `Site::pinned`.

- **Size:** ~24 s of 197 s.
- **Risk:** the highest on this list. It changes the model's windowing.
- **Verify:** byte-identical VCF and mosaic, and the chr20 switch rate unchanged.

### S4. Stop building a protobuf `Alignment` per read

~12% of thread time, against 2.1% for the scoring it feeds. Every read fetched -- 14.2 M on
chr20 -- is parsed from GAF text into an `Alignment` with a `Path` of `Mapping`s and `Edit`s,
cached, then destroyed. The allocator attribution puts 21.7% of its samples inside
`gaf_to_alignment` and 17.5% inside `get_next_record_from_gaf`, with the protobuf destructors
next.

`WindowedSiteReadSource` already keeps a compact `bounds` array beside the reads precisely so the
common rejection path never walks a `Path`. The rest of the object exists for `touches()` and for
scoring. A leaner record -- node-id run, offsets, sequence, quality -- would cut the parse, the
allocation and the destruction together.

- **Size:** potentially the largest working-set win on the list, but unquantified until tried.
- **Risk:** high, and it is a refactor rather than a fix. Listed because the profile says the
  ratio is 5:1 the wrong way, not because it should be done next.

### S5. Send the reads-only ask upstream

A third of every `gbz-base query` builds a GFA subgraph that goes to `/dev/null`.
`planning/gbz-base-c-api-request.md` already argues this and has been sitting unsent. Sending it
costs nothing and is the only lever that reduces the work inside a query rather than rescheduling
it.

### S6. `panel_alleles` runs twice for a respecified site

`record_site` computes `panel_alleles(graph, travs)` when the site is first recorded
(`graph_caller.cpp:4665`); the barrier computes it again from `pr.travs`
(`graph_caller.cpp:5044`) before `respecify`. The traversal list does not change between the two
-- only the ploidy and the genotype do -- so the second GBWT walk repeats the first. `panel_alleles`
is 5.8% of the thread budget.

The single-sweep design note already said the panel should be *retained* rather than recomputed.
Carrying the vector on the `PendingRecord` does that.

- **Size:** small -- a low single-digit percentage of the 5.8%, depending on how many records
  reach `respecify`. Worth instrumenting the call count before writing the fix.
- **Risk:** low. **Verify:** byte-identical output.

### S7. `deliver()` scans every read in the window

3.39 billion reads considered to hand over 24.5 million, 138 rejected per one kept. The bounds
array is compact and the loop is already the cheap form, so it is 0.4% of samples -- but sorting
a window's bounds once on fetch and binary-searching per query would remove most of it.

- **Size:** ~0.4%. **Risk:** the lowest on the list. Reads are already adjudicated by `touches()`
  afterwards, so this changes speed and not which reads a site sees -- exactly what the existing
  comment says about the bounds test.

### S8. The record key goes through a `stringstream`

`std::hash<string>{}(print_snarl(snarl, false))` appears at five sites, once per record or site.
`print_snarl` does two `to_string` allocations and builds a `stringstream`, and the result is
immediately hashed and thrown away. The key is two node IDs and two orientation bits; packing
them into a `size_t` is both faster and collision-free, where hashing a formatted string is not.

- **Size:** small, a few tenths of a second. **Risk:** low, but it changes every key at once, so
  it is a single commit with a byte-identity gate -- keys are internal, so the VCF must not move.

## Redundancies and simplifications

### R1. `record()` and `respecify()` share 42 lines

`linkage_model.cpp:1230-1352` and `1394-1508`: 42 identical code lines, 35 of them in two
contiguous runs -- the `compact_of` lambda, the ploidy derivation, the `gls` fill from
`genotype_ln_likelihood`, and the arena writes. `respecify` is "record this again at a different
ploidy", so one `fill_entry(Entry&, ...)` helper would carry both.

### R2. The PhaseCall pair rendering appears three times

`linkage_model.cpp:1958-1967`, `2208-2217`, `2305-2314` -- `trav_first`/`trav_second`,
`render_phase_pair`, the `allele_*` writeback and the `fell_back` count, verbatim in the diploid
chain, the nested-strand pass and the haploid pass. One `finish_phase_call(PhaseCall&, const
Entry&, size_t& fallbacks)` removes two copies of a block whose comment already warns what
happens when it is got wrong.

### R3. `stage_render_record` + `emit_variant` appear twice

`graph_caller.cpp:5599-5611` and `5713-5725`, same long argument list both times, inserted into a
structure that is duplicated upstream anyway. Worth folding only as far as the branch's own lines
go.

### R4. Four windowing wrappers with one loop between them

`posteriors`, `haploid_posteriors`, `phasing` and `haploid_phasing` each re-implement the same
overlapping-window loop over `params.window` and `params.margin`. The diploid and haploid decodes
genuinely differ, but the window arithmetic does not, and S3 would have to touch all four. Fold
the loop first, then do S3 once.

### R5. Two functions that are never called

- `WindowedSiteReadSource::get_cache_entries()` -- `site_read_source.hpp:235`,
  `site_read_source.cpp:163`.
- `AlleleReadLikelihoods::read_name(size_t)` -- `allele_likelihood.hpp:97`,
  `allele_likelihood.cpp:49`.

Both declared, defined, and referenced nowhere else in `src/`, tests included.

### R6. `--read-min-mapq` is referenced nowhere outside `call_main.cpp`

Of the 31 options the branch adds, this is the only one with no use in the eval harness, no use
in `test/t/*.t`, and no plumbing into any caller. `--no-atomize-blocks` is also unused by both,
but it is the escape hatch for a default that changed on 08-23 and wants a test rather than
deletion. `--atomize-blocks`, `--nested` and `--phased` look like no-ops now that each names a
default, but each sets an `*_explicit` flag that turns a silent decline into an error, which
earns their keep.

### R7. Five counters are zero on all 24 contigs

Checked across every log in `work/wgs-atomize`:

| counter | values seen genome-wide |
|---|---|
| chain steps spaced along the settled parent traversal | 0 |
| nested chains with no single carrying traversal | 0 |
| nested chains whose parent record was not found | 0 |
| nested sites carried on both parent strands | 0 |
| nested sites with no phased parent | 0 |

These are invariants, and the branch's own convention is that an invariant belongs in an
assertion rather than in a line of output. Two neighbours must stay counters: "the frame would
reorder" takes 0 and 4, and "on neither strand" runs to 76,808 -- entirely on chrX and chrY,
where a haploid parent has no second strand, with 0 of them carrying a VCF line. That one is
correct behaviour rather than a defect, but it is not constant and cannot be asserted away.

Demoting the five would shorten the `frames:` and `nested strands:` lines, which are currently
the longest in the log and mostly zeros.

## What upstream duplication is *not* ours

The clone detector's largest hits in `graph_caller.cpp` -- the reference-path lookup at
5265-5294 against 6084-6113, the reference traversal build at 5360-5370 against 6153-6163, the
VCF header FORMAT lines in four places -- are all present verbatim at the merge-base. They are
the pre-existing `FlowCaller` / `NestedFlowCaller` duplication, about 420 lines of the latter,
reachable only through `-A`. Not this branch's to fix, and the `child_support_map` bug already
filed against that path lives there.

## Suggested order

Cheap and safe first, so the byte-identity gate is exercised on small changes before it is
trusted on large ones:

1. **R5, R6, R7** -- one commit. Two dead functions, one dead option, five counters to
   assertions. No behaviour change; the gate is that the assertions do not fire on chr20 and
   chrX.
2. **S2** -- the `-r` measurement. An experiment, not a change, and it either buys 46 s per
   invocation or is closed with an answer.
3. **S6, S7, S8** -- three small speedups behind one byte-identity gate each.
4. **R1, R2, R4** -- the three extractions, R4 last because S3 depends on it.
5. **S1** -- prefetch. The largest win, and the first item that needs its own care.
6. **S5** -- send the letter. Independent of everything above and costs an afternoon.
7. **S3** -- chain cutting, once R4 has given it one window loop to change.
8. **S4** -- left explicitly unscheduled. The profile says the ratio is wrong by 5:1, but it is a
   read-path refactor and wants its own plan.


---

# Working the list: outcomes

Recorded as each item was attempted, including the ones that did not survive contact.

## Step 1 -- R5 done, R6 and R7 withdrawn

**R5 done.** `WindowedSiteReadSource::get_cache_entries` and `AlleleReadLikelihoods::read_name`
deleted. Their members stay: `cache_entries` sizes the per-thread cache and `read_names` is
written into the likelihood dump directly, so only the accessors were dead.

**R6 withdrawn.** `--read-min-mapq` is not dead. It has no use in the harness and none in the
test suite, which is how it reached the list, but `SiteReadFilter::min_mapq` is compared in both
read paths (`site_read_source.cpp:42` and `:172`). Deleting it would have removed a working
filter on the strength of nobody having tested it. It now has a test instead: the fixture's
mapping qualities run 0 to 60, so a floor of 61 must leave no record standing, and a floor of 30
takes total DP from 1077 to 290. `t/18_vg_call.t` 309 -> 311.

**R7 withdrawn.** The premise was that an invariant belongs in an assertion rather than a line of
output. It already is one: `t/18_vg_call.t:806` asserts the literal text
`0 carried on both parent strands (0 with a line), 0 on neither (0 with a line), 0 with no phased
parent`, and its comment records that three of those were real populations on chr20 -- 440, 0 and
19 -- each tracing to the same bug. The verbosity is the detector. Suppressing the zero clauses
would have broken the test that catches the regression.

Also wrong in the original item: "chain steps spaced along the settled parent traversal" is not
an invariant at all. It is zero because every nested chain on this graph has a reference position;
on a graph without a covering reference it is the path that fires.

## Step 2 -- S2 done, and it works

`vg snarls` reproduces the caller's decomposition exactly, but only with **both** flags:

- `-P <ref path>` -- `snarls_main.cpp:284` and `call_main.cpp:1424` apply the same
  `EXTRA_WEIGHT` to the same first and last node of each reference path.
- `-T` -- include trivial snarls. Without it the top-level count still matches at 165,408, but the
  nested structure does not: chr20 gained 236 records and 714 lines moved, because the symbolic
  projection keys chain symbols on child chain boundaries.

| chr20 | wall | peak RSS |
|---|---|---|
| decomposing in-process | 197.2 s | 3.93 GB |
| `-r` with `vg snarls -T -P` | **148.5 s** | 3.78 GB |

VCF body and mosaic byte-identical. `vg snarls -T -P` itself costs 44.6 s, so a single
whole-genome pass breaks even; the gain is every run after the first, which is what the A/B loop
is. Wired into `prep_wgs.sh`, and `call_wgs.sh` passes `-r` when the file is there -- optional,
not required, which the byte-identity is what licenses.

## Step 3 -- S6 and S7 dropped on arithmetic, S8 kept but inverted

**S6 dropped.** `panel_alleles` does run twice for a respecified site, but respecify fires only
where a chain's ploidy changes: 2,537 revised plus 518 gained on chr20, against about 222,000
sites. That is 1.4% of the panel calls, and `panel_alleles` is 5.8% of the thread budget, so the
item is worth **0.08% of the run**. Not worth threading the vector through `PendingRecord`.

**S7 dropped.** The 3.39 billion reads considered sounds large and is not: the profile puts
`deliver` and `touches` together at 61 samples, **1.2 thread-seconds of 319**. About one cycle an
iteration, which is what a streaming two-compare loop over a compact bounds array should cost.
Sorting and binary-searching would save on the order of a second.

**S8 kept, in the opposite direction.** Packing the boundary node IDs instead of hashing the
printed snarl is about 0.15 s faster, and the byte-identity gate rejected it: chr20's 19,472
linkage-quality patches all stopped landing.

The reason is a **seventh producer** that the review missed. `write_variants` recovers a buffered
line's key by re-hashing the line's own ID column (`graph_caller.cpp:990`), which `emit_variant`
set from `print_snarl`. That is what lets a compressed record find its linkage entry with nothing
carried alongside it, and what makes the key survive `--translation`, where both sides print the
translated form. A key in a file cannot be repacked.

So the form is fixed by that site, and the six copies are centralised into
`VCFOutputCaller::record_key_of` with the coupling written down once. No speed change; the point
is that six copies of an expression that must agree with a string in a file cannot drift any more.

## Step 4 -- R1, R2 and R4 all done, and R1 was worth 10 s

**R2.** The three copies of the PhaseCall pair split -- diploid chain, nested-strand pass, haploid
pass -- folded into `finish_phase_call`. They had drifted: only the diploid copy widened
`order_arbitrary` on a fallback. The other two are ploidy-1 populations where the pair is one allele
twice and the test cannot fire, so folding it in changes nothing. Byte-identical.

**R1.** `record()` and `respecify()` shared 42 lines building the compact space; `CompactSite` makes
their agreeing structural rather than maintained by hand. The unexpected part is the lookup: both
built a `map<int, int>` from candidate traversal to compact allele and threw it away -- one map and
up to 127 node allocations **per site, on every site of the contig**. `compact_allele_space` already
returns its space sorted, so a binary search answers the same question with no allocation.

**chr20 197.3 s -> 187.2 s.** A tidy-up worth 10 s, which is more than S6, S7 and S8 put together
were.

**R4.** The four windowing wrappers fold into two templates, not one -- which is the point. A
marginal is a per-site quantity so its windows are independent; a path's are not and each has to be
pinned to what the previous one decided. Two shapes, each written once, with the reason it is that
shape attached to it. Byte-identical.

## Step 5 -- S1 measured, and the measurement changes the recommendation

The ceiling on overlapping the subprocess wait is what concurrency alone can recover:

| threads | wall | CPU multiple | peak RSS |
|---|---|---|---|
| 5 | 187.6 s | 3.04 | 4.19 GB |
| 8 | **159.1 s** | 4.06 | 4.98 GB |
| 10 | 157.5 s | 4.69 | 5.09 GB |

Identical output at all three. The whole overlap is worth about 30 s, and it **saturates at 8** --
beyond that the limit is the read database's tolerance for concurrent random reads, not cores.

So prefetching is a **memory** optimisation, not a speed one: `-t 8` already buys the 28 s, for
0.8 GB of extra per-thread cache. Prefetch would buy the same for one extra in-flight query per
thread, which matters where memory binds -- `schedule_wgs.py` packing contigs against 24 GB -- and
not otherwise. **Not built.** The actionable finding is `-t 8` for single-contig work, and leave the
whole-genome scheduler at 5, where packing contigs fills the machine anyway.

## Step 6 -- S5 rewritten, NOT sent

`planning/gbz-base-c-api-request.md` said its own ask was wrong and needed rewriting. It has been
rewritten: the C-interface proposal is superseded by a request for a **reads-only query**, since a
third of every `gbz-base query` builds a GFA subgraph we send to `/dev/null`, plus a second small
ask to skip unknown node IDs rather than fail the query. The two questions that building it answered
-- MAPQ and base qualities both survive; `--alignments` defaults to a trap -- are now notes rather
than questions.

**Sending it is not mine to do.** It is a message to an upstream maintainer, and it needs a decision
on whether it goes as a GitHub issue or a direct message.

## Step 7 -- S3 sized, not built

The caller now reports how much of each generation's decode is live, because assuming was what got
the previous guard's reasoning wrong:

| generation | sites decoded | live | pinned |
|---|---|---|---|
| 1 | 209,244 | 17,199 | 192,045 |
| 2 | 211,480 | 2,236 | 209,244 |
| 3 | 211,930 | 450 | 211,480 |
| 4 | 211,952 | 22 | 211,930 |

Generations 1-4 cost 24.4 s of 187.6, and sub-chaining at the pins would take that to about 2 s.

The original framing was wrong in a way worth recording: it assumed only a handful of non-nested
sites arrive after generation 0, because nested sites are held out of the diploid runs. Generation 1
has 17,199 of them, so the sub-chains at that generation would be many and short rather than few.

Not built. Two things have to be right and neither is free: the phase set is the chain's first
position and has to keep coming from the whole chain rather than from a piece of it, and the
severing rests on every pin being *accepted*, where `window_phasing` declines one whose pair cannot
spell the site's constrained genotype. The gate available is one chromosome.

## Step 8 -- S4 left where it was

Unscheduled by the plan and still unscheduled. The profile says the caller spends five times longer
building protobuf `Alignment` objects than scoring them, which is worth acting on, but it is a
read-path refactor and wants its own plan rather than a slot at the end of a cleanup.

## Where chr20 ended up

| | wall |
|---|---|
| before any of this | 328.3 s |
| after the two profile fixes | 197.2 s |
| after the cleanup pass | 187.6 s |
| with cached snarls | **148.5 s** |
| with cached snarls, at `-t 8` | untested, expected near 120 s |

Every step byte-identical on the VCF body and the mosaic. `vg test` 12,547,761 assertions in 858
cases; `t/18_vg_call.t` 311/311.
