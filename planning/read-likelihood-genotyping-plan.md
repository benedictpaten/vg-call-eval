# Appendix: source reading, prior art, and settled decisions

> **This is not a plan any more.** It began as the pre-implementation plan for read-level genotype
> likelihoods; the feature is built (PR vgteam/vg#4990) and
> [vg-read-likelihood-design.md](vg-read-likelihood-design.md) is the authoritative description of
> what exists. What survives here is the part the design doc does not carry: **the reading of vg's
> source that constrained the design, the prior art already in the tree, and the reasoning behind
> each settled decision**.
>
> Six sections were removed as duplicative rather than edited, and where to find each:
>
> | Removed | Now in |
> |---|---|
> | §1 new components, interfaces | design §3 |
> | §2 the model, derived on the unnormalised score scale | design §4 — algebraically equivalent, in the row-normalised form that shipped |
> | §3 read retrieval, four backends | design §6 |
> | §4 performance projection | design §5.4, and measured in design §6.5 |
> | §6 staging table | design §8 |
> | §7 validation, including the scoring unit-test table | design §7 |
>
> The model section is the one worth knowing about: it derives the same result a different way, on
> the raw score scale. If the row-normalised form in design §4.2 is ever in doubt, that derivation
> is in this file's history.

## 0. The single most important structural fact

This change **hangs entirely off the `SnarlCaller` interface** and needs no modification to
the traversal finders, the snarl decomposition, `GraphCaller`, or `FlowCaller`'s control flow.

`SnarlCaller::genotype()` already receives everything needed to locate the site and its
candidate alleles:

```cpp
pair<vector<int>, unique_ptr<CallInfo>> genotype(
    const Snarl& snarl,
    const vector<SnarlTraversal>& traversals,   // <- the enumerated alleles
    int ref_trav_idx, int ploidy,
    const string& ref_path_name,
    pair<size_t,size_t> ref_range);
```

So the read×allele matrix can be built **inside** the new caller's `genotype()` as an
implementation detail. No interface change, no touching the 5 existing `genotype()` call
sites in `graph_caller.cpp`.

Two constraints on the new class, both verified in the source:

- `FlowCaller`, `NestedFlowCaller`, and `LegacyCaller` all `dynamic_cast` the caller to
  `SupportBasedSnarlCaller` ([call_main.cpp:896](vg/src/subcommand/call_main.cpp:896),
  [graph_caller.cpp:2060](vg/src/graph_caller.cpp:2060)) to reach `get_support_finder()` — used
  for the flow-traversal node/edge weights and the avg-vs-min flow decision. **The new caller
  must therefore subclass `SupportBasedSnarlCaller`**. `-k` (the pack file) is required only
  where *enumeration* needs support, which means only with `FlowTraversalFinder`. Reads are used for
  *genotyping* only, so with `-g`/`-z` — where `GBWTTraversalFinder` enumerates from recorded
  haplotypes rather than node/edge weights — nothing in the path consults support and the pack file
  is not needed at all. A `NullTraversalSupportFinder` satisfies the interface there. Note the caller
  must then also override `get_skip_allele_fn()`, whose support-based version would prune every
  allele at every site when support reads zero.
- `VCFGenotyper` takes a generic `SnarlCaller&` ([graph_caller.hpp:259](vg/src/graph_caller.hpp:259)),
  so `-v` VCF-regenotyping gets read-level likelihoods **for free**, with zero extra work.

## 5. Gotchas found in the existing code

**5.1 `update_vcf_info` receives *deduplicated* traversals, and two of them may not be traversals.**
`emit_variant` builds `site_traversals` / `site_genotype` by deduplicating alleles **by
sequence string** and dropping every uncalled allele, then calls
`snarl_caller.update_vcf_info(snarl, site_traversals, site_genotype, ...)`
([graph_caller.cpp:665](vg/src/graph_caller.cpp:665)). Three distinct hazards, not one:

- **The indices no longer match the matrix.** `PoissonSupportSnarlCaller` copes by simply
  recomputing support; recomputing the *matrix* would mean redoing every read.
- **⚠ `site_genotype` can contain negative sentinels.** `STAR_ALLELE_MARKER = -2` and
  `MISSING_ALLELE_MARKER = -1` ([graph_caller.hpp:29](vg/src/graph_caller.hpp:29)) are pushed into it
  in `--top-down` nested mode, and a default-constructed placeholder `SnarlTraversal()` is pushed
  into `site_traversals` to stand in for the star allele
  ([graph_caller.cpp:539-566](vg/src/graph_caller.cpp:539)). Indexing a GL table by `site_genotype`
  reads out of bounds. That is a crash, not a mismatch.
- **Dedup is by allele *string*, so two distinct scored traversals can collapse into one VCF allele.**
  GL for that allele needs an explicit aggregation rule — use `max` over the merged set, which keeps
  GL a likelihood of the emitted allele rather than of an arbitrary mixture.

Options: (a) memoize the matrix per-thread keyed by snarl bounds and remap `site_traversals`
back by equality; (b) carry the per-genotype LLs in `CallInfo` and emit `GL` only over
scored genotypes; (c) extend `emit_variant` to also return a
`site_index -> vector<scored_index>` mapping. **Recommend (c)**, with `GL` emitted only over the
emitted allele set (which is what VCF requires anyway) and the sentinels skipped explicitly — a
small, contained change to `VCFOutputCaller::emit_variant` that also fixes the same latent
awkwardness for the Poisson caller. Budget real time for this; it is the kind of mismatch that costs
a day if unplanned.

**5.2 Legacy prior art exists and is half-wired.** `src/genotyper.{hpp,cpp}` (the old
`vg genotype`) already contains this design:

- `Genotyper::Affinity` ([genotyper.hpp:36](vg/src/genotyper.hpp:36)) is a read×allele entry
  and **already has a `likelihood_ln` field**;
- `get_affinities` ([genotyper.hpp:228](vg/src/genotyper.hpp:228)) returns
  `map<const Alignment*, vector<Affinity>>` — literally the matrix — by realigning each read
  to a per-allele graph ([genotyper.cpp:707-753](vg/src/genotyper.cpp:707)). Its
  reads-at-a-snarl gather loop ([genotyper.cpp:663-745](vg/src/genotyper.cpp:663)) is the
  canonical pattern: for each node in the snarl contents, `aug.get_alignments(node_id)`, then
  collect the flanking node IDs those reads touch and build a surrounding graph. Note it
  **bails out entirely when `reads × paths > 1000`** ([genotyper.cpp:692](vg/src/genotyper.cpp:692))
  — independent confirmation that the R×K cost in §4 is the real constraint, and that the
  original authors hit it too;
- it populates `likelihood_ln` from `score_to_unnormalized_likelihood_ln`
  ([genotyper.cpp:748](vg/src/genotyper.cpp:748)).

**But `likelihood_ln` is written and never read.** Grep confirms the only references are the
two assignments. `get_genotype_log_likelihood`
([genotyper.cpp:1010](vg/src/genotyper.cpp:1010)) instead uses the **binary `consistent`
flag** plus per-strand counts and an interval-censored binomial — a count model, not the
continuous mixture. So the scaffolding is there and the actual likelihood was never
connected. This is worth reading closely as prior art (especially the informative-read
filtering and the fwd/rev handling), while **not** adopting its genotype model. It is built
on `VG`/`AugmentedGraph` protobuf structures, so it needs porting, not reuse.

**5.3 Option letters are nearly exhausted.** `vg call`'s getopt string leaves only
`jnquwxy` / `DEFHJKQUVWXZ` free. Use long-only options (`--gam`, `--gaf`,
`--read-likelihood`, `--max-realign-alleles`) rather than burning short letters.

**5.4 Two pre-existing bugs in the neighbourhood** (documented separately; fix or avoid):
`depth_err` is assigned from a malformed ternary at
[snarl_caller.cpp:602](vg/src/snarl_caller.cpp:602) — `depth_info.second ? !isnan(depth_info.second) : 0.`
— so it is always 0.0 or 1.0; and `GBWTTraversalFinder` has an index-mismatch in its backward-dedup
path ([traversal_finder.cpp:3519-3547](vg/src/traversal_finder.cpp:3519)).

**Correction: this bug is inert.** The only consumer of `depth_err` inside
`genotype_likelihood` is commented out, deliberately, with the rationale that the small bin sizes make
the binned-coverage error far too large to be useful; `depth_err` is otherwise only carried on the
`CallInfo` and printed in debug output, and never reaches the VCF. Fixing the ternary produces
byte-identical calls (verified over three simulated 400 kb replicates). So it is worth fixing as
hygiene — a latent trap for whoever re-enables that line — but it is **not** a prerequisite for stage 4
and the Poisson caller is **not** a mis-parameterised baseline. Earlier drafts of this document claimed
otherwise; that claim was wrong.

**5.5 The read source and the pack file must agree on which reads exist.** `vg pack -Q` sets *both*
`min_mapq` and `min_baseq` from the one value ([pack_main.cpp:148](vg/src/subcommand/pack_main.cpp:148)),
and `--trim-ends` drops read ends. Traversal enumeration is pack-driven while genotyping is
read-source-driven, so absent a matching eligibility filter the two disagree about what evidence
exists at a site. Stage 0 needs a stated policy: default to mirroring the pack file's filters, and
decide explicitly about secondary/supplementary alignments and duplicates, which nothing on vg's pack
path excludes either.

## 8. Decisions

**All settled** — the authoritative record is the design doc's §10 table, which also lists what each
decision leaves available as a later, purely additive addition. Retained here with the reasoning that
produced each one.

1. ~~**gbz-base integration mechanism**~~ — **decided, then revised by measurement** (design §6.4): the
   subprocess is what shipped, and it should stay. The `extern "C"` shim was framed as the target and
   the subprocess as an interim; timing showed `fork`/`execvp` is 3 ms of a 62 ms query, so the shim
   would recover ~5%. What is worth asking upstream for is a reads-only query that skips the 20 ms of
   subgraph construction we throw away. Direct SQLite reading is rejected — a
   400–600 line decoder against an explicitly unstable format (three breaking changes in ~4 months,
   exact-match version check, no C API) is a liability we can't absorb. Remaining sub-question for
   that conversation: **should gaf-base support be an optional compile-time feature in vg?** A Rust
   `staticlib` puts the Rust toolchain in every build environment including release CI, so this may
   matter more to vg than the ABI shape does.
2. ~~**Read source priority order**~~ — **decided** (design §6): in-memory prototype first (stage 0), then
   `.gai` GAM, then gaf-base, with tabix GAF only if a concrete user needs it. The decisive point, found
   later: the first two need **no new dependency of any kind**, so the whole first version and the scale
   path are buildable from code already in the tree, and every dependency question defers to stage 7.
3. ~~**Depth term**~~ — **decided** (design §4.6): pure `P(reads | G)`, no coverage term. Keeps the two
   evidence sources separable so a depth term added later is measurable against a clean baseline.
4. ~~**Genome-wide vs targeted / replace or option**~~ — **decided**: an **option**, not a replacement.
   `--read-likelihood` is opt-in, the Poisson caller stays the default, and a run without the flag must
   produce byte-identical output (design §1, §7). That constraint is what makes the feature safe to
   merge early, and it lowers the bar for the first version — the accepted limitations in design §4.6
   can ship documented rather than solved, because nobody gets them unless they ask. On affordability:
   genome-wide looks fine now that scoring is DP-free (design §5.4), with read retrieval rather than genotyping
   as the limiting factor.
5. ~~**Soft clips**~~ — **decided** (design §5.2): the scoring window is the maximum span of the read
   within the site, taken over all alleles. Soft clips lie outside that span and are **excluded by
   default**, with a flag to include them. Revisit if stage 4 shows SV recall suffering.
6. ~~**MAPQ calibration**~~ — **decided** (design §4.3): take MAPQ at **face value**. No recalibration table
   and no within-site correction — the suspected within-site-ambiguity problem was investigated against
   the source and does not apply. Ship only the `e_max` clamp (load-bearing because
   `phred_to_prob(0) = 1`) and a flag to switch the term off, then measure at stage 4. A 61-entry
   MAPQ→`e_r` table fitted on `vg sim -a` reads is the escape route if measurement demands it, not
   stage-1 work. Note any such table is specific to the mapper, graph, read length and error profile,
   and the clamps are scorer-specific — they will not transfer between the quality-adjusted and
   non-quality-adjusted scorers (design §5.2).
7. ~~**Prior**~~ — **decided**: **uniform** over genotypes. Matches the current Poisson caller's `GP`
   semantics and keeps likelihood and prior separable, so the likelihood model can be evaluated without
   a prior confounding the comparison. An HWE or het-biased prior — as in the legacy `Genotyper`
   ([genotyper.hpp:102-116](vg/src/genotyper.hpp:102)) — multiplies in later without touching anything
   else.

**One risk this document did not originally identify**, recorded here because it changes how the output
should be described: assuming reads independent (design §4.1) makes reported `GQ` **over-confident, and worse
as depth rises** — the product accumulates confidence like `R` rather than `√R`. Treat `GQ`/`GL` as
useful for ranking rather than as calibrated probabilities, and say so wherever they are documented.
The fix is a scalar `read_weight` discount fitted against binned `GQ`-vs-actual-error; it is additive
and does not change the model's shape (design §4.4, §4.6 risk 2).
