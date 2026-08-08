# Characterization of `vg call`

> **Status.** This describes `vg call` as it was *before* the read-likelihood work, plus the classes
> that work added, marked `[ADDED]`. It is the orientation document; the design and the as-built record
> live in [vg-read-likelihood-design.md](vg-read-likelihood-design.md). Two behaviours described below
> have since changed: nested `--top-down` now genotypes a child at the ploidy that actually traverses
> it (§7 of the design doc), and `-k` is no longer unconditionally required — but only for
> `--read-likelihood` with haplotype enumeration. Poisson *genotyping* still consumes support, so
> `vg call -z` without `-k` remains an error.

Source: `vgteam/vg` (cloned shallow into `./vg`), commit at clone time (2026-07-28).
Entry point: [src/subcommand/call_main.cpp](vg/src/subcommand/call_main.cpp)

## 1. What it does

`vg call` turns read-support evidence on a variation graph into variant calls. It has two
distinct modes selected by CLI flags:

- **De novo / genotyping-from-graph mode** (default): enumerate candidate alleles at each
  bubble ("snarl") in the graph, score them against read support (from `vg pack`), and emit
  a genotype (VCF or GAF).
- **VCF-regenotyping mode** (`-v`): given an existing VCF whose alleles are already embedded
  in the graph as alt-paths (from `vg construct -a`), just decide which alleles are present
  in this sample's read data.

It always needs:
1. A graph (plain `PathHandleGraph`, or GBZ) with one or more REFERENCE/GENERIC paths to
   call against.
2. A `vg pack` support index (`-k`), giving per-node/edge read coverage — this is what stands
   in for a BAM/pileup in linear variant calling.
3. A snarl decomposition (computed on the fly, or loaded via `-r`).

## 2. Top-level pipeline (`main_call`)

```
load graph (GBZ or PathHandleGraph)
  -> apply path-position / vectorizable overlays as needed
  -> resolve reference paths (-p/-P/-S, or all REFERENCE/GENERIC paths) + per-path ploidy (-d/-R)
  -> load or compute snarls (-r, else IntegratedSnarlFinder: Cactus-graph decomposition)
  -> build a SnarlCaller (genotyping model) wired to a TraversalSupportFinder over the pack file
  -> build a GraphCaller (orchestrator), choosing:
       -v          -> VCFGenotyper   (regenotype known VCF alleles)
       --legacy    -> LegacyCaller   (old RepresentativeTraversalFinder-based algorithm)
       else        -> FlowCaller     (Yen's/max-flow or GBWT traversal enumeration)
                        --top-down   -> FlowCaller in nested mode (parent->child propagation)
                        --bottom-up  -> NestedFlowCaller (deprecated bottom-up merge)
  -> graph_caller->call_top_level_snarls() or call_top_level_chains()
       (RecurseOnFail | RecurseAlways (-A) | RecurseNever (top-down))
  -> write VCF header + sorted VCF body, or stream GAF
```

Everything downstream is essentially: **decompose graph into sites -> enumerate allele
traversals per site -> score by support -> pick a genotype -> emit**.

## 3. Algorithmic stages in detail

### 3.1 Snarl/chain decomposition (site definition)

- `IntegratedSnarlFinder` ([integrated_snarl_finder.hpp](vg/src/integrated_snarl_finder.hpp)):
  builds a Cactus graph over the `HandleGraph`, finds bridge-tree paths / simple cycles, and
  recursively decomposes into a tree of **snarls** (bubbles, possibly containing cycles) and
  **chains** (linear runs of nested snarls). Rooted at whichever bridge path/cycle has the
  most "fixed sequence" (reference paths get extra weight via `extra_node_weight` so the
  decomposition roots along the reference).
- Can be precomputed once with `vg snarls` and reused (`-r`), since this is the expensive,
  purely-topological part and doesn't depend on sample data.
- `--chains` (`-I`) groups consecutive trivial snarls in a chain into one "fake snarl" so a
  single VCF record can span a longer stretch of reference (`GraphCaller::break_chain`).

### 3.2 Read support (evidence layer)

- `Packer` (from `vg pack`) stores per-node and per-edge coverage counts (and optionally base
  qualities/MAPQ) compressed/packed over the graph.
- `TraversalSupportFinder` hierarchy ([traversal_support.hpp](vg/src/traversal_support.hpp)):
  - `PackedTraversalSupportFinder` — raw node/edge support lookups from the `Packer`.
  - `CachedPackedTraversalSupportFinder` — adds an LRU cache per thread (support queries are
    the hot path).
  - `NestedCachedPackedTraversalSupportFinder` — additionally memoizes support already
    computed for child snarls, needed by the bottom-up nested caller.
  - Support of a traversal = combination of its nodes'/edges' support, with a switch between
    **minimum** support (small traversals) and **average** support (long traversals, so one
    low-coverage base doesn't kill an otherwise well-supported long allele) — threshold
    `avg_trav_threshold`/`avg_node_threshold` (50bp by default).
  - `get_traversal_set_support` implements support *splitting*: when several alleles share
    graph structure (e.g. two alleles sharing a flanking node), the shared support is divided
    among them rather than double counted, and can be computed "exclusively" (only structure
    private to one traversal) for tie-breaking.
- `algorithms::binned_packed_depth_index` / `coverage_depth.*` ([coverage_depth.hpp](vg/src/algorithms/coverage_depth.hpp)):
  precomputes expected background depth in bins along each reference path (multiple bin
  sizes, log-spaced), used as the Poisson model's expected-depth prior — this is the
  graph analogue of a "local coverage" estimate a linear caller gets from BAM depth.

### 3.3 Traversal (allele) enumeration

`TraversalFinder` implementations ([traversal_finder.hpp](vg/src/traversal_finder.hpp)) turn a
snarl into a set of candidate `SnarlTraversal` (start-to-end walks = candidate alleles):

- `FlowTraversalFinder` (default, non-GBWT path): finds the **K widest** (max-min edge/node
  weight, i.e. max-flow-flavored) traversals via a Yen's-algorithm-style k-shortest/widest
  path search, weighted by node/edge support callbacks. Falls back to average-flow instead of
  bottleneck-flow for long traversals (`greedy_avg_flow`). Always returns the reference
  traversal first.
- `GBWTTraversalFinder`: when `-g`/`-z` is given, only enumerates traversals that correspond
  to actual haplotypes recorded in a GBWT/GBZ haplotype index — a BFS that only branches
  where a real haplotype thread branches, so it can't hallucinate a haplotype combination
  that was never observed.
- `VCFTraversalFinder` (used by `VCFGenotyper`): maps VCF alleles to graph alt-paths (from
  `vg construct -a`) and brute-force / greedily assembles per-haplotype traversals, with
  support-based pruning (`skip_alt`) to keep combinatorics down at dense/multi-allelic sites.
- `RepresentativeTraversalFinder` (used only by `--legacy`/`LegacyCaller`): the original vg
  call algorithm — BFS out from each node/edge to the nearest points back on a reference
  `PathIndex`, forming "bubbles" one node/edge at a time. Superseded by `FlowCaller` because
  it's harder to guarantee it finds every well-supported off-reference traversal.
- `PathTraversalFinder` / `PathBasedTraversalFinder`: traversals taken directly from embedded
  paths (used as a component elsewhere, e.g. finding the reference traversal).

### 3.4 Genotyping (statistical model)

`SnarlCaller` hierarchy ([snarl_caller.hpp](vg/src/snarl_caller.hpp)) takes
`(snarl, traversals, ploidy, ref info)` and returns a genotype (indices into the traversal
list) + opaque `CallInfo` for VCF annotation:

- `RatioSupportSnarlCaller` (`-B`/`--bias-mode`, legacy default): picks the genotype by
  comparing the best/second-best allele's support against a **het-bias ratio**
  (`-b`, default 6): if best allele has >6x the support of the runner-up, call homozygous,
  otherwise heterozygous. Simple, threshold-based, no explicit error model.
- `PoissonSupportSnarlCaller` (new default): a proper probabilistic model, "inspired in part
  by Paragraph":
  1. Rank all traversals by support; consider the top `top_k` (20) as primary allele
     candidates and top `top_m` (100) as secondary/partner candidates, to bound the
     genotype search space combinatorially (candidate genotype pairs, not all C(n,2)).
  2. For each candidate genotype (1 or 2 alleles), compute a Poisson log-likelihood
     (`genotype_likelihood`): expected depth for the site comes from the binned depth index;
     observed support is split evenly across homozygous pairs; an error rate
     (`baseline_error_small`/`baseline_error_large`, chosen by variant size vs. a length
     threshold, adjustable via `-e`) models the chance that supporting reads for one allele
     are noise, and reads for *other* traversals are modeled as an error term too (so the
     genotype has to "explain" reads elsewhere in the site, not just maximize its own count).
  3. Best-likelihood genotype wins; **GQ** = phred-scaled difference between best and
     second-best log-likelihoods; **posterior** from log-likelihoods under a uniform
     genotype prior (`PoissonCallInfo`).
  4. Optional insertion-bias multiplier when the site looks like an insertion (longest allele
     ≫ reference length).
- Both share `SupportBasedSnarlCaller` base: minimum total support (`-m`, default 2/4),
  minimum site depth, and an `get_skip_allele_fn()` hook the traversal finders (esp.
  `VCFTraversalFinder`) use to prune poorly-supported alleles *before* combinatorial
  enumeration blows up.

### 3.5 Orchestration / recursion (`GraphCaller`)

- `GraphCaller::call_top_level_snarls` walks only the top-level (unnested) snarls in
  parallel, calling `call_snarl()` (pure virtual, implemented per caller type) on each.
  `RecurseType` controls what happens to children of a snarl:
  - `RecurseOnFail` (default): only recurse into children if the parent snarl's call failed
    (e.g. too big, no ref path) — keeps output flat/independent normally.
  - `RecurseAlways` (`-A`/`--all-snarls`): call every nested snarl independently regardless of
    parent outcome — useful for genotyping small variants nested inside a structural variant.
  - `RecurseNever` (`--top-down`): the caller (`FlowCaller` in nested mode) handles recursion
    itself.
- **`--top-down` nested mode** (current, actively-developed approach, in `FlowCaller`): calls
  a parent snarl, then for each child snarl restricts its candidate traversals to those
  consistent with the parent's called alleles (`ChildTraversalSets`/`find_child_traversal_set`)
  — i.e. genotype propagates from parent to child, and phase is inherited. Writes `LV`/`PS`
  VCF tags recording nesting level/parent site. Supports **star alleles** (`-Y`): a parent
  haplotype that doesn't traverse a given child snarl at all gets a `*` allele there rather
  than a missing/reference call, so ploidy stays consistent up and down the nesting.
- **`--bottom-up` mode** (`NestedFlowCaller`, marked deprecated in favor of top-down): calls
  children first into a `CallTable`, then flattens child alleles into the parent's allele
  strings recursively (`flatten_reference_allele`/`flatten_alt_allele`) when emitting the
  parent VCF record — the reverse data-flow direction from top-down.
- **Traversal clustering** (`-L`/`--cluster`, [traversal_clusters.hpp](vg/src/traversal_clusters.hpp)):
  optionally merges traversals whose Jaccard similarity (over node-set) exceeds a threshold
  before (default) or after (`--cluster-post`) genotyping, to avoid splitting support across
  near-identical traversals from minor graph noise.
- **Chains instead of snarls** (`-I`/`--chains`): `call_top_level_chains` merges runs of
  small/trivial nested snarls in a chain into one synthetic "fake snarl" (bounded by
  `max_chain_edges`/`max_chain_trivial_travs`) before calling, producing longer, less
  fragmented VCF records.

### 3.6 Output

- `VCFOutputCaller` (mixed into most caller classes): builds the VCF header (contigs from
  resolved reference paths/lengths), and buffers variant lines per-thread as
  `(sort_key, string)` pairs to avoid keeping the huge in-memory `vcflib::Variant` around;
  `write_variants` sorts and flushes. Handles allele-string flattening
  (`flatten_common_allele_ends`), path-based allele encoding (GFA W-line/GAF-style,
  `add_allele_path_to_info`), and `LV`/`PS` nesting tags.
- `GAFOutputCaller`: alternative output as GAF alignment records — either raw candidate
  traversals only (`-T`/`--traversals`, no genotyping, useful for debugging/inspecting what
  the traversal finder sees) or full genotype calls as GAF (`-G`/`--gaf`).
- Both can add reference-path padding to short traversals (`-M`/`--trav-padding`) so GAF
  alignments have some flanking context.

## 4. Class map

```
GraphCaller (abstract: call_top_level_snarls/chains, recursion, break_chain)
 ├─ VCFGenotyper       (+ VCFOutputCaller, GAFOutputCaller)   -- regenotype input VCF
 ├─ LegacyCaller       (+ VCFOutputCaller)                    -- old algorithm (RepresentativeTraversalFinder)
 ├─ FlowCaller         (+ VCFOutputCaller, GAFOutputCaller)   -- default de novo caller, optional nested/top-down mode
 └─ NestedFlowCaller   (+ VCFOutputCaller, GAFOutputCaller)   -- deprecated bottom-up nested caller

SnarlCaller (abstract: genotype(), update_vcf_info/header)
 └─ SupportBasedSnarlCaller (shared support-cutoff plumbing)
     ├─ RatioSupportSnarlCaller     (-B: het-bias ratio heuristic)
     ├─ PoissonSupportSnarlCaller   (default: Poisson depth-based likelihood model)
     └─ ReadLikelihoodSnarlCaller   (--read-likelihood: explicit P(reads|genotype))   [ADDED]

TraversalSupportFinder (abstract: support queries over graph)
 ├─ PackedTraversalSupportFinder (reads from vg-pack Packer)
 │   └─ CachedPackedTraversalSupportFinder (+ LRU caches)
 │       └─ NestedCachedPackedTraversalSupportFinder (+ child-snarl support memo, for bottom-up)
 └─ NullTraversalSupportFinder   (reports no support, needs no Packer, so -k can be     [ADDED]
                                  omitted when enumeration does not consult support)

SiteReadSource (abstract: random-access reads by graph locality; thread-safe)          [ADDED]
 ├─ InMemorySiteReadSource (one streaming pass over GAM/GAF, bucketed by node ID)
 └─ WindowedSiteReadSource (abstract: quantises fetches to node-ID windows and         [ADDED]
     │  caches the last few per thread; subclasses supply only fetch_span())
     ├─ IndexedGamSiteReadSource (sorted GAM + .gai, one cursor per thread;            [ADDED]
     │    --gam-index. This is the one that delivered the memory reduction)
     └─ GafBaseSiteReadSource    (GAF-Base via `gbz-base query` as a subprocess;       [ADDED]
          --gaf-base. Runtime dependency on an external binary, nothing linked)

AlleleLikelihoodCalculator (abstract: produce the reads x alleles matrix)              [ADDED]
 └─ GraphAlignedAlleleLikelihoodCalculator (scores from each read's existing alignment)

TraversalFinder (abstract: enumerate SnarlTraversals for a site)
 ├─ FlowTraversalFinder          (Yen's/max-flow k-widest paths; default)
 ├─ GBWTTraversalFinder          (haplotype-consistent BFS; -g/-z)
 ├─ VCFTraversalFinder           (VCF-alt-path-guided; used by VCFGenotyper)
 ├─ RepresentativeTraversalFinder(BFS-bubble; used only by --legacy)
 ├─ PathTraversalFinder / PathBasedTraversalFinder (embedded-path traversals)
 └─ ExhaustiveTraversalFinder / TrivialTraversalFinder / *RestrictedTraversalFinder (older/simple, mostly unused by call_main)

MCMCCaller (+ VCFOutputCaller) -- exists (mcmc_caller.{hpp,cpp}) but NOT wired into
                                   main_call's option parsing; effectively orphaned/experimental.
```

## 5. Key files

| File | Role |
|---|---|
| `src/subcommand/call_main.cpp` | CLI parsing, wiring, top-level driver |
| `src/graph_caller.{hpp,cpp}` | `GraphCaller`/`VCFOutputCaller`/`GAFOutputCaller`, `VCFGenotyper`, `LegacyCaller`, `FlowCaller`, `NestedFlowCaller` |
| `src/snarl_caller.{hpp,cpp}` | Genotyping models: `RatioSupportSnarlCaller`, `PoissonSupportSnarlCaller` |
| `src/traversal_finder.{hpp,cpp}` | Allele/traversal enumeration strategies |
| `src/traversal_support.{hpp,cpp}` | Read-support queries over `vg pack` data |
| `src/algorithms/coverage_depth.{hpp,cpp}` | Binned background depth estimation (Poisson prior) |
| `src/integrated_snarl_finder.{hpp,cpp}` | Cactus-graph snarl/chain decomposition |
| `src/traversal_clusters.hpp` | Jaccard-similarity traversal clustering (`-L`) |
| `src/gref.hpp` | "gref" reference-path-cover naming convention (pangenome multi-ref handling) |
| `src/mcmc_caller.{hpp,cpp}` | Unused/orphaned MCMC-based caller |
| `src/read_likelihood_caller.{hpp,cpp}` | **Added.** `ReadLikelihoodSnarlCaller` |
| `src/allele_likelihood.{hpp,cpp}` | **Added.** reads x alleles matrix, mixture model, graph-implied scoring |
| `src/site_read_source.{hpp,cpp}` | **Added.** random-access read retrieval by node-ID range |
| `src/unittest/allele_likelihood.cpp` | **Added.** the model, on hand-built matrices |
| `src/unittest/allele_likelihood_scoring.cpp` | **Added.** the scoring, on hand-built graphs/alignments |
| `test/t/18_vg_call.t`, `test/call/` | Bash-TAP integration tests / fixtures |

## 6. Notable adjustable behaviors (relevant for adaptation)

- Two independently swappable axes: **which traversals get enumerated** (`TraversalFinder`)
  and **how they get scored/genotyped** (`SnarlCaller`) — `FlowCaller` is generic over
  `TraversalFinder&`, so a new traversal source just needs to implement that interface.
- Poisson model's error rate is a simple constant-by-size-bucket, not learned from data — an
  obvious extension point.
- Nesting/recursion strategy is orthogonal to genotyping model and mostly lives in
  `GraphCaller`/`FlowCaller`, not `SnarlCaller` — nested-vs-flat calling and
  ratio-vs-Poisson scoring can be varied independently.
- `--top-down` is the actively maintained nested-calling path; `--bottom-up`
  (`NestedFlowCaller`) is explicitly marked deprecated in the code/help text — worth knowing
  before building on it.
- Everything downstream of "list of traversals + a way to score them" is graph-representation
  agnostic (works over `PathHandleGraph`/GBZ uniformly), so most of the machinery would carry
  over even if the site-decomposition or support source changed.
