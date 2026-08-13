# Read-level genotype likelihoods for `vg call`

**Audience:** an engineer comfortable with C++ and bioinformatics who has not worked in the `vg`
codebase. Sections 1–2 explain enough of the existing code to make the rest actionable.

**Status: built and measured.** Stages 0–8 are implemented (PR vgteam/vg#4990) and evaluated on real
data over two chromosomes and two graphs. This document describes **what exists**, not what is
intended; where it still says "we will", read it as a record of the reasoning that led to what was
built. The tier-2 summary at the end is the short version of how it performs.

Where to go for what:

| | |
|---|---|
| the model, the architecture, what is emitted | here — §3–§5 |
| how `vg call` worked before this | [vg-call-characterization.md](vg-call-characterization.md) |
| the source reading that constrained the design, prior art, settled decisions and why | [read-likelihood-genotyping-plan.md](read-likelihood-genotyping-plan.md) |
| the evaluation harness, and the log of every measurement | [vg-call-eval-plan.md](vg-call-eval-plan.md) — start at "Where this stands" |

---

## 1. Goal

`vg call` genotypes variants against a pangenome graph. Today it decides genotypes from **aggregate
read depth**: it asks "how many reads cover this allele's nodes and edges?" and fits counts to a
Poisson model. It never looks at an individual read.

We want to replace that with an explicit **`P(reads | genotype)`** model:

1. For each site, build a **reads × alleles matrix** of `ln P(read | allele)` — how well each
   overlapping read is explained by each candidate allele.
2. Enumerate **every allele combination** (all unordered pairs, for a diploid) and score each as a
   proper likelihood, marginalising over which haplotype produced each read.

### Why this is worth doing

- **Allele balance emerges instead of being imposed.** The current caller has an explicit
  `het_bias` knob to stop it over-calling homozygotes. In a mixture model the ½/½ weighting does
  that automatically.
- **Continuous evidence instead of binary support.** A read that matches an allele with one
  low-quality mismatch is currently either counted or not. Here it contributes a graded likelihood,
  with the mismatch charged at its own base quality.
- **Multi-allelic sites work unchanged.** No special-casing; a 5-allele site is 15 genotypes.
- **Mismapping becomes part of the model**, not a filter — a read that fits nothing is explained as
  "probably from elsewhere" rather than either forcing a bad call or being silently dropped.
- **Calibrated `GL`/`GQ`/`GP` fall out** of the likelihoods rather than being derived from counts.

### What we are not changing — this is an option, not a replacement

**The depth-based Poisson caller remains the default.** Read-level genotyping is opt-in behind a
`--read-likelihood` flag, and a `vg call` invocation that does not pass it must behave exactly as it
does today. That is a firm design constraint, not an interim state, and it shapes several things:

- Traversal enumeration, snarl decomposition, graph loading and VCF writing are untouched.
- The new caller is *added alongside* `PoissonSupportSnarlCaller` and selected by flag in
  `call_main.cpp`, not substituted for it.
- **Default-path output must stay byte-identical.** This is cheap to guarantee and cheap to test
  (§7, §8), and it is what makes the feature safe to merge early — it cannot regress anyone who does
  not opt in.
- `--read-likelihood` requires a read source (`--gam`/`--gaf`), so option validation must fail early
  and clearly when it is passed without one, rather than producing a silently read-free genotyping run.

Being an option also lowers the bar for the first version usefully: it can ship with the known
limitations in §4.6 documented rather than solved, because nobody gets them unless they ask.

---

## 2. How `vg call` works today

Enough background to place the changes. Concepts first, then the pipeline.

### 2.1 Concepts

**Pangenome graph.** Nodes hold sequence; edges join them. A haplotype is a *path* through the graph.
Variation appears as alternative routes between shared nodes.

**Snarl.** A "bubble" in the graph — a subgraph with exactly two boundary nodes, everything between
them being alternative routes. **A snarl is a site of variation.** Snarls nest (a SNP inside a larger
deletion), forming a decomposition tree managed by `SnarlManager`. This nesting matters later.

**Traversal (`SnarlTraversal`).** One path through a snarl, from one boundary node to the other.
**A traversal is a candidate allele.** The traversal following the reference path is the reference
allele.

**GAM / GAF.** vg's read alignment formats — the graph analogues of BAM/SAM. A GAM record
(`Alignment`, protobuf) holds the read sequence, base qualities, a MAPQ, and a `Path`: the list of
nodes the read traverses, each with `Edit`s recording matches, mismatches, and indels relative to
that node's sequence. **The read's alignment therefore already tells us which route through each
snarl the read took, and where it disagrees with the graph.** That fact is the basis of §5.

**`vg pack`.** A coverage index built by streaming a GAM. It stores per-base coverage, per-edge
coverage, and one *averaged* MAPQ per node. It does **not** store reads, read identities, or base
qualities.

### 2.2 The pipeline

```
graph (+ snarls)  ──►  GraphCaller::call_top_level_snarls()   [parallel over snarls, OpenMP]
                            │
                            ├─►  TraversalFinder::find_traversals(snarl)     ──► candidate alleles
                            │      (needs `vg pack` for flow weights)
                            │
                            ├─►  SnarlCaller::genotype(...)                  ──► the genotype
                            │      (currently: PoissonSupportSnarlCaller, from pack counts)
                            │
                            └─►  VCFOutputCaller::emit_variant(...)          ──► VCF record
                                     └─► SnarlCaller::update_vcf_info(...)      (fills GT/GQ/GL…)
```

Two facts about `SnarlCaller` drive the whole design.

**First: `genotype()` already receives everything we need.**

```cpp
pair<vector<int>, unique_ptr<CallInfo>> genotype(
    const Snarl& snarl,
    const vector<SnarlTraversal>& traversals,   // <-- the candidate alleles
    int ref_trav_idx, int ploidy,
    const string& ref_path_name,
    pair<size_t,size_t> ref_range);
```

So the matrix can be built **inside** a new `SnarlCaller` subclass as a private implementation
detail. **No interface change, and no modification to `GraphCaller`, the traversal finders, or the
snarl decomposition.** This is the single most important structural fact in the plan: it keeps the
change contained.

**Second: `vg call` has no read input at all today.** It takes a graph and a `vg pack` file. Getting
individual reads to the genotyper is genuinely new work, and it is the largest piece of this project
(§6) — not the modelling.

### 2.3 Two constraints discovered in the source

- `FlowCaller`, `NestedFlowCaller` and `LegacyCaller` `dynamic_cast` the caller to
  `SupportBasedSnarlCaller` ([call_main.cpp:896](vg/src/subcommand/call_main.cpp:896),
  [graph_caller.cpp:2060](vg/src/graph_caller.cpp:2060)) to reach `get_support_finder()`, which
  supplies the node/edge weights that `FlowTraversalFinder` uses to *enumerate* alleles. **The new
  caller must therefore subclass `SupportBasedSnarlCaller`** — but `-k` (the pack file) is
  required only for *enumeration*, never for genotyping. With `-g`/`-z`, `GBWTTraversalFinder`
  enumerates alleles from recorded haplotypes rather than from support, so nothing on that path
  consults the pack and it can be omitted entirely; a `NullTraversalSupportFinder` stands in.
  A caller doing that must override `get_skip_allele_fn()`, since the support-based version prunes
  any allele below a support threshold and would discard every allele when support reads zero.
- `VCFGenotyper` takes a generic `SnarlCaller&`
  ([graph_caller.hpp:259](vg/src/graph_caller.hpp:259)), so `vg call -v` (re-genotyping an existing
  VCF) picks up read-level likelihoods **for free**, with no extra work.

---

## 3. Architecture

Three new components plus a read-source abstraction. **Everything below is reached only via
`--read-likelihood`;** without the flag, `call_main` constructs `PoissonSupportSnarlCaller` exactly as
it does today and none of this code runs.

```
                          ┌─────────────────────────────────┐
   GAM / GAF / .gai       │  SiteReadSource                 │   reads overlapping
   / gaf-base       ─────►│  (pluggable backend, §6)        │──►  this snarl
                          └─────────────────────────────────┘        │
                                                                     ▼
                          ┌──────────────────────────────────────────────────┐
                          │  AlleleLikelihoodCalculator  (§5)                │
                          │  score each read against each allele from the    │
                          │  read's existing graph alignment — no DP         │
                          └──────────────────────────────────────────────────┘
                                                                     │
                                                                     ▼
                          ┌──────────────────────────────────────────────────┐
                          │  AlleleReadLikelihoods                           │
                          │  flat reads × alleles matrix of ln P(read|allele)│
                          └──────────────────────────────────────────────────┘
                                                                     │
                                                                     ▼
                          ┌──────────────────────────────────────────────────┐
                          │  ReadLikelihoodSnarlCaller : SupportBasedSnarlCaller │
                          │  enumerate all genotypes, score, emit GL/GQ/GP   │
                          └──────────────────────────────────────────────────┘
```

New files: `src/allele_likelihood.{hpp,cpp}` for the matrix and calculator;
`ReadLikelihoodSnarlCaller` added to the existing `src/snarl_caller.{hpp,cpp}`;
`src/site_read_source.{hpp,cpp}` for the read backends.

### 3.1 `AlleleReadLikelihoods` — the matrix

```cpp
/// Per-read relative likelihoods over the alleles at one site.
/// Rows are NORMALISED (§4.2): each is divided by that read's best entry, so every value is
/// in [0,1] and each row's max is exactly 1. The per-row divisor is kept in best_ll().
class AlleleReadLikelihoods {
public:
    size_t num_reads()   const;
    size_t num_alleles() const;
    /// Relative likelihood in [0,1]; 0 means "no valid placement" on this allele (§7).
    double rel(size_t r, size_t a) const;
    /// e_r: mismapping probability from MAPQ, already clamped (§4.3).
    double mismap_prob(size_t r) const;
    /// ln of the per-row divisor: the read's absolute best fit at this site.
    /// NOT used by the model -- diagnostics and the §5.3 escape-hatch trigger only.
    double best_ll(size_t r) const;
private:
    std::vector<double> matrix;                   // row-major, num_reads * num_alleles
    std::vector<double> read_mismap_prob;
    std::vector<double> read_best_ll;
    size_t n_reads = 0, n_alleles = 0;
};
```

Flat row-major, not a map — sites are small and this is the hot path.

### 3.2 `ReadLikelihoodSnarlCaller` — the genotyper

```cpp
class ReadLikelihoodSnarlCaller : public SupportBasedSnarlCaller {
public:
    ReadLikelihoodSnarlCaller(const PathHandleGraph& graph, SnarlManager& snarl_manager,
                              TraversalSupportFinder& support_finder,
                              AlleleLikelihoodCalculator& lik_calc);

    struct ReadLikelihoodCallInfo : public SnarlCaller::CallInfo {
        double gq = 0, posterior = 0;
        size_t n_reads = 0, n_informative = 0;
        std::vector<std::pair<std::vector<int>, double>> genotype_lls;   // for GL
    };

    pair<vector<int>, unique_ptr<CallInfo>> genotype(...) override;
    void update_vcf_info(...) override;
    void update_vcf_header(string& header) const override;
private:
    AlleleLikelihoodCalculator& lik_calc;
    // Circuit breakers for pathological snarls, not tractability requirements (§5.4).
    // Default off; log whenever one fires.
    size_t max_alleles_scored = 0;   // 0 = unlimited
    size_t max_reads_per_site = 0;   // 0 = unlimited
};
```

---

## 4. The model

### 4.0 The objective as built — superseding the derivation below

**Read this first.** §4.1–4.3 derive the model as originally specified and are kept because the
reasoning behind each piece is still the reasoning; §4.3a and this section record what the pieces
became. Where they disagree, this section is the code. Three things here are not in the derivation
at all: the mixture weights are length-aware, the observation is a *fractional* read count, and there
is a depth factor in `GQ`.

At each site `s`, with alleles `A = {0,…,K−1}`, reads `R(s)` — the rows of the site's likelihood
matrix — and ploidy `k`:

```
Ĝ(s) = argmax over G ∈ 𝒢(K,k) of  ℒ(G)

ℒ(G) = w_d · ln Poisson(N_eff ; λ_G)
     + Σ_{r ∈ R(s)} ln[ (1 − e_r) · Σ_{i=1..k} ω_i(G) · rel(r, h_i) + e_r ]
```

`𝒢(K,k)` is every multiset of size `k` from `A`, enumerated exhaustively. `R(s) = ∅` gives a
**no-call**, not hom-ref: absence of reads is absence of evidence to a read-conditioned model.

**The per-read term.** `rel(r,a) = exp( ℓ(r,a) − max_{a'} ℓ(r,a') )`, where `ℓ` is the edit-aware
alignment log-likelihood (quality-adjusted where base qualities exist) and `−∞` where the read
cannot be placed. The divisor runs over *all* alleles at the site, not those in `G`, which is what
makes it genotype-independent and lets it cancel from every comparison — and why the background
term below is exactly 1.

`e_r = clamp( phred_to_prob(MAPQ_r), e_min, e_max )`, with `e_min = 0.02`, `e_max = 0.7`. The
bracket therefore lies in `[e_r, 1]`: the log is always finite and no single read can penalise a
genotype without bound.

**The mixture weights** (§4.3a, sharpened). With `unique(a,b)` the length of nodes in allele `a`'s
content but not `b`'s, a node visited twice counted once, and `R̄` the mean read length in this
site's matrix:

```
U_i = min over j ≠ i of unique(h_i, h_j)
ω_i = max(U_i + R̄ − 1, 1) / Σ_j max(U_j + R̄ − 1, 1)
```

`U_i` is the sequence on which haplotype `i` can be *told apart* from the others; `U_i + R̄ − 1` the
read start positions that can see it. Equal-length alleles give exactly `1/2`, so SNVs and balanced
indels are unchanged bit for bit. A homozygote gives `unique(a,a) = 0`, hence `ω = (½,½)` and
`Σ ω_i rel = rel(r,a)` — hom and het need no special casing.

**The depth term.**

```
N_eff = Σ_{r ∈ R(s)} (1 − e_r)
λ_G   = c(s) · Σ_{i=1..k} max( T_{h_i} + R̄ − 1, 1 )
c(s)  = Σ_{r ∈ W(s)} (1 − e_r) / ( 2 · Σ_{v ∈ W(s)} |v| )
ln Poisson(n; λ) = n ln λ − λ − ln Γ(n+1)
```

Three details, each of which was got wrong once and is worth stating as a constraint:

- `N_eff` is **not** a read count. A read the mapper places at MAPQ 0 enters the per-read term at
  `1 − e_r` of its weight, so counting it as a whole read of depth contradicts the term one line
  above. `Γ` rather than a factorial follows: `N_eff` is fractional by construction.
- `c(s)` carries **the same weighting**, over the node-ID window the read source already fetched to
  answer this site's own query. Weighting one side only puts a constant factor between `N` and `λ`
  and pushes every ratio the same way — a bias, not a signal. Weighting both makes the correction
  *relative*: it cancels wherever a site's mapping quality matches its neighbourhood's.
- `T_a` is the allele's **interior** traversal length. A `SnarlTraversal` runs from the snarl's start
  visit to its end visit inclusive, but a read lying entirely within a boundary node cannot
  discriminate and is dropped before it becomes a row — so those bases recruit nothing. Counting
  them made `λ` too large by roughly the two anchors' length, a fixed overhead per site, which put
  the median `DR` at 0.59 instead of 1 and left typical sites being scored on the Poisson's steep
  low-count flank. A traversal with no interior is the deletion edge, where `max(·,1)` correctly
  leaves the `R − 1` junction positions.

**Reported quantities. None of these enter the argmax.** With `L₁ ≥ L₂` the two best `ℒ`:

```
share = min( 1, Σ_{a ∈ distinct(Ĝ)} AD_a / |R(s)| )
DR    = N_eff / λ_Ĝ
δ     = exp( −A · |ln DR| )  if A > 0 and max_a |T_a − T_ref| ≥ 50 bp, else 1

GQI  = (10 / ln 10) · (L₁ − L₂)
GQ   = clamp( GQI · share · δ, 0, 256 )
GL   = ℒ(G) / ln 10, per genotype
post = L₁ − ln Σ_G exp( ℒ(G) )
```

`share` and `δ` are ranking discounts and can only lower `GQ` — verified end to end: `GQ` was raised
at 0 of 780,356 records across four datasets. `GQI` keeps the raw ratio for anything needing a
posterior. `DR` is emitted whether or not `w_d > 0`, which is deliberate: the observable was measured
as a ranking signal before the model was allowed to act on it, the same order the share discount was
established in.

**Parameters, and where each number comes from.**

| symbol | flag | default | basis |
|---|---|---|---|
| `e_min` | `--mismap-min` | 0.02 | swept; interior optimum, and its old job is now the depth term's |
| `e_max` | `--mismap-max` | 0.7 | swept 0.5→0.7; 0.5 clamps every read at MAPQ ≤ 3 to a coin flip |
| `ω_i` | `--flat-mixture`, `--length-weight-whole-traversal` | unique content | §4.3a |
| `w_d` | `--depth-term` | 0.1 | 4-hap graphs prefer 0.25, 34-hap prefer 0.1; unimodal, turns over above 0.25 |
| `N_eff` | `--depth-count-raw` | effective | genotype-neutral; chosen for `DR`'s ranking power |
| `c(s)` window | none — source's own fetch window | 4096 in practice | saturated by 4096; a fixed fallback covers sources without one |
| `A` | `--depth-quality` | 0 | 0.5 when armed; 7 of 8 cells improve, one does not |

**Two invariants worth keeping in mind when changing any of this.** The weights sum to 1, so adding
an allele still costs and a clean homozygote still beats the heterozygote. And a global scalar on the
read term cannot change a genotype — `argmax_G w·ℒ(G) = argmax_G ℒ(G)` — which is why `--read-weight`
was removed rather than fitted, and why any effective-sample-size correction has to enter the
*per-read* term, as `e_r` and the mixture weights do.

### 4.1 Genotype likelihood

A genotype `G` is a multiset of allele indices with `|G| == ploidy`. Assuming each read is drawn from
one of `G`'s haplotypes uniformly at random, and reads are independent:

```
ln P(reads | G) = Σ_r  ln [ Σ_{h ∈ G} (1/|G|) · P(read r | allele h) ]
```

The inner sum is the marginalisation over which haplotype produced the read. For a homozygous
genotype the two ½ weights collapse to 1, so hom and het are handled by the same expression with no
special case. That is where allele balance comes from for free.

**On "reads are independent".** Mates overlapping the same site are not, and a read whose mate anchors
elsewhere carries placement evidence its own MAPQ does not reflect. Policy for the first version:
count a fragment once when both ends overlap the site, and state the assumption at the point in the
code where the product is taken. The broader consequence of assuming independence — over-confident
`GQ` — is §4.6 risk 2.

§4.2 then normalises the `P(read | allele)` values, and §4.3 adds a term for the read not having come
from this site at all. The result is a single expression evaluable in plain arithmetic.

### 4.2 Normalising each read's row

`score_to_unnormalized_likelihood_ln(score)` is just `log_base * score`
([mapping_quality_calculator.cpp:415](vg/src/mapping_quality_calculator.cpp:415)). So what §5 produces
is `ln P(read | allele)` **up to a per-read constant**: the true probability is
`(1/Z_r) · exp(log_base · score)`, and `Z_r` is neither known nor computable. It cancels in §4.1
because it multiplies every genotype equally — but that only holds as long as we do nothing but
compare genotypes.

Rather than rely on that, **normalise each row of the matrix up front** by its own best entry. With
`ell_h = log_base · score(r, h)` and `ell*(r) = max over ALL alleles at the site`:

```
lambda_h = ell_h − ell*(r)   ≤ 0            L_h = exp(lambda_h)  ∈ [0, 1]
```

Each read's likelihoods are now expressed **relative to that read's best available explanation**: the
best-fitting allele has `L = 1`, everything else is a fraction of it. `Σ_r ell*(r)` is
genotype-independent, so dropping it changes no comparison — and every downstream expression works
with proper numbers in `[0, 1]`.

Three things this buys, beyond readability:

- **`ell_h = -inf` ("no valid placement", §7) becomes `L_h = 0`** — no special case anywhere.
- **Underflow is harmless.** A hopeless allele gives `L_h = 0` rather than a large negative number
  needing careful handling; that is the correct answer, not a loss of precision.
- **It removes `logsumexp` from the hot path**, once combined with §4.3.

**Do not normalise the other way.** Making the alleles sum to 1 —
`exp(ell_h) / Σ_{h'} exp(ell_{h'})` — is a plausible-sounding mistake that must be avoided. That
quantity is a *posterior over alleles* under a uniform prior, not `P(read | allele)`. Substituting it
double-counts the prior, and it **destroys the absolute-fit information that §4.3 depends on**: a read
fitting one allele badly and all others worse becomes indistinguishable from a read fitting one allele
perfectly. Absolute fit is the evidence of mismapping.

Keep `ell*(r)` alongside the normalised row (`best_ll(r)` in §3.1). The model no longer needs it, but
it is the only remaining record of absolute fit — used by the §5.3 escape-hatch trigger ("this read
fits everything badly") and worth having in `--dump-likelihoods`.

Finally, these values remain valid for argmax, for `GQ` (a likelihood *difference*), and for a
posterior over genotypes, but **not** as calibrated absolute probabilities. Normalisation makes them
comparable within a read; it does not make them absolute.

### 4.3 Mismapping, folded into the same mixture

A read may not come from this site at all, and MAPQ estimates that. Rather than filtering, make
"came from somewhere else" **one more component of the same mixture**. With normalised rows this is
the textbook two-component form:

```
ln P(reads | G)  =  Σ_r  ln [ (1 − e_r) · Σ_{h∈G} (1/|G|) · L_h(r)   +   e_r · 1 ]
```

where `e_r` is the read's mismapping probability from MAPQ:

- `e_r = phred_to_prob(mapq_r)` — [statistics.hpp:215](vg/src/statistics.hpp:215)
- (log-space equivalents, if wanted: `phred_to_logprob` /
  `logprob_invert`, [statistics.hpp:225](vg/src/statistics.hpp:225) and
  [:207](vg/src/statistics.hpp:207))

**The background component is exactly 1**, which is what normalisation bought us. It reads directly:
*with probability `e_r` the read came from elsewhere, where it fits about as well as its best
explanation here.* That bounds how much any single read can penalise a genotype.

### 4.3a The mixture weights are not flat — superseding the `1/|G|` above

> Sharpened again since: the weight counts sequence *unique* to each allele rather than whole
> traversal length, and §4.0 states the shipped form. Whole traversal length is still available as
> `--length-weight-whole-traversal`, and one thing found later explains part of why it underperforms:
> a traversal includes the snarl's two boundary anchors, which every allele shares, so counting them
> inflates both alleles alike and pulls the weights back toward flat.

**`1/|G|` is wrong wherever the alleles differ in length, and the formula above is kept only because
section numbers are stable.** The shipped model uses

```
ln P(reads | G)  =  Σ_r  ln [ (1 − e_r) · Σ_{h∈G} w_h · L_h(r)   +   e_r · 1 ]

        U_h + R − 1
w_h = ─────────────────────
      Σ_{h'∈G} (U_{h'} + R − 1)
```

`w_h` is the mixture weight the flat version assumed away: the probability that a read observed *at
this site* came from haplotype `h`. A flat weight asserts each haplotype contributed half the reads,
which is false over an interval where one haplotype carries a deletion — that haplotype contributes
none. `U_h` is the sequence `h` visits that the genotype's other allele does not, and `R` is the mean
read length in the site's own matrix, so `U_h + R − 1` counts the read start positions that can yield
a read able to tell `h` apart from its partner.

Why unique content rather than whole traversal length: a traversal includes the site's shared
sequence, and reads landing there fit every allele equally, contribute the same factor to every
genotype and cancel. Counting them over-states the shorter allele — at one measured 2648 bp deletion
the traversals are 296 and 2945 bp, a ratio of 6.9, where the reads split about 14.6.

What it cost to get wrong, measured on HG002 chr6 against the GIAB draft benchmark: 94% of
heterozygous deletions above 1 kb were lost outright, and two thirds of heterozygous insertions above
1 kb were called homozygous. The insertion half was invisible to recall for a year because truvari
matches on locus, size and sequence but not genotype.

Three properties make this safe to have on by default:

- **equal-length alleles give exactly ½**, so every SNV and balanced indel is bit-for-bit unchanged —
  small-variant genotype F1 does not move to four decimal places on either graph tested;
- **the weights sum to 1**, so adding an allele still costs and a clean homozygote still wins. This is
  what a plain `max_{h∈G}` destroys: `max` is monotone in the allele set, so a heterozygote can never
  score below either homozygote and with any noise it always wins;
- **it is symmetric** in the direction of the imbalance, so one rule covers insertions and deletions.

`--flat-mixture` restores the pre-correction model exactly. Full measurements in the harness repo's
[docs/tier2-sv-errors.md](../docs/tier2-sv-errors.md).

**This does not remove the need for a depth term (§5.3.3).** It corrects the *relative* weight between
a genotype's haplotypes; it says nothing about whether the absolute number of reads at a site is
plausible, which is what the collapsed-repeat pile-ups need and what §4.4 below decided against.

Two properties worth understanding before touching this code:

**The background magnitude is not a free parameter.** A background of `c` with mismapping probability
`e` is nearly equivalent to a background of 1 with a rescaled `e`. So `c` and `e_r` are **not
separately identifiable** — fix the background at 1 and let the (clamped) MAPQ be the *only* knob.
Adding a tunable `P_bg` alongside `e_r` would look meaningful and do nothing independent.

**The whole expression is numerically bounded, so no `logsumexp` is needed.** Since
`L_h ∈ [0,1]`, the mixture lies in `[e_r, 1]` and its log lies in `[ln e_r, 0]` — **always finite**.
Compute it with plain arithmetic and one `log()`. No `add_log`, no stability care, no `-inf`, no
`NaN` (given the guard in §7). This is the concrete implementation payoff of normalising first.

#### Behaviour at the limits

| Situation | Result |
|---|---|
| high-MAPQ read fitting *a*, `G = {a,b}` | `≈ ln(½)` — full evidence, the familiar het factor of 2 |
| high-MAPQ read fitting *a*, `G = {a,a}` | `≈ ln(1) = 0` |
| high-MAPQ read fitting *a*, `G = {b,b}`, `L_b ≈ 0` | `≈ ln e_r`; ~13.8 nats against `{b,b}` at MAPQ 60 — strong but **bounded**, where without the background it would be unboundedly negative |
| MAPQ ≈ 3 (`e_r ≈ ½`) | `ln e_r ≈ −0.7`; the read barely discriminates. Correct shrinkage |
| read with **no valid placement** on *h* | `L_h = 0`, and the `+ e_r` term keeps the log finite |

#### Why MAPQ is the right signal here

One might worry that MAPQ measures the wrong ambiguity: it comes from the score gap to the next-best
alignment, and in a graph one could imagine that competitor being *another allele of the same snarl*
rather than another locus — which would mean a read spanning a heterozygous SNP got a low MAPQ for
being maximally informative, and `e_r` would then shrink exactly the reads we most want to hear from.

**vg's mappers are built so that this does not happen.** The set of alignments feeding the MAPQ
calculation is deliberately constructed to hold *distinct graph placements*, not competing alignments
of the same locus:

- Alignments are generated **per seed cluster / gapless-extension group** — one group per graph
  locus — and it is that per-locus set that becomes the score vector MAPQ is computed from
  ([minimizer_mapper.cpp:940-1005](vg/src/minimizer_mapper.cpp:940), and the score vector assembled
  at [:1086](vg/src/minimizer_mapper.cpp:1086)). Alternative traversals of one snarl fall in the same
  cluster and do not enter as separate competitors.
- The reported set is explicitly deduplicated
  (`Mapper::score_sort_and_deduplicate_alignments`, [mapper.cpp:4418](vg/src/mapper.cpp:4418)).
- Read-space overlap between primary and secondaries is measured explicitly rather than ignored —
  `sub_overlaps_of_first_aln` ([mapper.cpp:3597](vg/src/mapper.cpp:3597)) counts secondaries
  overlapping the primary beyond a fraction threshold, and that count feeds the MAPQ computation as a
  deliberate penalty term ([mapping_quality_calculator.cpp:250](vg/src/mapping_quality_calculator.cpp:250)).

So the score gap MAPQ is derived from reflects *"this read could be at a different locus"* — which is
precisely the quantity `e_r` is supposed to be. **No correction is needed, and informative reads at
heterozygous sites are not penalised for being informative.**

Two things to keep anyway, on their own merits rather than as a hedge:

- **Clamp `e_r` to `[e_min, e_max]`.** Ordinary numerical hygiene: it guarantees the mixture stays
  strictly inside `(0, 1)` so the `log()` is always finite, and it stops any single read dominating a
  site. The clamp values need choosing regardless.
- **Calibrate against simulated truth, and keep the term switchable off.** `vg sim -a` gives reads with
  known origin, so the MAPQ-to-mismapping relationship can be measured rather than assumed — which is
  how the clamp values get set, and how we confirm the term earns its place (§9).

For context if this is ever revisited: `MappingQualityCalculator::compute_group_mapping_quality`
([mapping_quality_calculator.hpp:77](vg/src/mapping_quality_calculator.hpp:77)) computes MAPQ over a
set of alignments treated as equivalent, should a future caller want explicitly group-level placement
confidence.

### 4.4 Calibration: assume it, and know the escape route

**Start simple.** Take MAPQ at face value:

```cpp
e_r = clamp(phred_to_prob(aln.mapping_quality()), e_min, e_max);   // that is the whole thing
```

No recalibration table, no reweighting. Two reasons this is a reasonable starting assumption rather
than wishful thinking:

- **vg's MAPQ is deliberately biased toward under-confidence.** Giraffe takes a `min` of the score-gap
  estimate and a separate cap derived from how much of the search space was actually explored, then
  clamps to 60:
  ```cpp
  mapq = round(min(mapq_explored_cap, min(mapq, 60.0)));   // minimizer_mapper.cpp:1171
  ```
  ([minimizer_mapper.cpp:1163-1182](vg/src/minimizer_mapper.cpp:1163); `mapq_explored_cap` is also set
  as an annotation on the alignment, so it is available downstream if ever wanted.) So `e_r` most
  likely *over*-states mismapping — erring toward shrinking a read's influence, which costs power
  rather than producing false calls.
- **`e_r` is genotype-independent**, so even a badly wrong `e_r` cannot favour one genotype over
  another. It can only compress a read's contribution. That is the benign direction for an error.

#### One part of this is not optional

**The `e_max` clamp is load-bearing, not hygiene.** Many mappers use MAPQ 0 to mean "multi-mapping"
rather than `P(wrong) = 1`, and `phred_to_prob(0) = 1` would set `e_r = 1`, collapsing the read's term
to `ln(1) = 0` for every genotype — the read silently contributes nothing at all. That may even be the
behaviour we want for a MAPQ-0 read, but it must be a deliberate clamped decision rather than an
accident of the phred conversion. Set `e_max` well below 1 and make the choice explicit.

Also worth having from the start because it is nearly free and makes the assumption testable: a flag to
**switch the mismapping term off entirely**, so its contribution can be measured rather than assumed.

#### The escape route, if the assumption fails

Recorded so nobody has to rediscover it. Miscalibration could enter in three places:

1. **MAPQ → `e_r`.** Fix: a 61-entry lookup table — MAPQ is an integer clamped to `[0,60]`, so a table
   is exact and assumes no functional form. It is a *measurement*, not a fit: `vg sim -a` emits reads
   annotated with their true origin, so the mismapping rate per MAPQ value can simply be counted. Note
   such a table is specific to the mapper, graph, read length and error profile — one fitted on
   Giraffe output will not transfer to a non-vg mapper's GAF.
2. **Score → likelihood.** Needs no new parameter: already expressible by overriding the scorer's
   `log_base`.
3. **The independence product over reads.** Fix: a scalar `read_weight ∈ (0,1]` multiplying the whole
   per-read sum — an effective-sample-size discount. Fitted against *outputs*, not inputs: bin sites by
   reported `GQ`, measure the actual error rate against the truth VCF, and adjust until the reliability
   curve is diagonal.

**None of these changes the shape of the model**, so all three are purely additive later — which is why
it is safe to defer them. (3) is the one most likely to be needed; it is recorded as an accepted risk
in §4.6.

#### The cap is not a spare knob — it is graph-dependent. **Default now 0.5, raised from 0.1**

`max_mismap_prob` was measured as inert in §9.14 of the harness plan, on a graph with 4
haplotypes. That measurement was right and the conclusion drawn from it was too narrow. On a
34-haplotype graph the cap becomes the single most important parameter in the model.

The mechanism is specific. Extra haplotypes do not mainly produce *unmappable* reads — MAPQ 0 reads get
slightly rarer. They produce **two-way ties**: reads that match the graph better than before (identity
0.921 → 0.965 at the failing sites) but cannot be placed among near-identical haplotypes. At the sites
that produce false calls, **MAPQ 1 alone is 23.3% of reads**, and MAPQ 1 means p(wrong) = 0.79. The cap
tells the model 0.1. Across the reads where it binds it **discards 8.1x of the mapper's stated doubt**,
and the model duly calls a heterozygote on evidence the mapper had already flagged as unreliable.

Raising the cap to 0.5 removes **94% of the excess false-positive SNVs** on that graph (1,597 → 443,
SNV precision 0.9776 → 0.9937) and is neutral on the 4-haplotype one (375 → 376 false SNVs). **This is
now the default.** The effect saturates above 0.5, so the value is not delicate; but the cap must stay
strictly below 1, because at `e_r = 1` the per-read term is `log(1) = 0` for every genotype and the read
silently contributes nothing at all. Plan §9.20 has the numbers.

Two consequences for this section. First, **the clamps are properties of the graph-and-mapper pair, not
constants**: the same code, the same reads and the same sample give opposite answers about which clamp
matters, purely because the graph changed what MAPQ reports. Second, they are **not interchangeable** —
the cap governs *placement* ambiguity and shows up in SNVs; the floor governs how hard one read may veto
an allele and shows up in indels. Tuning either against an aggregate F1 will hide what the other is
doing.

**And they interact, so a value fitted under one is not evidence about the other.** Re-sweeping the
floor once the cap was corrected moved its default from 0.01 to **0.02** (plan §9.21). At cap 0.1,
raising the floor cost about four times as much SNV precision as it does at cap 0.5, which is why §9.15
rejected it. With both defaults corrected the 34-haplotype graph gains **+0.0087 small-variant GT F1 and
+0.0006 SV F1** — neither clamp was wrong on the graph it was fitted to, and each was hiding part of the
other's effect.

**The floor is also where small variants and SVs stop agreeing.** Small-variant GT F1 keeps rising to
0.05 on the 34-haplotype graph while SV F1 peaks at 0.02 and falls away; 0.02 is the only point that
improves both against 0.01, which is why it is the default. Anything above it is a real trade, and a run
that cares only about small variants should set `--mismap-min 0.05`. That the two size classes prefer
different values of a clamp meant to bound one read's veto is itself unexplained, and worth a look.

### 4.5 Genotype enumeration and output fields

`K(K+1)/2` genotypes for `K` alleles, diploid — 15 genotypes at a 5-allele site, each costing `R`
two-term evaluations. This is cheap enough to enumerate exhaustively, so unlike the Poisson caller
(which prunes candidates with `top_k`/`top_m`) we do not need candidate pruning at all.

**The prior over genotypes is uniform.** So the posterior is proportional to the likelihood, and:

- `GL` — the full vector of `ln P(reads | G)`, one per genotype.
- `GQ` — the log-likelihood gap between best and second-best genotype.
- `GP` — normalised posterior, which under a uniform prior is just `GL` renormalised over the genotype
  set: `GP = LL_best − logsumexp(all LL)`.

**Do not copy the Poisson caller's posterior formula.** It computes
`posterior = LL_best − ln N − logsumexp(all LL)` ([snarl_caller.cpp:629](vg/src/snarl_caller.cpp:629)),
with `N` the candidate count. Under a uniform prior the `ln N` cancels analytically and does not belong
there. Nearly harmless at a fixed candidate count; not here, because exhaustive enumeration makes
`N = K(K+1)/2` rather than a small pruned set, so the term would deflate `GP` and make it **vary with
allele count** — defeating the comparability it was presumably meant to preserve. (That loop also uses
`0` as the "unset" sentinel for a quantity that is a log-probability; use `-inf` in new code.)

A uniform prior is the deliberate simple choice: it matches the current Poisson caller's `GP`
semantics, and it keeps likelihood and prior **separable**, so the likelihood model can be evaluated
without a prior confounding the comparison. An HWE or het-biased prior is a later refinement if wanted
— it multiplies in without touching anything above.

Note this is also why §4.1 needs no `het_bias` parameter: with a uniform prior, the allele-balance
behaviour comes entirely from the mixture weights rather than from a tuned correction. Those weights
are no longer flat — see §4.3a — which sharpens the point rather than weakening it: allele balance now
follows from the site's geometry instead of from either a tuned bias or an assumption of symmetry.

#### What is actually emitted, and one departure from the above

`DP`, `GL`, `GQ`, `GP`, plus three fields added after tier 2 measured what the model was throwing
away:

- `AD` — reads whose best-fitting allele is each emitted allele, ties split fractionally rather than
  awarded to the lowest index. **It does not sum to `DP`**, and at a busy site falls well below it,
  because only alleles that reached the record get a column while the genotyper scored every allele
  the site offered. That shortfall is the point: it is the share of reads the called genotype fails
  to explain.
- `BL` — mean over reads of `best_ln`, the row divisor of §4.2. It measures whether reads fit
  *anything* here, where `GQ` measures only the gap between the top two genotypes, so the two are
  nearly independent (r = +0.18).
- `GQI` — `GQ` from the likelihood ratio alone.

**The departure:** `GQ` is the likelihood-ratio quality **scaled by `sum(AD)/DP`**. The reason is
structural rather than empirical. Each read contributes to `GQ` only the *difference* its term makes
between the best genotype and the runner-up, so a read whose best-fitting allele is in neither has
both terms collapse to ≈ `e_r` and drops out entirely. `GQ` is therefore blind to whether the called
genotype explains the pile-up at all — a site where a third of the reads prefer an uncalled allele
scores the same as one where none do. Scaling restores them.

Measured over two chromosomes, two graphs and both benchmarks, this improved the ranking in all eight
combinations and at fifteen of sixteen operating points; the linear form was the only one of four
tried that made none of them worse. The cost is that `GQ` is now a quality score rather than a
calibrated posterior — anything needing the posterior must read `GQI`, which is emitted
unconditionally for exactly that reason. `--no-share-quality` restores the previous behaviour.

### 4.6 Accepted limitations and risks

Three things we are knowingly accepting in the first version. Each is recorded with the signal that
would tell us to act on it.

**1. The mixture is depth-agnostic — decided: pure `P(reads | G)`, no coverage term.** It uses only the
*relative* fit of reads to alleles, never the absolute read count against expected coverage. So
deletions are handled fine *given breakpoint-spanning reads*, but a coverage drop alone is not evidence;
and at a site no read overlaps, the model correctly reports flat likelihood — no evidence — where the
Poisson model would confidently call hom-ref from depth.

Going pure is the right first version because it keeps the two sources of evidence separable: if a
depth term is later added, we will be able to measure what it actually contributes rather than having
baked it in from the start. And since this is an opt-in mode (§1), users who want depth-based behaviour
still have it as the default caller. *Signal to act:* hom-ref or SV recall materially worse than the
Poisson caller at stage 4. *Fix:* add the depth term as an option — additive, since it is a separate
factor in the per-site likelihood.

> **Resolved. The signal fired and the fix is built, as specified.** Heterozygous deletion recall above
> 1 kb came in at 0.21–0.44 against the Poisson caller's 0.36–0.84, which is the "SV recall materially
> worse" condition. `--depth-term W` adds `W · ln Poisson(N ; λ_G)` — additive, optional, still off by
> default — and puts the read model ahead of both Poisson arms on all four datasets. Keeping the two
> sources separable paid off as argued: the mixture-weight defect (§4.3a) and the absolute-depth gap
> were separately diagnosable and separately fixable, and each was measured against a model that did not
> already contain the other. Numbers in `docs/tier2-depth-term.md`.
>
> One refinement the design did not anticipate. `N` is **not** a read count: it is `Σ_r (1 − e_r)`, and
> the local rate is measured the same way. Counting a MAPQ 0 read as a whole read of depth contradicts
> the per-read term, which already discounts that read to `1 − e_r`. Because both sides carry the same
> weighting the correction is relative and cancels wherever a site's mapping quality matches its
> neighbourhood's — it barely moves a genotype, but it takes `DR`'s power to rank false positives above
> true ones from 0.51–0.55 to 0.62–0.64 on every dataset.

**2. `GQ` will be over-confident, and this is the most likely of the three to bite.** §4.1 multiplies
per-read likelihoods as if reads were independent. They are not: reads at a site share PCR duplicates,
local misalignment, strand artefacts and graph misassembly. The product therefore accumulates
confidence like `R` rather than `√R`, so reported `GQ` will be too high — **and it gets worse with
depth**, which is exactly the regime we most want to be trusted in. This is the standard reason
likelihood-based callers emit implausible `GQ`s; we are not doing anything unusual, but we should not
pretend the numbers are calibrated.

*Signal to act:* `GQ` values that are implausibly high, or a `GQ` reliability curve that is not
diagonal — worth eyeballing the `GQ` distribution at stage 4, which is nearly free. *Fix:* the scalar
`read_weight` discount in §4.4, which is additive and does not change the model's shape.

Until that is done, treat `GQ`/`GL` as **useful for ranking, not as calibrated probabilities** — and say
so wherever they are documented.

**3. MAPQ is taken at face value** (§4.4). Mitigated by vg's MAPQ being conservative by construction and
by `e_r` being genotype-independent, so errors cost power rather than creating bias. *Signal to act:*
the mismap-term-off comparison showing the term is doing something unexpected. *Fix:* the 61-entry
recalibration table.

Recommendation for all three: build pure and simple, measure at stage 4, and add only what the
measurements justify.

---

## 5. Scoring: use the graph's implied alignment, not fresh DP

For each (read, allele) pair, score the read against the allele **from the read's existing alignment
in the graph** — walk the read's path against the allele's path over a fixed window (§7) and
account the differences. No dynamic programming.

### 5.1 Why this is right, not just fast

**In a variation graph, the snarl decomposition already is a multiple alignment.** Two traversals of
the same snarl share its boundary nodes by construction, so they are aligned *to each other* through
the graph topology. "How many differences between this read and this allele" is therefore already
well-defined structurally — we read off an alignment the graph asserts rather than approximating one.
DP could find a higher-scoring alignment corresponding to **no path in the graph**, which is not what
we want: our alleles *are* paths, and the quantity we need is `P(read | this path)`.

A second benefit: counting differences tells us **exactly which read bases mismatch**, so each
mismatch is charged its own base quality — more principled than DP's length-averaged handling. vg
already has the primitives ([alignment_scorer.hpp:95-125](vg/src/alignment_scorer.hpp:95)), all
overridden by `QualAdjAlignmentScorer`:

```cpp
int32_t score_exact_match(string::const_iterator seq_begin, string::const_iterator seq_end,
                          string::const_iterator base_qual_begin) const;   // matched runs
int32_t score_mismatch  (string::const_iterator seq_begin, string::const_iterator seq_end,
                          string::const_iterator base_qual_begin) const;   // mismatched runs
int32_t score_gap(size_t gap_length) const;                                // indels
double  get_log_base() const;                                              // -> likelihood scale
```

**Produce a `path_t`, not a number.** Rather than accumulating those primitives by hand, make the
output of the read-vs-allele walk the read's **implied alignment to the allele**, as a `path_t`, and
score it with the primitive that already exists — `score_partial_alignment`
([alignment_scorer.cpp:361](vg/src/alignment_scorer.cpp:361), `QualAdjAlignmentScorer` override at
[:570](vg/src/alignment_scorer.cpp:570)):

```cpp
int32_t score_partial_alignment(const Alignment&, const HandleGraph&, const path_t&,
                                string::const_iterator seq_begin,
                                bool no_read_end_scoring = false) const;
// ell = get_log_base() * score_partial_alignment(read, graph, implied_path, window_begin, true)
```

Three reasons this is the better shape:

- The §5.3 DP escape hatch also produces a `path_t`, scored by the same function — so the stage-5b A/B
  is *exactly* comparable rather than approximately.
- A `path_t` dumps as GAF, so `--dump-likelihoods` yields inspectable alignments rather than opaque
  numbers. That is worth a great deal during stage 1.
- `no_read_end_scoring` gives explicit control over full-length bonuses, which **must be off**:
  `score_alignment` includes them and they are not part of `P(read | allele)`.

> **As built: not yet done this way.** The implementation accumulates
> `score_exact_match` / `score_mismatch` / `score_gap` directly over the anchored walk rather than
> constructing a `path_t`. That works and is quality-adjusted, but it forfeits all three benefits
> above — in particular `--dump-likelihoods` emits bare numbers, which made diagnosing the two
> scoring bugs in §5.5 slower than it needed to be. Worth revisiting before the DP escape hatch is
> built, since that is where the comparability argument starts to bite.

Either way, use the scorer's own primitives, not a hand-rolled ±1 difference count. §4.3 only holds if
the score is on the same scale as the `log_base` converting it.

### 5.2 Procedure

1. Fetch reads overlapping the snarl (§6).
2. Extract each read's traversal of the snarl from its alignment path.
   `get_traversal_of_snarl` ([genotypekit.hpp:37](vg/src/genotypekit.hpp:37)) does exactly this — it
   filters the read's mappings to those on nodes in the snarl contents. It is ~15 lines and takes a
   legacy `VG&`, so port it to `HandleGraph` rather than reuse it.
3. Determine the read's **scoring window**: the maximum span of the read within the site, taken over
   all alleles, fixed once and shared by every allele (§7).
4. Drop reads with an empty window. **Not** reads that merely touch no interior node: a read
   traversing straight from one boundary to the other uses a deletion edge and is maximally
   informative (§5.5). A read is uninformative only if it touches no interior node *and* moves
   between no two distinct nodes inside the site — i.e. it sits within a single boundary node.
   Precedent for dropping the genuinely uninformative:
   [genotyper.cpp:690-703](vg/src/genotyper.cpp:690).
5. For each allele, walk the read's traversal against the allele's traversal over that window:
   - **where the paths agree** — score the read's own edits from the GAM (matches, mismatches with
     their base qualities, indels). Exact, not an approximation.
   - **where they diverge** — account the divergence as graph events: bases on the allele that the
     read's path skips (and vice versa) become gaps; an equal-length substituted node becomes
     mismatches at the differing bases.
   - **read bases in the window this allele cannot place** — charge as insertions. Never omit them;
     that is the calibration error §7 exists to prevent.
6. `ell = get_log_base() * score`, then normalise the row per §4.2.

Allele sequences still need materialising for base-level comparison — `VCFOutputCaller::trav_string`
([graph_caller.cpp:499](vg/src/graph_caller.cpp:499)) does this, and
`GAFOutputCaller::pad_traversal` ([graph_caller.hpp:235](vg/src/graph_caller.hpp:235)) adds
reference-path flanks so partially-overlapping reads have anchor sequence. Both are **per-allele**,
not per-(read, allele), so cache them per site.

**A scorer has to be built.** `vg call` constructs no `Aligner` or scorer at all today, so stage 1 must
instantiate one — and `log_base` follows from its scoring parameters, so this choice propagates into
every likelihood. One consequence worth deciding up front: `QualAdjAlignmentScorer` is inapplicable to
reads with no base qualities (GAF without a quality column, gaf-base `--no-quality`, §6.4). That needs
a stated fallback scorer and a loud warning, and note that the `e_min`/`e_max` clamps calibrated under
one scorer do not transfer to the other (§4.4).

### 5.3 Keep DP as a bounded escape hatch

The one real cost: scoring from the existing alignment **inherits that alignment**. If the mapper
placed a read wrongly — systematically likelier around indels, SVs and repeats — we score against the
wrong path and cannot recover. Independent DP realignment can.

vg's original authors reached the same conclusion. `Genotyper` ships **both** methods —
`get_affinities` (realignment) and `get_affinities_fast` (string comparison) — switched by a flag
whose comment says it all ([genotyper.hpp:78](vg/src/genotyper.hpp:78)):

```cpp
// Whould we do indel realignment, or should we use fast substring
// affinities for everything?
bool realign_indels = false;
```

Their considered default was the fast path, with realignment reserved for indels. Adopt the same
shape: a `RealigningAlleleLikelihoodCalculator` stays in the design as a **fallback for a bounded
subset** of pairs, triggered by low MAPQ, high edit count, or near-tied likelihoods. Build the cheap
path first; add this only if measurement (§9) shows it matters.

If it is built: `QualAdjAligner` ([aligner.hpp:218](vg/src/aligner.hpp:218)) needs a `HandleGraph`,
so build a single-node `bdsg::HashGraph` per allele. It is move-only and holds per-thread state, so
allocate **one per thread**. `align_global_banded` can throw `NoAlignmentInBandException` and
`BandMatricesTooBigException` — treat a throw as "no placement", not a fatal error. Note that
`SSWAligner`, vg's only linear-sequence aligner, has **no** quality-adjusted variant, which is why
the graph route is necessary.

### 5.3.1 MAPQ is location confidence, not haplotype confidence — and what to do short of realignment

The tier-2 insertion deficit (§9.12–9.14 of the harness plan) is best explained not by mismapping but by
**local misalignment**: MAPQ says the read is in the right *place*, and says nothing about whether its
path through *this site* is the right one. A read whose indel is placed a few bases off, or shifted
within a homopolymer, is still MAPQ 60 — so the mismapping term cannot discount it — yet it votes at full
strength for whichever allele its (wrong) path happens to match. That is consistent with everything
measured: the failing reads are MAPQ 60, the model is confident (median `GQ` 126), and the sites are
tandem-repeat-enriched, which is exactly where indel placement is ambiguous.

Full realignment of every read to every haplotype would fix it and is too slow to be the default. Four
cheaper options, in increasing cost:

**1. Floor the per-read veto (`--mismap-min`). Implemented.** The per-read term is
`log((1-e_r)·mixture + e_r)`, so a read fitting allele A perfectly and B not at all costs B exactly
`ln(e_r)`. At MAPQ 60 that is **−13.8 nats from one read**. Flooring `e_r` caps it: 0.01 gives −4.6, 0.05
gives −3.0. This does not detect misalignment, it bounds the damage any single read can do — and it
reinterprets the term as *P(this read's evidence here is unreliable)*, of which mismapping is one cause
and misalignment another. Note the *floor* governs the 90% of reads at MAPQ 60; the *cap*
(`--mismap-max`) only ever touched the 6.3% at MAPQ ≤ 9, which is why sweeping it was inert.

**2. Downweight reads that fit nothing well, using `best_ln`.** Free, and the infrastructure is already
there: `best_ln_likelihood()` retains each read's absolute best fit at the site and is documented as
kept "as a realignment trigger". A locally misaligned read has a characteristic signature — it fits its
*best* allele poorly, because no allele explains a misplaced indel. Reads whose absolute fit is far below
what a clean match would score are exactly the suspect ones, and can be discounted per-read rather than
globally. Unlike a uniform weight this *can* change calls, because it is non-uniform across reads.

**3. Score only the discriminating columns.** The likelihood *ratio* between two alleles depends only on
where they differ; scoring the shared flanks injects noise, and with a fixed inherited alignment those
flanks can score differently against each allele purely from coordinate shifts. Restricting the
comparison to the positions where alleles actually diverge removes that artifact at source. Cheap — it
needs the alleles' divergence set, not any read realignment.

**4. Bounded realignment, triggered rather than universal.** **Tried, twice, and abandoned — see
§9.28 of [vg-call-eval-plan.md](vg-call-eval-plan.md).** A bounded shift and full BiWFA realignment
were both built and reverted. The read-off is wrong constantly, but `AlleleReadLikelihoods`
normalises each row by its own maximum, so any improvement common across a site's alleles is divided
out; optimal alignment changes which allele a read prefers 40 times in 91,914. The paragraph below is
left as written because its reasoning about *where* to spend effort was sound — it was the premise,
that better alignment yields better genotypes, that failed.

§5.3's escape hatch, made affordable by the
concentration measured at tier 2: **the deficit lives in 0.7% of sites**. Realigning only where it
matters — sites whose top-two genotypes are close, or whose alleles differ in length, or where reads
show the poor-absolute-fit signature of (2) — costs ~1% of realigning everywhere. The objection that
realignment is "rather slow" is an objection to doing it *everywhere*; it does not apply to doing it
where the evidence says the inherited alignment is untrustworthy.

Not pursued: using the mapper's own secondary alignments as a direct measure of haplotype-level
ambiguity. It is the most principled signal — it is literally the quantity wanted — but it needs
remapping with secondaries retained and a much larger read store, so it is the expensive option, not the
cheap one.

**(1) was tried first and largely solved it, so (2)-(4) are no longer urgent.** Raising the floor from
1e-8 to 0.01 caps a single read's veto at −4.6 nats instead of −13.8, and on HG002 chr20 it moves 1,493
genotypes — **94% of them het → hom** — improving *every* class: SNV GT F1 0.9759 → 0.9766, insertion
0.7783 → **0.8231**, deletion 0.8231 → **0.8706**, overall 0.9370 → **0.9482**, SV 0.4991 → **0.5120**.
With it, the read-likelihood caller beats the Poisson caller on every class including the two it
previously lost. Full numbers in §9.15 of the harness plan.

**The default has since moved twice, and 0.01 is not it.** Those numbers were measured with
`--mismap-max` still at its original 0.1, which plan §9.20 later showed was itself wrong on
haplotype-rich graphs. Re-swept at the corrected cap of 0.5 and scored on *both* benchmarks, the
floor settled at **0.02** (§4.4, plan §9.21). The mechanism described here is unchanged; only the
value is.

The failure was therefore **spurious heterozygosity**: a few locally misaligned but confidently mapped
reads, each able to veto the homozygous hypothesis almost without bound, forcing a second allele that is
not there. The "preference for longer alleles" first observed was a symptom of that, not the cause.

The apparent residue — an insertion BASEPAIR precision gap that narrows from −0.186 to −0.139 but does
not close — turned out **not** to be a scoring defect at all. It is a benchmark-scope artefact, traced in
§5.3.2. (2), (3) and (4) are therefore not needed for insertion sequence, and remain on the table only as
general refinements.

### 5.3.2 The insertion BASEPAIR gap is a benchmark artefact, not a defect

The GIAB `smvar` truth set contains **no record >=50 bp** — that size class is in the separate `stvar`
benchmark — yet the two confident regions overlap almost completely (58.9 Mb vs 59.4 Mb). So a >=50 bp
insertion called inside the small-variant confident region has every one of its bases scored as a false
positive, however right the call is. It is unscoreable-as-correct by construction.

That is precisely where the gap lives. 246 `readlik-z` calls carry a >=200 bp insertion allele, and they
contribute **27,951 FP bases and zero TP bases** — the whole of the precision difference. The Poisson
caller scores better because it does not emit them; at the two largest sites it emits nothing at all.

Restricting **both** callers to the range the benchmark can adjudicate (dropping any record with a called
allele >=50 bp from REF, applied identically to each) settles it:

| arm | insertion BP recall | precision | F1 |
|---|---|---|---|
| `sm50-poisson-z` | 0.7637 | **0.8700** | 0.8134 |
| `sm50-readlik-z` | **0.8578** | 0.8624 | **0.8601** |

The precision gap collapses from **0.139 to 0.008**, and insertion BASEPAIR F1 flips from a 0.047 loss
into a 0.047 win. Overall BASEPAIR F1 goes 0.9184 → 0.9413 the same way. There is no insertion-sequence
defect in the likelihood model; the unrestricted number was measuring that one caller emits large
insertions and the other does not.

Whether those large calls are *right* is a different question, and `stvar` answers it: they are a net win
(SV insertion recall 0.4976 vs 0.4263), but of the 246, only **35 are confirmed true**, **73 confirmed
false**, and **138 fall outside the SV confident region** and cannot be judged. The genuinely bad ones are
characterised in §5.3.3.

### 5.3.3 What is actually wrong: no depth-plausibility term

`readlik-z` emits a handful of enormous homozygous insertions in and around the chr20 pericentromere at
impossible depths — 61,958 bp at DP 7,873; 57,716 bp at DP 5,337; 33,050 bp at DP 291; 23,450 bp at
DP 932 — all `1/1` at saturated GQ. Chromosome-median DP is 29, and the Poisson caller's expected depth
never exceeds 167 anywhere on chr20. Median DP rises monotonically with called insertion length (28 at
1 bp, 35 at 50–199 bp, **330 at >=1 kb**), so these are collapsed-repeat pile-ups, not haplotypes.

The model cannot reject them, and the reason is structural. It computes P(reads | genotype) **conditioned
on the reads it is handed**, and never asks whether that many reads should be there. The Poisson caller
gets this for free, because observed-versus-expected depth is the whole of its model. A depth-plausibility
guard is the obvious remedy, and expected depth is already reachable: `ReadLikelihoodSnarlCaller`
subclasses `SupportBasedSnarlCaller` and holds a `TraversalSupportFinder` for allele enumeration.

Filtering on depth is **not** that remedy, and the measurement says so. Dropping every call above DP 200
removes 195 records — including all the giants above — and moves insertion BASEPAIR precision by 0.0001
(0.6226 → 0.6227): no benchmark charges for them. Dropping above DP 58 does help (+0.087) but costs SV
insertion recall 0.4976 → 0.4167, because it is a blunt proxy for length that discards real SVs. These
calls should be fixed because they are wrong, not because they cost a score.

**The signature generalises, and it needs a second condition.** Measured across chr6 and chr20 on both
graphs (plan §9.26): on the 34-haplotype graph a false structural call sits at a site with a median
depth ~30% *above* a true one — 53.5 against 41.0 on chr6, 51.0 against 40.0 on chr20 — while only
~63% of its reads fit any called allele, against 100% for true calls. On the 4-haplotype graph neither
difference exists (41.0 vs 38.0, and 1.000 vs 1.000). So the pile-ups are not a chr20 pericentromere
curiosity; they are what a richer graph exposes, and depth alone does not identify them. The
explained-read fraction — now emitted as `AD`, whose shortfall against `DP` is the whole point of the
field — is the second condition, and the pair is far more specific than either alone.

**But the remedy is a size-conditional ranking term, not a guard.** Scored as a rule
(`DP > 1.3 ×` local median **and** `sum(AD)/DP < 0.8`), the pair flags a population that is 71–78%
false on the rich graph and lifts SV precision 0.4170 → 0.5061 on chr6 and 0.3986 → 0.5057 on chr20,
at a cost of 0.2% of small-variant true positives. That looks decisive and is not: the rule reaches
~0.72 recall, and simply thresholding GQ to the same recall gives about the same precision. It beats
*no* filter, not the ranking already available. Folding depth into the score does beat it — precision
at 70% SV recall goes 0.5564 → 0.6325 on chr6 34-hap.

**A depth threshold is not the answer either, and that is now measured rather than assumed.**
Sweeping a two-sided cut on DP over a rolling local median (plan §9.27), against the test any
hard filter has to pass — beat lowering the GQ threshold to the same recall:

- A **minimum** is not worth having. Every setting that removes a material number of calls
  loses F1 on both benchmarks, and GQ thresholding reaches higher precision at the same recall
  in all eight cells. Low depth already depresses GQ, because few reads means a small
  likelihood gap; a separate cut has nothing to add.
- A **maximum** survives only narrowly: at 5× the local median it beats GQ thresholding by
  about +0.025 SV precision on both 34-haplotype runs, and is dominated everywhere else — both
  4-haplotype runs at every threshold, and small variants in every cell. Tight cuts are far
  worse than they look: at 2× the local median GQ wins by 0.09 precision on chr6 4-hap SVs.

So the guard stays unimplemented, and `DP` is emitted for anyone who wants to cut on it.

**What did ship is the explained-share discount** (§9.27): `GQ = GQ_ratio × sum(AD)/DP`, with
`GQI` carrying the undiscounted value and `--no-share-quality` restoring the old behaviour.
This addresses the same blindness from the other side — rather than filtering the pile-up
sites, it stops GQ from being unable to see them — and it improved the ranking in all eight
dataset × benchmark combinations. It does **not** address the sign reversal below, which is a
separate signal and remains open.

The obstacle that has to be solved first is a sign reversal. Within one dataset, local depth ratio has
AUC 0.65 against small-variant labels — where *low* depth means false — and 0.37 against SV labels,
where *high* depth means false. A single linear depth term added to the likelihood would improve one
class and damage the other, so whatever is implemented must condition on called-allele size.
`best_ln` carries the mirror-image constraint: it is the strongest small-variant signal (AUC 0.79–0.84
alone) and is worthless to actively harmful for SVs.

### 5.4 Cost

This is the reason the project is feasible. DP realignment of every pair would cost ~60 reads × up to
50 alleles × 150bp vs 200bp ≈ **90M DP cells per snarl**, across millions of snarls — not viable
genome-wide, and it would force aggressive allele pruning and read subsampling just to run.

Walking the window instead costs ~`60 × 50 × 200` ≈ 600k integer operations for the same site,
roughly 10²–10³× less, scaling linearly in window length. Consequently:

- **No allele pruning needed**, so we can genotype against *all* enumerated traversals including the
  long tail a support-based prefilter would discard. This matters for `GBWTTraversalFinder`, which
  has no allele cap at all.
- **No read subsampling needed.**
- **`-k` is required only for enumeration** (§2.3), so with `-g`/`-z` it is not required at all.
- **The bottleneck moves to read retrieval (§6)** — SQLite queries, index scans, GAF parsing — not
  scoring. That is where optimisation effort should go.

Remaining safety valves, genuinely optional: the `max_alleles_scored` / `max_reads_per_site` circuit
breakers (default off, **and log whenever they fire** — a silent cap reads as full coverage when it
is not), plus an exact-match fast path when a read's traversal equals the allele's and the read has
no edits.

---

### 5.5 Two bugs this design walked into, and what they cost

Both were found by testing nested calling, and neither is nesting-specific. They are recorded because
each came from a plausible-sounding sentence in this document, and because of what they say about
which tests are worth writing.

**Reverse-strand reads were scored against the wrong allele.** Alleles all run from the snarl's start
to its end, so they impose a reading direction on the site. A read aligned to the other strand visits
the same nodes with the opposite orientation flag, so it anchored on no allele step, fell through to
the substitution path, and was compared against a different node's sequence without being
reverse-complemented. On a nested SNP fixture this split reads perfectly by strand: 81 forward reads
preferred the correct allele, 69 reverse reads preferred the wrong one. Roughly **half of all reads,
at every site.**

*Fix:* flip the read into the alleles' reading direction before comparing, with the direction decided
by vote over shared nodes so that a node different alleles visit in different orientations cannot
flip a whole read.

**Reads spanning a deletion were discarded as uninformative.** §5.2 step 4 used to say that reads
touching only the boundary nodes "contribute an identical constant to every allele". That is **false**,
and it is the kind of false that produces plausible output: a read going straight from one boundary to
the other uses the deletion edge, touches no interior node, and is the *only* direct evidence the
deletion allele ever gets. On the fixture all 150 deletion-carrying reads were dropped, so the site
saw reference-supporting reads only, called hom-ref, and `emit_variant` discarded the record — the het
deletion the Poisson caller finds simply vanished.

*Fix:* informativeness is edge-aware. **The discriminating signal for a deletion is in the edge, not
the node set.** A read is informative if it touches an interior node *or* moves between two distinct
nodes inside the site.

#### Why the tests did not catch either

This is the part worth internalising. Both bugs sat in the matrix *builder*, and:

- the model tests (§9) use hand-built matrices, so by design they never exercise scoring at all;
- the integration tests assert that `GL` is well-formed and that `GT` is the `GL` argmax — and **a
  uniformly strand-confused matrix satisfies both perfectly**.

Internal consistency cannot detect a consistently wrong matrix. Two things follow. First, scoring
needs its own unit tests on hand-built graphs and alignments, which now exist in
`src/unittest/allele_likelihood_scoring.cpp`: strand symmetry, deletion-edge reads, a read inside one
boundary node still being dropped, and the window invariant. Both regressions were confirmed to fail
those tests before the fixes — a regression test nobody has seen fail is not evidence of anything.
Second, the de novo truth-concordance harness (stage 3b, §8) is not a nicety to schedule after the
model works: it is the only check that would have caught either bug on the first run, and its absence
is why they survived to be found by accident.

After the fixes the caller reproduces the Poisson caller's answer on the nested star-allele fixture,
and on the HGSVC region emits 19 variants — the same count as the Poisson caller, against 26 before.

---

## 6. Read retrieval

The largest piece of new work, because `vg call` has never needed reads. All backends converge on
`vg::Alignment` behind one interface:

```cpp
/// Random-access source of read alignments by graph locality. Must be thread-safe:
/// GraphCaller visits snarls in parallel, in arbitrary order.
class SiteReadSource {
public:
    virtual ~SiteReadSource() = default;
    /// Visit every read touching any of the ranges, each at most once.
    virtual void for_each_read(const std::vector<std::pair<nid_t,nid_t>>& ranges,
                               const std::function<void(const Alignment&)>& iteratee) const = 0;
    /// Convenience wrapper that collects into a vector. Copies, so for tests only.
    std::vector<Alignment> get_reads(const std::vector<std::pair<nid_t,nid_t>>& ranges) const;
    virtual size_t get_read_count() const = 0;
};
```

Reads go to a callback rather than being returned by value: a site is visited many times over a run,
and an in-memory backend should not copy every read every time. `get_reads_in_snarl()` is not
implemented — no backend yet has a native snarl query, so it would be dead code; add it with the
gaf-base backend that motivates it. `InMemorySiteReadSource` also exposes `add()`, so a source can be
assembled without a file, which is what makes the scoring unit-testable at all (§5.5).

Convergence is cheap because the GAF plumbing already exists and is header-only:
`gafkluge::parse_gaf_record` ([gafkluge.hpp:100](vg/deps/libvgio/include/vg/io/gafkluge.hpp:100)) and
`vg::io::gaf_to_alignment` ([alignment_io.hpp:116](vg/deps/libvgio/include/vg/io/alignment_io.hpp:116)).
**Any backend that can produce GAF text gets an `Alignment` for free.** Caveat: that path requires
GAF in node-ID space; named-segment GAF is not yet supported.

**The whole first version, and even the scale path, needs no new dependencies.** This was worth
checking, and the answer is better than expected: three of the four backends are buildable entirely
from code already in the vg tree.

| Backend | Mechanism | New dependency? | When |
|---|---|---|---|
| `InMemorySiteReadSource` | one streaming pass, bucket reads by top-level snarl | **none** | **§6.1 — prototype** |
| `IndexedGamSiteReadSource` | `GAMIndex` (`.gai`) — *already in vg*, one cursor per thread | **none** | §6.2 — first extension |
| `TabixGafSiteReadSource` | `for_each_gaf_record_in_ranges` — *already in vg*, htslib already a dep | **none** | §6.3 — only if wanted |
| `GafBaseSiteReadSource` | gbz-base / gaf-base, via subprocess | **runtime only — an external binary; nothing links** | §6.4 — **done** |

Because everything sits behind `SiteReadSource`, **backend choice never touches the matrix or the
genotyper**, and each backend has an exact equivalence test available: the same reads must produce the
same calls. That is what lets us prototype on the simplest backend and defer every dependency question.

### 6.1 Prototype: in-memory, zero new dependencies

Build this first. One streaming pass over the GAM or GAF using the existing
`vg::io::for_each_parallel` / `gaf_unpaired_for_each_parallel` iterators, assigning each read to the
top-level snarls its path touches, held in a `unordered_map<const Snarl*, vector<Alignment>>`.

Why this is the right prototype:

- **No index, no new dependency, no new file format.** Reads whatever the user already has.
- **Correct by construction** — no over-fetch, no filtering, no cursor lifetime or thread-safety
  questions. If the matrix is wrong, it is the matrix's fault.
- **The test fixture fits comfortably.** `test/call/HGSVC_chr22_17119590_17880307.gam` covers ~760 kb;
  at 30× that is on the order of 10⁵ reads, tens of MB in memory. `vg call` already holds the whole
  graph, so this is not the dominant allocation.
- It stays useful permanently for regional and chunked workflows.

**The honest limit:** whole-genome at 30× is hundreds of millions of reads and will not fit. So this
backend is for prototyping, testing, and regional runs — genome-wide needs §6.2. Log the retained read
count so the limit is visible rather than discovered as an OOM.

### 6.2 First extension: sorted GAM + `.gai` (still no new dependencies)

Indexed read access already exists in vg and is used by nothing analytical.
`GAMIndex = StreamIndex<Alignment>` ([stream_index.hpp:403](vg/src/stream_index.hpp:403)) is a
BAI-like index over node IDs, built by `vg gamsort -i` or `vg index -l`, with exactly the query we
want:

```cpp
void find(cursor_t& cursor, const vector<pair<id_t,id_t>>& ranges,
          const function<void(const Message&)> handle_result,
          bool only_fully_contained = false) const;
```

Concurrent `find()`s are documented thread-safe, with one cursor per thread — the pattern `vg chunk`
uses. Helpers exist too: `sorted_id_ranges` builds the query ranges, `alignment_pieces_within` clips
a read to a node set.

**Important:** the index is an *over-fetch and filter* design — it yields group start offsets, so
each query scans groups until one is out of range. Querying per snarl across millions of snarls would
rescan the same GAM groups repeatedly. Instead keep a **small per-thread cache** of
`(node-ID range → reads)` populated at **top-level-snarl** granularity, and let nested child snarls
filter that cached set in memory. Since `call_top_level_snarls` visits a top-level snarl and its
descendants on one thread, hit rates should be high.

**Measured, and the premise only half holds.** The cache reports its hit rate under
`--progress`: 21% with `--top-down`, but **0% in the default flat mode** and 0% with `-A`.
`RecurseOnFail` visits only top-level snarls, whose node ranges are disjoint, so there is no
locality to exploit and every site is a fresh index scan. The cache therefore earns its keep
for nested calling and is inert otherwise. It was also wrong to call this "the main performance
lever": on the HGSVC fixture the indexed backend is *faster* than in-memory (0.35s vs 0.67s)
with no cache hits at all, because it avoids loading the whole GAM. The real lever is read-set
size, not the cache.

### 6.2.1 Align the visit order with the read windows

**The problem stage 6 left behind.** Memory is bounded, but `compute()` issues one index query per
genotyped site -- 920 on a 400 kb simulation -- and that scales with the number of sites rather than with
the data. Each query also over-fetches, because the index only bounds where a group starts. Visible in
both directions: on the HGSVC fixture (29 queries) the indexed backend *beats* in-memory, 0.35s vs
0.67s, by skipping a 14 MB parse; on the simulation (920 queries) it loses badly, 3.03s vs 0.36s.

**A cache alone is the wrong answer.** The first version of this section proposed keeping the arbitrary
visit order and caching fetches by node-ID window. That only converts *some* queries into hits, depends
on threads happening to revisit windows, and pays for it in memory: the ceiling becomes
`threads x windows_retained x reads_per_window`, which on 16 threads can exceed the in-memory backend
the whole exercise replaced. Mitigating a self-inflicted problem with a tunable is worse than not
creating it.

**Order the traversal instead.** If snarls are visited in node-ID order, grouped into windows, then each
window is fetched **exactly once**, processed, and released. No cache to size, no reliance on accidental
locality, and one window resident per thread rather than several retained speculatively. The memory
question disappears rather than being traded.

```
for each window of node-ID space, in parallel over windows:
    fetch the window's reads once          # one index query
    for each top-level snarl in the window, in order:
        call it                            # reads already resident
    release the window
```

#### Why this is cheap to do

Three facts from the source make it much smaller than it sounds:

- **`SnarlManager::roots` is already a materialised `vector<const Snarl*>`.** Grouping means collecting
  it via `for_each_top_level_snarl` and sorting by node ID. No new traversal, no index, no new data
  structure.
- **Output order is already independent of visit order.** `VCFOutputCaller` buffers `(sort_key, string)`
  per thread and `write_variants` sorts before emitting, so reordering the visit does not reorder the
  VCF. This has to be *verified* rather than assumed -- see the check below -- but the design intends it.
- **The windowed cache from the discarded design becomes trivial.** With access ordered, retaining a
  single window per thread is sufficient and every window is fetched once. The cache stops being a
  gamble and becomes bookkeeping.

Sketch, in `GraphCaller::call_top_level_snarls`:

```cpp
vector<const Snarl*> roots;
snarl_manager.for_each_top_level_snarl([&](const Snarl* s) { roots.push_back(s); });
sort(roots.begin(), roots.end(), [](const Snarl* a, const Snarl* b) {
    return min(a->start().node_id(), a->end().node_id())
         < min(b->start().node_id(), b->end().node_id());
});
// partition into contiguous windows, then:
#pragma omp parallel for schedule(dynamic, 1)
for (size_t w = 0; w < windows.size(); ++w) {
    for (const Snarl* snarl : windows[w]) process_snarl(snarl);
}
```

The recursion rounds get the same treatment: each round is already a `vector<const Snarl*>`, so sorting
it before the existing `parallel for` costs three lines and keeps queued children window-ordered too.

Node-ID space is the right axis and not arbitrarily: `StreamIndex` bins on **node-ID bit-string
prefixes** (a radix tree over node IDs) and a `.gai` GAM is sorted by node ID, so a node-ID window is a
contiguous stretch of the file. Reference-position windows would be converted to node ranges to query
the index anyway.

#### What this costs, honestly

**It changes `graph_caller.cpp`, which the default caller shares.** That is the real difference from the
cache-only design, which was confined to one class. Two consequences:

- **It must be gated.** A `bool` on `GraphCaller`, set by `call_main` only when an indexed read source is
  in use, so the default path keeps its current order and its byte-identical output by construction.
  Without the gate we would have to prove order-independence for every mode and option combination,
  which is a much larger claim than this change deserves.
- **Parallel granularity gets coarser** -- a task per window rather than per snarl. That is *better* for
  read locality and reduces task overhead, but it risks load imbalance if windows differ a lot in cost.
  `schedule(dynamic, 1)` plus windows large enough to hold many snarls should cover it; if it does not,
  the fix is smaller windows, which costs queries rather than correctness.

One correctness detail: a snarl whose contents extend beyond its boundary node IDs may need reads from
an adjacent window. The read source must still fetch whatever the requested ranges demand rather than
assuming the current window covers them -- so window alignment is an *optimisation* of the fetch
pattern, never a constraint on what a site is allowed to ask for.

Roughly 60-80 lines: collect and sort roots, a window-partitioned loop, three lines per recursion round,
the gate, and simplifying the existing cache to one entry.

#### Measured

Implemented and gated on `--gam-index`. On the 400 kb simulation, the case this was meant to fix:

| | queries | cache hits | time | peak RSS |
|---|---|---|---|---|
| in-memory | n/a | n/a | 0.40s | 84 MB |
| indexed, unordered | 920 | 0% | 3.03s | 47 MB |
| indexed, **ordered** | 920 asked, 17 fetched | **98%** | **0.42s** | **54 MB** |

Calls byte-identical to in-memory at `-t 1` and `-t 6`. That was the load-bearing check: output was
*expected* to be order-independent because `VCFOutputCaller` buffers with sort keys and sorts before
emitting, and it is now demonstrated rather than assumed.

**The window size is a genuine trade, and the default is a compromise rather than an optimum.** Cost
per window is the reads in its node span; benefit is the number of sites it serves. So density decides:

| dataset | w=32 | w=128 | w=512 | w=1024 | w=4096 |
|---|---|---|---|---|---|
| HGSVC, 29 sites / 24k nodes | 0.37s 56 MB | 0.36s 58 MB | 0.43s 77 MB | 0.54s 89 MB | 0.70s 195 MB |
| 400 kb sim, 920 sites | - | 0.75s 46 MB | - | 0.42s 53 MB | 0.42s 76 MB |

Sparse sites want small windows; dense sites want larger. `--read-window` exposes it and the default of
256 is within reach of the best on both without being right for either. Note the sparse case is *worse*
than unordered indexed at large windows -- windowing fetches reads no site wants, so it is not
unconditionally an improvement.

**Auto-sizing is the obvious follow-up and is deliberately not attempted yet**, because the natural
heuristic is wrong. "Aim for k sites per window" sizes windows *up* when sites are sparse, which is
exactly backwards: when sites are further apart than a window can usefully span, the right answer is not
to window at all. A correct rule has to compare the reads a window would pull against the sites it would
serve, and that wants measuring on real data rather than guessing.

#### Verification

1. **Default path byte-identical**, gate off -- the existing guarantee, unchanged.
2. **Gate on, calls byte-identical** to the unordered indexed run. This is the load-bearing check: it
   demonstrates the reorder is observationally neutral, which is the assumption the whole design rests
   on. If the VCF differs, something is order-dependent that the design assumed was not, and that is a
   finding worth having before building on it.
3. **Query count** falls to roughly `span / window`, visible in the existing `--progress` counter.
4. **Peak RSS** stays at or below the current indexed backend -- the point of ordering rather than
   caching is that this number should *not* regress.
5. **Wall clock** on the 400 kb simulation, where the unordered version is 8.4x slower than in-memory.
   That is the case this is meant to fix.

Doing this before gaf-base still matters, and more so now: gaf-base's per-query cost is higher (SQLite,
plus one query per node per alignment path to reconstruct targets), so it needs the ordered access
pattern more, and building it first tells us how much of gaf-base's benefit is the backend rather than
the access pattern.

### 6.3 Optional extension: tabix GAF (also no new dependencies)

Worth knowing this exists, though it is low priority. `for_each_gaf_record_in_ranges`
([alignment.hpp:81](vg/src/alignment.hpp:81), implemented at
[alignment.cpp:407-460](vg/src/alignment.cpp:407)) already does indexed GAF lookup via htslib tabix on
a bgzipped sorted GAF, querying pseudo-contigs named `{node}`. htslib is already a vg dependency, so
this needs no new one either — just wiring into `SiteReadSource`.

Three caveats, which are why it is not the recommended extension: it hands back **unparsed GAF line
strings**; it dedupes by retaining every returned line in an `unordered_set<std::string>` because
htslib lacks a multi-region iterator (a real memory hazard, flagged with a `TODO` at the call site);
and nothing in vg builds the `.tbi`, so users must produce it externally. There is also no in-tree test
for `vg chunk -F` / `vg find -F`, the only existing consumers. Implement only if users turn out to have
bgzipped+tabixed GAF and no appetite for converting.

### 6.4 gbz-base / gaf-base — implemented, and not the win it was expected to be

**Status: built and tested (stage 7).** `GafBaseSiteReadSource`, selected with `--gaf-base`. It is the
only backend with a dependency, and that dependency is *runtime only*: it runs `gbz-base query` as a
subprocess, so nothing links and no build environment changes.

**The headline result is that it works, meets the memory goal, and is still not the backend to reach
for.** Measured on the 400 kb / 80 k-read simulation, single-threaded, all four producing the same 889
variants:

| backend | wall | peak RSS | disk |
|---|---|---|---|
| in-memory GAM | 0.73 s | 309 MB | 7.6 MB GAM |
| in-memory GAF | 0.79 s | 263 MB | 22 MB GAF text |
| indexed GAM (`.gai`) | **0.92 s** | **56 MB** | 7.6 MB + index |
| GAF-Base (`--gaf-base`) | 1.90 s *(best window)* | **52 MB** | **2.8 MB db** |

So gaf-base does bound memory — 309 MB → 52 MB, a 5.9× reduction — but **the indexed GAM backend
already did that, at 56 MB, roughly twice as fast, and with no runtime dependency at all.** On memory,
which is what promoted this stage, the two are a tie. gaf-base's remaining real advantages are **disk
footprint** (2.8 MB against 7.6 MB, and against 22 MB of GAF text) and that it does not over-fetch the
way a `.gai` group scan does. Those are worth having, but they are not why this was prioritised, and on this
dataset the honest conclusion is that **stage 6 was the load-bearing memory work and stage 7 is an
alternative rather than an improvement.**

**That conclusion does not survive contact with a real chromosome, and it is worth saying why.** This
400 kb simulation has 80 k reads and a 2.8 MB database; the cost that dominates at scale is per-query
overhead paid thousands of times, which barely registers over 20 queries. After §6.5, GAF-Base calls
chr20 in 99 s against a 150 s floor for holding every read in memory — so at the scale the backend
exists for, it is faster *and* uses a fifth of the memory. Small-fixture timings ranked these backends
in the wrong order, which is an argument for tier-2 measurement rather than against the fixture.

Window size matters more here than for the GAM index, because a query is far more expensive:

| window | subprocess queries | wall | peak RSS |
|---|---|---|---|
| 256 | 69 | 2.99 s | 52 MB |
| 1024 | 20 | **1.90 s** | 82 MB |
| 4096 | 5 | 2.14 s | 202 MB |
| 16384 | 5 | 4.88 s | 262 MB |

The default of 256 is inherited from stage 6b and is *not* the right default for this backend; 1024 is
better on this dataset. **Now done** — §6.5 measured it at chromosome scale, where the gap is much
wider than here (256 → 180 s against 4096 → 99 s), and `--read-window` defaults per backend: 4096 for
`--gaf-base`, 256 for `--gam-index`.

#### Where the time actually goes — and why this changes the C-shim ask

The natural assumption is that the subprocess is the problem. It is not. Timing a single realistic
1024-node window query (5512 reads returned):

| component | cost | wanted? |
|---|---|---|
| process spawn | 3 ms | no — but nearly free |
| subgraph extraction (GFA to `/dev/null`) | 20 ms | **no — pure waste** |
| read decode | 41 ms | yes |
| **total** | **62 ms** | |

**The `fork`/`execvp` is 5% of the cost.** Replacing the subprocess with a library call would recover
almost nothing. What *would* help is not extracting a subgraph we immediately discard: `gbz-base query`
is a subgraph tool that can also return the overlapping reads, and we only ever want the reads, so a
third of every query is thrown away. That is 20 ms of 62 ms, and unlike the spawn it is not incidental.

This reverses the framing in [gbz-base-c-api-request.md](gbz-base-c-api-request.md), which treats the
subprocess as the compromise and a C ABI as the fix. The measurement says the subprocess is fine and the
useful ask is **a reads-only query** — by boundary nodes or node set, returning GAF text, doing no
subgraph work. Whether that arrives as a library symbol or as a CLI flag on the existing binary matters
much less than that it exists. The letter has been revised accordingly.

The remaining ~41 ms of decode is not obviously unreasonable: 7.4 µs per read, against a few µs per read
for vg's own GAF text parse. We also pay some of it twice, because reads spanning a window boundary are
returned by both windows (82,021 reads fetched for 80,000 reads at window 256, versus 80,072 at 4096).

#### Original rationale, retained

<https://github.com/jltsiren/gbz-base> — Rust, SQLite-backed. **GBZ-Base** stores the graph;
**GAF-Base** stores read alignments.

**Best conceptual fit of the four backends.** Its query surface includes **snarl-based extraction
between boundary nodes**, returning the reads overlapping the resulting subgraph. That is precisely
our access pattern — a snarl *is* two boundary nodes — and it is why `SiteReadSource` exposes an
overridable `get_reads_in_snarl()`.

**The route taken: subprocess now; an upstream reads-only query if upstream is willing.** (Written
before the measurement above, and correct in its conclusion for a reason it did not anticipate: the
subprocess turned out to be a fine permanent answer, not an interim one.)

Reading the SQLite directly is rejected. GAF-Base does not store GAF text: one `Alignments` row is a
*block* of many reads re-encoded column-wise into seven blobs (GBWT ByteCode varints, zstd, rANS
4x16, GBWT RLE, an ad-hoc bitvector), with no per-alignment offset table and no stored target path —
paths are recovered by LF-stepping a GBWT held in a second table. A C++ reader means ~400–600 lines
whose hard part is a conditional column layout that silently desynchronises if misordered. And the
format is explicitly unstable: the README says it can change without warning, the version tag is
exact-match enforced, and there were three format-breaking changes in about four months. There is
definitively no C API.

So:

1. **Ask upstream for a small `extern "C"` surface** over the existing Rust `ReadSet` — returning
   **GAF text** rather than a struct, so the ABI survives internal format changes and vg reuses
   `gafkluge` unchanged. A draft request is in
   [gbz-base-c-api-request.md](gbz-base-c-api-request.md) (not yet sent). Ask for the **snarl
   boundary-node query as the primary entry point**, not a node-id list: the boundary-node form is the
   entire reason this backend is the right fit, and a node-id list would just be the `GAMIndex` access
   pattern with extra steps. Two other things worth settling in the same conversation: whether the
   returned GAF carries **MAPQ and per-base qualities** (the model needs the first and wants the
   second, and `--no-quality` databases exist — if quality is commonly dropped, the fallback scorer of
   §5.2 becomes the normal case rather than an edge case); and whether gaf-base support should be an
   **optional compile-time feature** in vg, since a Rust `staticlib` would put the Rust toolchain in
   every build environment including release CI, and a `cdylib` would break vg's `-static` target.
2. **Subprocess — what was built.** `fork()`+`execvp()` of `gbz-base query`, writing GAF to a per-thread
   temp file via `--gaf-output` and parsing it with `vg::io::gaf_unpaired_for_each`. Upstream's binary
   does the decoding, so format churn costs us nothing. Precedent: vg already shells out to `kmc` this
   way ([index_registry.cpp:4469](vg/src/index_registry.cpp:4469)). Spawn cost is amortised over a
   **window of node IDs**, never per snarl, reusing the stage-6b machinery.

   The actual invocation, which the CLI's own defaults get wrong for this purpose in two ways:

   ```
   gbz-base query <graph.gbz.db> -n <id> [-n <id> ...] --context 0 \
       --gaf-base <reads.gaf.db> --gaf-output <tmp> --alignments overlapping
   ```

   `--context 0` because the default of 100 bp expands the subgraph past the nodes asked for;
   `--alignments overlapping` because the default is `clipped` (see the correctness note below). The
   subgraph GFA goes to `/dev/null`; only `--gaf-output` is read. Note it is `gbz-base query`, not
   `gaf-base` — `gaf-base` has only `construct`/`decompress`/`sort` and no random-access query at all,
   so the graph-side binary is the one that reaches the reads.

   Two implementation details worth recording. **Only node IDs that exist are passed**, since ID space is
   sparse and the tool is entitled to object to a node that is not there; vg holds the graph, so this is
   exact rather than a guess. And **a query too large for one `argv` is split into chunks and the results
   de-duplicated by read name** — sound because the filter keeps at most one alignment per name
   (`skip_secondary`, the default), which is also what the independence assumption in §4.1 requires.

This was not wasted work against a future shim: both paths consume GAF text behind the same class, so
they share the parsing, overlap mode, filtering, windowing, and caching. Swapping `execvp` for a library
call touches one function.

The subprocess route adds a *runtime* dependency on an external binary but **no build dependency** —
nothing links, no build environment changes, and a missing binary is an error only for users who ask for
this backend. That error is routed through `vg call`'s normal error path rather than allowed to escape as
an exception, so "install gbz-base" does not surface as a crash inviting a bug report. A failure *during*
calling reports and exits instead of throwing, because an exception must not propagate out of an OpenMP
parallel region, and a half-fetched site would silently produce a wrong genotype.

**Two correctness requirements for either path, both now enforced:**

- **Use the `Overlapping` overlap mode, not the CLI default `clipped`.** `Clipped` can split one read
  into multiple fragments, which would place a single read in several matrix rows and break the
  independence assumption in §4.1. Checked directly: on a 9-node query `overlapping` returned 290
  records with 290 distinct names, `clipped` 289, and `contained` **0** — reads are longer than a small
  window, so `contained` is useless here, which is worth knowing before anyone reaches for it.
- **Check graph consistency ourselves.** The reads' database may legitimately have been built against
  a *supergraph* of the graph passed to `vg call`, and upstream only checks compatibility softly via
  stable graph names, with no error if a name is absent. `check_setup()` runs one probe query on the main
  thread at startup, so a mismatched or unreadable database fails immediately with an actionable message
  rather than on the first snarl inside a worker thread.

#### A GAM converted to GAF does not give identical calls, and the reason is not quality loss

Worth stating because it looks like a bug in this backend and is not. Feeding `vg call` a GAM, versus
converting that GAM to GAF and building a gaf-base from it, gives the same variants and genotypes but
DP differing by ±1 at a few sites and correspondingly shifted `GL`s. The backend is exonerated
precisely: **calls from gaf-base are byte-identical to reading the same GAF in memory**, on both the
1 kb and 400 kb datasets. The difference is entirely in `vg convert -G`.

There are **two distinct causes**, and only one of them is arguably a vg defect. Both were found by
round-tripping GAM -> GAF -> GAM and diffing the alignment paths; 47 of 2000 reads differ.

**Cause 1 (42 of the 47): an insertion re-partitioned across a node boundary. Not a defect.** A GAM
mapping carries the insertion at the *start* of the following node's mapping; the GAF round trip attaches
it to the *end* of the preceding one:

```
GAM : node 135 [(1,1,"")]              node 137 [(0,1,"G"), (11,11,"")]
GAF : node 135 [(1,1,""), (0,1,"G")]   node 137 [(11,11,"")]
```

The node set is **unchanged**, and both describe the same alignment: an insertion between two nodes
belongs to neither, and GAF's difference string has no way to record -- and no reason to record -- which
side it came from. GAM simply has a redundant degree of freedom here that GAF does not.

It still perturbs our numbers, and that is a property of *this caller* rather than of the formats.
`step.read_length` is `mapping_to_length(mapping)`, which counts inserted bases, so moving an insertion
across a boundary moves one read base between two nodes' steps. If one of those nodes is in the site and
the other is not, the read's scoring window (§3.2) changes by a base and its likelihoods shift slightly.
This does **not** make a call wrong -- the window invariant still holds, so every allele is scored over the
same span and the comparison within a read is unaffected -- but it does mean the same alignment expressed
two ways yields slightly different `GL`s. Making the window insensitive to boundary attribution is a
deliberate design question, recorded as an open item rather than silently changed.

**Cause 2 (5 of the 47, i.e. 5 of 1872 mapped reads): a trailing path node dropped. This one is a real
round-trip asymmetry in `vg convert`.** When a read's final mapping consumes **zero reference bases** -- a
pure trailing insertion on a node the read does not actually align to -- `vg convert -G` omits that node
from the GAF path string:

```
GAM path : ... >132 >133 >134     with node 134 carrying only (0,6,"TTATCT")
GAF path : ... >132 >133          node 134 absent
```

Confirmed to be the *writer*, not the reader: the GAF `path` column itself has no `134`. Exactly
characterised -- of 1872 mapped reads, 5 have a final mapping consuming no reference, and those 5 are
precisely the 5 whose node set changes. Never a leading mapping, never mid-path. **This is what moves
DP**: the read no longer touches that node, so at a site whose only overlap with the read was that node it
stops being evidence, and DP falls by one at each of ~5 sites.

Whether the defect is in `vg convert -G` dropping the node or in `vg giraffe` emitting a zero-reference
trailing mapping at all is a fair question for upstream; either way the round trip is silently asymmetric
and vg has no test covering it. Worth a narrow issue; not worth blocking on.

The consequence for testing is the same in both cases: **the equivalence test compares gaf-base against
in-memory GAF, not against GAM**, or it would be testing `vg convert`.

#### Shared with the GAM index rather than duplicated

Stage 7 also refactored what stage 6 and 6b had built. The window arithmetic, the per-thread cache, the
boundary-straddling bypass, and narrowing a window's reads back down to the ranges a site asked about are
now in a `WindowedSiteReadSource` base class, and each backend supplies one primitive:

```cpp
virtual void fetch_span(nid_t min_id, nid_t max_id,
                        const function<void(const Alignment&)>& iteratee) const = 0;
```

This removed ~90 duplicated lines, made the second backend small, and stops the two paths drifting apart
in how they interpret a query — a filter applied in one and forgotten in the other would mean the two
backends genotyped different read sets. It also made that shared logic **unit-testable for the first
time**: a fake backend that records the spans it was asked for tests the windowing without a GAM index,
without a GAF-Base, and without the `gbz-base` binary (12 cases, 36 assertions in
`src/unittest/site_read_source.cpp`). The refactor is verified behaviour-preserving by the indexed GAM
backend still producing calls identical to the in-memory one.

---

### 6.5 Making it fast — measured, and 5.8x

**Status: done.** HG002 chr20 on the HPRC v2.1 MC CHM13 graph, 6 threads, GAF-Base
backend: **570 s → 99 s**, peak RSS 4.25 GB → 3.74 GB, and the VCF is **byte-identical**
at every step.

5.4 predicted that "the bottleneck moves to read retrieval, not scoring — that is where
optimisation effort should go." That was right, and the effort went there. Where it
guessed *within* retrieval — the database, the process spawn — it was wrong, and two of
the four things that looked expensive turned out not to be.

#### The baseline, and the floor

| | wall | user | sys | peak RSS |
|---|---|---|---|---|
| `poisson-z` (reads a pack, never touches a read) | 74 s | — | — | 2.5 GB |
| `readlik-z` GAF-Base, before | 570 s | 680 s | **520 s** | 4.25 GB |
| `readlik-z` in-memory GAF (`--gaf-reads`) | 150 s | 170 s | 30 s | 18.2 GB |
| **`readlik-z` GAF-Base, after** | **99 s** | 280 s | 102 s | **3.74 GB** |

The in-memory run was the first measurement taken and the most useful: same 105,936
variants with no windowing, no cache and no subprocess, so everything above 150 s was
retrieval and everything below it was not. sys being 88% of wall in the baseline said
the same thing from the other side — that is an I/O and allocation profile, not a
genotyper's. The final number is *below* that floor because GAF-Base parses reads
across six threads where `--gaf-reads` parses 4.9 GB of text on one.

Re-running the whole tier-2 arm set through the evaluation harness afterwards agrees: `readlik-z` 97 s,
`readlik` 115 s, `readlik-nomismap` 115 s, against `poisson` 156 s and `poisson-z` 72 s on the same
machine in the same session. **`readlik` is now faster than `poisson` at matched enumeration**, and
every accuracy number in `docs/tier2-chr20-results.md` is unchanged — the only lines that moved in that
file were the cost table.

#### What it actually was

**Fetching by span instead of by ranges. This is the whole win, and nothing suggested
it in advance.** `for_each_read` collapsed the caller's precise node ranges to
`[min_id, max_id]` before fetching. For a site inside one window that is exact. For a
snarl whose contents are *sparse in ID space* it is not: on chr20, **215 sites — 0.16%
of all sites — spanned 13.2 M node IDs while wanting 133 k of them**, a 99× over-fetch.
Those 215 sites accounted for most of the reads fetched and about half the subprocess
queries. Reads fetched fell **32.7 M → 13.4 M** against the 13.3 M actually on the
chromosome, so the 2.4× redundancy is gone rather than reduced. Both backends address
ranges natively — `GAMIndex::find` already took a range list — so `fetch_span` now takes
ranges, and the GAM backend got shorter.

**`posix_spawn` instead of `fork`/`exec`, worth 26%.** Not a style preference. Forking
from a process whose other threads are allocating hard makes libc take a fork lock
around `malloc`; a profile showed threads that never fork stalled in
`_xzm_fork_lock_wait`. `posix_spawn` never duplicates the address space, needs no child
branch, and needs no stdio flush — so it is also less code than what it replaced.

**Taking ownership of each read rather than deep-copying the protobuf, worth 14%.** The
cache did `push_back(aln)`, and `Alignment`'s copy constructor was 16% of the profile on
its own: 32 M reads cached over the chromosome, each copy allocating a `Path`, its
`Mapping`s and their `Edit`s. `fetch_span`'s iteratee now takes a mutable reference so
the cache can move; the backends hand over a per-record scratch alignment they `Clear()`
before reuse, so moving from it is safe. Everything above that still sees `const`.

**Rejecting cached reads on node-ID bounds before walking their mappings — nearly free
on its own.** The counters said 610 M reads were considered to deliver 7.5 M, a
selectivity of 1.2%, which looked like the obvious target. It was worth about 1%: the
scan was already cheap. It earns its place only because it is what makes a larger window
affordable, and the larger window is worth 45%.

#### Window size, per backend

`--read-window` now defaults to **4096 for `--gaf-base` and 256 for `--gam-index`**. One
default never suited both — a GAF-Base query is a process spawn where a `.gai` group
scan is a seek — which 6.4 noted as "a cheap future improvement with a real payoff".

| window | wall | peak RSS | subprocess queries |
|---|---|---|---|
| 256 (old default) | 180 s | 3.24 GB | 13,825 |
| 1024 | 114 s | 3.46 GB | 3,870 |
| **4096** | **99 s** | 3.74 GB | 1,130 |
| 16384 | 101 s | 6.03 GB | 806 |

#### What is left, and what is not worth doing

**The discarded subgraph is still discarded.** `gbz-base query` builds a GFA subgraph we
send to `/dev/null`; 6.4 measured that at 20 ms of a 62 ms query. `--reference-only`
looked like it made queries 15× faster and does not — it panics on node-based queries
with "Reference-only output is not supported", and the timing that suggested otherwise
was reading a stale output file. So the reads-only ask in
[gbz-base-c-api-request.md](gbz-base-c-api-request.md) stands, and is the main remaining
lever outside vg.

**Phase 3's streaming sweep is not needed.** It was reserved for the case where the
cheap work left a large gap. It did not: 99 s against a 150 s in-memory floor and 74 s
for a caller that reads a precomputed pack. Inverting `SiteReadSource`'s contract to buy
what is left would not pay for itself.

**Parallel efficiency is no longer the interesting number.** It was 2.2 threads of 6
before and is 3.9 of 6 now, but much of the remaining wait is `wait4` on a child that is
itself using a core, so the machine is busier than the vg-thread figure suggests.


---

## 7. Correctness requirements

These are the places where a wrong implementation still produces plausible-looking VCF. Each needs
its own unit test.

**Partial overlap: the scoring window belongs to the read and the site, never to the allele.** A read
often overlaps only part of a site. Define its window once, as the **maximum span of the read within
the site** taken over all alleles — every read base that any allele can place is informative, so the
window is a union, not a per-allele overlap — and then score *every* allele over that same window,
charging read bases a given allele cannot place as insertions rather than omitting them.

This is an invariant, not a preference. Score allele *a* over 150 read bases and allele *b* over only
the 100 it can place, and the difference contains 50 × match-reward: a fabricated likelihood ratio of
arbitrary magnitude that merely happens to point the right way. `P(read | allele)` may differ between
alleles only in how well the *same* observed bases are explained. **Assert it** — the read-base count
entering `score(r, a)` must be identical for every *a*.

- **Alleles identical over the window score identically.** A read whose window spans a region where
  alleles *a* and *b* agree must score the same on both — discriminating against *c* while staying
  neutral between *a* and *b*. Structural comparison over a fixed window gives this for free, and it is
  how a partial read contributes *partial* information.
- **Unequal window lengths *between reads* are fine and must not be "fixed".** Reads contribute unequal
  information, which is correct. It does not break genotype comparison, because normalisation (§4.2) is
  **per-read** and so is insensitive to how long that read's window was. Put a comment in the code so
  nobody normalises per base. Unequal spans *across alleles for one read* are the opposite: that is the
  calibration error the invariant above exists to prevent.

**"No valid placement" ≠ "uninformative".** A read covering only bases that a deletion allele deletes
places perfectly on the reference and **not at all** on the deletion. That is *strong evidence
against the deletion*, not missing information. Give it `rel = 0`; the `+ e_r` term in §4.3 bounds the
penalty and keeps the log finite. The distinction that matters: `rel = 0` must never be conflated with
dropping the read or giving it an equal entry across alleles. Getting this backwards would
systematically destroy deletion and SV genotyping — most of the point of the exercise.

In practice the window rule above usually gets there on its own: charging every unplaceable base as an
insertion yields a very small *non-zero* `rel`, which is the same evidence more gracefully expressed,
and §4.3 bounds it either way. Reserve exactly `rel = 0` for an allele admitting no alignment of the
window at all.

**Guard the read that places on nothing.** If a read has no valid placement on *any* allele, its row
max is zero and the §4.2 normalisation divides by zero (`ell*(r) = -inf`, giving `NaN`). This is
reachable in practice: the retrieval layer fetches by node-ID range, so a read can overlap the range
yet place on no traversal. **Drop such reads explicitly, and count them.** A rising count is a signal
in its own right — it means the read source is over-fetching, or the reads and graph do not match
(§6.4).

**Traversal indices shift before `update_vcf_info`, and two entries may not be traversals at all.**
`emit_variant` deduplicates alleles by sequence string and drops uncalled ones *before* calling
`update_vcf_info` ([graph_caller.cpp:665](vg/src/graph_caller.cpp:665)). Three distinct hazards, not
one:

- **The indices no longer match the matrix.** The Poisson caller copes by recomputing support;
  recomputing the *matrix* would mean redoing all the scoring.
- **⚠ `site_genotype` can contain negative sentinels.** `STAR_ALLELE_MARKER = -2` and
  `MISSING_ALLELE_MARKER = -1` ([graph_caller.hpp:29](vg/src/graph_caller.hpp:29)) are pushed into it in
  `--top-down` nested mode, and a default-constructed placeholder `SnarlTraversal()` is pushed into
  `site_traversals` to stand in for the star allele
  ([graph_caller.cpp:539-566](vg/src/graph_caller.cpp:539)). Indexing a `GL` table by `site_genotype`
  reads out of bounds — a crash, not a mismatch.
- **Dedup is by allele *string*, so two distinct scored traversals can collapse into one VCF allele.**
  `GL` for that allele needs an explicit aggregation rule; use `max` over the merged set, which keeps
  `GL` a likelihood of the emitted allele rather than of an arbitrary mixture.

**As built, `emit_variant` was left alone.** The caller keeps the traversals it scored in its
`CallInfo` and maps the deduplicated ones back by **structural comparison** of their visit lists —
exact regardless of how allele strings were flattened, and it needs no change to `graph_caller.cpp`
at all. An allele that matches nothing (the empty placeholder standing in for a star allele) simply
gets no `GL` entry. `GL` is additionally suppressed whenever the emitted genotype carries a star or
missing marker, because such a genotype was called at a lower ploidy than the record reports (see Nested calling in §7) and
the vector length would disagree with the ploidy `GT` implies. The option below is recorded as the
alternative if structural matching ever proves too fragile:
**extend `emit_variant` to return a `site_index -> vector<scored_index>` mapping**, emit `GL` only over
the emitted allele set (which is what VCF requires anyway), and skip the sentinels explicitly. Small
and contained, but it will cost a day if discovered late.

**Soft clips: excluded by default, with a flag to include them.** Soft-clipped bases lie outside the
read's aligned span, so excluding them means the window of the rule above is that aligned span. A read
clipped at the site boundary may be clipped *because* it disagrees with the graph, which is in principle
informative — but a clip is ambiguous evidence (it could equally be adapter, quality trim, or a
chimeric read), vg's existing precedent is to strip them (`augment_from_alignment_edits`), and widening
the window with bases *no* allele can place loads every allele equally and so buys nothing unless some
allele genuinely places them. Excluding is the conservative default: it loses signal rather than
inventing it. The flag exists so the alternative is measurable rather than theoretical.

Implementation note: this must be *deliberate* rather than incidental. Clip the read to its aligned span
before establishing the window, and do not let clipped bases silently become mismatches against an
allele — that would penalise exactly the reads carrying SV signal. *Signal to revisit:* SV recall
materially worse than the Poisson caller at stage 4, with soft-clipped reads over-represented at the
missed sites.

**Read eligibility must match the pack file's.** `vg pack -Q` sets *both* `min_mapq` and `min_baseq`
from the one value ([pack_main.cpp:148](vg/src/subcommand/pack_main.cpp:148)), and `--trim-ends` drops
read ends. Allele *enumeration* is pack-driven while *genotyping* is read-source-driven (§2.3), so
without a matching eligibility filter the two disagree about what evidence exists at a site. State a
policy at stage 0: default to mirroring the pack file's filters, and decide explicitly about
secondary/supplementary alignments and duplicates, which nothing on vg's pack path excludes either.

**Two pre-existing bugs in the neighbourhood** — fix or avoid, do not build on:
`depth_err` is assigned from a malformed ternary at
[snarl_caller.cpp:602](vg/src/snarl_caller.cpp:602) — `depth_info.second ? !isnan(depth_info.second) : 0.`
— so it is always 0.0 or 1.0; and `GBWTTraversalFinder` has an index mismatch in its backward-dedup path
([traversal_finder.cpp:3519-3547](vg/src/traversal_finder.cpp:3519)).

**Correction: this bug is inert.** The only consumer of `depth_err` inside
`genotype_likelihood` is commented out, deliberately, with the rationale that the small bin sizes make
the binned-coverage error far too large to be useful; `depth_err` is otherwise only carried on the
`CallInfo` and printed in debug output, and never reaches the VCF. Fixing the ternary produces
byte-identical calls (verified over three simulated 400 kb replicates). So it is worth fixing as
hygiene — a latent trap for whoever re-enables that line — but it is **not** a prerequisite for stage 4
and the Poisson caller is **not** a mis-parameterised baseline. Earlier drafts of this document claimed
otherwise; that claim was wrong.

**Use long-only options.** `vg call`'s getopt string has only `jnquwxy` / `DEFHJKQUVWXZ` free. Add
`--gam`, `--gaf`, `--read-likelihood` etc. rather than burning short letters.

**The default path must not change.** Since this is opt-in (§1), a `vg call` run without
`--read-likelihood` must produce byte-identical output to the current build. Concretely:

- Construct the new caller only inside the `--read-likelihood` branch of `call_main.cpp`; do not alter
  the existing construction of `PoissonSupportSnarlCaller`.
- Do not add `FORMAT`/`INFO` header lines or fields unconditionally — the new `GL`/`GQ` semantics and
  any new annotations belong inside `update_vcf_header` / `update_vcf_info` of the *new* caller only.
- Resist the temptation to "tidy" shared code in `SupportBasedSnarlCaller` while passing through it. The
  two pre-existing bugs above are real, but fixing them changes default output and so belongs in a
  separate, separately-reviewed change.
- **Validate option combinations early.** `--read-likelihood` without `--gam`/`--gaf-reads` must be a
  clear startup error, not a run that silently genotypes with no reads. (Note `--gaf` was already taken
  as the long form of `-G`, GAF *output*, hence `--gaf-reads` for the input.)

### Nested calling

This document was originally silent on nesting. Both modes are supported; the third is out of scope.

**Independent (`-A`) works with no changes.** `RecurseAlways` invokes the caller once per snarl, so each
site is genotyped against its own reads with no parent restriction. Every record is then a
self-contained `P(reads | G)` statement, which is the easiest thing to reason about. The cost is that
nothing ties the levels together: no phase, and a child call can contradict its parent. It emits no
`PS` tags, which is the honest signal that no phase is being claimed.

**Top-down (`--top-down`) needed a fix, and not in this caller.** When a parent allele does not traverse
a child snarl, the child is reached by fewer haplotypes than the parent's ploidy. `FlowCaller` asked the
genotyper for a full-ploidy genotype anyway and then overwrote the positions belonging to
non-traversing parent alleles with star or missing markers. Wrong twice over:

- **A site one haplotype reaches is not diploid.** Asking for a diploid call there lets a spurious
  heterozygote absorb minority reads for free, biasing which allele is picked before any marker is
  applied. This matters more for a likelihood model than for a count model, because the ½/½ mixture
  weight is doing real work.
- **`genotype()` returns a sorted allele multiset with no haplotype identity**, so overwriting position
  *i* bore no relation to which parent haplotype was empty. Which allele got discarded was effectively
  arbitrary. This affects `PoissonSupportSnarlCaller` identically.

The child is now genotyped at the ploidy that actually traverses it, and the called alleles are
scattered back onto the traversing haplotypes. The effect is measurable rather than theoretical: on a
site reached by one haplotype with 100 reads supporting the reference allele and 30 the alternative,
reported `GQ` goes from **38 to 256** — the diploid misspecification had been spending its confidence
explaining minority reads as a second haplotype that does not exist there.

⚠ **This part is not opt-in.** It changes `graph_caller.cpp`, which the default caller uses. It changed
no existing test outcome (`18_vg_call.t`'s failure set is identical before and after) but it can change
`--top-down` output on other data, deliberately. It is defensible on its own merits as a bug fix, so it
probably belongs in its own separately-reviewed change rather than inside an opt-in feature.

**Deliberately not done: a joint model over the snarl tree.** Two known consequences of leaving it out.
Reads are reused across levels — a read informative about a child is also informative about its parent,
so the two records are *not* independent evidence and their `GQ`s must not be combined. On a simple
nested SNP fixture the parent and child emitted *identical* `GL` vectors, because the parent's two
alleles differed only by the child SNP; correct, but it shows the same evidence being reported twice.
Properly, the parent likelihood would marginalise over child alleles. That is combinatorial and is a
research direction, not a fix.

**The opportunity not taken: read-backed phasing.** A read spanning both a parent and a child site
directly links a parent allele to a child allele. With matrices at both levels, joint (parent, child)
assignments could be scored and phase derived from evidence rather than inherited by construction —
something the support-based caller structurally cannot do, since it has no per-read identity. This is
the strongest argument for read-level genotyping in nested mode and is the obvious next extension:
keep the parent's matrix keyed by read name and intersect the read sets.

---

## 8. Staging

Each stage leaves the tree working and is independently testable. **Stages 0–5 add no dependency of
any kind** — no new submodule, no new library, no external binary, no changes to any external project.
The first thing needing a dependency decision is stage 7, which is why it is last.

**Off the critical path:** open the upstream conversation about gaf-base (§6.4). It has a lead time we do
not control and gates nothing. **What to ask for has changed** now that stage 7 is measured: the ask is a
reads-only query rather than a C ABI, because the subprocess is 3 ms of a 62 ms query while the subgraph
extraction we discard is 20 ms. The one-line `depth_err` fix (§7) is worth landing but is *not* a gate on
stage 4: it changes no output (§7).

Stages 0–3, 6, 6b and 7 are **built** (PR vgteam/vg#4990), together with the nested-calling support in §7
and the scoring fixes in §5.5.

**Sequencing: windowed batching (6b) before gaf-base (7).** This was the right order, and for the reason
given: stage 6 bounded memory but left query count scaling with the number of sites, and the same
batching gaf-base needs was cheaper to build and measure against the `.gai` backend first. It also
isolated how much of gaf-base's benefit is the backend and how much is merely batching — which turned out
to matter, because the answer is *almost all of it was the batching* (§6.4).

**Priority change: read retrieval comes before further accuracy work.** *(Resolved — retrieval is done
(stages 6/6b/7) and the accuracy work it unblocked has since run at tier 2. Retained for the reasoning.)*
The in-memory read source is
the binding constraint on everything — it caps runs at a small region, blocks tier 2 of the evaluation
harness, and means the accuracy comparison so far has only been run where memory is not a problem.

**Resolved, and not as expected.** Stage 6 was ordered before stage 7 because `.gai` indexed GAM bounds
memory too, is already in the vg tree, and needs no new dependency, whereas gaf-base needs a subprocess or
a shim. Both are now built and measured, and **stage 6 was the one that mattered**: 56 MB versus 52 MB
peak RSS, with the indexed GAM roughly twice as fast and dependency-free. gaf-base is the better *fit* and
much smaller on disk, but it is an alternative rather than an improvement, and the memory constraint that
promoted this work is fully addressed by stage 6 alone. **Tier 2 of the evaluation harness is therefore
unblocked, and should use `--gam-index`, not `--gaf-base`.**

**Confirmed against real data, and the requirement is harder than "unblocked" suggests.** Tier 2 is now
specified on real inputs — HG002 (~588 M reads, 28.8 GB GAF) on a 100 M-node HPRC v2.1 MC CHM13 graph, on
a 32 GB laptop; see [vg-call-eval-plan.md](vg-call-eval-plan.md) §9. Held reads cost ~3.2 KB each, so one
chromosome's ~12.5 M reads would need **~40 GB in memory** — more than the machine. A windowed backend is
therefore a *hard requirement* for tier 2, not an optimisation, which is the first time stages 6/7 have
been load-bearing rather than merely prudent.

A third thing measurement settled — and it substantially raises the value of stage 7. **On real data
the gaf-base route is ~4× cheaper to build than the indexed-GAM route, which reverses the §6.4 verdict
for this workload.** On the same 8 M-record sample of real reads: `vg convert -F` + `vg gamsort -i` took
273 s at 13.3 GB peak, while `gaf-base sort` + `gaf-base construct` took 88 s. The sort is where it
diverges — `gaf-base sort` does 8 M records in **17.2 s at 0.9 GB, single-threaded**, against
`vg gamsort`'s 249 s at 13.3 GB on 8 threads. Projected to 588 M records: **~1.5 h for gaf-base against
~6 h for GAM**, with `gaf-base construct`'s memory verified *fixed* (13.2 GB at 8 M records, 12.2 GB at
32 M) rather than proportional. Queries cost ~155 ms per 1024-node window on the real 6.27 GB GBZ-Base,
which is immaterial next to that.

So §6.4's conclusion — that the indexed GAM backend was the one that mattered — holds for *converting an
existing GAM*, and inverts when the input is already GAF, which is how mapped pangenome reads normally
arrive. Two corrections that follow: gaf-base is **not** smaller on disk than GAM on real data
(131 vs 112 bytes/record; the earlier 7.8× ratio came from a toy graph with 3-digit node IDs), and
`gbz-base construct` found top-level chains for only 2 of 46 components on this graph, so the
`--snarls`-based per-contig extraction of §6.4 is unavailable here without a distance index. See §9.7 of
the harness plan.

Two further things that measurement settled, both recorded in §9 of the harness plan. The **graph** side
does not need the per-chromosome gbz-base machinery here: this graph loads in 4.3 GB, so `vg chunk`
suffices and the genuine 32 GB risk is genome-wide *snarl* computation, avoided by working per
chromosome. And the graph carries only **4 haplotypes** (CHM13, GRCh38, two recombinants — it is a
haplotype-sampled personalised graph), so with `-z` the candidate allele set is narrow enough that
**enumeration, not genotyping, is likely the binding constraint on recall**. Both enumeration modes have
to be run or a graph limitation will be misreported as a caller limitation.

| Stage | Deliverable | Verification |
|---|---|---|
| 0 ✅ | `SiteReadSource` + `InMemorySiteReadSource` + `--gam`/`--gaf-reads` + a read-eligibility policy matching the pack file's filters (§7); log reads per snarl | reads-per-site counts sane on the existing `test/call` GAM fixture |
| 1 ✅ | `AlleleReadLikelihoods` + graph-implied-alignment calculator (§5); `--dump-likelihoods` TSV. *Built with accumulated primitives rather than a `path_t` (§5.1).* | the §9 unit tests plus `allele_likelihood_scoring.cpp` for the scoring itself (§5.5); matrix rows peak on the correct allele in simulated data |
| 2 ✅ | `ReadLikelihoodSnarlCaller`: exhaustive genotypes, uniform prior, MAPQ at face value with the `e_max` clamp and a term-off flag (§4.4); ties broken toward the reference | unit tests on synthetic matrices (clear het, clear hom, flat/no-evidence) |
| 3 ✅ | Wire into `call_main` behind `--read-likelihood`; `GL`/`GQ`/`GP` output (needs the §4.5 posterior fix; the §7 index mapping was solved by structural matching instead). `-k` optional with `-g`/`-z` (§2.3) | `t/18_vg_call.t` failure set identical to baseline; new flag produces valid VCF |
| 3b 🟡 | **De novo truth-concordance harness** — built as a separate repo, [benedictpaten/vg-call-eval](https://github.com/benedictpaten/vg-call-eval), so its dependencies stay out of vg's build. Uses `aardvark`, not truvari/vcfeval (see its plan). Tier-0 simulation, the five-arm matrix and comparison work; the `GQ` sweep and calibration remain. | positive and negative controls both pass: identical inputs score 1.0, corrupted inputs are detected |
| 4 ✅ | Accuracy comparison vs the Poisson caller. **Done on real data over two chromosomes and two graphs.** The model leads on genotype F1 on every class of every dataset, and its margin widens with graph richness. Structural variants remain the weak class for both callers. See the tier-2 summary at the end of this document | stage 3b's harness; `18_vg_call.t` for the `-v` path — plus the §4.6 risk triggers |
| 4b | *Only if stage 4 shows `GQ` inflation:* fit the `read_weight` discount (§4.4) | binned reported `GQ` vs actual error rate against truth — diagonal after fitting |
| 5 ✅ | Profile — expected read *retrieval* to dominate, not scoring. **The expectation was wrong, and then the gap closed anyway.** Retrieval was never the bottleneck at chromosome scale; the work in §6.5 removed 5.8× from the read path with byte-identical output, and `readlik-z` now runs at 99 s against `poisson-z`'s 75 s on chr20 | wall-clock on a whole chromosome, both graphs |
| 5b | *Only if stage 4 shows it matters:* bounded DP escape hatch (§5.3) | A/B on the HGSVC harness |
| 6 ✅ | `IndexedGamSiteReadSource` (`.gai`) + the per-thread cache (§6.2), selected with `--gam-index`. No new dependency | **done.** Identical calls asserted by test. On the HGSVC fixture: 354 MB to 58 MB peak RSS and 0.67s to 0.35s. Cache hit rate is reported under `--progress`: 21% with `--top-down`, 0% in flat mode |
| 6b ✅ | **Align the visit order with read windows (§6.2.1).** Sort top-level snarls by node ID and process them window by window, so each window is fetched exactly once and released. ~60-80 lines; touches `graph_caller.cpp`, so **gated** to keep the default path unchanged | default byte-identical with the gate off; **calls byte-identical with the gate on**, which is the load-bearing check; query count falls to ~span/window; peak RSS must **not** regress |
| 7 ✅ | **`GafBaseSiteReadSource` (§6.4), via `gbz-base query` as a subprocess.** Selected with `--gaf-base`, pointed at a graph with `--gbz-base`. Runtime dependency only; nothing links, and users who do not pass the flag are unaffected. Shares all windowing and caching with the GAM backend through a new `WindowedSiteReadSource` base | **done.** Calls byte-identical to the same reads read in memory, at `-t 1` and `-t 4`, and independent of window size. Peak RSS 309 MB → 53 MB. But **slower than the indexed GAM backend, which already met the memory goal** — see §6.4 |
| 7b | *If a C shim lands:* swap `execvp` for the library call behind the same class. **The measurement in §6.4 changes what to ask for**: the spawn is 3 ms of a 62 ms query, so the ask is a *reads-only* query that skips subgraph construction (20 ms), not merely a library entry point | stage 7's equivalence test unchanged, plus a schema-version guard |

**Stage 3 is the natural first merge point.** By then the feature works end to end on the in-memory
backend, the default path is provably unchanged, and nothing new has been added to the build — so it can
land and be iterated on in tree rather than growing on a long-lived branch.

**Stage 3b is what makes stage 4's verdict mean anything**, and stage 4's verdict is what gates 4b, 5b
and the whole question of promoting this past opt-in. It is the one stage whose cost is mostly tooling
rather than code, so it is the easiest to under-budget.

Stages 6–7 come *after* the model is validated deliberately: they are pure substitutions with an
exact equivalence test available — the same reads must give the same calls.

---

## 9. Validation

**A concordance harness already exists.** `test/t/18_vg_call.t:89-105` calls a real HGSVC chr22
region against a truth VCF and asserts fewer than 8 genotype differences, using
`test/call/HGSVC_chr22_17119590_17880307.gam` — **a real GAM fixture already present in the repo**,
already available at the point the pack file is built. So the `--gam` path needs no new test data.

**The de novo harness now exists** as
[benedictpaten/vg-call-eval](https://github.com/benedictpaten/vg-call-eval), outside this tree so its
dependencies stay out of vg's build. **It has since been run on real data — see the tier-2 outcome at
the end of this document, which supersedes the simulated signal described next.** First measured signal, simulated, 4x, 100 bp reads: the
read-likelihood caller reached recall 0.887 against the Poisson caller's 0.865, at marginally lower
precision. **That demonstrates the harness can measure; it is not a result** — one seed, 141 truth
variants, a 3-true-positive difference, on simulated reads that are unrealistically easy to map. At
20 kb and 20x every arm scores 1.0, so the easy regime cannot discriminate at all, and a
configuration that gives everything a perfect score is a statement about the configuration.

**⚠ The in-tree harness only covers the `-v` path.** That test runs `VCFGenotyper`, and it is the *only* truth-based
concordance assertion in the file — the de novo `FlowCaller` path has self-consistency assertions only
(line 121: the same 0/0 count for `-a` as for `-v`; lines 124-131: GBWT-vs-direct identity). So reusing
it validates precisely the path §2.3 says gets read likelihoods for free, and **not the default
caller**, which is what stage 4 is meant to judge. A de novo truth-concordance harness needs tooling not
currently in the test deps — truvari or vcfeval against the same truth VCF, restricted to the region.
That is stage 3b.

One integration consequence to expect at stage 3: the pure model reports no-call rather than hom-ref
where no informative read exists (§4.6 risk 1), so the new caller's 0/0 record count differs from the
Poisson caller's. That is harmless here only because the feature is opt-in — the line-121 assertion
exercises the default path, which is unchanged. It does mean the read-likelihood variant of the test
cannot reuse the same thresholds unmodified.

- Extend `18_vg_call.t` with a read-likelihood variant of the SV-genotyping test; assert concordance
  no worse than the support-based caller.
- Add `src/unittest/allele_likelihood.cpp` (Catch2, per `CLAUDE.md`) with hand-built matrices, so the
  genotyping maths is tested independently of scoring.
- **And `src/unittest/allele_likelihood_scoring.cpp` for the scoring itself**, on hand-built graphs and
  alignments. This layer is not optional: hand-built matrices cannot see a wrong matrix *builder*, and
  that is exactly where both bugs in §5.5 lived. Cases: a reverse-strand read must score identically to
  its forward equivalent; a deletion-spanning read must be kept and prefer the deletion allele; a read
  inside a single boundary node must still be dropped; every allele must be scored over the same span of
  read bases. **Confirm each regression test fails before its fix** — a regression test nobody has
  watched fail is not evidence.
- **Eyeball the `GQ` distribution** at stage 4. Nearly free, and it is the trigger for the §4.6 risk 2:
  if `GQ`s pile up implausibly high, especially as depth rises, that is the read-independence problem
  showing itself and the `read_weight` discount becomes worth fitting. Also compare against the
  mismapping term switched off, which is the trigger for risk 3.

**Scoring unit tests** — on hand-built graphs with hand-built alignments:

| Case | Expected |
|---|---|
| read path == allele path, no edits | max score; exact-match fast path |
| 1 mismatch at a **low-quality** base | small penalty, distinctly better than the same mismatch at high quality |
| read spans a SNP, scored against both alleles | mismatch on one, exact on the other |
| read spans a deletion junction | good on the deletion, gap-penalised on the reference |
| read placeable over 150 bases on the reference but only 100 on a deletion allele | **both scored over the same 150-base window** (§7); the 50 unplaceable bases charged as insertions on the deletion allele, never omitted |
| **read-base count entering `score(r, a)`** | identical for every *a* — invariant of §7, assert it per read |
| read covers **only deleted bases** | reference good; deletion gets `rel = 0` bounded by the mismap term — **not** dropped, **not** an equal entry |
| overlap where *a* == *b* but *c* differs | `rel[r][a] == rel[r][b]`, both ≠ `rel[r][c]` |
| read overlaps only boundary nodes | dropped as uninformative |
| two reads with very different window lengths | both usable; genotype ranking unaffected by the length difference |
| same read at MAPQ 60 vs MAPQ 3 | penalty shrinks with MAPQ; ranking never flips *direction* |
| **every row's max is exactly 1** | invariant of §4.2 — assert it on every matrix built |
| read placing on **no** allele | dropped, and counted; no `NaN` reaches the genotyper (§7) |
| per-read term for any genotype | lies in `[ln e_r, 0]` — finite, never `-inf`/`NaN` (§4.3) |

The last six are the ones most likely to be wrong in a way no end-to-end test would catch. The final
three, plus the read-base count row, are cheap invariants worth asserting in the code as well as
testing.

**Cross-check against DP.** Once the optional realigning calculator exists, run it over the same site
and confirm both agree on which allele each read prefers, even where absolute scores differ. That is
the cheapest available evidence that graph-implied scoring does what we think.

**Default-path regression test.** Because the feature is opt-in (§1, §7), the most valuable single test
is the cheapest: run the existing `vg call` invocations with no new flags and diff the VCF against the
pre-change build. It should be byte-identical. Worth adding to `18_vg_call.t` as an explicit assertion
rather than relying on the existing tests to notice — they assert concordance thresholds, which would
tolerate a small unintended drift.

---

## 10. Decisions

All settled. The consistent theme is **the simplest defensible option now, with the more complex
one kept additive** — every entry in the right-hand column can be added later without changing the
model's shape, the interfaces, or the default path.

| # | Question | Decision | Later, if measurements justify it |
|---|---|---|---|
| 1 | Read source | in-memory prototype (§6.1) → `.gai` GAM (§6.2), which is where the memory win came from; both dependency-free | gaf-base (§6.4), built and equivalent but not faster than `.gai`; tabix GAF (§6.3) if a user needs it |
| 2 | MAPQ calibration | assume reasonably calibrated; only the `e_max` clamp and a term-off flag (§4.4) | 61-entry MAPQ→`e_r` table; `read_weight` discount |
| 3 | Depth term | **pure `P(reads \| G)`**, no coverage term (§4.6) | optional depth factor |
| 4 | Soft clips | **exclude by default**, matching vg precedent, with a flag to include (§7) | model them as SV evidence |
| 5 | Prior | **uniform** over genotypes (§4.5) | HWE / het-biased prior |
| 6 | Scope | **opt-in `--read-likelihood`**; Poisson stays the default (§1, §7). One exception: the nested effective-ploidy fix touches shared code and is **not** opt-in (§7) | promote to default only if it measures better *and* the §4.6 risks are resolved |
| 7 | Nesting | **independent (`-A`) and top-down both supported** (§7); no joint model over the snarl tree | read-backed phasing across nesting levels |

Decisions 3–5 in particular are chosen to keep things *separable*: no depth term means the read evidence
can be evaluated on its own; a uniform prior means the likelihood is not confounded by a prior;
excluding soft clips means we lose signal rather than invent it. Each is a place where adding complexity
later will be measurable against a clean baseline — which is much harder if the complexity is there from
the start.

### Remaining sub-questions

Minor, none blocking:

- **How to describe `GQ`/`GL`** in the docs and VCF header, given decision 2. Recommend stating plainly
  that they are for ranking, not calibrated probabilities, and will be over-confident at high depth
  (§4.6 risk 2).
- **Whether the tabix GAF backend is worth wiring at all** (§6.3). Dependency-free and already in vg,
  but it returns line strings and carries a memory hazard. Recommend skipping unless a concrete user has
  that input.
- **Whether promoting this to the default is a goal**, or the two callers coexist indefinitely. Nothing
  in the design depends on the answer; it only affects how much of §4.6 needs resolving first.

---

## Tier-2 outcome (real data), summarised

Two chromosomes (chr20, chr6) × two graphs (4-haplotype, 34-haplotype) × five arms, HG002 against
the GIAB defrabb V0.019 draft benchmark on CHM13v2.0. Small variants scored by aardvark, structural
variants by truvari. Full tables in the harness repo's `docs/`; the derivation of every number below
is in [vg-call-eval-plan.md](vg-call-eval-plan.md) Appendix A.

**The model earns its place on real data, and by more on the richer graph.** Genotype F1:

| | chr20 4-hap | chr20 34-hap | chr6 4-hap | chr6 34-hap |
|---|---|---|---|---|
| `poisson-z` | 0.9359 | 0.9124 | 0.9466 | 0.9318 |
| `readlik-z` | **0.9490** | **0.9547** | **0.9586** | **0.9616** |

It leads on every class on every dataset, and its margin **widens** from +0.012 to +0.042 as the
graph goes from 4 to 34 haplotypes — the Poisson caller loses ground there while this one gains.
Tier 0 had the model *losing*; real reads reverse it, and the reversal is not marginal.

**The MAPQ mismapping term is vindicated.** `--no-mismap-term` is worse on every variant class and
both comparison types, and it does roughly ten times more work on the richer graph. The tier-0
anomaly that had `readlik-nomismap` ahead on BASEPAIR is best explained as an artefact of reads
simulated from the graph they are mapped back to, where MAPQ is uninformative.

**Structural variants are the weak class, and the deficit is real.** Against truvari, `readlik-z`
trails `poisson-z` on SV F1 on three of four datasets (chr20 4-hap 0.4824 vs 0.4930; chr6 4-hap
0.5349 vs 0.5478), winning narrowly only on chr20 34-hap. Both callers recover roughly half of true
SVs at roughly half precision. Two things are established about this: it is **not** a metric artefact
(truvari and aardvark agree on the direction), and it is largely an **exposure** effect rather than a
regression — the richer graph offers multi-allelic sites the sparse one cannot produce, and those
sites are harder. Enumeration matters far more than the genotyper here.

**Two defects were traced and closed as non-defects.** The insertion BASEPAIR gap is benchmark scope,
not a scoring bias (§5.3.2) — the `smvar` truth set holds no record ≥50 bp, so large insertions
inside its confident region are unscoreable-as-correct; size-matched, the gap closes and the sign
flips. And the scorer was cleared of favouring long alleles by three synthetic tests. What remains
genuinely wrong is the impossible-depth pericentromeric pile-ups of §5.3.3, which no benchmark
scores.

**Both clamp defaults were set by measurement, and one was actively wrong.** `--mismap-max` at its
original 0.1 was overriding the mapper on haplotype-rich graphs, where 23.3% of reads at the failing
sites sit at MAPQ 1 (p(wrong) = 0.79) and were being told 0.1; the default is now **0.5**. The floor
`--mismap-min` is now **0.02**, re-swept at the corrected cap and scored on both benchmarks. The
lesson worth carrying: a clamp that is inert on a sparse graph is not thereby harmless, and a sweep
that sets a default has to be scored on every benchmark the project runs.

**Cost is no longer a reason not to use this.** After the read-path work (§6.5), `readlik-z` on chr20
runs in 99 s against `poisson-z`'s 75 s, and `readlik` is *faster* than `poisson` (120 s vs 168 s) —
near parity at matched enumeration, against 4.5–9× slower before. Peak RSS 3.6–3.8 GB on chr20 and
9–10 GB on chr6. The windowed read source is what makes that possible: chr20's 13.3 M reads held in
memory would need ~40 GB.

**A pleasant surprise for §2.3:** haplotype enumeration (`-z`, no pack) is as accurate as support
enumeration and faster, so the **pack-free path is not just viable but preferable** here.

**Ranking, as distinct from accuracy.** Everything above is scored at every GQ. Separately, `GQ` is
now scaled by the fraction of reads the called genotype explains, with `GQI` carrying the raw
likelihood ratio — see §4.5. That changes no genotype and so moves none of the numbers on this page;
what it changes is the precision/recall curve, where it removes 44–58% of false small-variant calls
at unchanged recall.
