# Linkage between sites: a PanGenie-style HMM over `vg call`'s existing likelihood

**Status: designed and measured, not built.** The measurement came first deliberately, and it
changes what is worth building: the prize is not the elimination of impossible calls, it is a
modest reweighting of *undecided* sites. That argues for the cheap end of the design.

`vg call` genotypes each snarl independently. The emitted call set is the concatenation of
per-site argmaxes, which corresponds to a pair of haplotypes free to switch panel haplotype
at every site, at no cost. Nothing in the objective notices when consecutive called alleles
are carried by no single panel haplotype. PanGenie's model is the obvious remedy, so the
question is what it would buy here and how the transitions should be parameterised over a
*haplotype-sampled* panel, which is not a population sample.

---

## 1. What PanGenie does

[Ebler et al. 2022](https://www.nature.com/articles/s41588-022-01043-w)
([open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9005351/),
[preprint](https://www.biorxiv.org/content/10.1101/2020.11.11.378133v1.full)).

Hidden states at bubble `v` are **unordered pairs of panel haplotype paths**, `H_{v,i,j}` for
`i,j ≤ N`, giving `N(N+1)/2` states per site. Each state induces a copy-number vector over the
bubble's k-mers — 0, 1 or 2 according to how many of the two paths carry each k-mer.

Emissions are k-mer counts, independent across k-mers: copy number 2 → Poisson(λ), copy number
1 → Poisson(λ/2), copy number 0 → **geometric** (sequencing error), with λ the mean k-mer
coverage.

Transitions are Li–Stephens. With `x` the bp gap, genetic distance `d = x · (1/10⁶) · 4rN_e`,
and

```
p_r = (1 − exp(−d/N)) · (1/N)      switch to a specific other haplotype
q_r = exp(−d/N) + p_r              stay (a switch can land back where it started)
```

with the pair transition multiplying over the two strands. Genotypes come from forward–backward
posteriors summed over states implying the same genotype. Bubbles closer than k = 31 bp are
merged. No pruning; the authors note runtime depends on panel size and that larger panels would
need "additional engineering".

**What transfers and what does not.** The transition structure transfers. The emission does not
need to: `vg call` already computes `ln P(reads | G)` per site from per-read alignment
likelihoods with a mismapping term, a length-weighted mixture and a Poisson depth term — a
strictly richer emission than k-mer counts against a Poisson/geometric mixture. So this is
adding a transition model to a better emission, which is the opposite of the usual porting
direction and means the expected gain is smaller than PanGenie's headline and concentrated
where our emission is flat.

---

## 2. What vg already has

- **`GBWTTraversalFinder::find_gbwt_traversals`** returns, per traversal, the GBWT path
  identifiers taking it. That *is* the haplotype→allele matrix the HMM needs, and `vg call -z`
  computes it at every site already, then discards the identifiers. No new data structure.
- **The GBWT is an FM-index over haplotype paths**: `find`, `extend`, `bdExtendForward/Backward`,
  and `state.size()` in O(1). `locate()` — range to haplotype identities — is the expensive
  operation, and a good encoding never needs it.
- **`deps/sublinear-Li-Stephens`** is vendored and linked, with `penaltySet` referenced from
  `src/haplotypes.{hpp,cpp}`. An LS implementation exists in-tree.
- **`Haplotypes` / `Recombinator`** already partition the graph into top-level chains and
  subchains and score haplotype choices against k-mer counts — structurally a PanGenie-shaped
  DP built for haplotype *sampling*. Its subchain boundaries matter below.

---

## 3. The measurement

`apparent_recombination.py`. `vg deconstruct` on the same graph emits one column per panel
haplotype and uses the **same snarl IDs** as `vg call`, so the panel matrix joins to the call
set on the ID column with no coordinate matching.

For adjacent called sites, take the distinct called alleles at each and ask whether any single
panel haplotype carries one from each.

| | chr20-4hap | chr20-34hap |
|---|---|---|
| panel sites / called records | 130,895 / 104,165 | 165,404 / 105,251 |
| joined on snarl ID | 100,906 (96.9%) | 94,359 (89.7%) |
| panel haplotypes | 3 | 33 |
| adjacent called pairs | 100,905 | 94,358 |
| **no haplotype carries any called combination** | 654 (**0.6%**) | 534 (**0.6%**) |

**Eleven times the panel, the same rate.** So "the call set implies a combination the panel has
never seen" is not where the prize is.

### The statistic that matters is normalised against independence

The first version of this measurement reported co-carrying haplotype counts, and on the
34-haplotype graph the median called pair is co-carried by **17 of 33** — which reads as
"linkage constrains nothing". That reading is wrong: if two common alleles are in *perfect*
linkage and each is carried by half the panel, the intersection is still half the panel. Under
independence the expected intersection would be ~8.75. The count has to be normalised.

Lift = `P(allele at site 2 | allele at site 1) / P(allele at site 2)`:

| | chr20-4hap | chr20-34hap |
|---|---|---|
| median lift | **1.50×** | **1.83×** |
| p25 / p75 | 1.00× / 3.00× | 1.03× / 2.06× |
| above independence | 71.9% | 63.2% |
| at or below | 28.1% | 36.8% |

**The typical adjacent called pair is 1.8× more likely than independence predicts**, and the
34-haplotype panel carries more of that than the 4-haplotype one despite being synthetic —
consistent with within-block linkage being real biology.

### The prize is at undecided sites

chr20-34hap, apparent recombination by the lower `GQ` of the pair:

| min GQ | pairs | apparent recomb |
|---|---|---|
| <10 | 6,454 | **2.8%** |
| 10–40 | 10,575 | 1.3% |
| ≥40 | 77,329 | **0.3%** |

**Nine times more where the reads were undecided.** That is the right shape: a prior would act
where the emission is flat rather than fight it where it is confident. The addressable
population is ~6,500 low-`GQ` pairs per chromosome — small against 94,000, but these are sites
we currently decide close to arbitrarily.

### Distance profile

chr20-34hap: 0.6% (<200 bp), 0.5% (200 bp–1 kb), 0.6% (1–5 kb), 0.9% (5–20 kb), 2.6% (>20 kb).

It rises, so this is linkage decaying rather than noise. But **it rises far too slowly to be
distance-driven**: a 20× wider gap gives 1.8× the rate, where a pure `1 − exp(−x/L)` term would
give ~20×. Most apparent recombination is distance-*independent*. That is a direct argument
against leaning on genetic distance, and it is the main thing the measurement contributes to
the parameterisation below.

### Two caveats on the numbers

- `vg deconstruct` gives the reference sample (CHM13) no column, so a haplotype carrying the
  reference allele is invisible in the panel matrix. That undercounts support and therefore
  makes 0.6% an **upper** bound on apparent recombination.
- 10.3% of called records on the 34-haplotype graph did not join to a panel snarl, so this
  covers 89.7% of the call set.

---

## 4. The model

### Emission: unchanged

`ln P(reads | G)` exactly as it is today — §4.0 of the design doc. A state `(i,j)` implies the
genotype `(a_i, a_j)` at that site, and the emission is the existing per-site likelihood of that
genotype. Nothing is rescaled; this is why the whole thing is tractable.

### State: unordered pairs of panel haplotypes

`N(N+1)/2` states. For N = 33 that is 561 per site, ~100k sites, forward–backward ≈ 10⁸
operations per chromosome. Fine.

**Not GBWT search states, at least not first.** They are more elegant — the state is a BWT range
and a recombination is literally the operation of widening back to `find(node)`, so the model
conditions on the whole matched prefix rather than a first-order chain, and states merge
automatically. That is what makes LS sublinear in panel size. But the measured signal is a 1.8×
lift, which does not justify managing range merging and restart semantics. Hold it for
HPRC-scale panels where `N²` stops being free.

### Transitions, in three layers

Per strand, per site, a switch probability `ρ_t`; the diploid transition multiplies over the two
strands with PanGenie's convention, `q = (1−ρ) + ρ/N` to stay and `p = ρ/N` for a specific other.

```
ρ_t  =  ρ_min  +  (1 − ρ_min) · (1 − exp(−x_t / L))        distance
ρ_t  ←  1 − (1 − ρ_t) · (1 − β · b_t)                      sampling block boundary
```

`x_t` is the bp gap; `b_t` is 1 if a haplotype-sampling block boundary falls between the sites.

**Why the boundary term is separate rather than folded into `L`.** The panel in a sampled GBZ is
a set of recombinants of real assemblies, chosen in blocks. Within a block a panel haplotype
*is* a contiguous piece of a real assembly, so linkage inside a block is genuine. At a boundary
the sampler continues the same haplotype about **43%** of the time, so ~57% of boundaries carry
a switch that is an artefact of construction. With ~10 kb blocks and ~1 kb site spacing, roughly
one adjacent pair in ten spans a boundary, so about 6% of pairs carry an artefactual switch and
~94% carry real linkage.

Modelling that as a smooth function of genetic distance would fit a constant to a step function:
over-penalising switches inside blocks, where linkage is real, and under-penalising them at
boundaries, where a switch cost the panel nothing. **The boundaries are known** — they are the
subchain structure `Haplotypes`/`Recombinator` already builds and stores in the `.hapl` file. So
this replaces a fudge factor with a lookup.

### Parameters

| | source | value |
|---|---|---|
| `β` boundary switch probability | **known from construction**, 1 − 0.43 | ≈ 0.57 |
| `L` distance scale | fitted by grid search on accuracy | start 10–50 kb |
| `ρ_min` floor | fixed; a switch must never be impossible | ~10⁻³ |
| `ε` off-panel escape | fixed, see below | small |
| `w_t` transition weight | **searched, default 0** | — |

**No `4rN_e`.** Over a panel whose members are synthetic recombinants, effective population size
has no meaning; the number it produces would be doing `L`'s job while implying a provenance it
does not have.

**`L` fitted against accuracy, not against the measurement above.** The distance profile
conflates real switches with our own call errors, so fitting to it would partly fit our
mistakes. What it legitimately says is that `ρ_min` will do more work than `L`.

### `w_t`, and why it defaults to 0

Add `w_t · ln P(h_t | h_{t−1})` to the objective, exactly as `w_d` was added for depth, and for
the same three reasons: at 0 it reproduces the current caller bit for bit; it makes "how far do
we trust linkage against reads" an explicit searched quantity rather than an implicit
consequence of the transition formula; and it can be searched on chr20 and validated on chr6
like every other parameter here.

Expected useful range: small. Enough to break ties at undecided sites, not enough to overturn
confident ones. **The harm metric is the fraction of high-`GQ` genotypes that change**, and it
should be near zero at any `w_t` worth shipping.

### The off-panel escape, which is mandatory

A state `(i,j)` implies genotype `(a_i, a_j)`. If the sample carries an allele no panel
haplotype carries, **no state implies it and the model cannot call it.** The graph contains no
HG002, and the 34-haplotype analysis found that most of its extra false positives have no truth
candidate at all — so this is not a corner case.

So: a wildcard haplotype, entered with probability `ε`, able to carry any allele at any site at
a flat cost. Without it the HMM would systematically suppress novel alleles, which is the
opposite of what the SV recall work has been trying to achieve — and it would be easy to miss,
because it would present as a precision improvement.

### Chain structure

Transitions need a linear order with meaningful gaps: top-level chains from the snarl distance
index, resetting to the stationary prior at chain boundaries, contig boundaries, and any gap
beyond a few `L`. `x_t` comes from reference positions, which `vg call` already has.

---

## 5. What this costs downstream

Genotypes would come from forward–backward posteriors, so **`GQ`, `GL`, the explained-share
discount and the depth discount all need restating**: they are currently defined against the
ratio of the top two *per-site* genotype likelihoods. That is real bookkeeping and should be
costed before starting, not discovered halfway.

Against that, forward–backward yields a **phased** genotype per site as a by-product, which
`vg call` does not currently emit.

---

## 6. What would make this worth building

The measurement sets a modest ceiling, so the decision rule should be set in advance:

- **Go** if a small `w_t` improves genotype concordance at low-`GQ` sites on chr20 and holds on
  chr6, while changing under ~0.1% of high-`GQ` genotypes.
- **Stop** if the gain requires a `w_t` large enough to move confident calls, or if it is
  confined to the 4-haplotype graph — the sampled panel is where the linkage is supposed to be
  real, so a 34-haplotype-only failure would mean the model is fitting construction artefacts.
- **Stop** if the off-panel escape has to be tuned to avoid suppressing novel alleles. That
  would mean the prior is fighting recall, which is the thing this project has spent most of
  its effort recovering.

## 7. Not measured, and worth knowing before committing

- **How often HG002's true genotypes imply panel switching.** That bounds how much a linkage
  prior can help without hurting, and it is the single most useful missing number. It needs
  truth genotypes expressed in graph alleles, which is not a join we currently do.
- **The minimum number of switches needed to explain the whole call set** as two paths through
  the panel — a Viterbi parse rather than a pairwise count. The pairwise measure here is a
  lower bound on disagreement, because it cannot see a switch that is locally consistent but
  globally impossible.
- **Whether the 0.6% apparent recombinations are enriched for errors.** If they are truth-set
  false positives, the prior would fix real mistakes; if they are true novel haplotypes, it
  would suppress correct calls. Joining them to the truth would say which, and it decides
  whether this is a fix or a regression.
