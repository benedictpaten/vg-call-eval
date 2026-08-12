# Linkage between sites: a PanGenie-style HMM over `vg call`'s existing likelihood

**Status: built, searched on two chromosomes and two graphs, awaiting the full-matrix refresh
before the default flips.** Shipped as `--linkage-weight` (opt-in, 0 by default) and
`--linkage-freq-prior` (5). Read §7 for the staged results, and **Stage 2b first** — it is where
the conclusions changed most.

The measurement came first deliberately, and it changed what is worth building: the prize is not
the elimination of impossible calls, it is a modest reweighting of *undecided* sites. That
argues for the cheap end of the design — and, in the event, understated the prize by a factor of
two, because the parameter that turned out to dominate was excluded from the search by an
argument in this document that did not survive contact with the data.

Three things this document asserted and the implementation disproved, all recorded in place
rather than edited out: the block-switch term `β` was the distance scale under another name and
is deleted; the frequency prior `f` was capped at 1 by a guard that read as an optimisation, and
uncapped it is worth more than the transition model; and the forward–backward itself accounts for
about 12% of the gain, the rest being reachable per-site with no chain at all.

**Read §7 first if you are deciding whether to proceed.** Stage 0 tests the entire model offline
with no changes to vg, because `vg call` already emits per-genotype likelihoods (`GL`) and
`vg deconstruct` supplies the panel matrix on matching snarl IDs. Roughly a day, no calling runs,
and it can sweep every parameter. Everything after it is gated on that result.

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
operations per chromosome. Fine — though see §5 on windowing, which matters for parallelism
rather than for raw cost.

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

> The second line was implemented, measured and **removed**; only the first survives. Note that
> it is written here with `b_t` an *indicator* — the localised form, which would have been a real
> second mechanism. What shipped had no boundary positions to set `b_t` from, so it became
> `gap / L_block`, and at that point it is the first line again with a shorter `L`. See
> [the Stage 2 result](#stage-2-result-w--2-and-the-boundary-term-deleted).

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

> **`β` was removed after Stage 2 — see [that result](#stage-2-result-w--2-and-the-boundary-term-deleted).**
> Everything in this section about it is the design as written, kept for the record. Two things
> broke it: smeared over `gap / block_length` it is *algebraically* the distance scale rather than
> a second parameter, and the panel linkage break it exists to model is +0.008 NMI (z = 1.1) at
> the gaps that matter. `L` and `w_t` below are unaffected.

| | source | value |
|---|---|---|
| `β` boundary switch probability | **user-specifiable**, estimated per graph, default 0 | ≈0.57 for this sampled panel |
| `L` distance scale | estimated from the panel, overridable | ~10 kb here |
| `ρ_min` floor | fixed; a switch must never be impossible | ~10⁻³ |
| `ε` off-panel escape | fixed, see below | small |
| `w_t` transition weight | **searched, default 0** | — |

**`β` must be a user parameter with a per-graph default, not a constant.** It is a property of
how a particular GBZ was built: the block size and the sampler's continuation rate both vary
between graphs, and a full HPRC pangenome has no sampling blocks at all, so `β = 0` there. Any
fixed value would be wrong on some graph, silently. It should also be explorable, on the same
footing as `--depth-term` and `--mismap-max`, because there is no reason to believe the
construction's nominal rate is the one that genotypes best.

**No `4rN_e`.** Over a panel whose members are synthetic recombinants, effective population size
has no meaning; the number it produces would be doing `L`'s job while implying a provenance it
does not have.

### Both parameters are estimable from the panel alone

`vg deconstruct` output is enough — no reads, no truth, no metadata. Median lift against
distance:

| distance | chr20-4hap | chr20-34hap |
|---|---|---|
| 0–300 bp | 3.00 | 2.36 |
| 300 bp–1 kb | 3.00 | 1.94 |
| 1–3 kb | 1.50 | 1.94 |
| 3–10 kb | 1.50 | **1.81** |
| 10–30 kb | 1.50 | **1.50** |
| 30–100 kb | 1.50 | 1.27 |
| 100–300 kb | 1.50 | 1.24 |

The 34-haplotype curve is **flat from 300 bp to 10 kb and then breaks** — the knee sits at the
stated block size. That is the signature of a sampled panel: inside a block a haplotype is a
contiguous piece of one assembly, so linkage is near-constant; past the block scale, identity
has been reshuffled. So the block scale can be read off the panel rather than taken on trust,
and the same procedure applied to a real pangenome measures genuine LD decay instead. **The
fitting method does not need to know which kind of graph it has.**

The 4-haplotype curve is useless for this and shows why a thin panel cannot be calibrated: with
3 haplotypes the lift is quantised to 3.00 and 1.50 and never decays, even at 300 kb.

**`β` and block length are not separately identifiable** from an aggregate curve, since survival
goes as `(1 − β)^(x / L_block)`; only the product — the effective per-bp switch rate — is. That
is sufficient, because the effective rate is what the transition consumes. Knowing `β`
separately only helps when the boundary *positions* are known, in which case the switch mass can
be placed where it belongs instead of smeared uniformly.

> **This paragraph is where the error entered, and it is left standing because the tell is in it.**
> The non-identifiability was noticed and then filed as harmless. It was not harmless: it was the
> whole finding, one step short. If survival goes as `(1 − β)^(x/L_block)` and only the per-bp
> product is identifiable, then the term is `exp(−x · rate)` — which is what the distance term
> already is, so the two compose into a single scale and `β` adds no freedom whatever. "Only the
> product is identifiable" and "this is a reparameterisation of `L`" are the same statement; the
> first was written down and the second was not. The design then went on to build a grid crossing
> the two.

**Fit against accuracy too, not only against the panel.** The panel curve calibrates the
*panel's* structure; how far to trust it against read evidence is `w_t`'s job, and the
apparent-recombination-versus-distance profile in §3 conflates real switches with our own call
errors, so fitting to that alone would partly fit our mistakes.

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

### Three kinds of discontinuity, and only one of them is about sampling

Checking the GBWT metadata corrected an assumption. Path names are
`sample#phase#contig#fragment`, and on chr20 the 4-haplotype graph stores 4 haplotypes as **16
paths**, the 34-haplotype graph 34 as **137**:

```
recombination#1#chr20#0   recombination#1#chr20#1   recombination#1#chr20#2
GRCh38#0#chr20[60296]     GRCh38#0#chr20[28728974]  …   (9 fragments)
```

**GRCh38 — a real assembly, not a sampled recombinant — is in 9 fragments.** So continuity is
not a sampled-graph question, and the model has to distinguish:

1. **Fragment boundaries.** The haplotype path ends. Readable from `PathName`'s fragment field
   **in any graph**. There is no linkage across it whatsoever, so the haplotype must be treated
   as *absent* rather than as having switched — a different transition, not a more expensive
   one. This applies to full pangenomes too.
2. **Sampling block boundaries.** The path is contiguous but its identity changed. Sampled
   graphs only, and invisible in the path itself: detectable from the `.hapl` subchain
   structure, from the `recombination` sample name (vg's documented convention, see
   `recombinator.cpp`), or from the empirical knee above.
3. **Biological recombination** in the sample's own mosaic — what `L` and `ρ_min` model.

Conflating 1 with 2 would be the worst of the three errors: it would charge a switch penalty for
crossing a place where the panel simply has no information, and then let the HMM carry a
haplotype identity across a gap where none exists.

### Chain structure

Transitions need a linear order with meaningful gaps: top-level chains from the snarl distance
index, resetting to the stationary prior at chain boundaries, contig boundaries, fragment
boundaries (above), and any gap beyond a few `L`. `x_t` comes from reference positions, which
`vg call` already has.

---

## 5. Inference: forward–backward marginals, not Viterbi

**Forward–backward, and the argmax taken over genotypes rather than over states.** Those are
two different things and the distinction matters:

```
posterior over states     γ_v(H_{v,i,j})  =  α_v(H) · β_v(H) / Z
posterior over genotypes  P(g | reads)    =  Σ over states implying g of γ_v(H)
call                      ĝ_v             =  argmax_g P(g | reads)
```

Taking the most likely *state* and reading its genotype off is not the same as taking the most
likely *genotype*: a genotype reachable by many moderately-likely haplotype pairs can beat one
reachable by a single very likely pair. Since the panel is redundant — a called allele is
typically carried by 17 of 33 haplotypes — many states imply the same genotype, so the two
answers genuinely differ. Summing is the correct one, and it is what PanGenie does.

**Why not Viterbi.** Viterbi maximises the joint probability of the whole state path, which is
the wrong loss for how we are scored: aardvark and truvari score records, not haplotypes, so
what we want is to maximise expected per-site accuracy — exactly what marginal posteriors give.
Viterbi also yields no per-site confidence, and with an approximate transition model plus a
wildcard escape state, a single path is more brittle than a marginal. It would also be the
wrong tool for the population the measurement identified: undecided sites, where the point is
that the posterior is genuinely spread and we want to know it.

**Correcting an earlier claim in this note: forward–backward does *not* give phasing.** Marginal
posteriors are per-site and unphased; the per-site argmaxes need not lie on any single
high-probability path. Coherent phasing is what *Viterbi* gives, because a state path assigns
alleles to strands consistently across sites. So phasing is available, but as a **separate
Viterbi pass over the same model**, not as a by-product of the genotyping pass — and it should
be reported as such rather than as free.

### How this maps onto the existing quality fields

Better than I expected, and it follows the precedent already set by the explained-share
discount:

```
GQI  =  the per-site value, transitions off       (what GQ is today)
GQ   =  -10 log10 (1 - P(ĝ_v | reads))            from the posterior
```

`GQI` already means "the quality before the discount", and keeping it as the transition-free
per-site value preserves that meaning exactly while making the before/after pair directly
comparable — the same structure that let the share discount be measured end to end. It also
*improves* the semantics: the current `GQ` is the ratio of the top two genotypes and ignores
everything else at the site, whereas a posterior accounts for the whole distribution.

`GL` becomes the per-genotype posterior rather than the per-site likelihood. The explained-share
and depth discounts still multiply `GQ` and remain ranking-only.

### Windowed, not chain-wide

**`vg call` is parallel over snarls; a chain-wide forward–backward serialises within a chain**,
and top-level chains here are whole chromosome arms. That is a real cost and is an argument for
a bounded window rather than exact global inference.

The measurement says a window is nearly free of approximation error: lift decays to ~1.50 by
10–30 kb and ~1.25 by 100 kb, so sites beyond a few tens of kilobases contribute almost nothing
to the posterior at a given site. A window of ±50 kb, or a fixed number of neighbouring sites,
captures essentially all of the available signal while keeping memory bounded, preserving
parallelism over windows, and leaving the streaming structure intact. Exact chain-wide inference
should be the thing we compare against on one chromosome, not the thing we ship.

---

## 6. What this costs downstream

`GQ`, `GL` and the two ranking discounts all need restating as above. That is real bookkeeping
and should be costed before starting, not discovered halfway — and it means the tier-2 quality
analysis (`share_gq.py`, `depth_gq.py`) has to be re-run, since both are defined against the
current per-site ratio.

### Measured, after the fact — and not where this section expected

`resolve()` is now instrumented, so the estimates above can be replaced. chr20-34hap, 5 threads:

```
105251 sites, 7.87 MB retained, 5736 genotypes changed, 1.80 s
136.7 s with linkage against 125.1 s without   ->  +9.3%
```

Two corrections to what was assumed. The memory figure was arithmetic — sites times an assumed
per-site size — and `bytes()` was written so it could be an observation but was never called by
anything; at 7.87 MB against an estimated ~8 MB the arithmetic was good, and it is now checked.

The runtime split is the surprise. This document worried about the serial phase-two pass in an
otherwise-parallel caller, and that pass is **1.8 s of the 11.6 s** the feature costs. The other
~10 s is **phase one** — recording each site under a global mutex while calling is parallel over
snarls — which was dismissed in the class comment as obviously cheap. Cost lives where it was not
looked for, which is the usual place. That is the target for profiling, not the HMM.

---

## 7. Implementation plan

Staged so that the expensive work is gated on evidence, following the pattern that worked for
the depth term: predict offline, then build, then search, then refresh.

### Stage 0 — the whole model offline, with no changes to vg

**This is the stage that decides everything, and it needs no C++.** `vg call` already emits
everything the HMM consumes:

| needed | already available |
|---|---|
| per-site `ln P(reads | G)` for every genotype | **`GL`**, log10 P(reads \| genotype), per record |
| haplotype → allele matrix | `vg deconstruct`, joined on snarl ID |
| site order and gaps | `POS` |
| allele identity across the two files | `AT` (allele traversal) |

So Stage 0 is a Python forward–backward over the emitted `GL` values and the panel matrix,
re-deciding each genotype and writing a modified VCF, then scoring it with the existing
`score_vcf.py` on both benchmarks. It can sweep `w_t`, `β` and `L` freely because no run of
`vg call` is involved — the emission is fixed and cached.

Two approximations to state, both of which only *understate* the model:

- `GL` covers the alleles that reached the record, not every allele enumerated at the site. At a
  site where many alleles were enumerated and few emitted, the HMM sees fewer states than the
  real implementation would.
- Multi-allelic records that were split or merged downstream may not map cleanly onto panel
  allele indices; those sites get skipped and counted.

**Kill criterion.** If Stage 0 does not improve genotype concordance at low-`GQI` sites on
chr20 at any `(w_t, β, L)`, stop. No vg work happens. Cost: ~1 day, no runs beyond scoring.

### Stage 0 result: passed, and it found a defect in the design

`linkage_hmm_offline.py`, chr20-34hap. The baseline is the caller's own output scored through
`score_vcf.py`, which reproduces the refresh **exactly** (0.9546 / 0.4655), so these deltas are
not a scoring-path artefact:

| arm | small-var GT F1 | SV F1 | `GQI≥40` changed |
|---|---|---|---|
| current caller | 0.9546 | 0.4655 | — |
| **frequency prior only** (`w_t = 0`) | 0.9563 (+0.0017) | 0.4676 (+0.0021) | 0.029% |
| + linkage `w_t = 1` | 0.9570 (+0.0024) | 0.4685 (+0.0030) | 0.032% |
| + linkage `w_t = 2` | 0.9575 (**+0.0029**) | 0.4697 (**+0.0042**) | 0.062% |

Small but real, inside the 0.1% harm budget, and **roughly half the gain is a frequency prior and
half is linkage**.

**The inertness check failed, and that is how the frequency prior was found.** At `w_t = 0` the
transitions go uniform and the posterior should collapse to the emission, recovering the current
caller. Instead 40% of `GQI < 10` genotypes changed. Summing the state posterior over the states
implying each genotype weights that genotype by its **multiplicity** — how many haplotype pairs
spell it — and that multiplicity is a panel allele-frequency prior. The state space bundles two
effects the note had treated as one, so Stage 0 needs **three** arms, not two, and `--freq-prior`
exists to separate them.

**The frequency prior defaults to off because its size is unmeasured here, not because it is
wrong.** An earlier draft of this note called it invalid; that was too strong.

There is some double counting: haplotype sampling chose the panel using k-mers from *these reads*,
so panel allele frequency is already conditioned on the sample's data. But the two are not the
same statistic. Sampling works from k-mer counts aggregated over ~10 kb subchains; the genotyper
works from per-read alignment likelihoods at one site. At an undecided site inside a
well-determined block, the sampler's choice therefore carries information the site's own reads
genuinely lack — that is information transfer, and a coarse version of the linkage this whole note
is about. Double counting's usual penalty is also miscalibration rather than wrong answers, and
`GQ` is already documented as a ranking score rather than a posterior.

The caveat worth keeping is narrower, and it is about not knowing rather than about being wrong:
with a sampled panel the bias cannot be bounded, and the failure mode is **correlated** — where
sampling chose the wrong haplotypes the prior compounds the error, and the site's own weak reads
cannot overturn it, concentrated exactly where this model is meant to help. So: measure before
enabling. Separating it needs a panel chosen independently of these reads — a full pangenome, or
sampling from a held-out read set.

PanGenie does not face this, since its panel is the full assembly set. The question is created by
personalising the graph first, not by PanGenie's model.

**This cannot be settled with the graphs in this harness**, because both were sampled with these
reads. It needs a non-sampled pangenome for the same contig, or a panel sampled from a held-out
read set. Until then: `--freq-prior` defaults to 0, documented as defensible only on a full panel.
The linkage half is a statement about co-occurrence rather than frequency and stands on its own.

**A second defect, which unit test 4 was written for.** The wildcard state carries posterior mass
but implies no *specific* genotype, so it is excluded from the genotype sum — meaning a genotype no
panel haplotype pair can spell is unreachable, and 53 of 75,647 records could not be reproduced
even at `(w_freq = 0, w_t = 0)`. The wildcard must contribute per-allele mass. Stage 1 must fix
this, and it is the cheap place to have found it.

### Stage 1 — implement in the caller

Only if Stage 0 pays. Flags, all defaulting to off or to the measured value:

```
--linkage-weight W        w_t; 0 disables and must be bit-for-bit inert   [0]
--linkage-freq-prior F    exponent on the panel allele-frequency prior    [5]
--linkage-scale L         distance scale in bp                           [10000]
--linkage-window N        sites either side; 0 = exact chain-wide        [64]
--linkage-escape E        off-panel wildcard mass ε                      [small]
```

`--linkage-freq-prior` shipped defaulting to 0 and was raised to 5 once the axis was uncapped and
measured (Stage 2b). `--linkage-window` and `--linkage-escape` were never exposed as flags and
remain compile-time constants; of the two, `escape` is the one worth measuring, since it is the
recall floor for alleles no panel haplotype carries.

A `--linkage-block-switch B` flag was shipped alongside these and has since been **removed**. It
defaulted to 0 rather than 0.57 on the argument that a full pangenome has no blocks and that
defaulting to a sampled graph's value would silently mis-model every other graph. That reasoning
was right and the flag was still wrong; see the Stage 2 result.

Where the code goes: the haplotype→allele matrix comes from `find_gbwt_traversals`, which the `-z`
path already calls and whose path identifiers it currently discards. Fragment and block structure
come from `PathName` and, if present, the `.hapl` subchains. The HMM sits above the per-site
likelihood, so `AlleleReadLikelihoods` is untouched — this is a new layer in the caller, not a
change to the emission.

### Stage 2 result: `w = 2`, and the boundary term deleted

Searched on chr20 both graphs, validated on chr6 both graphs, `linkage_grid.py`. Deltas against no
linkage at `w = 2`:

| dataset | SV F1 | small-variant GT F1 |
|---|---|---|
| chr20-34hap | +0.0095 | +0.0047 |
| chr20-4hap | +0.0015 | +0.0022 |
| chr6-34hap | +0.0082 | +0.0036 |
| chr6-4hap | +0.0064 | +0.0019 |

Eight of eight cells positive — a stronger start than `--depth-term` had when it was defaulted on.

**`beta` has been removed from the model.** It lost at every weight, and the first reading here was
that it was starved of boundary positions — that smearing the switch mass as `gap / block_length`
made it a blunter copy of the distance scale, so wiring the `.hapl` subchains was a prerequisite
rather than a refinement. That reading was too generous by half, in both directions.

*It was not a blunter copy of the distance scale; it was exactly the distance scale.* Writing out
the smeared form:

```
1 − ρ' = (1−ρ_min)·exp(−g/scale)·(1−β)^(g/block_length)
       = (1−ρ_min)·exp(−g/scale_eff),   1/scale_eff = 1/scale + −ln(1−β)/block_length
```

`β = 0.57` over a 10 kb block **is** `--linkage-scale 5423`, to floating point. The Stage 2 grid
was crossing one axis with itself. Made `scale` a real axis (it had been pinned at 10 kb through
every experiment to date) and swept it: 10k → 20k → 40k moves chr20-34hap SV F1 by 0.0009 and
chr20-4hap by 0.0017, 20 kb weakly best. So the whole smeared-β family is worth ~0.001, in either
direction, and no value of β was ever going to be found.

*And the premise fails, so the localised form is not worth building either.* `subchain_linkage.py`
asks the prior question directly — does panel linkage actually break at subchain boundaries? —
measuring NMI between adjacent sites in the panel matrix, needing no genotyping. Boundaries from
`vg haplotypes -H` on the sampled chr20 GBZ: 5511 subchains, 12 kb median, 82 seconds.

Two corrections during, the second of which reversed the answer:

- **Coarse distance bins do not control for gap.** Crossing probability rises with gap *inside* a
  bin, so crossing pairs pile up at the top of each bin while NMI falls across it. Permuted
  boundaries produced a *larger* apparent drop below 200 bp (0.090) than real ones (0.059). The
  first verdict — "worth wiring" — was that artefact. Replaced with narrow gap-matched strata,
  after which the permutation control sits at zero, which is the check that the matching worked.
- The surviving `z = 3.9` is 253 pairs in the 5–20 kb bin, whose *non*-crossing members must sit
  inside unusually long subchains and so are a different population.

At the gaps where every adjacent call pair actually sits: **+0.0076 NMI, z = +1.1 under 5 kb**;
+0.0119, z = +1.6 under 2 kb. An offset scan peaks at zero shift, so the recomputed partition is
not badly misphased — which was the one caveat that would have left a null ambiguous.

Mechanism, in hindsight: at a boundary the sampler switches to **another assembly in the same
panel**, and human haplotypes agree at the overwhelming majority of sites. Switching source
changes the allele only where the two assemblies differ, so adjacent-site NMI — dominated by
shared background — barely moves. β assumed a switch destroys linkage; it mostly preserves it.
This is the 43% continuation figure holding up *better* than the model gave it credit for, and it
retires the per-haplotype ρ_h refinement for the same reason.

The flag is gone rather than pinned at 0: a knob that is secretly `--linkage-scale` invites
setting both and believing they are independent.

**The weight was pushed to 6 and structural-variant F1 never turned over**: 0.4801 on chr20-34hap
and 0.5047 on chr20-4hap, both the best SV numbers measured. What turns over is small variants on
the thin panel — GT peaks at `w = 2`, and SNV F1 goes *below baseline* by `w = 4` (0.9755 against
0.9756, and 0.9749 by `w = 6`). SNVs are the largest class, so past `w = 2` the model is trading
small-variant accuracy on thin panels for structural-variant precision. That is a judgement about
what the caller is for, of the same kind as the `--mismap-min` decision, not a maximisation.
`w = 2` is the largest weight where every cell on every dataset improves; `w = 6` is defensible if
structural variants are the only target.

**The gain is precision throughout.** Precision rises on all four datasets (chr20-34hap 0.4380 ->
0.4572, chr6-34hap 0.4716 -> 0.4886) while recall slips on three of four. The model suppresses
implausible calls, and where the panel cannot distinguish it suppresses some real ones too.

**One regression to watch**: chr6-34hap heterozygous deletion recall above 1 kb goes 0.656 -> 0.639.
Confined to that dataset, and small, but it is the class the depth term exists to recover, so it
should not be averaged away.

**On the harm metric, which was misused.** It was written here as a budget that disqualifies a point
regardless of F1. That was wrong: when both benchmarks improve on both graphs, changes at confident
sites are corrections more often than not, and a proxy should not overrule the outcome it stands in
for. Demoted to a diagnostic — and as one it earned its place, crossing 0.1% at `w = 3`, exactly
where the SNV regression begins.

### Stage 2b result: the frequency prior was capped, and it is the dominant parameter

`f` was not swept in Stage 2 because this document argued it should stay at 0 — panel allele
frequency over a read-sampled panel being "the same evidence twice." That argument is wrong.
Nothing connects the panel to the *truth set*: sampling reads k-mer counts, and the benchmark
informed neither the graph nor the selection. Reusing one's own reads is what mapping and calling
already do. What survives is a calibration concern — double counting can leave `GQ` overconfident
while the genotype improves — which is a separate question from whether to switch it on.

Worse, the axis had a **ceiling nobody had noticed**. Both application sites guarded on
`freq_prior < 1.0` before dividing by `multiplicity^(1-f)`. At exactly 1 the exponent is 0 and the
division is a no-op, so the guard looked free; it also made every larger value behave as 1,
silently, with no parse error and byte-identical output. Above 1 the exponent goes negative and
the prior is *amplified* past the state space's own multiplicity — a real setting, unreachable.
Removing the guard leaves `f ≤ 1` bit-identical, verified by md5 across the change.

Uncapped, `f` dominates `w_t`. Crossed on chr20-34hap, validated on chr6-34hap; best at
`(w_t = 2, f = 5)`, against no linkage:

| | chr20-34hap | chr6-34hap |
|---|---|---|
| small DEL F1 | **+0.0434** | **+0.0328** |
| small INS F1 | +0.0290 | +0.0255 |
| small ALL GT | +0.0099 | +0.0074 |
| SNV F1 | +0.0017 | +0.0010 |
| SV F1 | +0.0289 | +0.0209 |
| SV recall | −0.0065 | −0.0104 |

The transition model alone gave +0.0047 and +0.0036 of small-variant GT, so **most of the gain was
behind the cap**. The effect is almost entirely small indels — deletions most, insertions next,
SNVs flat to within 0.0002 across the whole axis. That is where the emission is flat and the panel
has something to add; SNVs are already settled by the reads. Past `f = 8` it inverts: by 12 the
prior overwhelms the reads, SNV F1 falls below the no-linkage baseline and SV recall collapses
0.4915 → 0.4706.

**A benchmark-resolution error, corrected.** The first reading of the `w_t` sweep led with SV F1
still climbing at `w = 13`. The SV truth set is 765 events on chr20, so one event is 0.0013 of F1
and the `w = 9` vs `w = 13` difference was *half an event*. The small-variant benchmark has 94,691
events and says the surface is flat from `w = 4` and then declines. Decide on the powered
benchmark; the SV column is directionally supportive and cannot locate a peak.

**Ablation: what the HMM itself is worth.** At `w → 0` the transitions go uniform, `γ ∝ E`, and the
posterior collapses to `μ(G)^f · L(G)` — a per-site prior with no chain at all. At `f = 5` that
reaches 0.9633 on chr20 and 0.9680 on chr6, against the full HMM's 0.9645 and 0.9689, from
baselines of 0.9546 and 0.9615. So **forward–backward contributes about 12% of the genotype-F1
gain, on both chromosomes independently**; a per-site prior would get the other 88% with no
windowing, no chain assembly and no retained state.

That is recorded as a fact about the present model, not an argument to delete it — the
forward–backward is kept because later work builds on it (§5's Viterbi phasing, for one). But it
means the honest description of this feature has changed: it began as a linkage model and measures
as an allele-frequency prior with a linkage correction.

*(One claim retracted: the ablation appeared on chr20 to show the HMM protecting recall against a
blunter prior, 0.4824 → 0.4902. On chr6 the same comparison is 0.5346 → 0.5352 — nothing. That was
one chromosome over-read. The 12% figure is what replicates.)*

### Stage 2 — search, on chr20, validate on chr6

`w_t` × `f`, crossed rather than swept — they are substitutes, and the substitution is measured:
`f = 1` versus `f = 0` is worth +0.0072 on small deletions at `w = 2` and only +0.0022 at `w = 13`,
because a tightly-linked chain carries some of the frequency information the prior would supply.
Sweeping either alone finds whichever corner it started nearest.

*(As first written this crossed `w_t` against `β`, on the same interaction argument. That was the
same axis twice — see the Stage 2 result — and `f`, the axis that mattered, was excluded by an
argument this document has since had to withdraw.)*

### Stage 3 — full matrix, then decide the default

`CANARY=1 JOBS=2 refresh_all.sh` with `READLIK_Z_EXTRA`, ~30 minutes. Default flips only if it
holds on all four datasets, as `--depth-term` had to. Two parameters to carry now, not one:
`--linkage-weight 2` on the command line, `--linkage-freq-prior 5` from the default.

### Stage 3 result: `readlik-z` is the best arm on every small-variant class, on all four datasets

Small-variant genotype F1, `readlik-z` against the best of the other four arms:

| dataset | best other | readlik-z | Δ |
|---|---|---|---|
| chr20-4hap | 0.9488 (nomismap) | **0.9507** | +0.0019 |
| chr20-34hap | 0.9513 (readlik) | **0.9645** | +0.0132 |
| chr6-4hap | 0.9583 (readlik) | **0.9602** | +0.0019 |
| chr6-34hap | 0.9588 (readlik) | **0.9689** | +0.0101 |

Clean on SNVs, insertions and deletions separately as well, on all four. Structural variants,
three of four:

| dataset | best other | readlik-z | Δ |
|---|---|---|---|
| chr20-34hap | 0.4592 | **0.4944** | +0.0352 |
| chr6-34hap | 0.4999 | **0.5268** | +0.0269 |
| chr6-4hap | 0.5547 | **0.5691** | +0.0144 |
| chr20-4hap | **0.5034** | 0.5016 | **−0.0018** |

**The chr20-4hap structural-variant loss is recorded rather than rounded away.** It is 0.0018 on a
765-event truth set — about 1.4 events, below what that benchmark resolves — and it is on the
3-haplotype panel, where the frequency prior is measurably inert and the transition model has
almost nothing to link against. That is the expected place for this to be neutral-to-slightly-
negative, which is not the same as it being fine. Anyone tempted to quote "best on all four"
should quote this row too.

The gap is much larger on the 34-haplotype graphs than the 4-haplotype ones, in both benchmarks
and by roughly sevenfold on small variants. That is the model working as designed rather than a
red flag: linkage and multiplicity are both panel-size effects, and a panel of three has neither
to offer.

### Three defects in the harness, found by running it

Worth recording because two of them made a broken run look like a successful one.

- **`refresh_all.sh`'s `JOBS>1` path had never worked on macOS.** It used `wait -n`, which needs
  bash 4.3; macOS ships 3.2, where it fails as an invalid option, and `set -e` turned that into
  an exit *after* the pool had filled — leaving orphaned `vg call` children under a dead parent.
  Replaced with PID collection and a `jobs -rp` poll, which also checks every dataset's status
  rather than whichever one `wait -n` happened to reap.
- **`--linkage-weight` through `READLIK_EXTRA` killed two arms.** It requires `-z`, and
  `READLIK_EXTRA` goes to all three read-likelihood arms. The two support-enumeration arms exited
  immediately. Added `READLIK_Z_EXTRA` for flags that need haplotype enumeration.
- **A dead arm scored as success.** `run_arms.py` logged `FAILED rc=…`, then continued, wrote the
  arm out as zero variants with empty metrics, and the run completed; the first visible symptom
  was a `KeyError` in the page build, forty minutes downstream. Zero variants is now fatal at the
  point of failure.

### Stage 4 — phasing, separately justified

A Viterbi pass over the same model, emitting `GT` with `|` and a `PS` tag. Deliberately last, and
not part of the case for stages 0–3 (§5).

---

## 8. Testing plan

### Unit tests (`src/unittest/`), mirroring the depth term's

1. **`w_t = 0` is inert** — the HMM layer must reproduce the per-site argmax exactly, not
   approximately. Same standard as `--depth-quality`: verified byte-identical on real data.
2. **Linkage decides a flat site.** Two sites; site 2's reads are uninformative between two
   alleles; the panel carries one of them on the same haplotype as site 1's call. Assert site 2
   is called from linkage.
3. **Linkage does not override decisive reads.** The same construction with site 2's reads
   strongly favouring the *unlinked* allele. Assert the call follows the reads at a reasonable
   `w_t`. This is the harm direction and the more important of the pair.
4. **Off-panel alleles remain callable.** An allele carried by no panel haplotype, with reads
   demanding it, must still be called. Guards the regression that would present as a precision
   improvement.
5. **Fragment boundaries block linkage.** A haplotype whose path ends between two sites must
   contribute no linkage across the gap — treated as absent, not as switched.
6. **A certain switch equals a chain reset.** Where the switch probability saturates, the
   transition must carry nothing across, giving the same posterior as starting a fresh chain.
   Written against `β = 1`; now reached by making the sites effectively infinitely far apart,
   which is the only route to `ρ = 1` since `β` was removed. The property is the valuable part
   and it survived the removal unchanged.
7. **Posteriors sum to 1** per site, and **argmax over genotypes ≠ argmax over states** on a
   constructed case, asserting we do the former (§5).
8. **Windowed equals exact** on a short chain, to tolerance, at a window wide enough to cover it.

### In-tree tests (`test/t/18_vg_call.t`)

- The existing invariant that **two read sources produce identical VCFs** must still hold.
- **New: thread count must not change output.** A windowed HMM with parallelism over windows is
  exactly the shape of change that introduces order dependence, and nothing currently asserts
  `-t 1` and `-t 5` agree. This should be added *before* Stage 1, so it is a standing invariant
  rather than a check invented to pass.

### Tier-2 evaluation

Fast tier first (`readlik-z`, four datasets, ~14 min), reporting on both benchmarks plus:

| metric | why |
|---|---|
| genotype concordance at `GQI < 10` | the addressable population; where the gain must appear |
| **fraction of `GQI ≥ 40` genotypes changed** | the harm metric; should be ≲0.1% |
| apparent recombination rate, before and after | direct confirmation the mechanism engaged |
| het fraction, SV recall by size, SNV F1 | the standing regression guards |
| runtime and peak RSS | forward–backward over 561 states is new work |

### Cost

| stage | cost | gated on |
|---|---|---|
| 0 — offline | ~1 day, no vg runs | — |
| 1 — implement | ~2–3 days | Stage 0 improving low-`GQI` concordance |
| 2 — search | ~2 h compute | Stage 1 passing tests |
| 3 — refresh | ~30 min | Stage 2 finding a stable operating point |
| 4 — phasing | ~2 days | its own justification |

---

## 9. What would make this worth building

The measurement sets a modest ceiling, so the decision rule should be set in advance:

- **Met.** `w_t = 2` improves both benchmarks on all four datasets, and holds on chr6, which was
  the validation set. The "under 0.1% of high-`GQ` genotypes" clause was a proxy and is not a
  criterion — see the Stage 2 result.
- **Stop** if the gain requires a `w_t` large enough to move confident calls, or if it is
  confined to the 4-haplotype graph — the sampled panel is where the linkage is supposed to be
  real, so a 34-haplotype-only failure would mean the model is fitting construction artefacts.

  *This clause fired, in the direction it was written to catch, and the result is worth stating
  plainly: `f = 5` is worth +0.0099 of small-variant genotype F1 on the 34-haplotype graphs and
  under 0.0004 on the 4-haplotype ones — safe there, useless there. It is not a construction
  artefact, though: the mechanism predicts exactly this, since three haplotypes spell a genotype
  with at most a couple of pairs and leave almost no multiplicity for an exponent to act on. A
  rule written to catch overfitting cannot distinguish that from a parameter that simply needs a
  panel to work with, so this was decided on the mechanism rather than on the rule.*
- **Do not count phasing as part of the case.** It needs a separate Viterbi pass (§5), so it is
  additional work rather than a by-product, and it should be justified on its own.
- **Stop** if the off-panel escape has to be tuned to avoid suppressing novel alleles. That
  would mean the prior is fighting recall, which is the thing this project has spent most of
  its effort recovering.

## 10. Not measured, and worth knowing before committing

- **How often HG002's true genotypes imply panel switching.** That bounds how much a linkage
  prior can help without hurting, and it is the single most useful missing number. It needs
  truth genotypes expressed in graph alleles, which is not a join we currently do.
- **The minimum number of switches needed to explain the whole call set** as two paths through
  the panel — a Viterbi parse rather than a pairwise count. The pairwise measure here is a
  lower bound on disagreement, because it cannot see a switch that is locally consistent but
  globally impossible.
- ~~**Whether `β` fitted for accuracy matches `β` read off the panel.**~~ **Answered, and the
  question was malformed.** No `β` could be fitted, because smeared it is not a free parameter at
  all — it is the distance scale. And the panel-side reading came out at +0.008 NMI (z = 1.1) at
  the gaps that matter, so there was nothing to match against. The instinct behind the question
  was right — a fitted value far from the construction's nominal rate would have meant the
  transition was absorbing something else — but it presumed the axis existed.
- **Whether the 0.6% apparent recombinations are enriched for errors.** If they are truth-set
  false positives, the prior would fix real mistakes; if they are true novel haplotypes, it
  would suppress correct calls. Joining them to the truth would say which, and it decides
  whether this is a fix or a regression.
