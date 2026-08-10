# A depth term: Stage 0, predicted offline before building it

The read-likelihood model scores `P(reads | G)` conditioned on the reads it was handed,
and never asks whether that many reads should be there. That blindness is behind two
known failures — collapsed-repeat pile-ups it cannot reject, and large heterozygous
deletions it still loses at 0.44 recall against the Poisson caller's 0.79–0.84 even after
the mixture fix.

Stage 0 predicts what a depth term would do, from artefacts already on disk, before any
of it is written. `depth_term_offline.py`.

**Verdict: build it, but it solves half the problem, and the other half is not what it
looked like.**

**Built, and Stage 1 beat the Stage 0 prediction.** `--depth-term W`, off by default,
with the rate measured over the read source's own fetch window. Structural-variant F1 at
`w = 0.25`, searched on chr20 and validated on chr6:

| dataset | `readlik-z` | **+ depth term** | `poisson` | `poisson-z` |
|---|---|---|---|---|
| chr20-4hap | 0.4890 | **0.5024** | 0.4954 | 0.4930 |
| chr6-4hap | 0.5542 | **0.5639** | 0.5490 | 0.5478 |
| chr6-34hap | 0.4963 | **0.5076** | 0.4944 | 0.4881 |

Heterozygous deletion recall above 1 kb — the class it was built for:

| dataset | flat mixture | + weighted mixture | **+ depth term** | `poisson-z` |
|---|---|---|---|---|
| chr20-4hap | 0.061 | 0.212 | **0.394** | 0.364 |
| chr6-4hap | 0.066 | 0.443 | **0.787** | 0.836 |
| chr6-34hap | 0.082 | 0.443 | **0.721** | 0.787 |

Small-variant genotype F1 does not move to four decimal places on any dataset, the
heterozygous fraction of calls shifts by 0.005, and chr20 runs in 111 s against a
120–170 s baseline. Stage 0 predicted 4 of 10 sites with a deliberately naive global
rate; a local rate does better, as Stage 0 said it would.

Details of what was built, and what the in-tree suite caught, are at the end of this
page.

---

## The formulation

A complete generative model factorises, so depth enters as an *additive* term rather
than a filter:

```
P(data | G) = P(N | G) · P(reads | N, G)

ln P(data | G) = ln Poisson(N ; λ_G) + Σ_r ln[ (1−e_r)·Σ_h w_h·rel(r,h) + e_r ]

λ_G = c · Σ_{h∈G} (L_h + R − 1)
```

**The footprint appears twice and means two different things.** The mixture weight `w_h`
wants sequence *unique* to an allele, because only reads over unique sequence can
separate genotypes — reads in shared sequence fit everything and cancel. `λ_G` wants the
*whole* traversal length, because every base generates reads whether or not those reads
discriminate. Conflating them would be an easy and invisible error.

`c` is calibrated as the median of `N / Σ_h(L_h + R − 1)` over the called genotype across
all dumped sites, which puts it in the same units as `N` — "rows in the likelihood
matrix", which is not coverage and not what a pack file would give.

## The signal is real and well-conditioned

The sharpest evidence is not the flip count but the implied read rate. For each missed
deletion, what per-position rate would each genotype need to explain the observed `N`?

| snarl | N | `c` implied by the **called** genotype | by the **het deletion** |
|---|---|---|---|
| `173428820` | 357 | 0.0577 | **0.1009** |
| `175751591` | 366 | 0.0556 | **0.1034** |
| `175895886` | 709 | 0.7447 | **0.0842** |
| `175991492` | 227 | 0.1998 | **0.0876** |
| `176377415` | 164 | 0.3306 | **0.0894** |
| `176391477` | 166 | 0.4368 | **0.0926** |
| `176749777` | 215 | 0.7072 | 0.0336 |
| `177587516` | 165 | 0.0448 | **0.0762** |
| `177703025` | 310 | 0.0491 | **0.0915** |
| `178573961` | 172 | 0.2671 | 0.0523 |

The global rate is **0.1022**. The heterozygous-deletion genotype implies 0.05–0.10 at
every site — tight, and right at the global rate. The genotype the caller actually
chooses implies anything from 0.045 to 0.745. Depth is not a weak signal here; the
correct answer is the one that makes the read count make sense.

## But the missed deletions are two different problems

Read density on the *called* genotype's footprint, against the global 0.1022:

| group | sites | called | reads per position | vs global |
|---|---|---|---|---|
| **A** | 4 | homozygous **long** | 0.045–0.058 | **0.4–0.6×** |
| **B** | 6 | homozygous **short** (the deletion) | 0.20–0.75 | **2–7×** |

**Group A is the problem we thought we had.** The caller picks a homozygous long
genotype that predicts about twice the reads observed. The depth term flips all four.

**Group B is a different failure wearing the same clothes.** The caller *does* call a
deletion — homozygous — and packs 2 to 7 times the expected reads onto it. These are
collapsed repeats. At `175895886`, 709 reads sit on a 326 bp allele: 0.745 reads per
position against a global 0.102, and `readlik-z` emits **no record at all** while
`poisson-z` calls the −7466 bp deletion from **DP 18**. The read caller is being handed
709 reads where the depth caller sees 18 informative ones.

The depth term detects group B correctly and emphatically — it penalises the called
genotype by 800 nats at that site — and still cannot fix it, because the read term
prefers the same genotype by 1,773 nats (−31.3 against −1803.9). No defensible weight
closes that gap.

| depth weight `w` | sites calling the heterozygous deletion |
|---|---|
| 0.00 (today) | 0/10 |
| 0.25 | **4/10** |
| 0.50 | **4/10** |
| 1.00 | **5/10** |

Exactly group A, plus one of group B at `w = 1`.

## The kill criterion was not triggered

The stated risk was that `λ_G` grows with allele length, so an anomalously large `N`
would mechanically favour whichever genotype presents the most sequence — a preference
for long alleles rather than a rejection of the site.

Measured on the 58 chr20 pile-up sites (`N` above three times the median of 293), the
number moved to a longer-footprint genotype is **5 of 58** at `w = 0.25` and `w = 0.5`.
The term churns some genotypes at these sites but does not systematically inflate them.
Proceed.

## Should low-MAPQ reads count less toward the observed depth?

Yes, and the model already carries the right quantity: `e_r` is "this read came from
somewhere else", so the expected number of reads genuinely from this locus is

```
N_eff = Σ_r (1 − e_r)
```

Comparing raw `N` against a locus-specific expectation is the inconsistency; `N_eff` is
the observable `λ_G` should be measured against. But measuring it splits group B again:

| group B sub-case | sites | reads at the `e_r` cap | median `e_r` | density `N` → `N_eff` |
|---|---|---|---|---|
| multi-mapping | 2 | 91%, 56% | 0.700 | 0.707 → **0.251**, 0.267 → **0.160** |
| confidently mapped | 4 | 0–2.4% | 0.020 (the floor) | 0.745 → 0.728 |

At four of six pile-ups the reads sit at the *floor* — MAPQ ≥ 17. The mapper is
confident, and it is confidently placing seven times too many reads there.

**That has a mechanism, and it is the argument for a depth term rather than against
one.** If a repeat is collapsed in the graph, reads from every copy pile onto the single
copy and there is no alternative placement to be ambiguous about, so MAPQ is high. MAPQ
measures ambiguity *within the graph*, not correctness relative to the real genome, so it
is blind to graph collapse by construction. Depth sees exactly the failure MAPQ cannot,
which is why the two are complementary rather than redundant.

**`N_eff` does not improve genotype flips** — it slightly reduces them, 4/10 to 3/10 at
`w = 0.5`. 38.7% of reads across these large sites sit at the cap, so recalibrating `c`
globally against `N_eff` mostly rescales it (0.1022 → 0.0795); group A sites have `e_r`
at the floor, so their `N_eff` is 98% of `N` while their `λ` falls 22%, which makes them
look better covered and weakens the heterozygote's advantage. The correction lands on the
sites that do not flip anyway.

So the two uses of depth separate cleanly, and should be implemented separately:
**`N_eff` for judging whether a site's depth is plausible** (the quality signal, and the
pile-ups), **raw `N` against a locally estimated rate for genotype selection** (group A).
This is also a second, independent argument for estimating `c` locally: a global
recalibration artefact is what cost the flip here.

## Where the local rate comes from: reads, not the pack — and not from within the site

The rate must be derived from the reads rather than a pack file, and not only to keep
`readlik-z` pack-free. `N` is *rows in the likelihood matrix*, which is neither coverage
nor what a pack reports: it depends on the read-fetch window and the placement filter.
Estimating the rate through the same fetch and placement path makes the units match by
construction, where a pack-derived rate would be a cross-unit comparison with a silent
scale error in it.

**The elegant version does not work.** Reads lying in a snarl's *shared* sequence fit
every allele equally, so they are genotype-independent and would give a local rate for
free, with no extra fetching, no global state and no second pass. Measured on the ten
missed deletions, the number of such reads is **0 at nine sites and 2 at the tenth**.
Snarl boundary nodes are too short to contain a read, and anything extending past them
reaches variable sequence and discriminates. There is no within-site control available at
exactly the sites that need one.

So the local rate has to come from outside the site. Three read-based options:

| | how | cost | risk |
|---|---|---|---|
| widen the fetch | extend the node-ID range to pull flanking reads; `--read-window` already exists | more reads per site, on the hot path | the flank may be repetitive too |
| coarse pre-pass | stream the read source once, bin read starts into ~10 kb windows | one extra pass over the read source | the cost is the open question |
| rolling estimate | accumulate `N_i / Σ e_i` from sites already processed | free | thread-unsafe under `-t`, plus a warm-up problem |

The pre-pass is the recommendation: it is computing the pack ourselves, coarsely, from the
same source that feeds the matrices, and it is read-only once built, so it is thread-safe
where a rolling estimate is not. Widening the fetch is the fallback if streaming proves
too expensive, and measuring that cost is the first task of Stage 1.

**A local rate does not replace the global one.** Any near-local control normalises away
exactly the whole-region anomaly the pile-up signal needs to see — if a repeat is
collapsed, its flank is frequently piled up as well. Group A wants *relative* depth within
a site; group B wants *absolute* depth against the genome. One term, two reference rates.

## Two requirements this puts on Stage 1

**A global `c` is too crude to ship.** The per-site ratio has an interquartile range of
0.062–0.221 — a 3.5× spread. At `N` in the hundreds, mis-estimating `λ` by 3.5× is worth
hundreds of nats, which is larger than the signal being extracted. `c` has to be
estimated locally, from flanking coverage, not globally. Stage 0's flip counts are
therefore a **lower bound**: they are what the term achieves with a deliberately naive
rate estimate.

**The depth term should also feed the quality field, not only genotype selection.** Group
B is the case where depth knows the site is wrong and cannot outvote the reads. That is
precisely the shape the explained-share discount already handles for a different
blindness: `GQ = GQI × share`. A depth-implausibility discount is its sibling, and it
would demote group B without needing to win an argument against 1,773 nats. It also
generalises to the pile-ups, which are the same sites.

## Caveats

- The read-side margins here use whole-traversal mixture weights rather than the shipped
  unique-content weights, because the dump records no allele identity and reconstructing
  unique node content needs per-node lengths the dump does not carry. The difference is
  single-digit nats against a depth term of tens to hundreds, so it changes no conclusion
  — but it means these are predictions, not measurements.
- Ten sites is a small set. It is the set that matters, but the flip counts should not be
  read to more than one significant figure.
- `--dump-likelihoods` still records no allele identity. Fixing that is a prerequisite
  for turning Stage 0 into a measurement, and is the same gap flagged during the
  heterozygous-deletion mechanism work.


---

# Stage 1: what was built

`--depth-term W` adds `W · ln Poisson(N ; λ_G)` to each genotype's log-likelihood, and
`DR` — observed reads over what the called genotype predicts — is emitted **whether or
not the term is armed**, so the observable can be measured as a ranking signal before the
model is trusted to act on it. That is the order the explained-share discount was
established in. On chr20 `DR` has a median of 0.60, a 99th percentile of 2.05 and a
maximum of 15.5.

## The rate comes from the fetch window

Stage 0 ruled out a within-site control and recommended a coarse pre-pass. The
implementation does something cheaper that was available all along: `WindowedSiteReadSource`
already fetches and caches a **4096-node window** around each site to answer the site's own
query, which is thousands of times wider than a snarl. Asking it for the window is a cache
hit. The per-window rate is memoised, so counting a window's reads happens once per window
rather than once per site — without that the diagnostic would cost more than the
genotyping. Measured cost: **111 s against a 120–170 s baseline**, i.e. inside the noise.

No pre-pass, no pack, no extra I/O.

## Two mistakes worth recording

**The unit test I wrote first asserted a false premise.** It set up a heterozygous
deletion and claimed the mixture alone would get it wrong. With alleles of 2000 and 50 bp
the length-weighted mixture already calls it correctly — that is the fix shipped earlier.
The test now uses a geometry where the mixture genuinely still fails (40 interior reads
against 1 junction read), which is the residual the depth term exists for.

**`DR` depended on which read source supplied the reads.** The neighbourhood width was
taken from the read source's own fetch window, and an in-memory source answers every range
exactly and so has none — it produced no rate, emitted no `DR`, and therefore a different
VCF from an indexed source over identical reads. `18_vg_call.t` asserts those two agree
and caught it immediately. The neighbourhood is caller policy, not a backend property, and
now falls back to `--depth-window` when the source has no window of its own.

## Still open

- **Off by default**, pending replication beyond two chromosomes.
- **Group B is untouched.** The term detects collapsed repeats emphatically and still
  cannot outvote the read evidence there; that needs `DR` wired into the quality field,
  which is why `DR` is emitted now.
- `w` was chosen from three values on one dataset. 0.25 beat 0.5 and 1.0 on chr20 and was
  not re-optimised on chr6, so it is a defensible setting rather than a fitted one.
