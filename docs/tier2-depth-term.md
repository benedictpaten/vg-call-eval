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
