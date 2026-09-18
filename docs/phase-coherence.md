# Phase coherence: the site statistic that finally moves switch error

Every previous intervention on switch error failed identically -- `--phase-break`, `--phase-cap`,
`--phase-min-gqn`, `--phase-confirm` -- and
[switch-error-causes.md](switch-error-causes.md) concluded that no re-weighting of the same read
evidence fixes them. That conclusion was right about re-weighting. It was wrong to stop there,
because the thing that works is not a re-weighting of the evidence but a different question asked
of each SITE.

## The two questions

`--phase-min-q` gates on **reliability**: can this site's reads tell its two alleles apart? It is a
purely local, within-site statistic, and it never looks at any other site.

**Coherence** asks instead: do this site's reads agree with the haplotype their OTHER sites imply?
A site can be perfectly discriminable and completely phase-incoherent, and that is precisely the
site a phase link must not be built from. Reliability cannot see it.

Held out by construction: a read's haplotype is recomputed with the site's own term removed, so a
site never votes on itself.

## As a predictor of switch position

Ranking chr20's 113,924 eligible sites worst-first and asking what share lie within 10 kb of a
whatshap switch (base rate 9.43%):

| worst N% of sites | coherence | reliability |
|---|---|---|
| 0.1% | **87.6% -- 9.3x** | 20.4% -- 2.2x |
| 0.5% | 70.8% -- **7.5x** | 10.9% -- 1.2x |
| 1.0% | 48.6% -- **5.1x** | 10.6% -- 1.1x |
| 5.0% | 23.0% -- 2.4x | 8.2% -- **0.9x, below the base rate** |

Pearson r between the two is **0.437**, so this is new information rather than a restatement.
Median coherence is 1.000 and 7.2% of sites fall below 0.8 -- the bad sites are a small, sharply
separated tail.

## As a gate, chr20 sweep

| `--phase-coherence` | sites demoted | assessed | true switches | flips | all_sw rate |
|---|---|---|---|---|---|
| 0 (off) | -- | 58,799 | **51** | 77 | 0.3486% |
| 0.60 | 944 | 58,880 | 32 | 80 | 0.3261% |
| **0.70** | 1,592 | 58,880 | **33** | 76 | **0.3142%** |
| 0.80 | 2,489 | 58,882 | 39 | 78 | 0.3312% |
| 0.90 | 4,300 | 58,863 | 36 | 79 | 0.3296% |

## The gate: both contigs, both metrics

| | chr20 switches | chr6 switches | chr20 ALL F1 | chr6 ALL F1 |
|---|---|---|---|---|
| off | 51 | 62 | 0.95848 | 0.96487 |
| 0.70 | **33** | **31** | **0.95900** | **0.96517** |

chr6 is the hold-out and did not choose 0.70; chr20's own all-switch rate does, and chr6 agrees
independently. Shipped on by default in vg `6772c5ae8`.

**Two artefacts ruled out.** The assessed denominator RISES on both contigs -- 58,799 -> 58,880 and
160,054 -> 160,137 -- so this is not whatshap dropping re-genotyped sites
([[switch-error-excludes-genotype-changes]]). And F1 moves with TP up AND FP down on both, chr20
+40/-58 and chr6 +98/-67, so it is not a precision/recall trade.

## What did not work, and why it is worth recording

The same read data supports a second statistic, `--phase-cp`: for each junction, the aggregate
log10 gain over reads SPANNING it of flipping everything downstream, using each read's full span.
S(j) > 0 implies the flip raises sum_r max_H L(H, r), so greedy flipping hill-climbs on the total
read likelihood and terminates. This is the transitive constraint the pairwise cascade discards.

It fires on **3 junctions in all of chr20**, each worth ~230 log10 units, and moves switches 51 ->
46. That near-silence is the finding: the cascade is already at a strong local optimum of read
likelihood, so the switches are not read-likelihood errors. It is the same answer the four earlier
interventions gave, reached for the first time from a global rather than a pairwise direction.

Which is exactly why coherence works where the rest did not. It does not ask the reads to overturn
a decision they collectively support; it removes the sites whose reads should never have been
allowed to vote.

## Method note

The per-read misfit that led here -- the delta between a read's likelihood under its assigned
haplotype and under its per-site optimum -- is itself a usable switch detector when ranked by how
well a single CHANGEPOINT explains the read's disagreement pattern: 70% of the top 50 reads sit
within 10 kb of a switch, 8.5x enriched, against 1.8x for the unranked set and 2.1x for raw
disagreement density. Shape beats magnitude. That detector is what suggested looking at sites
rather than junctions.
