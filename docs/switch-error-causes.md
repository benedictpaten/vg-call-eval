# Where the switches in the long phase blocks come from

chr20 ONT, vg `e31c264c9`, whatshap against T2T-Q100. The block is one 66.2 Mb chain carrying **52
true switches** (plus 71 flips, which are a separate and genotype-driven error -- see
[switch-error-in-long-blocks.md](switch-error-in-long-blocks.md)).

## The algorithm uses ZERO flanking variants

`src/read_phasing.cpp:118-137`. Stage 1 decides 59,118 of 60,202 links -- **98.2%** -- like this:

```cpp
d[m] = phase_link(sites[begin + rel[m]], sites[begin + rel[m + 1]], params.cap);
...
o[rel[m + 1]] = o[rel[m]] ^ (d[m] < 0.0 ? 1 : 0);
```

`phase_link` merge-intersects the two sites' read lists and sums `log10(cis/trans)` over the reads
they SHARE. `PhaseSite` holds no cross-site state, so no third site can contribute. The orientation
is then a **sign-only XOR cascade**: `|d[m]|` is read once, for the break test at `:125`, and never
again. A link decided at `|d| = 10.001` propagates with exactly the authority of one at `|d| = 400`,
and one wrong sign inverts every reliable site to the end of the segment -- mean 55.5 sites -- with
stage 3 hanging the unreliable ones off the inverted chain so they invert too.

A read spanning sites m, m+1, m+2 contributes to `d[m]` and to `d[m+1]` independently, and **nothing
ever asks whether the implied m->m+2 relation agrees with a direct link**. The transitive constraint
the read supplies is discarded. Only the 1,084 break junctions get a wider view, and that is still a
sum of pairwise links compared against 0.0 with no threshold.

So the answer to "is it using enough flanking variants" is **no, it uses none**. That is a real
architectural weakness. It is also, on this data, **not what produces these 52 switches.**

## Four interventions, all measured, all negative

| arm | chain breaks | true switches | flips |
|---|---|---|---|
| **default** | 1,084 | **55** | 85 |
| `--phase-break 20` | 4,113 | 54 | 88 |
| `--phase-break 40` | 56,857 | 61 | 80 |
| `--phase-break 80` | 60,202 | 61 | 80 |
| `--phase-cap 1 / 2 / 3` | 60,202 | 71 / 60 / 59 | 78 / 82 / 81 |
| `--phase-min-gqn 0` | -- | 58 | **261** |

- **Raising the break threshold does nothing.** At 20 the weak links become breaks and are relinked
  over 3 reliable sites either side; switches go 55 -> 54. At 40, 94% of ALL links are relinked, and
  switches get *worse*. **The switches are therefore not at weak links** -- they are at links with
  large `|d|` that are confidently wrong, and stage 2 cannot rescue them because it re-sums the same
  pairwise evidence.
- **Capping a read-pair's contribution does nothing**, so the independence inflation that lets `|d|`
  reach 10^40 is not the proximate cause either.
- **Gating on GQN makes it worse.** Removing the 3,192 panel-overridden sites costs their links; the
  sites are then hung off the chain by the weaker stage 3 and flips triple, 85 -> 261.

## What the switches actually are: a few bad loci

| | observed | expected if uniform |
|---|---|---|
| switches within 100 kb of another | **23 of 52 (44%)** | 7.6 (15%) |
| busiest 1 Mb bin | **11 switches at 26 Mb** | ~0.8 |
| minimum inter-switch distance | **15 bp** | -- |

The 26 Mb cluster, listed: `26,459,358` `26,539,627` `26,539,819` `26,540,024` `26,540,300`
`26,540,357` `26,540,372` `26,645,320` `26,645,401` `26,744,760` `26,745,948`.

**Six of them fall within 745 bp**, two more 81 bp apart. A chain making independent per-link errors
does not do that; a locus whose reads are mismapped does, because the phase thrashes across it. The
cluster sits in pericentromeric sequence, immediately proximal to chr20's 1.99 Mb variant-free
centromere gap at 27,053,143-29,039,480.

The GQN signal fits the same story rather than competing with it. Sites the panel overrode are 10x
enriched at switch junctions (25.0% against 2.5%, surviving stratification by gap width at 23.4x,
7.3x, 10.8x, 10.1x in the four narrow bands) -- but gating on them hurts, because the enrichment is a
**symptom** of the underlying bad region, not an independent lever. In a mismapped locus the reads
are wrong, so the panel overrides the genotype AND the phase goes wrong, from one cause.

## Resolution

**The 52 switches are dominated by a handful of hard loci, chiefly pericentromeric, where read
mapping is unreliable. No re-weighting of the same read evidence fixes them, and all four available
knobs were measured and do not.** Excluding the clustered ones leaves 29 isolated switches over
66.2 Mb -- one per 2.3 Mb -- which is the rate the algorithm achieves where the data are sound.

Two things are worth doing, and neither is a parameter change:

1. **Retain `d[m]`.** It is computed at `src/read_phasing.cpp:118`, used twice, and dies at end of
   scope. It is the per-link margin and the only quantity that ranks junctions by how narrowly they
   were decided. Reporting it -- or the count of links in a 10-15 band per segment -- costs nothing
   and would let a consumer distrust a region rather than a whole chromosome. Read it as a RANK, not
   a probability: `src/read_phasing.hpp:119-121` records that 40 reads give 10^40 against a measured
   read-only error of 10^-2.4.
2. **Treat the pairwise cascade as the known architectural debt it is.** It is genuinely fragile --
   sign-only, no transitive constraint, unbounded downstream propagation -- and a locus with mismapped
   reads is exactly the input it handles worst. Fixing it means a joint solve over the read-site
   graph, not a wider window in the same framework, which `--phase-break 40` already shows buys
   nothing.

`--phase-min-gqn` ships defaulted to -1, which admits everything and is inert. It is kept because the
10x enrichment is real and re-measurable, not because any other value is recommended.
