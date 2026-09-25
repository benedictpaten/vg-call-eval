# Dropping low-MAPQ reads, and why the mismap clamp is not the lever

> **A record of 2026-09-16; current state (2026-09-23).** The switch-error case for a floor of
> 20-30 below does not hold, and `--preset ont` sets `--read-min-mapq 5` instead, on reliable het
> count ([mapq-floor-is-five.md](mapq-floor-is-five.md)). The cap default is now 0.95, not 0.7. On
> ONT both levers are inert at the benchmark: on the current preset, dropping the floor, at a cap of
> 0.7 or 0.95, moves no F1 on chr20 or chr6 by more than 2e-5 and no truvari SV call, and whatshap
> switch counts stay within noise (floor on vs off at 0.95: chr20 172 vs 182, chr6 340 vs 340;
> `work/wgs-run5`). The saturation argument
> below is about *phasing*. For short-read genotyping the cap is a live lever: 0.7 → 0.95 gains
> whole-genome SV F1 +0.0021 over chr1-22+X with small variants unchanged within noise. Short
> reads keep `--read-min-mapq 0`, because MQ5 there adds SV F1 +0.0067 over chr1-22+X at 3,290 more
> net SNV errors, both arms at the then-default cap of 0.7; against the current default the MQ5 arm
> still leads on SV F1 by 0.0047 and trails on SNV F1 by 0.0005, both significant
> (`work/wgs-run3/compare_mq5_vs_mm095.out`; [tier2-parameters.md](tier2-parameters.md)).

chr20 ONT, `--preset ont`, vg `54e1b31c9`. Prompted by a switch-error regression that turned out to
be a mapping problem, and by the observation that the reads causing it are not the ones the mapper
flags.

## The loci that go wrong are a mapping artefact, not a coverage one

Sites within 25 kb of a switch that admitting off-reference chains newly created:

| window +/- 25 kb | DP | DR | BL | GQ |
|---|---|---|---|---|
| new switch loci | 41.0 | **1.574** | **12.0** | 38.0 |
| the 65.5 Mb cluster | 42.0 | **1.594** | **7.0** | 33.0 |
| kept switch loci | 33.0 | 1.074 | 33.7 | 43.0 |
| all sites | 42.0 | 1.014 | 111.2 | 84.0 |

**Depth is normal.** What is abnormal is 57% more reads than the called genotype explains (DR) and a
nine-fold collapse in how well any allele fits them (BL). That is a paralogous pileup: reads from a
locus the graph does not contain, mapping here at full depth and matching nothing.

## The reads are NOT low-MAPQ, but the low-MAPQ minority is what does the damage

MAPQ of the reads actually placed at those loci, from the GAF:

| reads at | median | mean | MAPQ 0 | MAPQ < 10 | MAPQ 60 |
|---|---|---|---|---|---|
| new switch loci (n=1,045) | 60 | 47.4 | **3.1%** | **12.2%** | 69.6% |
| random loci (n=3,695) | 60 | 59.6 | 0.1% | 0.4% | 99.0% |

MAPQ 0 is 31x enriched and MAPQ<10 is 30x -- but 70% of the reads there are MAPQ 60, because MAPQ
asks "could this read be somewhere else IN THIS GRAPH", and when the true source locus is absent the
honest answer is no.

**I concluded from this that a MAPQ filter could not help, and that was wrong.** The 12% minority
turns out to carry the damage. The experiment settles it where the inference did not.

## The scan

| `--read-min-mapq` | het sites | reliable hets | true switches | flips | ALL F1 | SNV F1 | SV F1 |
|---|---|---|---|---|---|---|---|
| **0 (default)** | 76,135 | 60,203 | **55** | 85 | 0.95848 | 0.98541 | 0.5589 |
| 5 | 76,056 | 61,261 | 51 | 77 | -- | -- | -- |
| 11 | 76,036 | 61,776 | 49 | 80 | 0.95840 | 0.98545 | 0.5589 |
| **20** | 76,064 | 61,814 | **45** | 80 | 0.95846 | 0.98541 | 0.5580 |
| **30** | 76,095 | 61,842 | **45** | 81 | -- | -- | -- |
| 60 | **77,251** | **62,771** | 48 | 79 | **0.95862** | **0.98555** | 0.5589 |

Three things at once, and the third is the surprise:

- **Switches fall 55 -> 45**, flattening at 20-30.
- **F1 does not move** -- ALL, SNV, insertion and deletion all flat to the fourth decimal, and MAPQ 60
  is fractionally the best of them. Discarding these reads costs nothing in calling.
- **Reliable hets RISE monotonically**, 60,203 -> 62,771, and at MAPQ 60 the het count rises too.
  Throwing reads away makes MORE sites phaseable, because site reliability is a MEAN and the
  discarded reads were dragging it below `--phase-min-q`.

## The clamp is already saturated, so raising it does nothing

| `--mismap-max` | reliable hets | true switches | flips |
|---|---|---|---|
| **0.7 (default)** | 60,203 | **55** | 85 |
| 0.9 | 60,202 | 57 | 83 |
| 0.99 | 60,201 | 57 | 83 |

The reasoning that a MAPQ 0 read is worse than 70% likely mismapped is sound, and the lever still
does nothing -- because of the escape mixture in `phase_link`:

```
cis = pr*same + (1-pr)/2      trans = pr*diff + (1-pr)/2
```

As `pr` -> 0 the escape term dominates and `log10(cis/trans)` -> `log10(1)` = **0**. A heavily
down-weighted read already contributes exactly nothing to a phase link, so more down-weighting has
nowhere to go. The clamp sits deep in the saturated regime.

**Removal and down-weighting are therefore not the same lever.** Down-weighting neutralises the
read's phase VOTE, which is already achieved. Removal additionally takes the read out of:

1. **the genotype** -- its responsibilities still shape AD, GL and which traversals settle, and every
   other read's `q0` is measured against that settled pair;
2. **site reliability** -- a MEAN over the site's reads, so one bad read drags the site under
   `--phase-min-q` and removes it from the chain entirely. This is the direct mechanism for the
   +2,568 reliable hets, and nothing in the mismap term can reach it.

So the clamp is not what wants redesigning. The harm flows through genotyping and through the
reliability mean, and only exclusion touches either.

## chr6, held out

| `--read-min-mapq` | het sites | reliable hets | raw rate | true switches | flips |
|---|---|---|---|---|---|
| **0 (default)** | 197,662 | 162,622 | 0.2866% | **63** | 198 |
| 11 | 197,563 | **163,527** | 0.2860% | 60 | 199 |
| 20 | 197,426 | 163,404 | 0.2835% | 60 | 197 |
| 30 | 197,425 | 163,413 | 0.2817% | 59 | 196 |
| 60 | 197,333 | 163,406 | **0.2761%** | **58** | 192 |

**It replicates, and that matters because the last chr20 optimum did not.** Switches fall monotonically
63 -> 58, flips 198 -> 192, the raw rate 0.2866% -> 0.2761%, and reliable hets rise. **No arm is worse
than baseline on either contig.**

What does NOT transfer is the precise optimum: chr20 bottoms at 20-30 and ticks back up at 60, while
chr6 keeps improving to 60. So the shape of the curve is contig-specific and only the direction is
shared -- which argues for a conservative value rather than chr20's argmin.

The magnitudes differ too: chr20 -10 switches (18%), chr6 -5 (8%). chr20 carries the pericentromeric
and subtelomeric clusters this began with; chr6's are presumably elsewhere.

## Status and recommendation

`--read-min-mapq` already exists and defaults to 0. On the evidence -- both contigs improved, monotone,
no arm worse anywhere, F1 flat to the fourth decimal in every class on chr20, and reliable het counts
up on both -- a default in the **20-30** band is defensible: chr20's optimum, and within one switch of
chr6's best.

The conservative alternative is 11, which improves both (chr20 55 -> 49, chr6 63 -> 60) while
discarding the least. Either is a default change to a shipped preset and should be a deliberate call,
not a consequence of this measurement.

---

**CORRECTION (vg `05acbfc6a`).** The switch-error argument on this page does not
hold: it is one pericentromeric megabase on chr20, does not replicate on chr6, and the shipped
value was never the measured arm. The ONT floor is now 5, set on reliable het count instead. See
[mapq-floor-is-five.md](mapq-floor-is-five.md).
