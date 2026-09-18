# Phasing parameters, measured on both contigs

Follow-on from [backbone-consistency-results.md](backbone-consistency-results.md). Baseline is the
shipped `--phase-coherence 0.70`. **The metric is TRUE SWITCHES**: a switch inverts everything
downstream -- mean 55.5 sites -- while a flip corrupts one. `S + 2F` is whatshap's raw boundary
count and prices a flip at TWICE a switch, which is backwards for long-range use; it is reported
only as a guard against an arm reducing switches by converting them into flips.

**Chain breaks do not fragment the output.** Measured directly: 707 breaks and 7,565 breaks both
give ONE `PS` block for the whole contig. A break cuts stage 1's cascade and stage 2 relinks it, so
more breaks means more of the chain is decided by the 9-pair magnitude-weighted relink instead of
the single-link sign cascade. An earlier reading of breaks as fragmentation was wrong, and it
inverted the recommended direction for `--phase-break`.

## Single parameters

| arm | chr20 switches | chr6 switches | note |
|---|---|---|---|
| default | 33 | 31 | |
| `--phase-coh-rounds 2` | **31** | **30** | helps slightly on both |
| `--phase-break 20` | **31** | **30** | helps slightly on both; 5-7x more breaks |
| `--phase-break 5` | 33 | **33** | WORSE on chr6 -- less relinking |
| `--phase-relink 10` | **29** | 31 | chr20 only |
| `--phase-lookback 3 / 8` | 33 | -- | overrules the adjacent link at 1-5 sites in a whole contig |
| `--phase-backbone 2/3/5` | 33/33/35 | -- | negative |
| `--phase-triangle 0.70-0.95` | 37/38/38 | -- | negative; 99.5% of triangles already close |

## Combinations, and the ablation

| arm | chr20 switches | chr6 switches | chr20 ALL F1 | chr6 ALL F1 |
|---|---|---|---|---|
| default | 33 | 31 | 0.95900 | 0.96517 |
| `rounds2 + break20` | **27** | 30 | 0.95900 | 0.96516 |
| `rounds2 + break20 + relink10` | **27** | **28** | 0.95900 | **0.96548** |

**The ablation says relink 10 adds nothing on chr20** (27 either way) and 2 on chr6 -- while relink
10 ALONE helps chr20 and not chr6. No consistent direction on switches.

But F1 separates it, on a far larger denominator than a ~30-event switch count: chr6
`+relink10` is the only arm that moves F1 at all, and it moves it UP by 0.00031 with **TP +80 and
FP -97**, both directions right. Either relink 10 is genuinely doing something on chr6 that a
30-event switch count cannot resolve, or it is one draw from the genotype perturbation that any
phasing change causes through re-genotyping. One observation on one contig does not separate those.

**`rounds2 + break20` passes the gate**: chr20 33 -> 27, chr6 31 -> 30, F1 flat on both. The chr20
gain is 6 switches and the chr6 gain is 1, so the replication is in sign, not magnitude.

## Where the survivors are

| | switches | clustered (<100 kb) | scattered |
|---|---|---|---|
| chr20 default | 33 | 26 (79%) | 7 |
| chr20 combination | 27 | 22 (81%) | 5 |
| chr6 default | 31 | 14 (45%) | 17 |
| chr6 combination | 28 | 13 (46%) | 15 |

The combination removes clustered and scattered switches in equal proportion; it is not
preferentially cleaning up either. **The two contigs have different shapes** -- chr20 is 79%
clustered into two megabases, chr6 only 45% with 15 isolated switches spread from 3 to 172 Mb --
which is part of why chr20-tuned parameters keep failing to replicate there. chr6's scattered set
also partly TURNS OVER between arms rather than shrinking, which is what noise looks like.

## Computational cost: all of it is free

chr20, `/usr/bin/time -l`, 6 threads:

| arm | wall (s) | user CPU (s) | max RSS |
|---|---|---|---|
| coherence off | 234.9 | 1287.1 | 7.01 GB |
| coherence 0.70 (shipped) | 233.1 | 1286.3 | 6.54 GB |
| `relink 3` | 243.9 | 1300.2 | 6.37 GB |
| `relink 10` | 230.6 | **1300.4** | 6.21 GB |
| `break 20 + relink 3` | 237.4 | 1312.2 | 6.16 GB |
| `break 20 + relink 10` | 234.2 | **1291.1** | 6.35 GB |
| `--phase-backbone 3` | 258.3 | 1323.0 | 6.12 GB |

Stage 2 costs `breaks x K^2` calls to `phase_link`, so `break 20 + relink 10` is **57x the default's
relink work and 12x all of stage 1** -- roughly 757,000 calls. User CPU across every arm spans
1287-1323 s, a 2.8% range with no ordering, and the most expensive configuration on paper measured
the cheapest. The reason is that a long-read call is its read fetch: 76% of the time is in the
`gbz-base` subprocess and 63% is blocked in `waitpid`, and `phase_link` is a merge over two
pre-sorted lists with no allocation. **Cost is not a reason to prefer any of these settings.**

The `K^2` term only bites if the workload changes shape -- far denser sites, or a much lower
`--phase-break`. It is the term to watch if the preset moves.

## Status

Nothing here is defaulted on beyond `--phase-coherence 0.70`. `rounds2 + break20` is the candidate
with a clean gate; relink 10 is unresolved and would need chr6 repeated before it could ship.
