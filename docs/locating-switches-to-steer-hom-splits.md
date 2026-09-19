# Can the switches be located well enough to steer hom splitting?

The concern: splitting homozygous anchors next to a switch lengthens a haplotype-labelled chain
that *contains a switch error* -- contiguity bought against correctness. The remedy would be to
find the probable switches from the run's own evidence and suppress splits there.

**Short answer: yes, with real and replicated precision -- 4x to 13x enrichment for 3-6% of the
splits -- but the intervention does not do what it is meant to, because the chain does not run
through the splits.** Both halves are measured below.

## Where the switches are: stage 2, not stage 1

Stage 1 cascades over reliable sites and BREAKS where `|d| < --phase-break`; stage 2 relinks the
blocks across each break with a K x K pairwise vote. A switch is a wrong orientation, so it is
either a wrong sign in the cascade or a wrong relink.

**18 of chr20's 24 isolated switches sit at a junction stage 1 already broke** -- 5.7% of
junctions, so **13x enrichment before any new statistic is computed**. The phasing already knows
where it had no cascade evidence. What it does not know is which of its 6,892 relinks it got
wrong: the per-relink error rate is 0.26%.

## How precisely they can be ranked

Every signal below is available at run time and uses no truth. `phase_link` returns a SUM of
per-read log10 odds, so a junction can look confident while its readers are split -- two confident
disagreeing reads cancel. The minority *vote fraction* is therefore the statistic to beat `|d|`
with, and it does.

Suppressing a junction means dropping every split inside its interval plus the 2 nearest each side.

| rule | chr20 | chr6 (hold-out) |
|---|---|---|
| worst 1% of relinks, minority vote fraction | 7/28 switches, **6.2%** of splits, **4.1x** | 12/28, **3.4%** of splits, **12.8x** |
| worst 5% of relinks, same | 7/28, 6.9%, 3.6x | 12/28, 5.4%, 7.9x |
| worst 25% of relinks | 13/28, 28.0%, 1.7x | 13/28, 19.5%, 2.4x |
| all 6,892 / 9,251 relinks | 22/28, 40.8%, 1.9x | 18/28, 34.4%, 1.9x |
| worst 1% of *all* stage-1 junctions, by `\|d\|` | 12/28, 23.3%, 1.8x | -- |

Region-level does about as well and is worth recording because the mechanism is different --
**switches cluster**. On chr20, 82% of them fall in 5.4% of the contig: 6 at 26.5 Mb
(pericentromeric) and 10 at 65.3-65.5 Mb (subtelomeric). On chr6 the clustering is weaker (25% in
2.1%), so this is not a contig-independent fact, but the windowed statistic still ranks:

| 100 kb windows by mean phase coherence | chr20 | chr6 (hold-out) |
|---|---|---|
| worst 2% | 0/28, 0.9% of splits | 6/28, 3.5%, **6.2x** |
| worst 5% | 9/28, 4.5%, **7.2x** | 7/28, 7.2%, 3.5x |
| worst 20% | 21/28, 25.0%, 3.0x | 14/28, 23.1%, 2.2x |

Coherence is the best window statistic on both, which is consistent -- it is the same quantity
that works for site demotion. Break rate, median `|d|`, shared-read count, site depth and mean
`|lambda|` were all weaker, several below 1x.

**So: 4-13x at the sharp end, covering a quarter to a little under half the switches for 3-6% of
the splits.** That is genuine precision, and affordable.

## Why it will not lengthen anything, and what would

Suppressing splits near a switch **does not shorten the chain that contains it**. At every switch
junction on both contigs -- 28/28 chr20, 26/26 chr6 -- the two flanking het sites are joined
**directly** by reads pinned at both (see `hom-split-does-not-bridge-switches.md`; the library
spans a median 21 kb while assessed hets sit 260 bp apart). Remove the splits and the chain is
exactly as long, with the same wrong join asserted by the hets, minus some homozygous labels.

What would shorten it is a **phase-set boundary** -- refusing the weak relink instead of making it.
Priced against block N50:

| chr20 | blocks | block N50 | switches left inside |
|---|---|---|---|
| today (always relink) | 1 | 66.21 Mb | 28/28 |
| **oracle: cut at every true switch** | 29 | **9.11 Mb** | **0/28** |
| refuse worst 0.2% of relinks | 14 | 14.11 Mb | 25/28 |
| refuse worst 1% | 69 | 6.77 Mb | 21/28 |
| refuse ALL relinks | 6,893 | 0.09 Mb | 6/28 |

| chr6 (hold-out) | blocks | block N50 | switches left inside |
|---|---|---|---|
| today | 1 | 172.09 Mb | 28/28 |
| **oracle** | 29 | **16.03 Mb** | **0/28** |
| refuse worst 1% | 93 | 6.61 Mb | 16/28 |
| refuse ALL | 9,252 | 0.12 Mb | 10/28 |

The oracle is the row that settles it: cutting at exactly the true switches gives **better N50 and
zero switches** than any achievable rule. Refusing every relink -- the maximum any relink-based
rule can do -- costs 700x the N50 and still leaves 6 of 28. A 4-13x signal cannot reach that
frontier.

And the quantity that actually measures long-range accuracy is invariant to all of it:

    mean correct-phase run  =  span / (switches + 1)  =  2.28 Mb chr20,  5.93 Mb chr6

Cutting at a true switch does not lengthen that run. It relabels a silent error as an honest
boundary -- which is worth something, because a 66 Mb block with 28 switches in it is a claim no
consumer can use -- but it is honesty, not contiguity. Cutting anywhere else shortens the block and
touches nothing.

## Recommendation

Do not build the suppression as a way to shorten chains; it cannot. If it is built, build it for
what it does do: **remove wrong haplotype labels from the 1% of junctions that carry a quarter to
a half of the switches, for 3-6% of the splits.** That is a defensible, cheap, replicated
insurance policy, and it needs no new statistic -- the minority vote fraction falls out of the
relink loop that already runs.

It will not move F1 or the whatshap switch count, because it is anchors-only. The honest way to
report it is as a reduction in wrongly-labelled anchor sequence, not as a phasing improvement.

The larger opportunity is the one the oracle points at: switches are 13x enriched at stage-2
relinks, which is where the phasing itself admits it had no evidence. Sharpening the relink -- so
that 0.26% of them stop being wrong -- would beat any downstream masking of their consequences.
