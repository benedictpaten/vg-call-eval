# Building a consistent backbone: three designs, one works, and the reason is the same each time

Follow-on from [coherent-backbone-plan.md](coherent-backbone-plan.md). chr20 ONT, `--preset ont`,
switch error by whatshap. The shipped baseline is `--phase-coherence 0.70`, at 33 true switches
against 51 with no gate.

## The three things tried

| design | question it asks | frame-free? | when |
|---|---|---|---|
| `--phase-coherence` (shipped) | do a site's reads agree with the haplotype their OTHER sites imply? | no | after phasing |
| `--phase-backbone` | reconsider each site against K neighbours weighted by \|d\| instead of one adjacent sign | no | after stage 1 |
| `--phase-triangle` | does sign(d_ij)·sign(d_jk)·sign(d_ik) > 0 close? | **yes** | **before phasing** |

## Iterating the coherence gate to a fixed point: converges, and does not help

| `--phase-coh-rounds` | demoted | chain breaks | switches | S+2F | converged |
|---|---|---|---|---|---|
| 1 | 1,592 | 1,494 | **33** | 185 | no |
| 2 | 1,630 | 1,805 | 31 | 181 | no |
| 3 | 1,636 | 2,114 | 31 | 183 | no (cap) |
| 5 / 10 | 1,636 | 2,114 | **33** | 183 | yes, in 3 rounds |

It reaches a fixed point in three rounds and the fixed point scores exactly what one round scores,
for **41% more chain breaks**. The 31 at two rounds is not a principled stopping point and 31 against
33 on ~30 events is noise. The predicted fragmentation cost is real and cancels the gain.

## Backbone local search: negative

| `--phase-backbone` | sites moved | switches |
|---|---|---|
| off | -- | **33** |
| 2 | 28 | 33 |
| 3 | 41 | 33 |
| 5 | 87 | **35** |

## Triangle pre-screen: sound, and there is almost nothing to screen

| arm | sites excluded | chain breaks | switches | S+2F |
|---|---|---|---|---|
| coherence 0.70 (baseline) | -- | 1,494 | **33** | 185 |
| triangle 0.70 alone | 142 | 1,003 | 37 | 193 |
| triangle 0.85 alone | 278 | 880 | 38 | 190 |
| triangle 0.95 alone | 530 | 743 | 38 | 190 |
| triangle 0.85 + coherence | 278 | 1,190 | 35 | 185 |
| triangle 0.95 + coherence | 530 | 1,036 | 35 | 187 |

Better than no gate at all (37-38 against 51), worse than coherence (33), and no better in
combination. The measurement that explains it is in the same log line:

> **1,771 of 367,218 triangles did not close -- 0.482%.**

The reads' pairwise phase relations are **99.5% transitive**. There is no large inconsistent subset
to exclude, so the largest set of mutually consistent sites is very nearly all of them, and a
criterion that selects for mutual consistency has almost nothing to select on.

## Why all three point the same way

Three independent probes of the same question, each from a different direction:

- **locally**: 99.5% of triangles close;
- **globally**: `--phase-cp`, the aggregate log10 gain of flipping everything downstream of a
  junction over reads spanning it, fires on **3 junctions in all of chr20**;
- **by refinement**: reconsidering every site against 2-5 weighted neighbours moves 28-87 sites and
  changes nothing.

**The read evidence is both locally and globally self-consistent.** The switches are not
inconsistencies in it. They are places where the reads agree with each other and are wrong together,
which is what [switch-error-causes.md](switch-error-causes.md) concluded from four failed
interventions and what these three confirm from three new angles.

That is also why coherence works where these do not, and the distinction is worth stating precisely.
Coherence does not look for inconsistency in the site-level pairwise structure -- there is almost
none. It finds sites at which individual READS disagree with their own other sites. That is a
statement about read reliability at a site, not about the site's phase relations being incoherent
with its neighbours, and it is a different quantity from anything a triangle or a local search sees.

## And the deserts are not empty

The other candidate explanation was that long junctions assert phase across intervals with no
evidence. Counting reads that carry both flanking sites:

| gap | junctions | mean reads spanning | with zero |
|---|---|---|---|
| 0-500 bp | 50,015 | 39.4 | 0.0% |
| 5-20 kb | 1,760 | 28.0 | 0.1% |
| 20-50 kb | 240 | **17.6** | 1.2% |
| >50 kb | 16 | **13.4** | 12.5% |

Even past the 33 kb read N50 a mean of 17.6 reads span, and only 1.2% of those junctions have none.
So breaking the block at long gaps would discard real evidence rather than acknowledge its absence.
The evidence is there, at roughly half depth, carried by the longest reads in the distribution.

## Where this leaves it

`--phase-coherence 0.70`, one round, is the configuration to keep: 51 -> 33 on chr20 and 62 -> 31 on
chr6 with F1 up on both. `--phase-backbone`, `--phase-cp`, `--phase-triangle` and
`--phase-coh-rounds > 1` are all measured, all negative or neutral, and all default to off. They are
kept because each is a re-usable probe of a question that will be asked again, and because the
0.482% and the 3-junction figures are the most informative numbers produced by this whole line of
work.

The residual 33 are three regions -- 26-35 Mb at 0.64x depth, 50-52 Mb, and 65.3-65.5 Mb at 6.7x het
density and normal depth, independently the most fragmented 100 kb in the anchor graph. Those need a
better graph or better mapping, not a better phasing algorithm.
