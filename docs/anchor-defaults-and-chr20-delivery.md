# Anchors imply nesting, ONT drops MAPQ<10, and what chr20 looks like now

Two defaults changed in vg `7adcff51f`, and the chr20 deliverables were rebuilt on top of them.
Both changes are measured here against the PR base, `2646dfe0a`, on chr20 and on the chr6 hold-out.

## The two defaults

**`--anchors-out` implies off-reference nesting.** Nesting off the reference is the thing anchors
are *for*: a nested chain no reference path crosses has no VCF line to express it, so before this
it was genotyped, phased and then dropped on the floor of every output that could have carried it.
Turning it on only for the VCF's benefit had it backwards. `--no-off-ref-nesting` turns it back off.

**`--read-min-mapq` is 10 under `--preset ont`.** Measured in
[mapq-filter-and-the-clamp.md](mapq-filter-and-the-clamp.md): monotone on both contigs, no arm worse
than baseline, F1 flat to the fourth decimal. The global default stays **0**, not 10, because a
global 10 emits no variants at all against the MAPQ-0 reads in the test fixtures -- 65 of 434 TAP
checks fail with "got 0, expected 1". Keying it to the preset is what the evidence supports anyway:
the filter was fitted on ONT and nothing says an Illumina arm wants the same number.

## Cost: none that the F1 can see

| | ALL F1 | SNV F1 | INS F1 | DEL F1 | SV F1 |
|---|---|---|---|---|---|
| chr20 base `2646dfe0a` | 0.95848 | 0.98541 | 0.85417 | 0.87574 | 0.55888 |
| chr20 new `7adcff51f` | 0.95830 | 0.98535 | 0.85366 | 0.87514 | 0.55737 |
| chr6 base (hold-out) | 0.96518 | 0.98813 | 0.87137 | 0.89158 | 0.60973 |
| chr6 new (hold-out) | 0.96519 | 0.98812 | 0.87128 | 0.89175 | 0.61010 |

chr20 moves **down** in the fourth decimal -- ALL -0.00018, SV -0.0015 -- and chr6 moves fractionally
up. Worth stating as a small negative rather than rounding to "no change": it is not identically
zero, and the two contigs disagree on its sign, which is what a null looks like.

## Benefit: switches down on both contigs

whatshap, decomposed. `all_switches` counts a flip as two switches; the identity `S + 2F =
all_switches` holds on every arm below, which is the check that the decomposition is being read right.

| | assessed | `all_switches` | raw rate | true switches | flips |
|---|---|---|---|---|---|
| chr20 base | 58,803 | 225 | 0.3826% | **55** | 85 |
| chr20 new | 58,689 | 206 | 0.3510% | **46** | 80 |
| chr6 base | 160,138 | 459 | 0.2866% | **63** | 198 |
| chr6 new | 160,127 | 455 | 0.2841% | **59** | 198 |

-16% and -6% in true switches. The denominators move by ~0.1%, far too little for the shrinking-
denominator artefact of [switch-error-excludes-genotype-changes]; this is a real reduction.

The phasing chain also gets bigger, which is the nesting change rather than the MAPQ one:

| | het sites on the chain | of which reliable |
|---|---|---|
| chr20 base | 76,135 | 60,203 |
| chr20 new | 77,293 | **62,408** |
| chr6 base | 197,662 | 162,622 |
| chr6 new | 197,568 | **163,524** |

chr6 trades 94 het sites for 902 more reliable ones. More sites *and* fewer switches is the
outcome that was not guaranteed -- adding sites to a chain usually costs accuracy on it.

## chr20 deliverables, vg `7adcff51f`

`--preset ont -s HG002 -d 2 -t 6 --read-likelihood --phased`, anchors built twice.

| file | size | contents |
|---|---|---|
| `HG002.chr20.ont.anchors.homsplit.7adcff51f.tsv` | 286M | 654,830 anchors, 12,901,972 placements |
| `HG002.chr20.ont.anchors.7adcff51f.tsv` | 278M | 489,116 anchors, 12,903,624 placements |
| `HG002.chr20.ont.mosaic.7adcff51f.tsv` | 2.4M | 34,109 rows, `#mosaic-version 5` |
| `HG002.chr20.ont.7adcff51f.vcf` | 35M | 115,180 records |

Against the previous `ec0a6ecee` set: **+18,346** anchors with hom-split and **+14,664** without,
all of it off-reference nesting. Read counts fall ~49k, which is the MAPQ filter. `check_anchors.py`
passes on both: every read placement sits at a distinct pinned position.

## What hom splitting buys, measured on the delivered files

A site is **phased** when its anchors name a haplotype: two slots (a het, or a homozygote split by
cross-site phase), or one slot that *is* a strand (haploid). An unsplit homozygote holds one slot
carrying no haplotype claim, so it breaks a run. See [[one-slot-anchor-is-not-unphased]].

| chr20, 179,448 sites | hom-split ON | hom-split OFF |
|---|---|---|
| two slots | 163,562 (91.1%) | 76,424 (42.6%) |
| haploid, one slot | 7,450 (4.2%) | 7,165 (4.0%) |
| **phased** | **171,012 (95.3%)** | **83,589 (46.6%)** |
| unsplit hom (positioned) | 2,211 (1.2%) | 36,340 (20.3%) |
| one slot, no VCF line | 6,225 (3.5%) | 59,521 (33.2%) |

The last row is the honest residual: a lone slot-0 anchor on an off-reference site has no GT to join
to, so haploid-strand-0 and failed-to-split-homozygote are indistinguishable there. With splitting on
it is 3.5% of sites; with it off the class swells to 33.2% because it absorbs every unsplit
homozygote, so the 46.6% is a floor and the true gap is if anything wider.

Runs over the 114,328 sites that have a reference position -- 65,120 off-reference sites cannot be
placed in a linear order and are excluded rather than guessed at:

| phased runs, chr20 | hom-split ON | hom-split OFF |
|---|---|---|
| runs | 391 | 10,988 |
| span | 63.86 Mb | 34.77 Mb |
| **N50** | **1,420.1 kb** | **18.9 kb** |
| longest | 4,130.4 kb | 132.4 kb |
| sites in longest | 6,077 | 576 |

**75x on N50, 31x on the longest run**, and the span nearly doubles because a homozygote that stays
collapsed ends a run even when everything around it is phased. Re-breaking at gaps wider than the ONT
read N50 (33 kb) barely dents it -- 1,420 -> 1,272 kb -- so these runs are walkable, not an artefact
of naming haplotypes across holes.
