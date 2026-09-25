# Phasing and the mosaic genome

`vg call --phased` emits phased genotypes and a phase set; `--mosaic-out` writes the inferred
genome as a mosaic of panel haplotypes. By default both come from the same thing: the linkage
layer's most probable path of haplotype pairs through a chain. Under `--preset ont` the reads
re-orient that phase, and the mosaic then describes the read-derived phase as panel segments.

## Where the phase comes from, and why that matters for reading these numbers

For short reads, not from reads. A read-based phaser links two heterozygous sites when one read,
or one fragment, spans both; its blocks are therefore read- or fragment-length, and it phases
nothing across a gap no read crosses. The default phases from the **panel**: the Li–Stephens layer
asks which pair of panel haplotypes best explains the calls, and the order of that pair is the
phase. Linkage carries across any distance the transition model allows, so a phase block is a whole
chain. `--read-phasing` is off by default, so short-read phase is panel-only.

Under `--preset ont` the reads come back in. `--read-phasing` links het sites through the reads
that span them and re-orients the panel phase to agree, and where the read chain breaks the pieces
are relinked against the panel phase rather than left as separate blocks, so the ONT blocks are
whole chromosomes too (one block on chr20, one on chr6). `--regenotype`, also on under the preset,
lets that read phase move genotypes.

Block length is the first thing to hold on to when reading a switch error. **A phaser can make
switch error arbitrarily small by emitting shorter blocks** — in the limit, one block per site has
no switch errors at all and says nothing. These blocks are chromosome-length, which is the hardest
case, so the two numbers are only meaningful together.

## Results

HG002, truth from T2T-Q100 v1.1 (GIAB defrabb V0.019 draft benchmark), which is fully phased.
Short reads run on the 32-haplotype hap32 graph (34 panel haplotypes with the CHM13 and GRCh38
paths), as the tier-2 readlik arm on the current binary; ONT runs `--preset ont` on the
16-haplotype E821 graph (18 panel). HG002 is excluded from the hap32 graph, so the short-read rows measure imputation
against a panel that does not contain the sample.

**Switch error** is switches over assessed pairs, from `whatshap compare` (`all_switches` /
`all_assessed_pairs`), where a pair is two consecutive heterozygous sites genotyped identically in
the calls and the truth. The switchflip rate is whatshap's variant that counts an isolated
single-site flip as one event instead of two switches; it is shown so other figures can be placed,
and switch error is the headline.

| dataset | graph (panel) | phase from | assessed pairs | switches | **switch error** | switchflip rate | block N50 |
|---|---|---|---|---|---|---|---|
| chr20, short reads | hap32 (34) | panel | 59,033 | 1,753 | **2.97%** | 2.44% | 66,209,624 |
| chr6, short reads | hap32 (34) | panel | 164,094 | 3,686 | **2.25%** | 1.86% | 172,123,900 |
| chr20, ONT preset | E821 (18) | reads + panel | 58,743 | 172 | **0.293%** | 0.169% | 66,206,855 |
| chr6, ONT preset | E821 (18) | reads + panel | 160,078 | 340 | **0.212%** | 0.114% | 172,091,966 |

One block per chromosome in every case — the whole contig. The short-read and ONT rows differ in
graph, panel and phase source at once, so the roughly tenfold gap is not a test of the reads alone.
Raising the `--mismap-max` default from 0.7 to 0.95 is neutral-to-better here: chr20 3.01% → 2.97%
(28 fewer switches), chr6 2.249% → 2.246%.

(2026-08-13: this page gave 2.30% and 1.74%, whatshap switchflip rates printed beside switch
counts; as switches over pairs those older-binary runs were 2.80% and 2.10%, not comparable with
the table.) Against the oldest harness readlik arm still on disk (2026-08-23) there is no detectable
regression: 2.97% against 3.04% over each arm's own pairs (59,033 and 59,249), and on the 58,754
pairs the two share, level within noise (2.95% against 2.92%, 18 more switches for the current
binary).

**The panel-size contrast is the result that says the metric measures what it claims.** Same
chromosome, same reads, same caller, on the 2026-08-13 binary: going from 4 panel haplotypes to 34
took switch error from 4.02% (2,303 / 57,350) to 2.80% (1,580 / 56,330), about 30% fewer. That is
the mechanism showing through. A mosaic built from more haplotypes fits the sample with fewer
switches, and if the number had *not* moved, it would have been measuring something else. No
4-haplotype run exists on the current binary, so this pair stands on its own and is not comparable
with the table above.

For scale: statistical phasers working from thousands of reference haplotypes reach roughly 0.5–2%
chromosome-wide. Short reads reach 2.25–2.97% from 34 panel haplotypes, and for them the panel is
the obvious lever. For long reads the reads are the lever: the ONT preset reaches 0.293% on chr20
and 0.212% on chr6 from 18 panel haplotypes, with each contig still one block.

## What the switch error does not include

`whatshap compare` assesses only variants that are heterozygous **and identically genotyped** in
both files. So sites we call wrongly do not count against the phasing — they leave the denominator
entirely.

This was not the intent. The benchmark was written to report two numbers, all het sites and
correctly-genotyped ones, on the reasoning that the phasing is constrained to our own calls and so
a mis-genotyped site could force a switch that is really a calling error. Applying that filter by
hand changes nothing: on every run in the table the assessed pair count and switch count are
identical to the digit (chr20 short reads 59,033 and 1,753 both ways, ONT 58,743 and 172). The two
rows are kept in the output so that the equality is visible rather than assumed, and so a scorer
that intersects differently would show up.

The honest reading is that this is phase error **given** a correct genotype, and a user calling and
phasing in one pass experiences both error sources.

It also makes the denominator depend on the arm. Under `--preset ont`, `--regenotype` moves
genotypes using the read phase, and whatshap assesses only sites whose genotype matches the truth,
so the assessed set moves with the genotypes: a re-genotyping arm can lower its switch error by
shrinking the denominator as well as by phasing better.

## Hamming distance is not a quality here

The tables above omit it deliberately. Over a single chromosome-length block, every switch flips
the relative phase of everything downstream, so blockwise Hamming says more about where the
switches fall than how many there are. Short reads report 45.15% on chr20 (26,656 of 59,033) and
48.77% on chr6; the ONT preset, at about a tenth of the switch rate, still reports 34.25% and
30.54%. It is measuring block length, not phasing. Switch error is the metric; Hamming would only
be informative over short blocks, which is the regime this deliberately is not in.

## The mosaic file

The phasing is piecewise, which is what makes a compact format possible. Measured over the emitted
likelihoods before any of this was built:

| | sites | switch points | share of sites |
|---|---|---|---|
| chr20, 34-hap | 105,251 | 2,064 | 1.96% |
| chr6, 34-hap | 284,529 | 4,434 | 1.56% |

Run in the caller itself, on the same chr20 short-read run as the Results row (hap32, 34 panel),
the mosaic is **6,517 segments over 111,066 sites, in a 482,977-byte file**, with a median run of
16 sites and a longest of 689.

**The first in-caller build (2026-08-13) already produced 1.8x more segments than the offline
measurement predicted** (3,673 over 105,251 sites), and the gap is worth recording rather than
quietly replacing the estimate. The offline harness reconstructs the panel by joining
`vg deconstruct` output to the call set on allele traversals, and it failed to map 10.4% of
chr20 records; an unmapped record has no panel row, so every state is free and the path has no
reason to move at it. The real implementation has a panel row at every site. The offline number was
therefore a lower bound on switching, not an estimate of it — which does not change the conclusion
(a file under half a megabyte is still nothing) but does mean the offline harness should not be
quoted as a predictor of segment counts.

```
#mosaic-version	5
#graph	…/work/wgs-tt/chr20/chr20.gbz
#sample	HG002
#reference	CHM13#0#chr20
#decoding	constrained-viterbi
#patch	reference
#nested	kept
#unexplained	connected
#haplotype	0	CHM13#0
#haplotype	1	GRCh38#0
#haplotype	2	recombination#17
#H	contig	strand	fragment	ref_start	ref_end	start_node	end_node	hap_index	haplotype	sites	gbwt_offset
H	chr20	0	0	22	72	229637742	229637788	6	recombination#22	3	22
H	chr20	0	0	116	439	229637788	229638012	7	recombination#25	13	28
```

(The `#note` lines, which explain the columns in the file itself, and the remaining `#haplotype`
lines are omitted.)

Anchored on **node IDs** rather than reference positions: a node ID is intrinsic to the graph, while
a position is a statement about one reference path, so `ref_start`/`ref_end` are advisory, in the
`#reference` path's coordinates. `start_node` and `end_node` are oriented (id × 2 + is_reverse),
and `(start_node, gbwt_offset)` is a GBWT position: a consumer reconstructs a haplotype by
extracting it there and following the named sequence to `end_node`. A segment's `end_node` is the
next segment's `start_node`, so the rows of one `(contig, strand, fragment)` are one walk. `*` in
the haplotype column means the panel cannot name a haplotype for that stretch; `ref` in `hap_index`
marks a stretch filled with the reference because no panel haplotype could be carried across it.

### Why this is small, measured rather than estimated

The mosaic and an explicit path list describe the *same two walks* through the graph. They differ in
how they say it: the mosaic names a panel haplotype and a node range and lets the graph supply the
steps, while an explicit list enumerates every node.

Extracting one chr20 haplotype as GAF (`vg paths -A`) gives the enumeration exactly: **2,031,992
steps in a 20.3 MB record**, so the two strands are ~4.06 M node references and **~40.6 MB** of
text. Against the mosaic's 482,977 bytes that is a factor of about **84**.

(The reference path averages 66,210,255 / 2,031,992 = 32.6 bp per node, close to the graph-wide 26.6
bp, so the earlier back-of-envelope estimate of ~45 MB was sound — this replaces it with the
measurement.)

The trade is that the mosaic is written **by reference**: it cannot be read without the GBZ it names,
where an explicit list is self-contained. That is the whole reason the header carries the graph name.

No segment of the chr20 short-read mosaic carries `*`, although the run reports 8 sites with a
strand the panel cannot explain: by default (`#unexplained connected`) the flanking haplotype is
carried straight through such a stretch, so the strand stays one walk, and
`--mosaic-break-unexplained` leaves a `*` row there instead. Separately, where no panel haplotype
can be carried across a stretch, the caller fills it with the reference (`#patch reference`;
`--no-mosaic-patch-gaps` leaves the gap instead). 128 rows are `ref`: 6 with `sites` = `.` that
fill a boundary between two segments and cover no called site, and 122 with a site count that
replace a haplotype the graph does not carry across that segment.

## Constraining to the calls is free

The phasing is a *constrained* Viterbi: at each site the states are restricted to those spelling the
called genotype, so the decoded path agrees with the emitted VCF by construction. The cost of that
constraint was the number this was designed around, on the reasoning that forcing the path through a
mis-called site would buy switches.

It does not. Measured on chr20, 2,115 switches unconstrained against 2,064 constrained, and on chr6
4,617 against 4,434 — consistency costs 0.98× and 0.96×, which is to say nothing at all. The
reasoning was wrong in a way worth recording: restricting the state space necessarily lowers the
path *probability*, but switch count is not monotone in it, because removing states takes away
opportunities to switch as readily as it forces them.

The written mosaic departs from the decoded path in two places: it substitutes the reference for a
haplotype the graph does not carry across a segment (122 segments on chr20), and it carries the
flanking haplotype through a site the panel cannot explain (8 on chr20), spelling that haplotype's
alleles there rather than the called ones.

This is a panel-phasing result. Under `--preset ont` the reads re-orient the phase, and
`--regenotype` can move genotypes, after the constrained decode, so it says nothing about read
phasing.

## Caveats

- **One sample, two chromosomes, three panel sizes (4, 18, 34), two read types.** The same scope
  limit as every other number in this harness, and the panel-dependence is sharper here than
  anywhere else: for panel phasing, phasing quality *is* panel quality, as the 4-vs-34 contrast
  shows directly. Here the 18-haplotype panel appears only with ONT reads and read phasing, so
  panel and reads cannot be separated in the ONT rows.
- **The 8 in 111,066.** That is how many chr20 short-read sites have a strand the panel cannot
  explain (26 of 294,746 on chr6; under the ONT preset, 202 of 114,084 on chr20 and 343 of 301,168
  on chr6). Where that happens the panel phase either side rests on the transition model alone. It
  is rare enough not to qualify the result, and it is reported per run so it cannot become
  rare-by-assumption.
- **Phase sets are per chain.** Phase is not comparable across chains, and `PS` says which is which.
