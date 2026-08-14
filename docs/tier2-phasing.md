# Phasing and the mosaic genome

`vg call --phased` emits phased genotypes and a phase set; `--mosaic-out` writes the inferred
genome as a mosaic of panel haplotypes. Both come from the same thing: the linkage layer's most
probable path of haplotype pairs through a chain.

## Where the phase comes from, and why that matters for reading these numbers

Not from reads. A read-based phaser links two heterozygous sites when one read, or one fragment,
spans both; its blocks are therefore read- or fragment-length, and it phases nothing across a gap
no read crosses. This phases from the **panel**: the Li–Stephens layer asks which pair of panel
haplotypes best explains the calls, and the order of that pair is the phase. Linkage carries across
any distance the transition model allows, so a phase block is a whole chain.

That difference is the first thing to hold on to when reading a switch error. **A phaser can make
switch error arbitrarily small by emitting shorter blocks** — in the limit, one block per site has
no switch errors at all and says nothing. These blocks are chromosome-length, which is the hardest
case, so the two numbers are only meaningful together.

## Results

HG002 against HPRC graphs, truth from the T2T-Q100 assembly-derived benchmark, which is fully
phased. HG002 is excluded from the graph, so this measures imputation against a panel that does not
contain the sample.

| dataset | panel | assessed pairs | switches | **switch error** | block N50 |
|---|---|---|---|---|---|
| chr20, 4-hap | 4 | 57,350 | 2,303 | **3.43%** | 66,205,242 |
| chr20, 34-hap | 34 | 56,330 | 1,580 | **2.30%** | 66,208,830 |
| chr6, 34-hap | 34 | 160,026 | 3,365 | **1.74%** | 172,123,900 |

One block per chromosome in every case — the whole contig.

**The panel-size contrast is the result that says the metric measures what it claims.** Same
chromosome, same reads, same caller: going from 4 panel haplotypes to 34 takes switch error from
3.43% to 2.30%, a third fewer. That is the mechanism showing through. A mosaic built from more
haplotypes fits the sample with fewer switches, and if the number had *not* moved, it would have
been measuring something else.

For scale: statistical phasers working from thousands of reference haplotypes reach roughly 0.5–2%
chromosome-wide. Reaching 1.7–2.3% from 34 haplotypes is a reasonable place to be, and the panel is
the obvious lever.

## What the switch error does not include

`whatshap compare` assesses only variants that are heterozygous **and identically genotyped** in
both files. So sites we call wrongly do not count against the phasing — they leave the denominator
entirely.

This was not the intent. The benchmark was written to report two numbers, all het sites and
correctly-genotyped ones, on the reasoning that the phasing is constrained to our own calls and so
a mis-genotyped site could force a switch that is really a calling error. Applying that filter by
hand changes nothing: on chr20 it drops 815 het-in-both sites whose genotypes disagree, and the
assessed pair count and switch count are identical to the digit. The two rows are kept in the
output so that the equality is visible rather than assumed, and so a scorer that intersects
differently would show up.

The honest reading is that this is phase error **given** a correct genotype, and a user calling and
phasing in one pass experiences both error sources.

## Hamming distance is not a quality here

The tables above omit it deliberately. Over a single chromosome-length block, every switch flips
the relative phase of everything downstream, so blockwise Hamming approaches 50% at any non-zero
switch rate — chr20-34hap reports 27,780 of 56,330, which is 49%. It is measuring block length, not
phasing. Switch error is the metric; Hamming would only be informative over short blocks, which is
the regime this deliberately is not in.

## The mosaic file

The phasing is piecewise, which is what makes a compact format possible. Measured over the emitted
likelihoods before any of this was built:

| | sites | switch points | share of sites |
|---|---|---|---|
| chr20, 34-hap | 105,251 | 2,064 | 1.96% |
| chr6, 34-hap | 284,529 | 4,434 | 1.56% |

Run in the caller itself, chr20-34hap produces **3,673 segments over 105,251 sites, in a 255 KB
file** — a factor of 29 against a per-site format, with a median run of 33 sites and a longest of
696.

**That is 1.8x more segments than the offline measurement predicted**, and the gap is worth
recording rather than quietly replacing the estimate. The offline harness reconstructs the panel by
joining `vg deconstruct` output to the call set on allele traversals, and it failed to map 10.4% of
chr20 records; an unmapped record has no panel row, so every state is free and the path has no
reason to move at it. The real implementation has a panel row at every site. The offline number was
therefore a lower bound on switching, not an estimate of it — which does not change the conclusion
(255 KB is still nothing) but does mean the offline harness should not be quoted as a predictor of
segment counts.

```
#mosaic-version	1
#graph	chr20_0_chr20.gbz
#sample	HG002
#decoding	constrained-viterbi
#H	contig	strand	ref_start	ref_end	start_node	end_node	hap_index	haplotype	sites
H	chr20	0	603	3605	114819056	114842605	9	recombination#30	60
```

Anchored on **node IDs** rather than reference positions: a node ID is intrinsic to the graph, while
a position is a statement about one reference path. A consumer reconstructs a haplotype by walking
the named GBWT sequence from `start_node` to `end_node`. `*` in the haplotype column means the panel
does not explain that strand there.

### Why this is small, measured rather than estimated

The mosaic and an explicit path list describe the *same two walks* through the graph. They differ in
how they say it: the mosaic names a panel haplotype and a node range and lets the graph supply the
steps, while an explicit list enumerates every node.

Extracting one chr20 haplotype as GAF (`vg paths -A`) gives the enumeration exactly: **2,031,992
steps in a 20.3 MB record**, so the two strands are ~4.06 M node references and **~40.6 MB** of
text. Against the mosaic's 255 KB that is a factor of **159**.

(The reference path averages 66,210,255 / 2,031,992 = 32.6 bp per node, close to the graph-wide 26.6
bp, so the earlier back-of-envelope estimate of ~45 MB was sound — this replaces it with the
measurement.)

The trade is that the mosaic is written **by reference**: it cannot be read without the GBZ it names,
where an explicit list is self-contained. That is the whole reason the header carries the graph name.

Only one segment on chr20-34hap carries `*`, matching the two sites the run reports as having a
strand the panel cannot explain.

## Constraining to the calls is free

The phasing is a *constrained* Viterbi: at each site the states are restricted to those spelling the
called genotype, so the emitted genome agrees with the emitted VCF by construction. The cost of that
constraint was the number this was designed around, on the reasoning that forcing the path through a
mis-called site would buy switches.

It does not. Measured on chr20, 2,115 switches unconstrained against 2,064 constrained, and on chr6
4,617 against 4,434 — consistency costs 0.98× and 0.96×, which is to say nothing at all. The
reasoning was wrong in a way worth recording: restricting the state space necessarily lowers the
path *probability*, but switch count is not monotone in it, because removing states takes away
opportunities to switch as readily as it forces them.

## Caveats

- **One sample, two chromosomes, two panel sizes.** The same scope limit as every other number in
  this harness, and the panel-dependence is sharper here than anywhere else: phasing quality *is*
  panel quality, as the 4-vs-34 contrast shows directly.
- **The 2 in 105,251.** That is how many chr20-34hap sites have a strand the panel cannot explain
  (58 on the 4-haplotype graph, 7 on chr6-34hap). Where that happens the phase either side rests on
  the transition model alone. It is rare enough not to qualify the result, and it is reported per
  run so it cannot become rare-by-assumption.
- **Phase sets are per chain.** Phase is not comparable across chains, and `PS` says which is which.
