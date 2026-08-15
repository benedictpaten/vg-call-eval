# Whole-genome results: HG002 against T2T-Q100

Called per contig on the 34-haplotype HPRC graph, `--read-likelihood` with panel
enumeration, phasing and mosaic on. chrY haploid; chrX haploid outside the
pseudoautosomal regions and diploid inside them, spliced from two runs.

**chrY is called but excluded from every total below.** The graph's CHM13 chrY path
is 57,686,750 bp where the truth's chrY runs past 62,111,784, and the two do not
correspond at any constant offset -- REF alleles match the graph's own FASTA at
chance level whether shifted by a PAR length or not at all. CHM13v2.0's chrY is
HG002-derived at 62.46 Mb and this graph's matches neither it nor GRCh38's 57.23 Mb.
Scored anyway it returns recall 0.09 at precision 0.000, which measures the
coordinate mismatch and not the caller. The calls remain in the VCF and the mosaic.

## Small variants (aardvark, GT)

- **ALL**: TP 4,039,861  FP 88,259  FN 226,958  recall 0.9468  precision 0.9786  **F1 0.9625**
- **SNV**: TP 3,239,492  FP 18,590  FN 146,695  recall 0.9567  precision 0.9943  **F1 0.9751**
- **Indel**: TP 800,369  FP 69,669  FN 80,263  recall 0.9089  precision 0.9199  **F1 0.9144**

## Structural variants (truvari, >=50 bp)

- TP 12,694  FP 12,696  FN 11,423  **F1 0.5128**

## Phasing (whatshap, autosomes)

| | |
|---|---|
| assessed het pairs | 2,442,552 |
| switches | 58,885 |
| **switch error** | **2.41%** |
| blocks | 22 — one per chromosome |
| longest block | 248,384,435 bp (chr1, end to end) |

One phase block per chromosome, so this is switch error at chromosome scale rather than over short
islands. It sits between the two tier-2 chromosomes measured on the same graph (chr20 2.30%, chr6
1.74%), which is the consistency check worth having. Per chromosome it spans 2.10% (chr6) to 3.22%
(chr17) — no single contig carries the total.

**Autosomes only, and that is a constraint rather than a choice.** whatshap requires uniform
ploidy in a file and rejects the genome-wide VCF outright -- "Inconsistent ploidy (2 and 1)" --
because chrY and non-pseudoautosomal chrX are haploid. Haploid records carry no phase to score, so
excluding them loses nothing measurable, but it does mean this number covers 22 of 24 contigs.

As on the tier-2 pages, the switch error already excludes genotype errors: whatshap assesses only
variants het and identically genotyped in both files, so the two rows above agree by construction.
Hamming is reported by the tool and is not a quality here -- over chromosome-length blocks every
switch flips everything downstream, so it approaches 50% at any non-zero switch rate.

## Why chrX scores lower

chrX is the worst scoring contig in the table below -- small-variant F1 0.9364 against an autosomal
range of 0.941-0.972, and SV F1 0.4056 against 0.46-0.55. It is worth saying what that is and what
it is not, because the obvious suspect turns out to be the opposite of the cause.

**It is entirely a precision deficit, and almost entirely in SNVs.** Recall is close to normal
(0.9335 against 0.9485 on chr7); SNV precision is 0.9569 against chr7's 0.9943, a ninefold higher
false-positive rate. So the caller is not missing chrX variants, it is emitting extra ones.

**The extra calls are ones the caller already marks as worthless.** Median GQ over chrX FPs is 3;
over chrX TPs it is 211. Filtering on GQ makes chrX's precision converge on the autosomes' and then
overtake them:

| GQ >= | chrX | chr7 | chr20 |
|---|---|---|---|
| 0 | 0.9334 | 0.9514 | 0.9523 |
| 10 | **0.9504** | 0.9459 | 0.9470 |
| 20 | **0.9490** | 0.9394 | 0.9407 |

(F1 recomputed as `2TP/(2TP+FP+FN)` for all three so the columns are comparable; the headline 0.9364
is aardvark's own asymmetric F1.)

The asymmetry in that table is the finding. A GQ>=10 filter *raises* chrX by 0.017 and *lowers* both
autosomes. Nothing is wrong with the ranking of chrX calls -- the caller separates its good chrX
calls from its bad ones at least as well as it does on an autosome. What is different is how much
mass sits in the bad bin: 13.6% of chrX calls have GQ<10, against 5.8% on chr7 and 9.3% on chr20.
chrX's *median* GQ is in fact the higher one, 184 against 92.

### The mechanism: a balanced pileup has no haploid genotype

That bimodality is the whole story, and it follows from ploidy directly. chrX FPs look like this:

```
chrX  47700009  GG>AT  GT=1  DP=13  AD=6,7  GQ=0
chrX  47700013  C>A    GT=1  DP=13  AD=6,7  GQ=0
chrX  47700026  T>A    GT=1  DP=12  AD=6,6  GQ=0
```

Half the reads support the reference, half support the alternate. On a diploid contig that is not an
error at all -- it is the definition of a heterozygote, and the model calls it confidently. Under
ploidy 1 there is no genotype that explains it: the model must pick one allele, the pileup gives it
almost no reason to prefer either, and the likelihood gap between the top two genotypes collapses.
The result is a near coin-flip call carrying GQ 0.

The numbers say this is what chrX FPs are, not an anecdote. Minor-allele fraction across chrX calls:

| | n | median minor-AF | share with minor-AF > 0.3 |
|---|---|---|---|
| TP | 83,193 | 0.000 | 1.4% |
| FP | 3,863 | 0.333 | 54.1% |
| FP, in the two hotspots below | 1,559 | 0.375 | 72.2% |

True hemizygous calls are unanimous. False ones are split down the middle.

### Where the split pileups come from

Two loci of about 200 kb each -- chrX:47.6-47.8 Mb and chrX:48.7-48.9 Mb -- are 0.26% of the
chromosome by length and produce **29% of its false positives** (1,559 of 5,422). They are 90%+ FP
internally and overcall the truth 8.5x and 4.2x respectively. Excluding just those two windows moves
chrX from 0.9334 to 0.9449.

Both are multi-copy. Counting 31-mers against chrX and chrY together, the fraction of each region's
k-mers with a copy outside it:

| region | distinct k-mers per 200 kb | k-mers with an external copy |
|---|---|---|
| chrX:47.6-47.8 Mb | 174,559 | 43.4% |
| chrX:48.7-48.9 Mb | 43,839 | 34.7% |
| control chrX:80.0-80.2 Mb | 194,791 | 11.3% |
| control chrX:140.0-140.2 Mb | 196,800 | 13.6% |

Three to four times the control paralogy rate, and the second region is additionally a tandem array
-- only 22% of its 200 kb of positions are distinct k-mers, against 97% for the controls. So reads
from a paralogous copy land on these loci, carry that copy's divergent bases, and manufacture a
balanced pileup where the sample's actual X is homozygous. Haploid calling then converts each into a
hom-alt FP.

### What it is not

**Not the haploid path being wrong.** The opposite: chrX was also called end-to-end at ploidy 2, and
scoring that pass against the same truth gives F1 0.8957 against the shipped 0.9362, with 4,746
genotype errors that the haploid pass does not make. Haploid handling is worth about +0.04 F1 here.

**Not a thinner panel.** The chrX subgraph carries 194 haplotype-sense paths, against 190 for chr7
and 127 for chr20.

**Not chrX being intrinsically harder to genotype, once depth is accounted for.** chrX's median DP
is 14 against 29 on the autosomes -- exactly half, which is what a male X should be. Stratifying
precision by depth, chrX is not the weaker contig at any matched depth:

| DP | chrX precision | chr7 precision |
|---|---|---|
| 6-10 | 0.9236 (n=12,022) | 0.7624 (n=404) |
| 10-14 | 0.9502 (n=27,346) | 0.8666 (n=1,642) |
| 14-18 | 0.9512 (n=26,875) | 0.9400 (n=6,867) |
| 24-32 | 0.9140 (n=2,896) | 0.9861 (n=97,431) |

Read that comparison with care rather than as a claim that chrX is easier. The `n` columns show why:
on chr7 a call at DP 12 is *abnormal* -- coverage collapsed there for a reason, usually a repeat --
so chr7's low-depth calls are a self-selected hard population. On chrX DP 12 is simply Tuesday. The
honest reading is the weaker one: chrX is not anomalous once depth is controlled for, and its whole
distribution sits in the depth range where every contig does worse.

### Summary

chrX's deficit is roughly two effects and no defect. About 29% of its excess FPs come from two
paralogous loci where stray reads build balanced pileups that ploidy 1 cannot represent; the rest is
half coverage moving the whole GQ distribution down into the range where marginal calls survive
emission. Both are visible in GQ, and a GQ>=10 filter recovers chrX to better-than-autosomal F1.

The one actionable item: the default emission threshold is tuned for diploid contigs at full depth,
and is too permissive for a haploid contig at half depth. A ploidy-aware or depth-aware emission
floor would remove most of these without costing recall. Not changed here, since it would alter
autosomal output too, and this run is the baseline everything else is measured against.

The SV gap (0.4056 against 0.4865/0.4944) rests on 496 truth SVs for the whole chromosome, against
1,614 for chr7, so it is a much noisier number than the small-variant one and is not decomposed here.

## Per contig

| contig | small F1 | SV F1 | notes |
|---|---|---|---|
| chr1 | 0.9649 | 0.5366 |  |
| chr2 | 0.9598 | 0.5285 |  |
| chr3 | 0.9675 | 0.5509 |  |
| chr4 | 0.9650 | 0.5243 |  |
| chr5 | 0.9639 | 0.5117 |  |
| chr6 | 0.9690 | 0.5268 |  |
| chr7 | 0.9638 | 0.4865 |  |
| chr8 | 0.9656 | 0.5297 |  |
| chr9 | 0.9631 | 0.5187 |  |
| chr10 | 0.9592 | 0.4807 |  |
| chr11 | 0.9641 | 0.5148 |  |
| chr12 | 0.9654 | 0.5227 |  |
| chr13 | 0.9688 | 0.5077 |  |
| chr14 | 0.9636 | 0.5495 |  |
| chr15 | 0.9523 | 0.5324 |  |
| chr16 | 0.9501 | 0.4612 |  |
| chr17 | 0.9593 | 0.5018 |  |
| chr18 | 0.9716 | 0.4708 |  |
| chr19 | 0.9410 | 0.4915 |  |
| chr20 | 0.9646 | 0.4944 |  |
| chr21 | 0.9592 | 0.5213 |  |
| chr22 | 0.9602 | 0.4815 |  |
| chrX | 0.9364 | 0.4056 |  |
| chrY | 0.0039 | - | excluded: reference mismatch with truth; truvari_error |
