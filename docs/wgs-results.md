# Whole-genome results: HG002 against T2T-Q100

Called per contig on the 34-haplotype HPRC graph, `--read-likelihood` with panel
enumeration, phasing and mosaic on. chrY haploid; chrX haploid outside the
pseudoautosomal regions and diploid inside them, in one run via --ploidy-bed.

**chrY is called but excluded from every total below.** The graph's CHM13 chrY path
is 57,686,750 bp where the truth's chrY runs past 62,111,784, and the two do not
correspond at any constant offset -- REF alleles match the graph's own FASTA at
chance level whether shifted by a PAR length or not at all. CHM13v2.0's chrY is
HG002-derived at 62.46 Mb and this graph's matches neither it nor GRCh38's 57.23 Mb.
Scored anyway it returns recall 0.09 at precision 0.000, which measures the
coordinate mismatch and not the caller. The calls remain in the VCF and the mosaic.

**How to run this, and how long it takes**: [wgs-performance.md](wgs-performance.md).
**Behaviour across coverage and ploidy**: [coverage.md](coverage.md).

**Compared against PanGenie on the same graph and reads**: see
[pangenie-comparison.md](pangenie-comparison.md). Briefly, on the autosomes vg is ahead on every
small-variant class on both recall and precision (ALL F1 0.9703 against 0.9505) and PanGenie is
ahead on structural variants (0.5739 against 0.5486). What is inside that SV gap, and whether
nested calling reached it: [sv-residual-errors.md](sv-residual-errors.md).

**The mosaic** this run also emits: 182,950 segments over 5,037,820 sites, 14 MB.
See wgs-performance.md for why assembling it is not `cat`.

**Nested calling and phasing are the defaults** as of this run, which is why these
numbers moved: SNV F1 0.9752 -> 0.9833, ALL F1 0.9626 -> 0.9699, SV F1 0.5134 -> 0.5468,
with 58,289 SNV false negatives recovered, at no runtime or memory cost. `--no-nested`
and `--no-phased` restore the old behaviour. See
[nested-calling-design.md](nested-calling-design.md).

**One caveat that belongs with these numbers.** The gain is a rich-panel effect: on the
4-haplotype tier-2 graphs nested calling is flat to 0.0005 *down* on ALL F1, because a
small panel enumerates few of the long collapsing ALTs it exists to break up while the
extra-records cost still applies. Parent/child ploidy incoherence, which cost 0.15% of
records a FILTER in earlier arms, is gone by construction: a nested chain is genotyped
at the ploidy its parent's settled genotype implies, so the three coherence FILTERs are
now an invariant check that fires zero times genome-wide.

## Small variants (aardvark, GT)

- **ALL**: TP 4,110,836  FP 99,548  FN 155,983  recall 0.9634  precision 0.9764  **F1 0.9699**
- **SNV**: TP 3,297,690  FP 23,714  FN 88,497  recall 0.9739  precision 0.9929  **F1 0.9833**
- **Indel**: TP 813,146  FP 75,834  FN 67,486  recall 0.9234  precision 0.9147  **F1 0.9190**

## Structural variants (truvari, >=50 bp)

- TP 13,765  FP 12,467  FN 10,352  **F1 0.5468**

## Per contig

| contig | small F1 | SV F1 | notes |
|---|---|---|---|
| chr1 | 0.9700 | 0.5645 |  |
| chr2 | 0.9665 | 0.5591 |  |
| chr3 | 0.9744 | 0.5859 |  |
| chr4 | 0.9745 | 0.5686 |  |
| chr5 | 0.9739 | 0.5494 |  |
| chr6 | 0.9750 | 0.5672 |  |
| chr7 | 0.9707 | 0.5159 |  |
| chr8 | 0.9742 | 0.5639 |  |
| chr9 | 0.9724 | 0.5515 |  |
| chr10 | 0.9624 | 0.4973 |  |
| chr11 | 0.9712 | 0.5560 |  |
| chr12 | 0.9723 | 0.5647 |  |
| chr13 | 0.9761 | 0.5511 |  |
| chr14 | 0.9726 | 0.5787 |  |
| chr15 | 0.9598 | 0.5559 |  |
| chr16 | 0.9581 | 0.4927 |  |
| chr17 | 0.9666 | 0.5267 |  |
| chr18 | 0.9736 | 0.5097 |  |
| chr19 | 0.9524 | 0.5382 |  |
| chr20 | 0.9699 | 0.5133 |  |
| chr21 | 0.9733 | 0.5481 |  |
| chr22 | 0.9675 | 0.5030 |  |
| chrX | 0.9493 | 0.4601 |  |
| chrY | 0.0034 | - | excluded: reference mismatch with truth; truvari_error |
