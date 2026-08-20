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
ahead on structural variants (0.5739 against 0.5488). What is inside that SV gap, and whether
nested calling reached it: [sv-residual-errors.md](sv-residual-errors.md).

**The mosaic** this run also emits: 180,858 segments over 5,037,872 sites, 14 MB.
See wgs-performance.md for why assembling it is not `cat`.

**Nested calling and phasing are the defaults** as of this run, which is why these
numbers moved: SNV F1 0.9752 -> 0.9833, ALL F1 0.9626 -> 0.9699, SV F1 0.5134 -> 0.5470,
with 58,336 SNV false negatives recovered, at no runtime or memory cost. `--no-nested`
and `--no-phased` restore the old behaviour. See
[nested-calling-design.md](nested-calling-design.md).

**One caveat that belongs with these numbers.** The gain is a rich-panel effect: on the
4-haplotype tier-2 graphs nested calling is flat to 0.0005 *down* on ALL F1, because a
small panel enumerates few of the long collapsing ALTs it exists to break up while the
extra-records cost still applies. Parent/child ploidy incoherence, which cost 0.15% of
records a FILTER in earlier arms, is now structural rather than flagged: a nested chain
is genotyped at the ploidy its parent's settled genotype implies. The guarantee holds
wherever the parent's crossing mask can be computed, which is not everywhere -- where it
cannot, the chain keeps its sweep-time ploidy and the coherence FILTERs stay live to say
so. They fire on no record in this run, and a handful in some single-contig runs, so
treat a nonzero count as a pointer at those chains rather than as a regression.

## Small variants (aardvark, GT)

- **ALL**: TP 4,110,906  FP 99,341  FN 155,913  recall 0.9635  precision 0.9764  **F1 0.9699**
- **SNV**: TP 3,297,737  FP 23,625  FN 88,450  recall 0.9739  precision 0.9929  **F1 0.9833**
- **Indel**: TP 813,169  FP 75,716  FN 67,463  recall 0.9234  precision 0.9148  **F1 0.9191**

## Structural variants (truvari, >=50 bp)

- TP 13,756  FP 12,424  FN 10,361  **F1 0.5470**

## Per contig

| contig | small F1 | SV F1 | notes |
|---|---|---|---|
| chr1 | 0.9700 | 0.5637 |  |
| chr2 | 0.9665 | 0.5590 |  |
| chr3 | 0.9744 | 0.5877 |  |
| chr4 | 0.9745 | 0.5688 |  |
| chr5 | 0.9739 | 0.5488 |  |
| chr6 | 0.9750 | 0.5666 |  |
| chr7 | 0.9707 | 0.5172 |  |
| chr8 | 0.9742 | 0.5636 |  |
| chr9 | 0.9725 | 0.5518 |  |
| chr10 | 0.9624 | 0.4975 |  |
| chr11 | 0.9713 | 0.5563 |  |
| chr12 | 0.9724 | 0.5643 |  |
| chr13 | 0.9761 | 0.5508 |  |
| chr14 | 0.9726 | 0.5796 |  |
| chr15 | 0.9598 | 0.5557 |  |
| chr16 | 0.9581 | 0.4948 |  |
| chr17 | 0.9667 | 0.5274 |  |
| chr18 | 0.9737 | 0.5105 |  |
| chr19 | 0.9524 | 0.5382 |  |
| chr20 | 0.9700 | 0.5140 |  |
| chr21 | 0.9734 | 0.5467 |  |
| chr22 | 0.9676 | 0.5051 |  |
| chrX | 0.9494 | 0.4598 |  |
| chrY | 0.0034 | - | excluded: reference mismatch with truth; truvari_error |
