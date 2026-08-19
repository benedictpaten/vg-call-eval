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
[pangenie-comparison.md](pangenie-comparison.md). Briefly, on the autosomes vg is ahead on small
variants (ALL F1 0.9630 against 0.9505, driven by indels) and PanGenie is ahead on structural
variants (0.5739 against 0.5152) -- though 42% of that gap is a representation artefact
rather than an evidence one, and the recall remainder sits at 50-300 bp; see
[sv-delta.md](sv-delta.md).

**The mosaic** this run also emits: 182,328 segments over 5,041,066 sites, 14 MB.
See wgs-performance.md for why assembling it is not `cat`.

**Nested calling and phasing are the defaults** as of this run, which is why these
numbers moved: SNV F1 0.9752 -> 0.9833, ALL F1 0.9626 -> 0.9699, SV F1 0.5134 -> 0.5467,
with 59,413 SNV false negatives recovered, at no runtime or memory cost. `--no-nested`
and `--no-phased` restore the old behaviour. See
[nested-calling-design.md](nested-calling-design.md).

**Two caveats that belong with these numbers.** The gain is a rich-panel effect: on the
4-haplotype tier-2 graphs nested calling is flat to 0.0005 *down* on ALL F1, because a
small panel enumerates few of the long collapsing ALTs it exists to break up while the
extra-records cost still applies. And 0.15% of records carry a ploidy-coherence FILTER
(`nested_diploid` 2,458, `nested_unreachable` 5,000, `nested_haploid` 0), meaning the
child's ploidy and its parent's final genotype disagree; those calls are flagged rather
than corrected.

## Small variants (aardvark, GT)

- **ALL**: TP 4,111,896  FP 100,316  FN 154,923  recall 0.9637  precision 0.9762  **F1 0.9699**
- **SNV**: TP 3,298,256  FP 24,108  FN 87,931  recall 0.9740  precision 0.9927  **F1 0.9833**
- **Indel**: TP 813,640  FP 76,208  FN 66,992  recall 0.9239  precision 0.9144  **F1 0.9191**

## Structural variants (truvari, >=50 bp)

- TP 13,781  FP 12,513  FN 10,336  **F1 0.5467**

## Per contig

| contig | small F1 | SV F1 | notes |
|---|---|---|---|
| chr1 | 0.9700 | 0.5666 |  |
| chr2 | 0.9665 | 0.5615 |  |
| chr3 | 0.9744 | 0.5879 |  |
| chr4 | 0.9745 | 0.5707 |  |
| chr5 | 0.9739 | 0.5502 |  |
| chr6 | 0.9749 | 0.5662 |  |
| chr7 | 0.9708 | 0.5152 |  |
| chr8 | 0.9744 | 0.5636 |  |
| chr9 | 0.9724 | 0.5447 |  |
| chr10 | 0.9624 | 0.4982 |  |
| chr11 | 0.9713 | 0.5562 |  |
| chr12 | 0.9725 | 0.5678 |  |
| chr13 | 0.9761 | 0.5461 |  |
| chr14 | 0.9726 | 0.5750 |  |
| chr15 | 0.9598 | 0.5538 |  |
| chr16 | 0.9583 | 0.4937 |  |
| chr17 | 0.9667 | 0.5294 |  |
| chr18 | 0.9737 | 0.5115 |  |
| chr19 | 0.9525 | 0.5332 |  |
| chr20 | 0.9698 | 0.5131 |  |
| chr21 | 0.9733 | 0.5426 |  |
| chr22 | 0.9674 | 0.5025 |  |
| chrX | 0.9494 | 0.4617 |  |
| chrY | 0.0034 | - | excluded: reference mismatch with truth; truvari_error |
