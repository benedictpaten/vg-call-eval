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

**Compared against PanGenie on the same graph and reads**: see
[pangenie-comparison.md](pangenie-comparison.md). Briefly, on the autosomes vg is ahead on small
variants (ALL F1 0.9630 against 0.9505, driven by indels) and PanGenie is ahead on structural
variants (0.5739 against 0.5152).

## Small variants (aardvark, GT)

- **ALL**: TP 4,039,700  FP 86,970  FN 227,119  recall 0.9468  precision 0.9789  **F1 0.9626**
- **SNV**: TP 3,239,401  FP 17,958  FN 146,786  recall 0.9567  precision 0.9945  **F1 0.9752**
- **Indel**: TP 800,299  FP 69,012  FN 80,333  recall 0.9088  precision 0.9206  **F1 0.9147**

## Structural variants (truvari, >=50 bp)

- TP 12,686  FP 12,615  FN 11,431  **F1 0.5134**

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
| chrX | 0.9422 | 0.4260 |  |
| chrY | 0.0047 | - | excluded: reference mismatch with truth; truvari_error |
