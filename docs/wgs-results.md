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

**Nested calling and phasing are the defaults**, and **decide-then-render** is how records are now
built: a site's genotype is settled by the linkage barrier before its record exists, so nothing is
patched after the fact. Cumulatively over the two changes: SNV F1 0.9752 -> 0.9846, ALL F1
0.9626 -> 0.9725, SV F1 0.5134 -> 0.5577. `--no-nested` and `--no-phased` restore the older
behaviour -- note `--no-phased` also disables nested calling, since a nested site's ploidy comes from
its parent's phased genotype, so it is not a control for phasing alone. See
[nested-calling-design.md](nested-calling-design.md).

**One caveat that belongs with these numbers.** The gain is a rich-panel effect: on the
4-haplotype tier-2 graphs nested calling is flat to 0.0005 *down* on ALL F1, because a
small panel enumerates few of the long collapsing ALTs it exists to break up while the
extra-records cost still applies. Parent/child ploidy incoherence, which cost 0.15% of
records a FILTER in earlier arms, is now structural rather than flagged: a nested chain
is genotyped at the ploidy its parent's settled genotype implies. The guarantee holds
wherever the parent's crossing mask can be computed, which is not everywhere -- where it
cannot, the chain is dropped rather than emitted at a ploidy its parent contradicts. The three
coherence FILTERs this paragraph used to describe are gone: a record is built from the settled
genotype, so a genotype naming an allele the record has no ALT for, and a record carrying a hom-ref
genotype, are both impossible by construction rather than flagged. Both counts are zero on all 24
contigs of this run, and they are asserted rather than reported.

## Small variants (aardvark, GT)

Current arm is **decide-then-render**: every site's genotype is settled before its record is built,
so no record is patched after the fact. The arm it replaced is kept alongside because the whole
comparison rests on it being the same binary, reads and scoring path.

| | decide-then-render | previous (inline) |
|---|---|---|
| **ALL** | TP 4,126,222  FP 92,497  FN 140,597  recall 0.9670  precision 0.9781  **F1 0.9725** | F1 0.9699 |
| **SNV** | TP 3,304,342  FP 21,547  FN 81,845  recall 0.9758  precision 0.9935  **F1 0.9846** | F1 0.9833 |
| **Indel** | TP 821,880  FP 70,950  FN 58,752  recall 0.9333  precision 0.9205  **F1 0.9269** | F1 0.9191 |
| Insertion | recall 0.9224  precision 0.9134  **F1 0.9179** | F1 0.9102 |
| Deletion | recall 0.9440  precision 0.9379  **F1 0.9409** | F1 0.9333 |

**Both precision and recall improve in every class**, which the chr20 development runs did not show --
there the gain was recall-only with false positives nearly flat. Genome-wide FP falls 6.9% and FN 9.8%.
All 23 scoreable contigs improve; none regresses.

## Structural variants (truvari, >=50 bp)

| | decide-then-render | previous (inline) |
|---|---|---|
| SV >= 50 bp | TP 14,401  FP 13,123  FN 9,716  **F1 0.5577** | F1 0.5470 |

SVs gain from recall: FN 10,361 -> 9,716 with FP rising 12,424 -> 13,123. That narrows the PanGenie gap
quoted above from 0.5488-vs-0.5739 to 0.5577-vs-0.5739 without any SV-specific work.

## Per contig

| contig | small F1 | SV F1 | small, previous | notes |
|---|---|---|---|---|
| chr1 | 0.9722 | 0.5871 | 0.9700 |  |
| chr2 | 0.9705 | 0.5594 | 0.9665 |  |
| chr3 | 0.9765 | 0.5994 | 0.9744 |  |
| chr4 | 0.9763 | 0.5841 | 0.9745 |  |
| chr5 | 0.9758 | 0.5528 | 0.9739 |  |
| chr6 | 0.9773 | 0.5820 | 0.9750 |  |
| chr7 | 0.9734 | 0.5235 | 0.9707 |  |
| chr8 | 0.9761 | 0.5789 | 0.9742 |  |
| chr9 | 0.9750 | 0.5748 | 0.9725 |  |
| chr10 | 0.9647 | 0.5174 | 0.9624 |  |
| chr11 | 0.9736 | 0.5756 | 0.9713 |  |
| chr12 | 0.9745 | 0.5785 | 0.9724 |  |
| chr13 | 0.9777 | 0.5571 | 0.9761 |  |
| chr14 | 0.9745 | 0.5897 | 0.9726 |  |
| chr15 | 0.9648 | 0.5580 | 0.9598 |  |
| chr16 | 0.9635 | 0.5119 | 0.9581 |  |
| chr17 | 0.9694 | 0.5364 | 0.9667 |  |
| chr18 | 0.9757 | 0.5413 | 0.9737 |  |
| chr19 | 0.9563 | 0.5472 | 0.9524 |  |
| chr20 | 0.9722 | 0.5258 | 0.9700 |  |
| chr21 | 0.9762 | 0.5561 | 0.9734 |  |
| chr22 | 0.9710 | 0.5260 | 0.9676 |  |
| chrX | 0.9567 | 0.4699 | 0.9494 | haploid outside PAR |
