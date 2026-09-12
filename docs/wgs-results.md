# Whole-genome results: HG002 against T2T-Q100

> **Re-measured 2026-09-12** against `ae68ffd08` — 24 contigs in 61.4 min under the scheduler, then
> scored per contig. **The short-read whole-genome result is unchanged**: ALL F1 0.9726 and SNV
> 0.9846 to four decimals, Indel 0.9272 -> 0.9273, SV 0.5620 -> 0.5625. That is the intended
> reading, not a null result — a long-read preset, read-backed phasing and phase-driven
> re-genotyping all landed since the previous measurement, and every one of them is **off by default**
> on short reads. This run is the check that they are.
>
> For what those features do when they *are* on, see the long-read sections of
> [tier2-chr20-results.md](tier2-chr20-results.md) and [tier2-chr6-results.md](tier2-chr6-results.md).
> There are **no whole-genome long-read numbers**: see the note at the end of this page.

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
small-variant class on both recall and precision (ALL F1 0.9729 against 0.9505) and PanGenie is
ahead on structural variants (0.5739 against 0.5643). What is inside that SV gap, and whether
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

| | block emission (current) | decide-then-render | previous (inline) |
|---|---|---|---|
| **ALL** | TP 4,127,291  FP 92,779  FN 139,528  recall 0.9673  precision 0.9780  **F1 0.9726** | F1 0.9726 | F1 0.9699 |
| **SNV** | TP 3,305,130  FP 22,351  FN 81,057  recall 0.9761  precision 0.9933  **F1 0.9846** | F1 0.9846 | F1 0.9833 |
| **Indel** | TP 822,161  FP 70,428  FN 58,471  recall 0.9336  precision 0.9211  **F1 0.9273** | F1 0.9272 | F1 0.9191 |
| Insertion | recall 0.9228  precision 0.9122  **F1 0.9174** | F1 0.9179 | F1 0.9102 |
| Deletion | recall 0.9443  precision 0.9370  **F1 0.9406** | F1 0.9409 | F1 0.9333 |

**Small variants are unmoved by block emission**: ALL +0.0001, SNV flat, Indel +0.0003, with
insertions and deletions marginally down. The gain over inline is the earlier work's, not this one's.
Against decide-then-render, the current arm recovers 1,064 more true small variants and adds 592
false ones.

The current arm carries two changes over decide-then-render: block emission became the default, and
`resolve_site` stopped rejecting reversed snarls. The second is measured separately at 10 false
positives removed on chr20 and essentially nothing on chr6, so it is a small part of the FP movement
and none of the SV movement.

## Structural variants (truvari, >=50 bp)

| | block emission (current) | decide-then-render | previous (inline) |
|---|---|---|---|
| SV >= 50 bp | TP 14,452  FP 12,814  FN 9,665  **F1 0.5625** | F1 0.5620 | F1 0.5470 |

**This is where block emission pays, and it is the only place it does.** F1 0.5577 -> 0.5620,
**+0.0043**, from 48 more true SVs and 266 fewer false ones -- so unlike the previous step's
recall-only gain, this one improves both sides. Autosomes alone give 0.5596 -> 0.5642. The PanGenie
gap quoted above narrows from 0.0143 to 0.0097.

The per-contig spread is wide and worth knowing before quoting a single figure. Measured arm against
arm on one binary, chr20 gives **+0.0099** and chr6 **+0.0017** -- a 6x range, and with opposite
mechanisms: chr20 gained 11 true SVs at unchanged FP, chr6 removed 14 false ones and lost 2 true.
The genome-wide +0.0043 is the aggregate over 22 autosomes and is the number to use; chr20's figure
is the favourable tail, not the typical case.

`truvari refine` puts the same two comparisons at +0.0233 and +0.0121, so the record-matching metric
understates this change -- expected, since it penalises splitting one record into several and drops
any resulting block under the 50 bp size floor. Reported unrefined regardless, because every other SV
figure in this repository is unrefined and a refined number would compare to none of them.

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


## Whole-genome long reads: what it would take

Not measured, and the reason is disk rather than time. ONT alignments exist genome-wide —
`data/alignments-combined.processed.gaf.gz`, 68.5 GB, 4,017,467 reads against the 16-haplotype
`E821-16-sampled.gbz` — but `gaf-base sort` is an external merge sort and the whole file is roughly
**480 GB sorted, against 169 GB free**. That is why chr20 and chr6 were each filtered to their own
component first and built separately, and it is the only way this machine can do it.

So a whole-genome long-read run needs a **per-contig build**: one pass over the 68.5 GB GAF
splitting it 24 ways on the first node of each alignment's path (exact here — a Minigraph-Cactus
graph numbers each component contiguously), then 24 `gaf-base sort | construct` pairs, then the
calls. Rough cost: a few hours for the split, and the final databases total roughly 47 GB by
extrapolation from chr20's 1.0 GB for 85,373 reads.

Until then the long-read evidence is two contigs, chr20 and chr6, and both are in the tier-2 pages.
