# Whole-genome results: HG002 against T2T-Q100

> **Measured 2026-09-23** on vg `0cab3fbd4`: 24 contigs in 86.5 min, packed about two at a time.
> Short reads, chr1-22+X: **ALL 0.9726, SNV 0.9844, Indel 0.9313, SV 0.5627**. Autosomes: ALL
> 0.9729, SNV 0.9847, Indel 0.9315, SV 0.5643. `0cab3fbd4` changes one default, `--mismap-max`
> 0.7 -> 0.95; the run passed `--mismap-max 0.95` to the build before it, and the two write a
> byte-identical chr20 VCF. That step buys SV F1 +0.0021 on chr1-22+X (+0.0020 on the autosomes)
> and leaves small variants unchanged within noise -- see *The mismap ceiling* below.
>
> Earlier short-read measurements: 2026-09-13 (`28f5b88e2`, the anchor-walk fix), within 0.0002 of
> the 2026-09-12 build in every class; 2026-09-12 (`ae68ffd08`), run to check that the long-read
> preset, read phasing and phase-driven re-genotyping stay off by default on short reads.

Called per contig on the 32-haplotype hap32 graph (34 panel haplotypes with the CHM13 and GRCh38
paths), `--read-likelihood` with panel enumeration, phasing and mosaic on, at vg's defaults. Scored
against T2T-Q100 v1.1 (GIAB defrabb V0.019 draft benchmark). chrY haploid; chrX haploid outside the
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
**The same 24 contigs on ONT**: *Whole-genome long reads*, at the end of this page.

**Compared against PanGenie on the same graph and reads**: see
[pangenie-comparison.md](pangenie-comparison.md). Briefly, on the autosomes vg is ahead on every
small-variant class on both recall and precision (ALL F1 0.9729 against 0.9505) and PanGenie is
ahead on structural variants (0.5722 against 0.5643, a gap of 0.0079). What is inside that SV gap,
and whether nested calling reached it: [sv-residual-errors.md](sv-residual-errors.md).

**The mosaic** is written per contig (`chr*.mosaic.tsv`). No genome-wide file has been assembled
since the format moved to mosaic-version 5, which `concat_mosaic.sh` rejects; wgs-performance.md
explains why assembling it is not `cat`.

**Nested calling and phasing are the defaults**, and **decide-then-render** is how records are now
built: a site's genotype is settled by the linkage barrier before its record exists, so nothing is
patched after the fact. `--no-nested` and `--no-phased` restore the older behaviour -- note
`--no-phased` also disables nested calling, since a nested site's ploidy comes from its parent's
phased genotype, so it is not a control for phasing alone. What each change bought when it landed is
under *How we got here*. See [nested-calling-design.md](nested-calling-design.md).

**One caveat that belongs with these numbers.** Nested calling's gain is a rich-panel effect: on the
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
contigs of the 2026-09-23 run, and they are asserted rather than reported.

## Small variants (aardvark, GT)

Records are built by **decide-then-render**, every site's genotype settled before its record is
built, with block emission on by default. Totals are over chr1-22+X, with the autosomes alongside
because the PanGenie and long-read comparisons are autosomal.

**How these are scored.** Recall is over truth records and precision over calls, each with its own
side's true-positive count, which is how aardvark computes its own F1. The two counts differ because
one truth record can be matched by several calls and the reverse, so both are shown. A genome figure
sums the four counts over contigs and then takes the rates; every per-contig rate is checked against
aardvark's and truvari's own (`scripts/bench_metrics.py`, `tests/test_bench_metrics.py`).

| | chr1-22+X | autosomes |
|---|---|---|
| **ALL** | truth TP 4,125,379  FN 141,440  query TP 4,139,226  FP 91,106  recall 0.9669  precision 0.9785  **F1 0.9726** | F1 0.9729 |
| **SNV** | truth TP 3,303,445  FN 82,742  query TP 3,218,808  FP 21,066  recall 0.9756  precision 0.9935  **F1 0.9844** | F1 0.9847 |
| **Indel** | truth TP 821,934  FN 58,698  query TP 920,418  FP 70,040  recall 0.9333  precision 0.9293  **F1 0.9313** | F1 0.9315 |
| Insertion | recall 0.9223  precision 0.9183  **F1 0.9203** | F1 0.9204 |
| Deletion | recall 0.9442  precision 0.9411  **F1 0.9427** | F1 0.9427 |

## Structural variants (truvari, >=50 bp)

| | chr1-22+X | autosomes |
|---|---|---|
| SV >= 50 bp | TP-base 14,406  FN 9,711  TP-comp 14,258  FP 12,553  recall 0.5973  precision 0.5318  **F1 0.5627** | TP-base 14,163  FN 9,458  TP-comp 14,015  FP 12,287  recall 0.5996  precision 0.5328  **F1 0.5643** |

Against PanGenie the SV deficit is precision, not recall: on the autosomes vg matches more true SVs
(TP-base 14,163 against 13,749) and makes more false ones (FP 12,287 against 10,544). Reported
unrefined, because every other SV figure in this repository is unrefined and a refined number would
compare to none of them.

## The mismap ceiling: 0.7 to 0.95

`--mismap-max` caps the mismapping probability a read's MAPQ implies. At 0.95 it binds only at
MAPQ 0, 4.97% of chr20's short-read alignments (a 5x sample); at 0.7 it bound MAPQ 0 and 1, 7.11%.
The two arms below are the same binary, reads and scoring, differing only in the ceiling.

| 0.7 -> 0.95 | chr1-22+X: F1 change | net FN+FP | autosomes: F1 change | net FN+FP |
|---|---|---|---|---|
| ALL | +0.000015 | −205 | −0.000042 | +282 |
| SNV | +0.000006 | −43 | −0.000055 | +349 |
| Indel | +0.000057 | −162 | +0.000012 | −67 |
| SV >= 50 bp | **+0.002073** | −239 | **+0.001951** | −222 |

**Structural variants gain, in both scopes.** SV F1 +0.0021 on chr1-22+X (95% CI +0.0008 to
+0.0033, paired bootstrap over 1 Mb blocks) and +0.0020 on the autosomes; held out of chr20, the
contig the ceiling was fitted on, it is +0.0020 again. The gain is precision: 274 fewer false SVs for
35 fewer true ones on chr1-22+X. It is broad rather than one contig's: SV F1 rises on 17 of 23
contigs and falls on chr1, chr8, chr9, chr13, chr15 and chr18.

**Small variants are unchanged within noise.** ALL, SNV and Indel F1 each move by less than 0.0001,
inside their bootstrap intervals (ALL −0.00023 to +0.00024 on chr1-22+X); deletions move by less
than 0.0001 too, and insertions by +0.00013 on chr1-22+X. The sign depends on scope:
on chr1-22+X net errors fall by 205, on the autosomes they rise by 282 (SNV +349). chrX is the
difference -- its ALL F1 rises 0.0026 with the ceiling, the largest move of any contig.

**Why 0.95 and not higher.** 0.99 was run on the same binary. Against 0.95 it costs small variants
significantly (ALL −0.00009, 95% CI −0.00019 to −0.00001, chr1-22+X) and does not move SVs
(−0.00016, n.s.).

## Per contig

chr1-22+X. Small F1 is aardvark ALL and SV F1 truvari >=50 bp, both four-count; the 0.7 columns are
the same binary at the old ceiling.

| contig | small F1 | SV F1 | small, 0.7 | SV, 0.7 | notes |
|---|---|---|---|---|---|
| chr1 | 0.9722 | 0.5869 | 0.9724 | 0.5871 |  |
| chr2 | 0.9716 | 0.5671 | 0.9706 | 0.5647 |  |
| chr3 | 0.9765 | 0.6004 | 0.9765 | 0.5985 |  |
| chr4 | 0.9766 | 0.5923 | 0.9765 | 0.5922 |  |
| chr5 | 0.9758 | 0.5557 | 0.9758 | 0.5556 |  |
| chr6 | 0.9774 | 0.5860 | 0.9775 | 0.5834 | same records as the tier-2 short-read arm |
| chr7 | 0.9733 | 0.5186 | 0.9733 | 0.5170 |  |
| chr8 | 0.9761 | 0.5837 | 0.9762 | 0.5839 |  |
| chr9 | 0.9750 | 0.5753 | 0.9750 | 0.5778 |  |
| chr10 | 0.9647 | 0.5326 | 0.9644 | 0.5287 |  |
| chr11 | 0.9717 | 0.5758 | 0.9734 | 0.5718 | largest small-variant drop with the ceiling |
| chr12 | 0.9745 | 0.5881 | 0.9745 | 0.5823 |  |
| chr13 | 0.9779 | 0.5686 | 0.9780 | 0.5689 |  |
| chr14 | 0.9747 | 0.5941 | 0.9748 | 0.5932 |  |
| chr15 | 0.9636 | 0.5638 | 0.9641 | 0.5666 |  |
| chr16 | 0.9640 | 0.5139 | 0.9636 | 0.5083 |  |
| chr17 | 0.9690 | 0.5390 | 0.9694 | 0.5339 |  |
| chr18 | 0.9758 | 0.5420 | 0.9758 | 0.5430 |  |
| chr19 | 0.9560 | 0.5501 | 0.9563 | 0.5469 |  |
| chr20 | 0.9725 | 0.5365 | 0.9724 | 0.5321 | same records as the tier-2 short-read arm (112,207); the ceiling was fitted here |
| chr21 | 0.9761 | 0.5638 | 0.9761 | 0.5612 |  |
| chr22 | 0.9704 | 0.5375 | 0.9711 | 0.5317 |  |
| chrX | 0.9609 | 0.4836 | 0.9583 | 0.4766 | haploid outside PAR; largest small-variant gain with the ceiling |

## How we got here

Dated records of the steps before the ceiling, each measured when it landed. Their deltas use the
older single-TP F1 (one true-positive count for both rates), so they do not add to the four-count
figures above; recomputed that way, a difference between any two short-read arms still on disk moves
by at most 0.0004 for ALL, SNV and Indel and 0.0008 for insertions, deletions and SVs. The
block-emission and decide-then-render arms are no longer on disk, and nothing here is recomputed
against the current run. The inline arm that decide-then-render replaced (2026-08-20) scores,
four-count on chr1-22+X: ALL 0.9699, SNV 0.9832, Indel 0.9231, insertions 0.9127, deletions 0.9352,
SV 0.5448.

**Nested calling and decide-then-render (to 2026-08-23).** Cumulatively over the two changes,
against the `--no-nested` arm: SNV F1 +0.0094, ALL F1 +0.0099, SV F1 +0.044.

**Block emission (2026-08-24).** Small variants were unmoved: ALL +0.0001, SNV flat, Indel +0.0003,
with insertions and deletions marginally down. The gain over inline is the earlier work's, not this
one's. Against decide-then-render, block emission recovered 1,064 more true small variants and added
592 false ones.

The block-emission arm carried two changes over decide-then-render: block emission became the
default, and `resolve_site` stopped rejecting reversed snarls. The second is measured separately at
10 false positives removed on chr20 and essentially nothing on chr6, so it is a small part of the FP
movement and none of the SV movement.

SVs were where block emission paid, and the only place it did. F1 **+0.0043**, from 48 more true SVs
and 266 fewer false ones -- so unlike the previous step's recall-only gain, it improved both sides.
Autosomes alone gave +0.0046, and the PanGenie gap narrowed by the same amount.

The per-contig spread was wide. Measured arm against arm on one binary, chr20 gave **+0.0099** and
chr6 **+0.0017** -- a 6x range, and with opposite mechanisms: chr20 gained 11 true SVs at unchanged
FP, chr6 removed 14 false ones and lost 2 true. The genome-wide +0.0043 is the number to use; chr20's
figure is the favourable tail, not the typical case.

`truvari refine` puts the same two comparisons at +0.0233 and +0.0121, so the record-matching metric
understates this change -- expected, since it penalises splitting one record into several and drops
any resulting block under the 50 bp size floor.

## Whole-genome long reads

> **Measured 2026-09-24** on the same build as the short-read run above (vg `0cab3fbd4`), with
> `--preset ont`: 24 contigs, `work/wgs-ont`.
>
> **Superseded for indels by `--hp-prior`** (a stronger panel prior at homopolymer-run indels, which
> `--preset ont` now sets; [ont-hp-prior.md](ont-hp-prior.md)). Genome-wide on chr1-22+X it takes ONT
> indel F1 from 0.8684 to 0.8918 and ALL from 0.9591 to 0.9648, every autosome up, with SNVs and SVs
> unchanged and no added cost. The indel gap to short reads narrows from 0.063 to 0.040.

24 contigs of ONT (about 44x, from chr20's 43.1x and chr6's 45.4x) against the same T2T-Q100 truth and the same confident regions as the
short-read run, on the **16-haplotype E821 graph** (18 panel haplotypes, `E821-16-sampled`) where
the short-read arm uses the 32-haplotype hap32 graph (34 panel). chrY called, excluded from the
totals for the same coordinate reason.

| chr1-22+X | short reads | ONT | ONT − short, 95% CI |
|---|---|---|---|
| ALL F1 | **0.9726** | 0.9591 | −0.0135 [−0.0142, −0.0127] |
| SNV F1 | 0.9844 | **0.9853** | **+0.0009** [+0.0002, +0.0016] |
| Indel F1 | **0.9313** | 0.8684 | **−0.0629** [−0.0641, −0.0617] |
| SV ≥50 bp F1 | 0.5627 | **0.5820** | **+0.0193** [+0.0136, +0.0252] |

Paired 1 Mb block bootstrap, 10,000 replicates (`work/wgs-run/wgs_compare.py`). Autosomes alone:
ONT ALL 0.9593, SNV 0.9854, Indel 0.8686, SV 0.5834, against 0.9729, 0.9847, 0.9315 and 0.5643. All
four differences hold with chr20 and chr6, the contigs the preset was fitted and checked on, both
left out: SNV +0.0009 [+0.0002, +0.0017], SV +0.0193 [+0.0132, +0.0255].

**ONT wins SNVs and structural variants, and loses indels by enough to lose overall.**

*Structural variants*, +0.0193, are mostly recall. truvari matches 15,327 truth SVs against the
short-read arm's 14,406, recall 0.6355 against 0.5973, at precision 0.5368 against 0.5318. ONT is
ahead on 20 of the 23 scored contigs and significantly so on five alone (chr7, chr8, chr10, chr13,
chr17); it trails on chr1, chr19 and chr21, none significantly. It is also ahead of **PanGenie**,
0.5820 against 0.5701 (+0.0119 [+0.0058, +0.0181]), in the one class where PanGenie leads short-read
vg; see [pangenie-comparison.md](pangenie-comparison.md).

*SNVs*, +0.0009, are a narrow but significant win: ONT at about 44x edges out Illumina at about 30x
(both measured on chr20 and chr6) on the class Illumina is supposed to own, on a graph with about half
the panel.

*Indels*, −0.0629, are ONT's loss: precision 0.8598 against 0.9293 and recall 0.8772 against 0.9333.
This is the homopolymer weakness the tier-2 pages measure on chr20 and chr6, where the preset trails
short reads by 0.066 and 0.059, holding at genome scale. Against the 2026-09-12 ONT build, indel F1
is up 0.024 on the autosomes (0.8445 to 0.8686), most of it precision (0.8261 to 0.8600), so the gap
to short reads has narrowed from 0.087 to 0.063. The other classes moved less: ALL 0.9531 to 0.9593,
SNV 0.9851 to 0.9854, SV 0.5826 to 0.5834.

**The caveat.** The two arms use different graphs, because alignments are graph-specific and the ONT
reads are aligned to E821. Panel size is what the linkage layer and the frequency prior feed on, so
the ONT figures are a floor on what ONT does with this caller, not a like-for-like.

### What it cost

| | short reads | ONT |
|---|---|---|
| CPU, 24 contigs | 8.58 h | **14.31 h** |
| peak RSS, worst contig | 10.0 GiB | 12.9 GiB (chrY); 9.0 GiB over the scored contigs (chr1) |
| wall clock, calls | 86.5 min (packed about two at a time, `-t 5`) | 2.18 h (two at a time, `-t 5`) |

**1.67x the CPU** for the same 24 contigs, down from 3.48x at the 2026-09-12 build: ONT's own CPU
fell from 30.59 h to 14.31 h while the short-read run's barely moved. chrY, which is not scored, is
the outlier: 6.7x the short-read CPU and the run's peak memory. Wall clock is not comparable between
the two runs, which is why the table leads with CPU; `scripts/wgs/runtimes.py` recomputes both from
the runs' own `/usr/bin/time -l` blocks.

### How it was built, and the constraint that shaped it

`gaf-base sort` is an external merge sort, and the whole 63.8 GiB genome-wide GAF is far too large
to sort in one pass on this laptop. It is built **per contig** (`work/run-ontg/ontg.py`):

1. one pass over the GAF (4,017,467 reads) splitting it on the first node of each alignment's path,
   compressed on the fly with `gzip -1`. The 22 contigs not already built take 64.6 GiB, a little more
   than the source because the split compresses less hard; chr20 and chr6 reuse their tier-2
   databases;
2. per contig, `gaf-base sort --preset long | gaf-base construct -r <the whole graph>`, deleting each
   split once its database lands. The sort spills to anonymous files in `$TMPDIR`, on the same disk,
   at about 2.4x the compressed split, and a database is about 0.7x its split. With ~21 GiB free after
   the split, the builds run **one at a time, smallest first**, so each landed database frees space
   before chr1 and chr2 (about 24 GiB each) come up. 22 builds took 50 minutes;
3. 24 calls, two at a time.

Node ranges are each graph chunk's minimum and maximum node ID, merged within a contig: chr13's small
component chunks fall inside its main chunk's span, which is harmless because both name chr13. An
overlap between two contigs would stop the build. 20,363 reads (0.5%) start on a node in no contig's
range and are not used.

**The partition check is the read count.** chr20 got 85,373 reads from the split, the same count as
its separately built tier-2 database, and each of the 22 rebuilt databases holds exactly the reads
its split was given. The genome run's chr20 VCF is also byte-identical to the tier-2 ONT arm
(`o20base`, 114,861 records), but chr20 reuses the tier-2 database, so that confirms the binary and
settings rather than the partition.
