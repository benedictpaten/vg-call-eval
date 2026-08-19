# vg call against PanGenie, same graph and same reads

An apples-to-apples comparison of two ways of genotyping the same panel from the same data:
`vg call --read-likelihood`, which scores **read alignments** against graph traversals, and
PanGenie v4.2.1, which scores **k-mer counts**.

Both were run on the HPRC v2.1 MC CHM13 graph with HG002 held out of the panel, from the same
30x NovaSeq PCR-free reads, and scored against T2T-Q100. The vg numbers are the shipped default:
symbolic-allele nested calling and panel phasing both on ([nested-calling-design.md](nested-calling-design.md)).

## The result in one table

Autosomes, summed counts, rates recomputed from them:

| autosomes | vg call | PanGenie |
|---|---|---|
| ALL F1 | **0.9703** | 0.9505 |
| SNV F1 | **0.9837** | 0.9722 |
| SNV recall | **0.9742** | 0.9659 |
| Indel F1 | **0.9195** | 0.8687 |
| SV ≥50 bp F1 | 0.5485 | **0.5739** |

vg leads every small-variant class on both recall and precision; PanGenie leads structural variants
on both, by 0.0254.

These are autosome-only, which is the scope this comparison is drawn at because chrX measures a
ploidy-handling difference rather than an evidence one (below). [wgs-results.md](wgs-results.md)
quotes the same run including chrX, so its figures run a few ten-thousandths lower — ALL F1 0.9699,
SNV 0.9833, SV 0.5467.

## What makes it like for like

Both call sets go through the *same code path*: `scripts/wgs/prep_external_vcf.sh` lays the
external VCF out in the shape `bench_wgs.py` already expects, with the truth VCFs, confident-region
BEDs and reference FASTAs **symlinked from the vg run**. So both see byte-identical truth inputs,
the same aardvark and truvari invocations, the same per-contig-then-summed aggregation, and the
same chrY exclusion. Writing a second scorer and hoping it matched was the alternative, and it is
worse.

Two adjustments to the input, both necessary and both stated rather than buried:

- **Hom-ref records dropped.** PanGenie emits a record for every panel site including `0/0`;
  `vg call` emits only non-reference calls. Keeping them would have handed the comparison 30M
  records that are not calls. 4,733,256 non-reference records remained.
- **Contig headers added**, since the file had none and could not otherwise be indexed.

Not adjusted: the allele representation. PanGenie's output is biallelic-split where ours is
multiallelic. aardvark compares by local haplotype rather than by record, so this is exactly the
difference it exists to absorb; normalising by hand would have been the riskier move.

## Autosomes — the like-for-like result

|  | vg call | | | | PanGenie | | | |
|---|---|---|---|---|---|---|---|---|
| | TP | FP | FN | **F1** | TP | FP | FN | **F1** |
| ALL | 4,026,067 | 95,173 | 150,921 | **0.9703** | 3,960,421 | 195,585 | 216,567 | 0.9505 |
| SNV | 3,230,752 | 21,325 | 85,542 | **0.9837** | 3,203,093 | 69,979 | 113,201 | 0.9722 |
| Indel | 795,315 | 73,848 | 65,379 | **0.9195** | 757,328 | 125,606 | 103,366 | 0.8687 |
| SV ≥50 bp | 13,540 | 12,206 | 10,081 | 0.5485 | 13,749 | 10,544 | 9,872 | **0.5739** |

Recall and precision behind those:

| | vg recall | vg precision | PanGenie recall | PanGenie precision |
|---|---|---|---|---|
| ALL | **0.9639** | **0.9769** | 0.9482 | 0.9529 |
| SNV | **0.9742** | **0.9934** | 0.9659 | 0.9786 |
| Indel | **0.9240** | **0.9150** | 0.8799 | 0.8577 |

**The result is a clean split by variant class.** vg leads every small-variant class on *both*
axes; PanGenie leads structural variants on both.

- **SNVs**: vg finds 27,659 more true SNVs and emits 48,654 fewer false ones — a 3.3x lower
  false-positive count at higher recall. This is the one place the ranking has actually changed
  rather than merely widened: recall used to be PanGenie's, and the alleles that took it back are
  the ones nested calling stopped burying inside longer records.
- **Indels**: vg leads by 0.051 F1, the largest small-variant margin. PanGenie emits 125,606 indel
  false positives against 73,848.
- **Structural variants**: PanGenie leads by 0.0254 F1, with slightly more true calls (13,749
  against 13,540) and fewer false ones (10,544 against 12,206). It is the one class where one tool
  is better in both directions, and it should be taken at face value rather than explained away.

## chrX, reported apart

| | vg call | PanGenie |
|---|---|---|
| ALL | **0.9494** | 0.8467 |
| SNV | **0.9631** | 0.8766 |
| Indel | **0.9022** | 0.7449 |
| SV ≥50 bp | 0.4617 | **0.4768** |

**This is a ploidy-handling difference, not an evidence one, and folding it into a genome-wide F1
would misreport it.** HG002 is male, so chrX outside the pseudoautosomal regions carries one copy;
the truth is haploid there (119,112 bare `1` genotypes). PanGenie as run here calls chrX diploid
throughout — 94,691 `1/1` and 29,286 heterozygous calls, the latter wrong by construction. `vg
call` expresses the real ploidy with `--ploidy-bed`.

That is a genuine capability difference worth recording, and it is *not* evidence that one model
reads the data better. Note that PanGenie still wins on chrX SVs despite it.

## How much of the recall ceiling is shared

Both tools draw alleles from the same panel, so neither can call what the panel does not carry.
Comparing which truth variants each misses:

| contig | vg FN | PanGenie FN | missed by both | vg only | PanGenie only |
|---|---|---|---|---|---|
| chr1 | 11,940 | 17,475 | 8,552 | 3,388 | 8,923 |
| chr20 | 3,535 | 5,263 | 2,432 | 1,103 | 2,831 |

Read down the vg column rather than across: **72% of what vg still misses on chr1 is also missed by
PanGenie** (8,552 of 11,940), and 69% on chr20. vg's residual recall deficit is now mostly the
shared, panel-limited floor rather than anything specific to the read model — which is a different
statement from the one these two contigs supported before nested calling, when the shared share was
58%.

This bounds the shared limitation; it does not fully separate *not offered* from *offered and not
called*. A variant missed by vg but found by PanGenie proves the panel carried it, but the converse
inference is not available from these files alone.

## What the split means

The two tools see the same panel and the same reads and disagree in a patterned way, which makes
the pattern more informative than the ranking.

**Where alignment evidence wins: small variants, on both axes.** An indel changes k-mer content
over a short window, and distinguishing a real short indel from a homopolymer miscount is exactly
the case where counting k-mers is weakest and where aligning a read across the site and asking how
well it fits is strongest. The SNV lead is newer and has a specific cause: a SNV inside a long
alternative allele is invisible to a caller that only emits the long allele, and descending into
those nested bubbles recovered 59,413 SNV false negatives without costing precision.

**Where k-mer evidence wins: structural variants.** PanGenie leads by 0.0254 F1 with both more true
calls and fewer false ones. The numbers behind that, since one F1 hides which side it comes from:

| autosomal SVs ≥50 bp | vg call | PanGenie |
|---|---|---|
| TP | 13,540 | 13,749 |
| FP | 12,206 | 10,544 |
| FN | 10,081 | 9,872 |
| F1 | 0.5485 | **0.5739** |
| F1 requiring the right genotype | 0.4805 | **0.5213** |
| distinct truth SVs missed | 10,053 | 9,841 |
| of those, missed by the other tool too | 8,172 | 8,172 |
| missed by this tool alone | 1,881 | 1,669 |

The first four rows are truvari's own row counts, kept so the F1s match the published figures; the
last three are over *distinct* truth variants, which is lower because truvari emits a row per match
and a multi-matched variant repeats. Set arithmetic needs the distinct form.

Requiring a correct genotype widens the gap from 0.0254 to 0.0408, so part of it is genotyping
rather than detection. The recall floor is mostly shared: 8,172 truth SVs are missed by both, 69.7%
of everything either tool misses, and 89.4% of those sit in a tandem repeat.

Of vg's 12,206 false positives, 9,748 (79.9%) were compared against a real nearby truth SV and
rejected on sequence or size similarity — they are near-misses, not inventions — and the 2,458 with
no truth SV within reach account for 63% of the whole 1,662 false-positive excess. 80.8% are under
300 bp. No threshold on GQ, GQN or DR raises SV F1, because false negatives already outnumber that
excess six to one. Full anatomy, including what nested calling did and did not reach:
[sv-residual-errors.md](sv-residual-errors.md).

**A caution on reading the FP counts.** vg's lower small-variant false-positive counts are partly a
property of what each tool emits: PanGenie genotypes every panel site and reports what it decides,
while `vg call` emits a record only where it calls non-reference. Both were reduced to
non-reference records before scoring, so the comparison is fair, but the two are not making the
same *number* of decisions and a per-decision error rate would differ from a per-record one.

## Caveats

- One sample, one graph, one read technology. Nothing here is replicated.
- PanGenie was run by someone else with its own defaults; no attempt was made to tune either tool
  for this comparison, and both are at their shipped settings.
- The SV numbers rest on ~24k truth SVs genome-wide and are much noisier than the small-variant
  ones.
- chrY is excluded from both, identically, for the reference mismatch documented in
  [wgs-results.md](wgs-results.md).

Reproduce with `scripts/wgs/prep_external_vcf.sh`, then `bench_wgs.py --work work/pangenie`, then
`scripts/wgs/compare_callsets.py`.
