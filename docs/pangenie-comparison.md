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
| ALL F1 | **0.9729** | 0.9505 |
| SNV F1 | **0.9849** | 0.9722 |
| SNV recall | **0.9762** | 0.9659 |
| Indel F1 | **0.9275** | 0.8687 |
| SV ≥50 bp F1 | 0.5642 | **0.5739** |

vg leads every small-variant class on both recall and precision; PanGenie leads structural variants
on both, now by 0.0097 — down from 0.0143, because block emission became the default and took vg's
autosomal SV F1 from 0.5596 to 0.5642.

**The vg column moved with decide-then-render** -- genotypes settled before records are built. It was
ALL 0.9703, SNV 0.9837, Indel 0.9195, SV 0.5488. PanGenie's column is unchanged: same run, same
scoring path, nothing about it re-measured. So the SV gap has nearly halved, from 0.0254 to 0.0143,
with no SV-specific work -- the change was about when a record is built, not about how SVs are
scored.

These are autosome-only, which is the scope this comparison is drawn at because chrX measures a
ploidy-handling difference rather than an evidence one (below). [wgs-results.md](wgs-results.md)
quotes the same run including chrX, so its figures run a few ten-thousandths lower — ALL F1 0.9725,
SNV 0.9846, SV 0.5577.

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
| ALL | 4,040,932 | 88,680 | 136,056 | **0.9729** | 3,960,421 | 195,585 | 216,567 | 0.9505 |
| SNV | 3,237,294 | 20,186 | 79,000 | **0.9849** | 3,203,093 | 69,979 | 113,201 | 0.9722 |
| Indel | 803,638 | 68,494 | 57,056 | **0.9275** | 757,328 | 125,606 | 103,366 | 0.8687 |
| SV ≥50 bp | 14,208 | 12,539 | 9,413 | 0.5642 | 13,749 | 10,544 | 9,872 | **0.5739** |

Recall and precision behind those:

| | vg recall | vg precision | PanGenie recall | PanGenie precision |
|---|---|---|---|---|
| ALL | **0.9674** | **0.9785** | 0.9482 | 0.9529 |
| SNV | **0.9762** | **0.9938** | 0.9659 | 0.9786 |
| Indel | **0.9337** | **0.9215** | 0.8799 | 0.8577 |

**The result is a clean split by variant class.** vg leads every small-variant class on *both*
axes; PanGenie leads structural variants on both.

- **SNVs**: vg finds 27,659 more true SNVs and emits 48,654 fewer false ones — a 3.3x lower
  false-positive count at higher recall. This is the one place the ranking has actually changed
  rather than merely widened: recall used to be PanGenie's, and the alleles that took it back are
  the ones nested calling stopped burying inside longer records.
- **Indels**: vg leads by 0.059 F1, the largest small-variant margin. PanGenie emits 125,606 indel
  false positives against 68,816.
- **Structural variants**: PanGenie still leads, by 0.0143 F1, but no longer in both directions: vg
  now makes MORE true calls (14,151 against 13,749) and still more false ones (12,805 against
  10,544). Under the previous arm PanGenie was ahead on both axes at once, and that is what changed
  -- the remaining gap is precision, not sensitivity.

## chrX, reported apart

| | vg call | PanGenie |
|---|---|---|
| ALL | **0.9563** | 0.8467 |
| SNV | **0.9691** | 0.8766 |
| Indel | **0.9172** | 0.7449 |
| SV ≥50 bp | 0.4569 | **0.4768** |

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

**This table is from an earlier vg arm and is the one thing on this page not re-measured.** vg's
small-variant FN has since fallen to **10,791 on chr1 and 3,200 on chr20**, so the shared-floor
percentages below are lower bounds on today's. It is left whole rather than part-updated because the
decomposition has to sum: 8,552 + 3,388 is the old 11,940, and replacing only the total would leave
a table that does not add up. Refreshing it needs a small-variant FN intersection, which does not
exist as a script -- `sv_delta.py` does this for structural variants only.

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

**Where k-mer evidence wins: structural variants.** PanGenie leads by 0.0097 F1 -- on precision only,
now that vg makes more true calls. The numbers behind that, since one F1 hides which side it comes
from:

| autosomal SVs ≥50 bp | vg call | PanGenie |
|---|---|---|
| TP | **14,208** | 13,749 |
| FP | 12,539 | **10,544** |
| FN | **9,413** | 9,872 |
| F1 | 0.5642 | **0.5739** |
| distinct truth SVs missed | **9,385** | 9,841 |
| of those, missed by the other tool too | 8,184 | 8,184 |
| missed by this tool alone | **1,201** | 1,657 |

**The daggers are gone: the overlap analysis has been re-run** against this arm, so every row above
comes from the same callset. `sv_delta.py` over both callsets now reports 8,184 truth SVs missed by
both, 1,201 by vg alone and 1,657 by PanGenie alone.

The change in "missed by vg alone" is large -- the previous, un-remeasured figure was 1,881 -- but do
not read all of it as this revision's doing. That number came from an arm several revisions back, so
the 680 spans everything between, not just block emission. What *is* attributable to block emission
is the row-count movement measured directly, arm against arm on one binary: TP 14,151 → 14,208 and
FP 12,805 → 12,539.

The first four rows are truvari's own row counts, kept so the F1s match the published figures; the
last three are over *distinct* truth variants, which is lower because truvari emits a row per match
and a multi-matched variant repeats. Set arithmetic needs the distinct form.

Requiring a correct genotype widens the gap from 0.0254 to 0.0408, so part of it is genotyping
rather than detection. The recall floor is mostly shared: 8,172 truth SVs are missed by both, 69.7%
of everything either tool misses, and 89.4% of those sit in a tandem repeat.

Of vg's 12,116 false positives, 9,672 (79.8%) were compared against a real nearby truth SV and
rejected on sequence or size similarity — they are near-misses, not inventions — and the 2,444 with
no truth SV within reach account for 65% of the whole 1,572 false-positive excess. 80.1% are under
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
