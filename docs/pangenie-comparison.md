# vg call against PanGenie, same graph and same reads

An apples-to-apples comparison of two ways of genotyping the same panel from the same data:
`vg call --read-likelihood`, which scores **read alignments** against graph traversals, and
PanGenie v4.2.1, which scores **k-mer counts**.

Both were run on the HPRC v2.1 MC CHM13 graph with HG002 held out of the panel, from the same
30x NovaSeq PCR-free reads, and scored against T2T-Q100.

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
| ALL | 3,956,007 | 82,837 | 220,981 | **0.9630** | 3,960,421 | 195,585 | 216,567 | 0.9505 |
| SNV | 3,173,587 | 15,723 | 142,707 | **0.9756** | 3,203,093 | 69,979 | 113,201 | 0.9722 |
| Indel | 782,420 | 67,114 | 78,274 | **0.9150** | 757,328 | 125,606 | 103,366 | 0.8687 |
| SV ≥50 bp | 12,476 | 12,335 | 11,145 | 0.5152 | 13,749 | 10,544 | 9,872 | **0.5739** |

Recall and precision behind those:

| | vg recall | vg precision | PanGenie recall | PanGenie precision |
|---|---|---|---|---|
| ALL | 0.9471 | **0.9795** | 0.9482 | 0.9529 |
| SNV | 0.9570 | **0.9951** | **0.9659** | 0.9786 |
| Indel | **0.9091** | **0.9210** | 0.8799 | 0.8577 |

**The result is a clean split by variant class, not a winner.**

- **SNVs are close, and the two err differently.** PanGenie finds more of them — recall 0.9659
  against 0.9570 — and vg rejects far more false ones, precision 0.9951 against 0.9786. Net F1
  favours vg by 0.0034, which is small enough that the interesting statement is the trade, not the
  ranking: k-mer evidence is more sensitive, alignment evidence more specific.
- **Indels favour vg on both axes**, by 0.046 F1. PanGenie emits 125,606 indel false positives
  against 67,114. This is the largest small-variant difference and it is not a trade-off.
- **Structural variants favour PanGenie on both axes**, by 0.059 F1 — more true calls (13,749
  against 12,476) *and* fewer false ones (10,544 against 12,335). It is the one class where one
  tool is better in both directions, and it should be taken at face value rather than explained
  away.

## chrX, reported apart

| | vg call | PanGenie |
|---|---|---|
| ALL | **0.9422** | 0.8467 |
| SNV | **0.9542** | 0.8766 |
| Indel | **0.9004** | 0.7449 |
| SV ≥50 bp | 0.4260 | **0.4768** |

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
| chr1 | 16,217 | 17,475 | 9,330 (58%) | 6,887 | 8,145 |
| chr20 | 4,782 | 5,263 | 2,680 (56%) | 2,102 | 2,583 |

So roughly **57% of missed variants are missed by both** — a shared, panel-limited recall floor —
and the remaining 43% is tool-specific in both directions.

This bounds the shared limitation; it does not fully separate *not offered* from *offered and not
called*. A variant missed by vg but found by PanGenie proves the panel carried it, but the converse
inference is not available from these files alone.

## What the split means

The two tools see the same panel and the same reads and disagree in a patterned way, which makes
the pattern more informative than the ranking.

**Where alignment evidence wins: indels.** PanGenie emits 125,606 indel false positives against
vg's 67,114, and finds fewer true ones. An indel changes k-mer content over a short window, and
distinguishing a real short indel from a homopolymer miscount is exactly the case where counting
k-mers is weakest and where aligning a read across the site and asking how well it fits is
strongest. This is the largest small-variant difference and it is not a trade-off — vg is ahead on
both precision and recall.

**Where they trade: SNVs.** PanGenie finds 29,506 more true SNVs; vg emits 54,256 fewer false
ones (a 4.5x lower FP count). A SNV is one k-mer-length window's worth of signal, which k-mer
counting handles well and is why its recall is higher; the alignment model's advantage is that a
read spanning the site can be judged against every candidate allele, which is why its precision is.
Neither is obviously the better error to make, so the near-tie in F1 is the honest summary.

**Where k-mer evidence wins: structural variants.** PanGenie has 1,273 more true SVs *and* 1,791
fewer false ones. This is the one class where a tool is better in both directions. The plausible
reason is that a large event's k-mer signature is distinctive and does not require a read to align
across the breakpoint at all, whereas the read model needs traversals that reads can be scored
against — and `vg call`'s SV performance has a known weakness here, documented separately in
[tier2-sv-errors.md](tier2-sv-errors.md). Treat this as the strongest single result in the
comparison, and as the clearest direction for work on the read model.

**A caution on reading the FP counts.** vg's lower false-positive counts are partly a property of
what each tool emits: PanGenie genotypes every panel site and reports what it decides, while
`vg call` emits a record only where it calls non-reference. Both were reduced to non-reference
records before scoring, so the comparison is fair, but the two are not making the same *number* of
decisions and a per-decision error rate would differ from a per-record one.

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
