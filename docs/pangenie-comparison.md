# vg call against PanGenie, same graph and same reads

An apples-to-apples comparison of two ways of genotyping the same panel from the same data:
`vg call --read-likelihood`, which scores **read alignments** against graph traversals, and
PanGenie v4.2.1, which scores **k-mer counts**.

Both were run on the 32-haplotype hap32 graph (34 panel haplotypes with the CHM13 and GRCh38
paths), sampled from HPRC v2.1 MC CHM13 with HG002 held out of the panel, from the same 30x NovaSeq
PCR-free reads, and scored against T2T-Q100 v1.1 (GIAB defrabb V0.019 draft benchmark). The vg
numbers are the shipped default, vg `2a6a228a5` (vgteam/vg#4990): symbolic-allele nested calling
and panel phasing both on
([nested-calling-design.md](nested-calling-design.md)), `--mismap-max 0.95`, and no
`--max-snarl-edges` cap under `--read-likelihood`; long reads use `--preset ont`, which includes
`--hp-prior 20` ([ont-hp-prior.md](ont-hp-prior.md)). The whole-genome vg runs are
`work/wgs-mm095` (short reads, 2026-09-23, on `0cab3fbd4`; the later change is ONT-only and leaves
short-read output byte-identical) and `work/wgs-ont-hp` (ONT, 2026-09-24).

## Long reads

PanGenie is a short-read k-mer genotyper, so there is no long-read PanGenie arm to compare against.
What can be asked is the other question: **where does `vg call` on ONT sit relative to both?** That
is answerable on chr20 and chr6, the two tier-2 contigs, and all three rows below are scored against
the same T2T-Q100 truth, the same confident-region BED and the same truvari settings. The short-read
rows are the shipped default, record for record the chr20 and chr6 of the whole-genome run; the ONT
rows are `--preset ont`, record for record the chr20 and chr6 of the ONT whole-genome run. The ONT rows run on a smaller graph, the 16-haplotype
E821 graph (18 panel), against the 32-haplotype hap32 graph (34 panel haplotypes with the CHM13 and
GRCh38 paths) under the other two, because alignments are graph-specific; the caveat under the
whole-genome table applies here too.

| chr20 | ALL F1 | SNV F1 | Indel F1 | SV ≥50 bp F1 |
|---|---|---|---|---|
| PanGenie, 30x Illumina | 0.9492 | 0.9726 | 0.8672 | 0.5167 |
| `vg call`, 30x Illumina | **0.9725** | 0.9851 | **0.9288** | 0.5365 |
| `vg call`, ONT 43x, `--preset ont` | 0.9646 | **0.9860** | 0.8892 | **0.5577** |

| chr6 | ALL F1 | SNV F1 | Indel F1 | SV ≥50 bp F1 |
|---|---|---|---|---|
| PanGenie, 30x Illumina | 0.9572 | 0.9764 | 0.8872 | 0.6016 |
| `vg call`, 30x Illumina | **0.9774** | 0.9879 | **0.9396** | 0.5860 |
| `vg call`, ONT 45x, `--preset ont` | 0.9700 | **0.9882** | 0.9037 | **0.6046** |

**Four things worth reading off this.**

*vg on ONT has the best SNV F1 of the three*, on both contigs — 0.9860 and 0.9882, narrowly ahead of
vg's own short-read arm (0.9851, 0.9879) and well ahead of PanGenie. Long reads place SNVs at least
as well as either short-read method here.

*Indels are where short reads still win, and ONT sits above PanGenie.* 0.8892 against PanGenie's
0.8672 on chr20 and 0.9037 against 0.8872 on chr6, in both recall (0.8868 against 0.8735 on chr20)
and precision (0.8917 against 0.8610), while vg on short reads sits 0.040 and 0.036 above ONT. The
residual is still homopolymer: 1 bp indels in a reference homopolymer of 5 bp or more are 34% of
ONT's indel false positives on both contigs, against 17% and 16% on short reads, and that is where
both the reads and the benchmark are least reliable. Before `--hp-prior` those were 44% and 43%, and
ONT trailed PanGenie on indels here (0.8632 and 0.8807).

*ONT also has the best SV F1*, 0.5577 and 0.6046, and it is a recall lead: it matches the most truth
SVs on both contigs (recall 0.6288 and 0.6535, against 0.5765 and 0.6063 for vg on short reads and
0.5268 and 0.6057 for PanGenie), at precision level with short-read vg and below PanGenie.
Pooled over the two contigs, only the recall lead over short-read vg is significant, +0.0489
[+0.0213, +0.0747]; the F1 lead is +0.0194 [-0.0001, +0.0379] (paired 1 Mb block bootstrap over
the truvari records). Short-read vg is ahead of PanGenie on chr20 SVs and behind on chr6. These
rest on 765 and 1,547 truth SVs, so read the SV column as direction, not margin.

*Overall vg-on-ONT leads PanGenie* (0.9646 vs 0.9492; 0.9700 vs 0.9572) but stays behind vg on
short reads (0.9725, 0.9774). Nothing here says long reads beat short reads for this caller — on
these two contigs they win SNVs narrowly and SVs, and lose enough on indels to lose overall.

### Whole genome, all three

The ONT arm also runs genome-wide (`work/wgs-ont-hp`). Autosomes, same truth and regions:

| autosomes | PanGenie | vg, short reads | vg, ONT |
|---|---|---|---|
| ALL F1 | 0.9505 | **0.9729** | 0.9650 |
| SNV F1 | 0.9719 | 0.9847 | **0.9854** |
| Indel F1 | 0.8744 | **0.9315** | 0.8920 |
| SV ≥50 bp F1 | 0.5722 | 0.5643 | **0.5834** |

**ONT leads PanGenie in every class.** Over chr1-22+X, by paired block bootstrap: SV +0.0119
[+0.0058, +0.0181] (+0.0115 [+0.0052, +0.0179] with chr20 and chr6 held out), indel +0.0201
[+0.0185, +0.0217], SNV +0.0155 [+0.0140, +0.0171] and ALL +0.0166 [+0.0151, +0.0181]. The indel lead
is `--hp-prior`'s: before it ONT's autosomal indel F1 was 0.8686, 0.0058 behind PanGenie. Against
short-read vg, ONT wins SNVs (+0.0009) and SVs (+0.0193) and loses indels by 0.0396, so overall
short-read vg stays ahead of both; [wgs-results.md](wgs-results.md) has that comparison and what the
ONT run costs.

**The caveat that matters.** The ONT arm runs on the **16-haplotype E821 graph**
(`E821-16-sampled`, 18 panel haplotypes), the other two on the 32-haplotype hap32 graph (34 panel
haplotypes with the CHM13 and GRCh38 paths), because the ONT reads are aligned to the former and
alignments are graph-specific. Panel size is exactly what the linkage layer and the frequency prior
feed on, so the ONT column is handicapped by an amount this comparison cannot separate out. Read it
as a floor on what ONT does here, not as a like-for-like.

## The result in one table

Autosomes, summed counts, rates recomputed from them. Recall is over truth records and precision
over calls, each from its own side's true-positive count, as aardvark's and truvari's own F1s are.

| autosomes | vg call | PanGenie |
|---|---|---|
| ALL F1 | **0.9729** | 0.9505 |
| SNV F1 | **0.9847** | 0.9719 |
| SNV recall | **0.9757** | 0.9659 |
| Indel F1 | **0.9315** | 0.8744 |
| SV ≥50 bp F1 | 0.5643 | **0.5722** |

vg leads every small-variant class on both recall and precision; PanGenie leads structural variants
by 0.0079, on precision -- vg's SV recall is the higher, 0.5996 against 0.5821.

PanGenie's column is one run, scored once; only vg's moves. Before decide-then-render and block
emission (the `wgs-current` arm) vg stood at ALL 0.9703, SNV 0.9836, Indel 0.9235, SV 0.5466, a
0.0256 SV gap; those two changes are about how records are built, not how SVs are scored, and they
took the gap to 0.0094 (2026-09-12, vg SV 0.5628). The step to 0.0079 (2026-09-23, vg
`0cab3fbd4`) comes from the mismap ceiling (`--mismap-max`), moved from 0.7 to 0.95 for its SV
gain: on one binary it is worth +0.0020 autosomal SV F1 (0.5623 → 0.5643), and +0.0020 over
chr1-22+X with chr20 held out. The net step is smaller, +0.0015 (0.5628 → 0.5643), because the
other caller changes between the two builds, the snarl-size cap coming off under
`--read-likelihood` among them, cost 0.0005 at the old ceiling (0.5628 → 0.5623). The ceiling
leaves small variants unchanged within noise: on these autosomes it adds 282 net small-variant
errors (+349 on SNVs), where on chr1-22+X it removes 205.

These are autosome-only, which is the scope this comparison is drawn at because chrX measures a
ploidy-handling difference rather than an evidence one (below). [wgs-results.md](wgs-results.md)
quotes the same run including chrX, so its figures run lower — ALL F1 0.9726, SNV 0.9844, Indel
0.9313, SV 0.5627, against PanGenie's SV 0.5701 there.

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

Truth TP and FN give recall, query TP and FP give precision; for SVs these are truvari's TP-base,
FN, TP-comp and FP. The two TP counts differ because one truth record can be matched by several
calls and the reverse. FN counts missed truth variants and FP false calls. vg on ONT is added for
reference: it runs on a different graph (16 sampled haplotypes) from different reads, so its rows
are not like-for-like with the other two (see *Long reads* above). Bold marks the best of the three.

| | | truth TP | FN | query TP | FP | **F1** |
|---|---|---|---|---|---|---|
| **ALL** | vg call, short reads | 4,038,968 | 138,020 | 4,052,250 | 87,472 | **0.9729** |
| | PanGenie | 3,960,421 | 216,567 | 3,959,329 | 195,585 | 0.9505 |
| | vg call, ONT | 4,001,019 | 175,969 | 4,033,215 | 115,532 | 0.9650 |
| **SNV** | vg call, short reads | 3,235,583 | 80,711 | 3,152,411 | 19,260 | 0.9847 |
| | PanGenie | 3,203,093 | 113,201 | 3,125,714 | 69,979 | 0.9719 |
| | vg call, ONT | 3,236,879 | 79,415 | 3,173,834 | 15,936 | **0.9854** |
| **Indel** | vg call, short reads | 803,385 | 57,309 | 899,839 | 68,212 | **0.9315** |
| | PanGenie | 757,328 | 103,366 | 833,615 | 125,606 | 0.8744 |
| | vg call, ONT | 764,140 | 96,554 | 859,381 | 99,596 | 0.8920 |
| **SV ≥50 bp** | vg call, short reads | 14,163 | 9,458 | 14,015 | 12,287 | 0.5643 |
| | PanGenie | 13,749 | 9,872 | 13,562 | 10,544 | 0.5722 |
| | vg call, ONT | 15,050 | 8,571 | 14,816 | 12,723 | **0.5834** |

Recall and precision behind those:

| | short recall | short precision | PanGenie recall | PanGenie precision | ONT recall | ONT precision |
|---|---|---|---|---|---|---|
| ALL | **0.9670** | **0.9789** | 0.9482 | 0.9529 | 0.9579 | 0.9722 |
| SNV | 0.9757 | 0.9939 | 0.9659 | 0.9781 | **0.9761** | **0.9950** |
| Indel | **0.9334** | **0.9295** | 0.8799 | 0.8691 | 0.8878 | 0.8961 |
| SV ≥50 bp | 0.5996 | 0.5328 | 0.5821 | **0.5626** | **0.6371** | 0.5380 |

**The result is a clean split by variant class.** Short-read vg leads every small-variant class on
*both* axes; PanGenie leads structural variants on precision alone.

- **SNVs**: vg finds 32,490 more true SNVs and emits 50,719 fewer false ones — a 3.6x lower
  false-positive count at higher recall. This is the one place the ranking has actually changed
  rather than merely widened: recall used to be PanGenie's, and the alleles that took it back are
  the ones nested calling stopped burying inside longer records.
- **Indels**: vg leads by 0.057 F1, the largest small-variant margin. PanGenie emits 125,606 indel
  false positives against 68,212.
- **Structural variants**: PanGenie leads by 0.0079 F1, on precision alone: vg matches more truth
  SVs (14,163 against 13,749) and makes more false calls (12,287 against 10,544). The gap is
  precision, not sensitivity.
- **vg on ONT**, for reference, makes the fewest SNV errors of the three -- 79,415 missed and 15,936
  false, against 80,711 and 19,260 for short-read vg and 113,201 and 69,979 for PanGenie -- and
  misses the fewest SVs (8,571, against 9,458 and 9,872), at the cost of more false SV calls than
  either (12,723). On indels it sits between the two: 96,554 missed and 99,596 false, against
  short-read vg's 57,309 and 68,212 and PanGenie's 103,366 and 125,606.

## chrX, reported apart

| | vg call | PanGenie |
|---|---|---|
| ALL | **0.9609** | 0.8509 |
| SNV | **0.9722** | 0.8780 |
| Indel | **0.9243** | 0.7640 |
| SV ≥50 bp | 0.4836 | 0.4730 |

**This is a ploidy-handling difference, not an evidence one, and folding it into a genome-wide F1
would misreport it.** HG002 is male, so chrX outside the pseudoautosomal regions carries one copy;
the truth is haploid there (119,112 bare `1` genotypes). PanGenie as run here calls chrX diploid
throughout — 94,691 `1/1` and 29,286 heterozygous calls, the latter wrong by construction. `vg
call` expresses the real ploidy with `--ploidy-bed`.

That is a genuine capability difference worth recording, and it is *not* evidence that one model
reads the data better. Note that on chrX SVs the two are within noise of each other despite it:
vg's 0.4836 against 0.4730 is a 0.0106 lead over 496 truth SVs.

## How much of the recall ceiling is shared

Both tools draw alleles from the same panel, so neither can call what the panel does not carry.
Comparing which truth variants each misses, joined record by record on (CHROM, POS, REF, ALT) across
the two tools' aardvark truth labels. "Missed" is aardvark's FN, which includes a truth variant
called with the wrong genotype.

| scope | vg FN | PanGenie FN | missed by both | vg only | PanGenie only |
|---|---|---|---|---|---|
| autosomes | 138,020 | 216,567 | 105,977 | 32,043 | 110,590 |
| autosomes, SNV | 80,711 | 113,201 | 65,811 | 14,900 | 47,390 |
| autosomes, indel | 57,309 | 103,366 | 40,166 | 17,143 | 63,200 |
| chr1 | 10,976 | 17,475 | 8,526 | 2,450 | 8,949 |
| chr20 | 3,199 | 5,263 | 2,422 | 777 | 2,841 |

Read down the vg column rather than across: **76.8% of what vg still misses on the autosomes is also
missed by PanGenie** (105,977 of 138,020) — 81.5% of its SNV misses and 70.1% of its indel misses —
and 77.7% on chr1, 75.7% on chr20. vg's residual recall deficit is mostly the shared, panel-limited
floor rather than anything specific to the read model; before nested calling the shared share on
these two contigs was 58%.

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

**Where k-mer evidence wins: structural variants.** PanGenie leads by 0.0079 F1 -- on precision only,
since vg makes more true calls. The numbers behind that, since one F1 hides which side it comes
from:

| autosomal SVs ≥50 bp | vg call | PanGenie |
|---|---|---|
| TP-base | **14,163** | 13,749 |
| FN | **9,458** | 9,872 |
| TP-comp | **14,015** | 13,562 |
| FP | 12,287 | **10,544** |
| F1 | 0.5643 | **0.5722** |
| of the FN, missed by the other tool too | 8,222 | 8,222 |
| of the FN, missed by this tool alone | **1,236** | 1,650 |

The first five rows are truvari's own counts, so the F1s match the published figures. The last two
split each tool's FN over the same 23,621 truth SVs, joined record by record on (CHROM, POS, REF,
ALT); every truth SV is TP-base or FN exactly once, so the two rows sum to the FN row.
`sv_delta.py` keys truth on CHROM:POS:SVTYPE:SVLEN instead, which merges compound-heterozygous
alleles at one site and slightly undercounts these sets. Block emission's own effect, measured arm
against arm on one binary (2026-08-24), was TP-base 14,151 → 14,208 and FP 12,805 → 12,539.

Requiring a correct genotype widens the gap from 0.0079 to 0.0107 (vg 0.5145, PanGenie 0.5252), so
part of it is genotyping rather than detection. That is an autosome figure: on chrX PanGenie's
diploid genotypes fail 215 of its 242 matched haploid SVs, which measures ploidy, not evidence.

The recall floor is mostly shared: 8,222 truth SVs are missed by both, 74.0% of everything either
tool misses and 34.8% of all 23,621, so the two together could reach at most 0.652 SV recall; 89.6%
of the shared misses sit in a tandem repeat.

Most of vg's 12,287 false positives are near-misses, not inventions: 10,274 (83.6%) were compared
against a real nearby truth SV and rejected on sequence or size similarity, and 2,013 had no truth
SV within reach. Of the 1,743-call excess over PanGenie, 1,146 are near-misses and 597 (34%) have
no truth SV nearby. 85.2% of vg's SV false positives are under 300 bp. No threshold on GQ, GQN or DR
raises SV F1 -- the best in a scan of all three, DR ≥ 0.05, gives 0.5639 against 0.5643 ungated --
because false negatives outnumber that excess about 5.4 to one.

A quarter of the F1 gap is a representation artefact: same-length substitutions, which vg's
multiallelic output carries, PanGenie's biallelic-split output essentially does not, and truvari
sizes by allele length and so scores as structural. Excluding them from both sides narrows the gap
from 0.0079 to 0.0060. The rest is extra insertion and deletion false positives (12,086 against
10,515), not missed truth. Full anatomy, measured on an earlier arm, including what nested calling
did and did not reach: [sv-residual-errors.md](sv-residual-errors.md).

**A caution on reading the FP counts.** vg's lower small-variant false-positive counts are partly a
property of what each tool emits: PanGenie genotypes every panel site and reports what it decides,
while `vg call` emits a record only where it calls non-reference. Both were reduced to
non-reference records before scoring, so the comparison is fair, but the two are not making the
same *number* of decisions and a per-decision error rate would differ from a per-record one.

## Caveats

- One sample. The PanGenie comparison proper is one graph and one read technology; the ONT rows add
  a second of each, on the smaller 16-haplotype E821 graph (18 panel). Nothing here is replicated.
- PanGenie was run by someone else with its own defaults and was not tuned here. vg is at its
  shipped settings, but those were fitted on HG002 itself, and `--mismap-max 0.95` was chosen for
  the SV gain it shows on this whole-genome run (+0.0020 on the autosomes, and +0.0020 over
  chr1-22+X with chr20 held out), so the vg column is not out-of-sample the way PanGenie's is.
- The SV numbers rest on ~24k truth SVs genome-wide and are much noisier than the small-variant
  ones.
- chrY is excluded from both, identically, for the reference mismatch documented in
  [wgs-results.md](wgs-results.md).

Reproduce with `scripts/wgs/prep_external_vcf.sh`, then `bench_wgs.py --work work/pangenie`, then
`scripts/wgs/compare_callsets.py --a 'vg call=work/wgs-mm095/score/per-contig.json' --b
'PanGenie=work/pangenie/score/per-contig.json'`. The SV false-positive split is
`scripts/wgs/sv_quality_gates.py --score work/wgs-mm095/score`, and the genotype-required F1 and
substitution split come from `scripts/wgs/sv_delta.py`, whose defaults regenerate
[sv-delta.md](sv-delta.md) from `work/wgs-mm095/score`. The finer gate scan and the two record-level joins, small-variant FN
overlap and SV set partition, are
`work/run-docs/analysis/sv-anatomy/fp_supplement.py`,
`work/run-docs/analysis/fn-intersection/fn_intersection.py` and
`work/run-docs/analysis/sv-anatomy/recount_sv_sets.py`.
