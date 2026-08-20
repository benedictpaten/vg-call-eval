# The structural variants vg still gets wrong

PanGenie leads SV F1 on the autosomes 0.5739 to 0.5488 ([pangenie-comparison.md](pangenie-comparison.md)),
with both slightly more true calls (13,749 against 13,516) and fewer false ones (10,544 against
12,116). This is what is actually inside that gap, and whether the nested-calling work reached it.

Measured on the single-sweep nested arm with the code-review fixes, `work/wgs-fixed` (vg a27149728). Every count and rate on this page is
from the scripts named at the foot of it, re-run against that arm -- except the record-level tracing
of the calls with no truth SV in reach (the section on what those 2,444 actually are, and the
offsetting insertion/deletion pair analysis that follows it). That tracing was done by hand two arms
ago, on a population of 2,458 that has since moved by 14 records, so its rates carry over; it has
not been repeated. Where a figure comes from that tracing rather than from a re-run, it is the
mechanism that is being asserted and not the third digit.

Two questions are answered separately, because they have different answers:

1. Is the caller still writing one long insertion or deletion where the sample carries several
   small changes?
2. What are the remaining errors made of -- how long, what type, and what distinguishes them from
   the calls that were right?

Everything below is autosomes only, truth T2T-Q100, from the same run reported in
[wgs-results.md](wgs-results.md). Where a comparison against the pre-nested-calling behaviour is
useful it is against `--no-nested` on the same graph, same reads and same scoring.

## 1. The long-record pathology is largely gone

Three independent measurements agree, and none of them is the F1 that motivated the work.

**Same-length substitutions.** A record whose REF and ALT are the same length has no size change at
all, so a long one is by construction several small changes written together -- and truvari sizes
it by allele length and scores it as structural. This is the pathology in its purest observable
form.

| | same-length substitution false positives |
|---|---|
| `--no-nested` | 2,219 |
| **current default** | **294** |
| PanGenie | 29 |

An 87% reduction, and it is most of the way to PanGenie's near-zero rate, which comes for free from
biallelic-split output rather than from better calling.

**What vg writes where it misses an SV that PanGenie finds.** These are the 1,881 truth SVs whose
alleles the panel demonstrably carried, since PanGenie called them. Asking what vg put within 100 bp
of each:

| what vg emitted at the locus | `--no-nested` | current default |
|---|---|---|
| no record within 100 bp | 1,230 (46.7%) | **258 (13.7%)** |
| a record of comparable size, which truvari rejected | 890 (33.8%) | 869 (46.1%) |
| only records far too small, not accounting for the event | 415 (15.8%) | 601 (31.9%) |
| several small records that *do* sum to the event | 98 (3.7%) | 156 (8.3%) |
| total | 2,633 | 1,884 |

**The "no record at all" population fell by 79%.** That row was the one nested calling was aimed at:
a bubble the default scope never descended into leaves no record anywhere near the event. It is now
the smallest category rather than the largest.

**Per-variant, following the same truth SVs across both arms**, which is the only way to read the
table above without being misled by two differently-sized denominators:

- 1,174 truth SVs recovered
- 119 lost
- net **+1,055**, of which 81% are 50-300 bp

| size | recovered | newly lost | net |
|---|---|---|---|
| 50-100 | 537 | 33 | +504 |
| 100-300 | 413 | 62 | +351 |
| 300-700 | 156 | 16 | +140 |
| 700-2000 | 49 | 6 | +43 |
| 2000-10000 | 17 | 2 | +15 |
| 10000+ | 2 | 0 | +2 |

### Taking each false positive apart, allele by allele

The three measurements above are all indirect: they count records in categories rather than looking
at what any individual record contains. The direct version strips each false positive of its common
prefix and suffix with the reference, aligns what is left, and treats a run of 30 bp or more of exact
match as separating one change from the next. A record resolving to several changes is several
variants written as one; where *none* of those changes reaches 50 bp, the record is scored as
structural only because of the bundling. Full tables in [sv-fp-anatomy.md](sv-fp-anatomy.md).

| | `--no-nested` | current default | PanGenie |
|---|---|---|---|
| false positives decomposed | 12,135 | 12,129 | 10,541 |
| resolve to one change | 60.4% | **79.5%** | 96.0% |
| resolve to several changes | 39.6% | **20.5%** | 4.0% |
| several changes, **none reaching 50 bp** | 2,240 (18.5%) | **310 (2.6%)** | 29 (0.3%) |

**The bundling pathology is down 86%**, and what is left is 205 same-length substitutions, 61
insertions and 44 deletions out of 12,129. In the `--no-nested` arm the same population was 2,099
substitutions. (Per-record classes from `sv_fp_anatomy.py --json`, which is where the breakdowns
below come from too.)

**And it was never a long-insertion or long-deletion phenomenon.** Of the 2,253 INS/DEL false
positives of 300 bp or more in the current call set, **zero** decompose entirely into changes under
50 bp -- as do zero of the 2,035 in the `--no-nested` arm and zero of PanGenie's 1,267. A long
insertion vg calls wrongly is a genuine long insertion in the wrong place or of the wrong length; it
has never been a handful of SNPs wearing an insertion's clothes. That shape only ever appeared in
same-length substitutions, which is exactly where symbolic collapsing was aimed.

Long records *can* still carry more than one change -- 896 of those 2,253 do -- but the largest
single change is a median **86%** of the record's length, so they are one real structural event with
hitchhikers rather than a bundle of equals. PanGenie's bundled long records sit at 82% by the same
measure, so this is a property of writing multi-kilobase alleles at all, not of this caller.

The remaining 20.5%-versus-4.0% difference against PanGenie in multi-change records is mostly
representation: PanGenie's output is biallelic-split, so its alleles are simpler by construction.

### The inverse pathology exists, and is smaller

Decomposition can go too far: the last row of the table above is loci where vg wrote several records
that individually fall under the 50 bp structural threshold but together account for the truth
event. That population grew from 98 to 156.

It is worth naming and not worth alarm. It is 8.3% of the residual actionable false negatives
against the 79% reduction in the category nested calling was built to remove, and it is **not
straightforwardly caused by nested calling**: the clearest example, a 245 bp insertion at
chr1:2,718,630 written as four records of about 61 bp, appears identically in the `--no-nested` arm.
Writing a tandem-repeat expansion as its repeat units is a separate representation choice that
predates this work.

## 2. What the remaining errors are made of

### The false positives are small, and most are near-misses

**80.8% are under 300 bp.** Median 102 bp.

| type | size | n | share |
|---|---|---|---|
| DEL | 50-100 | 2,932 | 24.0% |
| INS | 50-100 | 2,908 | 23.8% |
| INS | 100-300 | 1,904 | 15.6% |
| DEL | 100-300 | 1,899 | 15.6% |
| INS | 300-700 | 599 | 4.9% |
| DEL | 300-700 | 469 | 3.8% |
| INS | 700-2000 | 394 | 3.2% |
| INS | 2000-10000 | 338 | 2.8% |
| DEL | 700-2000 | 217 | 1.8% |
| INS | 10000+ | 135 | 1.1% |

The first split that matters is whether truvari had any truth SV to compare the call against:

| | vg call | PanGenie |
|---|---|---|
| no truth SV near enough to compare with | 2,444 (20.2%) | 1,416 (13.4%) |
| compared against a nearby truth SV and rejected | 9,672 (79.8%) | 9,128 (86.6%) |

**65% of the entire false-positive excess is the first row** -- +1,028 of the +1,572. Taken apart
below, it is not one thing, and it is not mostly calls out of nowhere: at 71% of these loci the truth
carries variation, just not variation truvari scores as structural.

### What the 2,444 with nothing to compare against actually are

This is the population worth a mechanism rather than a rate, so each one was traced back to what the
truth carries at that locus and what vg wrote there.

| | n | share |
|---|---|---|
| truth has a 20-49 bp indel within 300 bp -- just under truvari's threshold | 792 | 33.6% |
| truth has only variants under 20 bp within 300 bp | 891 | 37.8% |
| truth has no variant of any size within 300 bp | 676 | 28.7% |

So at 71% of them the truth is not quiet; it carries something smaller than 50 bp, and vg has written
it as something larger. A truth indel of 20-49 bp that vg calls at 50 bp or more becomes a structural
false positive that no truth SV can match, because the truth variant is a small variant. Selecting on
the *truth* indel's size, so nothing is conditioned on our own call, that is how often it happens:

| truth indel | type | n | median truth | median ours | our call reaches 50 bp |
|---|---|---|---|---|---|
| 10-19 | INS | 5,210 | 12 | 14 | 4.8% |
| 10-19 | DEL | 5,154 | 12 | 13 | 3.6% |
| 20-49 | INS | 2,936 | 26 | 26 | **12.6%** |
| 20-49 | DEL | 2,962 | 27 | 26 | **10.7%** |
| 100-299 | INS | 671 | 147 | 128 | 83.0% |
| 300+ | INS | 723 | 573 | 406 | 92.8% |

(chr1, chr6 and chr20.) Called sizes track the truth well in the median -- if anything vg *under*-calls
the largest -- so this is a tail, not a bias: about 10% of calls come out at twice the truth's size or
more, slightly more often for insertions than deletions. It is enough. 12.6% of a 20-49 bp population
that large is comparable to the whole 2,444.

### Offsetting insertion/deletion pairs: two real errors, and a population that is mostly not

Reading individual records turned up something a rate could not. At chr8:1,769,212 vg calls a 132 bp
deletion, and 137 bp later at chr8:1,769,349 it calls an insertion of **the same 132 bp of sequence**,
on the same haplotype. Reconstructing that haplotype from the calls and diffing it against the
reference over the span, the two records together express a net change of **+5 bp**, and the truth
carries no variant within 300 bp. At chr1:206,041,378 a 177 bp deletion pairs with a 176 bp insertion
241 bp later, together expressing **-1 bp**, where the truth carries exactly a 1 bp deletion.

Both of those are real errors and the reconstruction is not in doubt. What does not survive is the
generalisation from them.

Counting pairs genome-wide -- opposite direction, sizes within 20%, same genotype, within 500 bp --
gives 1,100 records in 550 pairs, 3.00% of calls of 50 bp or more. **Against a null that keeps every
position and shuffles the size and direction labels within each contig, chance alone produces
1.31%.** So the effect is real but 2.3x over chance, and roughly 480 of those 1,100 records are
coincidence -- two unrelated calls in the same repeat-dense neighbourhood.

Three mechanisms were tested against the population and none of them survived:

| hypothesis | prediction | measured |
|---|---|---|
| nested calling creates them | absent, or much rarer, under `--no-nested` | present at 2.58% there against 3.00% now |
| the same sequence is being placed at the wrong copy | the two members carry identical sequence | identical in 0.4%; alignment-aware similarity ≥0.90 in 12.3% |
| the likelihood is tied and the tie-break is arbitrary | paired records have much lower GQ | median GQ 8 paired against 10 unpaired |
| consecutive snarls genotyped independently | the pair shares a snarl boundary node | 6.4%, against 5.8% for any nearby ≥50 bp call |

So the chr8 locus -- identical sequence, consecutive snarls sharing boundary node 189964199 -- is a
genuine instance of a comprehensible failure mode and is *not* representative of the 1,100.

And on those individual loci the framing "bad representation" is itself too kind. If the caller
believes the repeat unit sits at copy B rather than copy A, then a deletion at A and an insertion at
B is the *correct* VCF for that belief. The belief is what is wrong. It is a genotyping error in an
ambiguous tandem array, which then necessarily costs two structural false positives instead of one
error -- so the doubling is a consequence of the mistake, not a separate defect in how records are
written.

What remains true and useful: the pair rate is 9.3% among false positives with nothing to compare
against, against 1.6% among calls that matched. Whatever produces the enrichment is concentrated in
the population this section is about. It has not been identified.

### The large-insertion excess

One anomaly is not accounted for by either mechanism. Among calls changing length by 2 kb or more:

| | INS | DEL | ratio |
|---|---|---|---|
| truth SVs | 579 | 526 | 1.10 |
| vg calls that matched | 461 | 380 | 1.21 |
| PanGenie calls scored false | 118 | 82 | 1.44 |
| **vg calls scored false** | **473** | **117** | **4.04** |

The truth is balanced and vg's *correct* large calls are balanced, so a 4:1 insertion skew among its
wrong ones is not biology and is not shared by the other tool. These sit at very high depth -- median
DP 563 against 42 for SV calls generally -- so they are collapsed repeats, and their median DR is
**0.38** against 0.95 for large insertions that matched: they claim about 2.6x more sequence than the
reads support. The depth model is already saying so; nothing acts on it.

### The near-misses are near

Of the 9,672 compared and rejected: median sequence similarity to the truth variant is 0.75, and

- 5,419 clear truvari's 0.70 sequence bar
- 3,867 clear its 0.70 size bar
- **3,140 clear both and are scored false anyway** -- 1,562 because the truth variant they best match
  was already taken by a different call
- 2,279 clear the sequence bar and fail on size, at a median size similarity of 0.62

**That last group is not a vg problem.** PanGenie's equivalent share is *higher*: 3,573 of its false
positives (33.9%) clear both bars, against vg's 3,140 (25.9%). Truvari matches one-to-one, so where
truth SVs cluster some good calls are scored false whatever the caller does, and this measures the
metric rather than either tool. Reading vg's 25.9% on its own would have made it look like a defect.

This is an upper bound on the assignment artefact rather than an exact count: it reads truvari's
annotation of each call's nearest candidate, not a re-run of its assignment.

### They are low-confidence, and that does not help

| | matched a truth SV | scored false |
|---|---|---|
| median GQ | 34 | 5 |
| median GQN | 0.275 | 0.040 |
| median DR | 0.81 | 0.42 |
| DR < 0.75 | 44.3% | 80.9% |
| heterozygous | 72.3% | 79.7% |

Both quality signals separate the populations. **Neither separates them enough to be worth using:
no threshold on GQ, on the depth- and ploidy-invariant GQN, or on DR raises SV F1.**

| gate | TP kept | FP kept | SV F1 |
|---|---|---|---|
| none | 13,291 | 12,116 | **0.5488** |
| DR >= 0.3 | 11,869 | 8,457 | 0.5454 |
| GQ >= 3 | 11,349 | 7,594 | 0.5387 |
| DR >= 0.5 | 10,009 | 4,916 | 0.5235 |
| GQ >= 10 | 9,650 | 4,845 | 0.5104 |
| GQN >= 0.02 | 9,788 | 5,962 | 0.5026 |
| GQN >= 0.05 | 8,806 | 4,688 | 0.4801 |
| GQ >= 20 | 8,207 | 3,385 | 0.4695 |

Recall in that table is accounted on the truth side, as truvari scores it, so the ungated row
reproduces truvari's own SV F1 exactly; the generator asserts that it does. Full working:
[sv-quality-gates.md](sv-quality-gates.md).

The reason is arithmetic, and it can be stated exactly. Removing a true call costs twice -- once in
the numerator and once as a new false negative -- so a gate improves F1 only if it removes more than
(TP + FP + FN) / TP = **2.64** false positives per true call it discards. Nothing available reaches
that. The closest is a depth-ratio gate restricted to the large insertions in the section above,
where the signal is strongest: dropping calls of 2 kb or more with DR below 0.5 removes 345 false
positives for 142 true ones, a ratio of 2.43, and still lands slightly *below* the ungated F1 at
0.5485. Confidence filtering is a lever for precision-limited
call sets, and on SVs this one is recall-limited. **Zygosity is not a discriminator either** --
79.7% of false positives are heterozygous against 72.3% of the calls that matched, a real
difference in the wrong direction to act on: there is no gate that keeps homozygotes only.

The two false-positive populations are not alike, and the difference points away from confidence as
the explanation:

| | no truth SV nearby | compared and rejected |
|---|---|---|
| n | 2,444 | 9,672 |
| median GQ | 13 | 4 |
| median DR | 0.46 | 0.41 |
| heterozygous | 87.2% | 77.6% |
| 700 bp and over | 14.4% | 8.7% |

The calls with nothing in the truth nearby are the *more* confident group, are 87.2% heterozygous,
and carry a long tail -- 148 insertions of 2-10 kb and 90 over 10 kb. A wrong call that the model
is confident about is a different problem from a call that was nearly right.

That long tail has a consistent shape. All 135 insertion false positives of 10 kb or more are
heterozygous, they sit at very high read depth (DP in the thousands, so these are collapsed
repetitive regions), nearly every read is assigned to the alternative allele, and their median DR is
0.37 -- a third of the reads a 40 kb insertion should produce. The depth ratio is already saying the
call is implausible; nothing acts on it, and per the gate table above nothing profitably could at
this operating point.

Records at an identical position are mildly enriched among the false positives -- 13.9% share a
position with another structural call against 10.2% of the calls that matched -- so the occasional
"heterozygous insertion called twice on opposite strands" pattern is real but is not a major
contributor.

**None of these carry a ploidy-coherence FILTER, because none exists any more.** Under single-sweep
nested calling a nested chain is genotyped at the ploidy its parent's settled genotype implies, so
parent and child cannot disagree; all 5,037,820 records genome-wide are `PASS`. The earlier arm
flagged 198 of its structural false positives this way, which was already too few to be the
explanation -- so closing the coherence defect was never going to move this gap, and it did not:
SV F1 went 0.54854 inline to 0.54861 here. See
[nested-calling-design.md](nested-calling-design.md).

### The false negatives are 50-300 bp, heterozygous, and in tandem repeats

Of 10,053 distinct missed truth SVs, 8,172 (69.7%) are also missed by PanGenie -- a shared,
panel-limited floor. 1,881 are missed by vg alone, and those are the ones whose alleles the panel
demonstrably carried.

76% of that actionable set is 50-300 bp and heterozygous:

| type | size | zygosity | vg-only FN | PanGenie-only FN |
|---|---|---|---|---|
| INS | 50-100 | het | 411 | 320 |
| INS | 100-300 | het | 387 | 264 |
| DEL | 50-100 | het | 337 | 282 |
| DEL | 100-300 | het | 301 | 291 |

Repeat context explains the shared floor and nothing about the difference between the two tools:
89.5% of vg-only misses are in a tandem repeat, against 86.6% of PanGenie-only misses and 89.4% of
those missed by both -- while only 56.0% of vg's true calls are. Being in a tandem repeat predicts
that an SV is hard; it does not predict which tool will miss it.

Genotyping is a real but second-order loss: among locus-matched SVs vg gets the genotype right 88.0%
of the time against PanGenie's 91.2%, and requiring a correct genotype widens the F1 gap from 0.0254
to 0.0408.

## Where this leaves the work

- **The recall side is the binding constraint**, and it is no longer about records that were never
  written. 13.7% of the actionable misses have no record nearby, down from 46.7%; the rest have a
  record that is the wrong size or the wrong shape.
- **The false-positive excess is concentrated in 2,444 calls with no truth SV within reach**, and it
  is mostly a representation problem rather than an evidence one. At 71% of those loci the truth
  carries variation under 50 bp that vg has written as something larger; 221 are offsetting
  insertion/deletion pairs whose net effect on the haplotype is a few bases. Only 28.7% sit where the
  truth carries nothing at all within 300 bp.
- **The offsetting pairs are 2.3x over chance, not a wholesale defect.** 550 pairs against a null of
  roughly 240, six times commoner in this population than among calls that matched. Two loci were
  traced to a genotyping error in an ambiguous tandem array, which costs two structural false
  positives for one wrong choice; no mechanism yet explains the population, and four hypotheses have
  been ruled out.
- **The 4:1 insertion skew among large false calls is unexplained**, against 1.10 in the truth and
  1.21 among vg's own correct large calls. Their median DR of 0.38 says the reads do not support the
  sequence claimed.
- **Nothing here is fixed by filtering.** Every confidence gate measured costs more F1 than it
  saves.
- **The long-record pathology is closed as a line of work.** No insertion or deletion false positive
  of 300 bp or more, in either arm or in PanGenie, is a bundle of sub-structural changes. Further
  effort on decomposition would be spent on the 310 records that still bundle, two thirds of them
  same-length substitutions, against a 12,163-record false-positive count.

Full working: [sv-delta.md](sv-delta.md), [sv-fn-mechanism.md](sv-fn-mechanism.md),
[sv-unmatched.md](sv-unmatched.md), [sv-fp-anatomy.md](sv-fp-anatomy.md).
Regenerate with `scripts/wgs/sv_delta.py`, `scripts/wgs/sv_fp_anatomy.py`,
`scripts/wgs/sv_fn_mechanism.py`, `scripts/wgs/sv_unmatched_why.py` and
`scripts/wgs/sv_quality_gates.py`, all against `--score work/wgs-fixed/score`.
