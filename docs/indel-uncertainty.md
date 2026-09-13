# Indel uncertainty is modelled as a penalty, not a probability

Measured 2026-09-13 on chr20 ONT 43x. Four probes, each adversarially re-derived by a second
independent analysis; **three of the four original headlines were refuted and replaced**. What
follows is the corrected set.

The question: the ONT preset raises `--mismap-min` to 0.05, and a homopolymer-conditional floor
buys a further +0.018 indel F1 held out (`tier2-parameters.md`). Is the mismapping clamp the right
instrument for indel uncertainty, or a proxy for something else?

**It is a proxy.** `e_r` is the mixture weight on "this read did not come from this locus". A
homopolymer miscount is a read that *did* come from here and is correctly placed, merely misspelled.

## 1. The gap term is not normalised by anything

`rel(r,a) = exp(log_base * (score(r,a) - max_a' score(r,a')))`, `allele_likelihood.cpp:1068`, with
`log_base = 1.3833252687`. That constant is obtained (`alignment_scorer.cpp:30-99`) by bisecting a
partition function that sums **only over the 4x4 substitution matrix** weighted by nucleotide
background frequencies (`:118-132`). There is no gap state in it. Nothing normalises
`exp(lambda * gap_score)`.

A gap is charged as a bare affine integer `-(gap_open + (len-1)*gap_extend)`
(`alignment_scorer.cpp:24-26`), at five call sites in `allele_likelihood.cpp` -- `:611` and `:614`
inside `score_shared_node` are the dominant path, `:701/:729/:740` the node-mismatch path. Its
entire signature is `(size_t gap_length)`: no sequence, no read offset, no quality is passed or
available, and neither `MatrixAlignmentScorer` nor `QualAdjAlignmentScorer` overrides it.

A substitution, by contrast, is genuinely probabilistic: `qual_adjusted_matrix`
(`alignment_scorer.cpp:438-492`) builds a per-quality log-odds table from `err = 10^(-q/10)`. So the
model is **half-calibrated** -- substitutions are probabilities, indels are penalties, and both are
exponentiated with the same constant as though they were the same kind of object.

Empirical confirmation of context-freeness: across 3,701,501 paired `(site, read, allele)` rel cells
from two runs differing only in `gap_open` (1 vs 6), 82.57% are bit-identical and **every** changed
cell moves by an exact integer multiple of the gap-open delta. A context-sensitive cost could not
produce that.

**A consequence worth stating separately.** The substitution matrix is hard-coded (`call_main.cpp`
passes `default_score_matrix`; there is no `--match` or `--mismatch` flag) while `gap_open` and
`gap_extend` are free integers fitted by F1. So every time the gap scale is tuned, the whole
score-to-nats conversion moves against a normalisation that does not. `--gap-open 1` is not merely
"cheaper gaps": it rescales the implied probability of every indel relative to a substitution
likelihood that stayed put.

## 2. The defect is DIRECTIONAL, and roughly constant in run length

The first analysis claimed a log-linear error rate rising x1.60 per base and a flat penalty "82x too
permissive at L=2". **Refuted.** Counting every deletion against every run it destroys (the original
discarded the 32% of deletions spanning a run boundary while keeping those runs in the denominator),
the clean-stratum per-run indel error rate is 0.55% at L=1, 0.96% at L=2, 4.11% at L=5, 12.6% at
L=8, 47.9% at L=12. A rise of 87x, not 238x, averaging x1.50 per base -- and **not log-linear**: a
quadratic cuts the SSE by 72% and the fitted line predicts a probability above 1 by L=13.

The real defect is not magnitude but direction, and it is roughly flat in L. In a two-allele vote
the flat penalty demands **66.7% of reads to call a contraction** where the data demand 55-61%, and
only **33.3% to call an expansion** where the data demand 40-46%. The model forces
`A/|B| = 2.00` and `A - |B| = match = 1` score unit; empirically `A/|B|` is 1.24-1.51 and `A - |B|`
falls from 1.19 units at L=2 to ~0.4 at L=10-11. **At long runs the mis-set parameter is `match`,
not `gap_open`** -- and `match` has no flag.

**Confirmed independently in the calls.** chr20 1 bp indels at reference HP >= 5:

| | INS | DEL | INS share |
|---|---|---|---|
| TP | 3007 | 3268 | 47.9% |
| FP | **1501** | 726 | **67.4%** |
| FN | 578 | 646 | 47.2% |

Insertion precision 0.667 against deletion 0.818; recall identical (0.839 / 0.835). The bias shows
up **entirely as insertion over-calling, not deletion under-calling**, because the basecaller's own
error is directional the same way (43.6% insertion against 25.5% deletion miscount at HP>=13): the
two compound for insertions and cancel for deletions. Matching insertion precision to deletion
precision would remove 832 FPs -- **17.5% of all chr20 ONT small-variant FPs**.

## 3. It is per-read mis-specification, NOT a read-independence violation

The first analysis claimed overdispersion proved correlated reads, so ~40 reads carry ~8-11 reads'
worth of information. **Refuted on a point of principle**: by de Finetti, an exchangeable binary
sequence is exactly a mixture of iid Bernoulli, so "reads within a site are correlated" and "reads
are independent but per-read alt probability varies between sites" are the *same distribution*. One
binomial observation `(y_i, n_i)` per site cannot separate them.

The decisive measurement is a control the first analysis had in hand and did not use. At
**truth-verified homozygous-alt** sites -- where there is no second haplotype for reads to come
from -- the per-read miscount rate rises from 2.7% at HP 1-2 to **35.6% at HP >= 13**, and **39.3%
of those hom-alt sites show a het-looking read split** against 0.4% at HP 1-2. Correlated error
cannot produce that; only per-read accuracy can.

Conditioning on insertion-vs-deletion alone cuts phi from 4.70 to 3.42 -- because the true het
support fraction is **0.403 for insertions and 0.597 for deletions**, where the model assumes 0.5.
Four coarse VCF-visible covariates take it to 2.69, removing 54-62% of the run-length-attributable
excess. Residual phi 2.69 against a 1.46 baseline, so an effective-sample-size penalty **may** still
be wanted -- but it is an open question, not an established need.

**This connects to the mixture weights.** `w_h ~ U_h + R - 1` handles the *sampling* asymmetry and
correctly flattens to 0.5 at ONT read lengths. The 0.403/0.597 split is a second source of
support-fraction asymmetry -- directional basecaller miscount -- worth ~0.1 in read fraction, which
the model represents nowhere. The weight is right about sampling and silent about spelling.

Also: GQ scores a 30% error rate identically to a 3% one (4.4x at GQ 40-60).

## 4. "More faith in the prior": frequency, not presence -- and it cannot help recall

- **Presence is not the deficit.** The panel spells the truth 1 bp indel allele at 95.8% of HP>=5
  sites, against 95.8% for SNVs and 94.3% at HP 1-2.
- **Presence can never veto a false positive**, near-tautologically: the caller can only emit an
  allele the graph holds, and at a panel-monomorphic long run no +-1 traversal exists at all in
  94.7% of cases. So the called allele is panel-present at 95.9% of HP>=5 FPs. Asking "is this in
  the panel?" re-asks a question the graph construction already settled.
- **The exploitable signal is panel allele FREQUENCY**: mean panel AF 0.208 at FPs against 0.559 at
  TPs, and `r(panel AF, GQ) = 0.11` -- orthogonal to what the caller already reports. An `AF < 0.25`
  cut removes **2.50 FPs per TP lost** where a TP-loss-matched GQ cut removes 0.45. Worth about
  +0.03 F1 on this class, additional to the panel-frequency prior that is already active
  (`graph_caller.cpp:2154`, exponent default 5) and currently slightly FP-enriching here.
- **The panel IS a bottleneck at false negatives.** Truth alleles the caller missed are panel-carried
  at 84.4% against 98.3% for those it found -- a 10.6x enrichment of panel absence at misses. About
  **2.8 of the 17.9 points of recall deficit are unrecoverable by any read-side change.**
- At FP-relevant run lengths the panel is not crisp anyway: 25.8% monomorphic, mean 2.20 alleles,
  modal length == reference only 67.7%.

So leaning on the prior is a precision instrument only, and a modest one.

## 5. A separate defect, found in passing

`score_read_against_allele`'s greedy single-pass anchor desynchronises when the read visits a node
the allele lacks: the read's extra node is matched against the allele's *next* node
(`allele_likelihood.cpp:713-731`) and the remainder falls through to the allele-exhausted branch
(`:733-741`), double-charging the downstream boundary. At snarl `>96335708>96335710` the same 1 bp
event costs 1 score unit in one direction and 67-68 in the other. 22.18% of paired read-rows have a
best allele carrying at least one gap.

This is an asymmetry **bug**, not a calibration choice, and it plausibly contributes to the
insertion/deletion direction bias in section 2. It should be fixed before the calibration question
is measured again, because it contaminates exactly that measurement.

## Phase 0, done: the anchor desynchronisation was most of the direction bias

Fixed 2026-09-13. The walk searched ahead in the ALLELE for each read node, but on finding
none assumed the read node SUBSTITUTED for the allele's current node and consumed it. When
the read node was a pure insertion the consumption burned an anchor the read still needed;
its own visit then found the allele exhausted and was charged again, so the FLANKING node's
whole length was billed twice for a one-base event.

Measured on the existing test fixture, at the default gap_open 6: a one-base **deletion**
costs one gap_open, `rel = 2.485e-4`; the identical one-base **insertion** gave
`rel = 0.0` -- the header's value for "this allele cannot place the read at all". The fix
adds the mirror of the anchor search, "is the allele's current node still to come in the
read?", and charges an insertion. `rel` becomes 6.2317e-5, and the residual against the
deletion is exactly one score unit -- the forgone match credit, which belongs to the score
model and not the walk.

| | chr20 ONT | chr6 ONT (held out) | chr20 short reads |
|---|---|---|---|
| SNV F1 | 0.98581 -> 0.98573 | 0.98800 -> **0.98816** | 0.98517 -> 0.98518 |
| Indel F1 | 0.83725 -> **0.85452** | 0.86023 -> **0.87483** | 0.92832 -> 0.92858 |
| ALL F1 | 0.95152 -> **0.95574** | 0.95962 -> **0.96318** | 0.97240 -> 0.97246 |

**+0.0173 indel F1 on chr20 and +0.0146 on the hold-out**, improving precision AND recall,
at no SNV cost -- chr6's SNV F1 rises. chr6 carries 85% of chr20's magnitude, better transfer
than any fitted parameter in this project, which is what fixing a bug rather than tuning
should look like.

And it was most of the asymmetry section 2 measured:

| HP>=5 1 bp indels | chr20 pre -> post | chr6 pre -> post |
|---|---|---|
| INS/DEL precision gap | -0.1512 -> **-0.0517** | -0.1132 -> **-0.0226** |
| FP insertion excess | +19.5% -> **+7.6%** | +16.5% -> **+3.8%** |
| FPs recoverable by closing it | 833 -> **257** | 1606 -> **295** |

**66% of the direction bias on chr20, 80% on chr6.** So the "17.5% of all ONT FPs" figure in
section 2 was measured on contaminated data and the real remaining target is about a third of
it. That is the whole argument for fixing the bug before measuring the model.

**Insertion coalescing: tried and dropped.** The deletion side sums every skipped allele node
into one affine gap; the insertion side charged one gap_open per node, so two adjacent
inserted nodes cost 2*gap_open where the mirrored deletion costs gap_open + gap_extend.
Making them symmetric is **byte-identical on ONT** (0 genotype changes on chr20 AND chr6) and
on short reads moves 53 genotypes plus 2,020 phase swaps for -0.00014 indel F1. No measured
benefit, a small cost, so it was removed. The asymmetry is real but it dissolves under a
calibrated gap cost, which has no gap_open/gap_extend structure to coalesce.

## What this implies, in order

1. ~~Fix the anchor desynchronisation.~~ **Done, above.**
2. **Give the gap a probability** -- either a gap state in the partition function so `lambda`
   normalises over it, or bypass the score path for indels and charge `ln P(L' | L)` from the
   measured spectrum. Note that at long runs the mis-set parameter is `match`, which has no flag, so
   this cannot be reached by tuning the two knobs that exist.
3. **Then test the falsification**: does the floor's SNV/indel conflict collapse? If it does, the
   floor was a proxy, confirmed, and it can return to 0.02 for every class.
4. **Panel AF as a separate orthogonal term** -- explicitly a precision instrument, not a recall fix.

## Phase 2, done: --insertion-nats, and why no existing knob could do it

The flat affine gap charges one constant whichever side carries the extra bases. Two
reasons neither existing flag can correct that:

**`--gap-open` couples the ratio to the scale.** The flat model forces
`A/|B| = (gap_open + match)/gap_open`, so the preset's `gap_open 1` gives exactly 2.00
against an empirical 1.24-1.51, and the two-allele decision boundary sits at 1/3 where the
data want 0.40-0.46. Moving gap_open to fix the ratio also moves the absolute weight of
indel evidence, and the scale dominates. Re-swept on the FIXED walk, chr20 indel F1:

| gap_open | 1 | 2 | 3 | 4 | 6 (vg default) |
|---|---|---|---|---|---|
| indel F1 | **0.85452** | 0.81411 | 0.79297 | 0.78663 | 0.78490 |

`gap_open 1` survives the bug fix unchanged, and is worth **+0.0696** against the default
on the fixed walk -- more than the +0.048 recorded pre-fix, so the bug was partly masking
the preset's own value.

**The score is an int32.** One score unit is 1.3833 nats and the wanted correction is about
0.4 units, which the integer path cannot represent. This is a plausible part of why
`--qual-gap` failed: its absolute form moved 2 units where the preset had fitted 1.

So `--insertion-nats` adds real-valued nats **after** the score-to-nats conversion, at the
four sites where the read carries bases the allele lacks. Default 0.0 and byte-identical
there (115,267 records, cmp clean).

The value was predicted BEFORE the sweep, at 0.7-1.05 nats, from the measured `A/|B|`
target. The measured chr20 optimum is interior at **0.9** (1.3 is worse):

| nats | 0 | 0.5 | **0.9** | 1.3 |
|---|---|---|---|---|
| chr20 indel F1 | 0.85452 | 0.86076 | **0.86237** | 0.86006 |
| chr20 ALL F1 | 0.95574 | 0.95739 | **0.95787** | 0.95729 |
| chr6 indel F1 | 0.87483 | 0.87818 | **0.88005** | -- |
| chr6 ALL F1 | 0.96318 | 0.96397 | **0.96444** | -- |

**It is context-blind but effectively targeted.** 92.5% of the FP reduction lands in HP>=5
without the term knowing homopolymers exist, because that is where the insertion gaps are;
HP 1-2 improves too (FP 101 -> 87) rather than paying for it. That is the difference from
the mismap floor, which touches 96.6% of reads to fix a 40% cell.

**Direction-neutrality is NOT the objective, which corrected an earlier recommendation.**
0.5 nats is direction-neutral on chr20 (precision gap +0.0053, FP excess -0.8%) and 0.9
over-corrects into a mirror asymmetry (+0.0612, -9.9%). Neutrality looked like the
principled setting -- but it is a symptom measure, it is not even a single value (~0.5 on
chr20, ~0.3 on chr6), and both contigs prefer 0.9 on F1. The error surface is asymmetric
because the insertion FP pool is far larger than the insertion TP pool at risk: on chr6,
0.9 removes **913 insertion FPs for 98 insertion TPs**.

**What it costs, which the aggregate hides.** It biases toward shorter alleles generally,
not against insertions specifically -- a REF-supporting read scored against a deletion
allele also carries extra bases. Deletion FPs rise 1,894 -> 2,338 on chr6 and 778 -> 930 on
chr20. Net is clearly positive (chr6: FP -469, FN -139, TP -58) but it is a redistribution.
**Not added to `--preset ont`** for that reason. The recall asymmetry also worsens
monotonically with the flag (+0.0038 pre-fix, -0.0291 after it, -0.0550 at 0.5, -0.0773 at
0.9), and no single scalar can correct precision and recall asymmetries at once. The
run-length-conditional table remains the end state.

## Phase 3: the floor was a proxy, and the second floor is retired

Re-sweeping `--mismap-min` on the fixed walk, chr20:

| floor | 0.02 | **0.05** | 0.10 | 0.15 | 0.20 | 0.30 |
|---|---|---|---|---|---|---|
| SNV F1 | 0.98549 | **0.98573** | 0.98443 | 0.98303 | 0.98120 | 0.97524 |
| indel F1 | 0.84111 | 0.85452 | 0.84791 | 0.85446 | 0.86025 | **0.86249** |
| ALL F1 | 0.95207 | **0.95574** | 0.95299 | 0.95369 | 0.95378 | 0.94988 |

The SNV/indel conflict **persists in direction** -- indel F1 still peaks at 0.30 -- but:

| raising the floor 0.05 -> 0.20 | pre-fix | post-fix |
|---|---|---|
| indel F1 gain | +0.0211 | **+0.0080** |
| ALL F1 | +0.0014 (prefers 0.20) | **-0.0196 (prefers 0.05)** |

The prize shrank by 62-73%, the aggregate optimum moved back to the shipped 0.05, and the
curve is now noisy (0.10 dips below both neighbours) where it used to be cleanly monotone.

**The decisive comparison.** Floor 0.30 reaches indel F1 0.86249 by spending 0.0105 of SNV
F1 and 0.0059 of ALL F1. `--insertion-nats 0.9` reaches 0.86237 -- the same number -- while
SNV F1 *rises*. The floor was a proxy for two things it does not model: the anchor-walk bug,
and the direction miscalibration. Fix both properly and it has nothing left to buy.

So **`--mismap-min` stays a single value at 0.05**, and the context-conditional two-floor
scheme in `tier2-parameters.md` is unnecessary -- it was correcting the bug. That measurement
stands as a record of what the bug was worth, not as a proposal.
