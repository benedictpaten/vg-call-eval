# Quality signals: what ranks a call, and what does not

Everything on the accuracy pages is scored at every GQ, so none of it depends on the quality
field. This page is about the other question: given the calls, how well are they *ordered*,
and can the caller order them better using evidence it already computes and throws away.

Measured on four datasets — chr6 and chr20, each on a 4-haplotype and a 34-haplotype graph —
against two benchmarks, small variants (aardvark, GIAB `smvar`) and structural variants
(truvari, GIAB `stvar`). Eight dataset × benchmark combinations. A signal is only reported as
useful if it helps in all eight.

Scripts: `coverage_model.py`, `share_gq.py`, `depth_filter_sweep.py`, `pileup_guard.py`,
`allele_balance_by_length.py`, sharing `filter_lib.py`. Full derivations in
[planning/vg-call-eval-plan.md](../planning/vg-call-eval-plan.md) §9.25–§9.27; the model these
signals sit on top of is in [planning/vg-read-likelihood-design.md](../planning/vg-read-likelihood-design.md) §4.

---

## The blind spot

The genotype likelihood is

```
ln P(reads | G) = Σ_r  w · ln [ (1 − e_r) · Σ_{h∈G} (1/|G|) · rel(r,h)  +  e_r ]
```

and GQ is the phred-scaled gap between the best genotype and the runner-up:

```
GQ_ratio = (10/ln10) · Σ_r  w · [ ln term_r(G₁) − ln term_r(G₂) ]
```

Each read contributes only the *difference* its term makes between the two genotypes. Take a
read whose best-fitting allele is in **neither** `G₁` nor `G₂`: its mixture is small under
both, both terms collapse to ≈ `e_r`, and the read drops out of GQ entirely.

So GQ measures how far ahead the winner is and never whether the winner accounts for the
data. A site where a third of the reads prefer an uncalled allele scores the same as a site
where none do. That is not hypothetical — on a 34-haplotype graph, false SV calls have a
median `AD/DP` of 0.63 against 1.00 for true ones, on both chromosomes.

## The fields

| field | definition |
|---|---|
| `DP` | rows in the site's likelihood matrix: reads that placed on **at least one** allele. Reads placing on nothing have no row maximum to normalise by and are dropped and counted separately, so `DP` is already a post-filter count |
| `AD` | per emitted allele, the number of reads whose best-fitting allele it is. Ties **split** fractionally rather than going to the lowest index — at multi-allelic sites many reads are genuinely undecided |
| `BL` | mean over reads of `best_ln`, the row divisor: how well reads fit their best explanation here in absolute terms. Nearly independent of GQ (r = +0.18) |
| `GQI` | GQ from the likelihood ratio alone, undiscounted |
| `GQ` | `GQI × sum(AD)/DP`, using the exact fractional support rather than the rounded `AD` |

**`AD` does not sum to `DP`**, and at a busy site falls a long way below it. Two reasons: split
ties round, and — the large one — only alleles that reached the VCF record get a column, while
the genotyper scored every allele the site offered. That shortfall is the useful part. It is
the share of reads the called genotype fails to explain.

---

## What works: the explained-share discount

`share` was the only signal positive in all eight combinations, which made it the one
candidate for changing the emitted quality rather than leaving it to a downstream filter. The
supporting evidence used per-dataset logistic weights, which cannot ship, so the question was
re-asked in shippable form: **is there a single fixed formula, nothing fitted, that improves
the ranking everywhere?**

| form | small variants | SVs | verdict |
|---|---|---|---|
| `GQ · share` | +0.005 to +0.008 AUC | +0.004 to +0.031 | **all 8 up, and up at nearly every operating point** |
| `GQ · share²` | +0.008 to +0.012 | +0.006 to +0.046 | better AUC, worse at high recall in 3 cells |
| `GQ · share⁴` | +0.012 to +0.015 | +0.008 to +0.056 | better AUC still, worse at high recall in 3 cells |
| `min(GQ, −10log₁₀(1−share))` | **−0.009, −0.006** on chr20 | +0.021 to +0.087 | breaks small variants |

The phred form is the one with a probabilistic reading, and the arithmetic shows why it fails:
at `DP` 30 a single stray read gives share 0.967, capping a perfectly good SNV at Q14.8. A
tolerance before the cap engages repairs the AUC but still loses at the operating points.
Linear was the only form that made none of the eight worse, so linear shipped.

Measured end to end on the caller's own output — `GQI` against `GQ` — which beats the offline
prediction in every cell because it uses the unrounded support:

| dataset | benchmark | AUC `GQI` → `GQ` | FP at moderate recall | FP at high recall |
|---|---|---|---|---|
| chr6 4-hap | small variants | 0.8128 → 0.8217 | 3,494 → 3,354 | 4,754 → 4,694 |
| chr6 4-hap | SVs | 0.6326 → 0.6422 | 165 → 160 | 224 → 224 |
| chr6 34-hap | small variants | 0.9171 → 0.9262 | 1,911 → **1,702** | 3,180 → 3,041 |
| chr6 34-hap | SVs | 0.6835 → **0.7207** | 239 → **199** | 351 → 315 |
| chr20 4-hap | small variants | 0.8151 → 0.8241 | 1,546 → 1,525 | 2,594 → 2,593 |
| chr20 4-hap | SVs | 0.6276 → 0.6396 | 103 → 103 | 148 → 150 |
| chr20 34-hap | small variants | 0.9214 → 0.9287 | 847 → **798** | 1,618 → 1,586 |
| chr20 34-hap | SVs | 0.6524 → **0.6907** | 134 → **121** | 200 → 190 |

Eight of eight improve on AUC; fifteen of sixteen operating points improve or tie. The single
exception is +2 false calls out of 148, on the smallest set in the table (n = 368).

**Unfiltered F1 is unchanged on every benchmark** — 286,557 variants and SV F1 0.5349 on chr6
4-hap before and after — as it must be, since the discount rescales a quality and does not
change a genotype.

### The cost

A discounted GQ is no longer the phred-scaled odds of the top two genotypes. It is a quality
score, not a calibrated probability. Anything that assumes GQ is a posterior — an HWE prior, a
Bayesian merge across samples — must read `GQI`, which is why `GQI` is emitted unconditionally,
including when the discount is off and the two are equal. `--no-share-quality` restores the old
behaviour exactly.

The pre-existing caveat is unchanged: reads are treated as independent, so confidence
accumulates like `R` rather than `√R`. GQ was already over-confident at depth.

---

## What does not work

### Refuted at the assumption, not just the outcome

**Length-aware prior on expected reads per allele.** The geometry is compelling: a 151 bp read
cannot span a 120 bp insertion with anchors either side, so `w_i = max(1, R − len_i − 2F)`.
`allele_balance_by_length.py` measures what that predicts, over true ref/het calls:

| ‖len ALT − REF‖ | n (chr20 34-hap) | mean ALT share | mean AD/DP |
|---|---|---|---|
| 0 | 45,564 | 0.499 | 1.002 |
| 4–7 | 1,539 | 0.502 | 0.958 |
| 16–31 | 277 | 0.499 | 0.938 |
| 32–63 (truvari) | 19 | 0.531 | 0.938 |
| 128–255 (truvari) | 9 | 0.529 | 0.859 |
| 256–1000 (truvari) | 19 | 0.605 | 0.834 |

The ALT share is 0.50 in **every** class from 0 bp to over 1 kb, on both chromosomes under both
benchmarks. The geometry never applied — the caller never required a read to span an allele, so
a read covering one breakpoint still fits one allele better than the other. Scored against this
prior, true het SV calls look impossible and the statistic ranks **worse than chance**
(AUC 0.46).

What the same table *does* show is `AD/DP` falling with length, 1.00 → 0.83. That is goodness of
fit, not skew, and it is what the share discount uses.

**Binomial coverage skew, even with the correct flat null.** The weakest of five signals in all
eight cells: +0.003 to +0.006 AUC over GQ. Allele balance carries real signal — an earlier claim
that it "does not separate at all" was wrong and was retracted — but it is not where the value
is.

### Real signal, dominated by GQ at matched recall

The test any hard filter must pass: it discards calls to buy precision, and so does lowering the
GQ threshold. The only question is which buys more.

**Minimum depth**, absolute or relative to a local median. Fails in all eight cells — GQ
thresholding reaches higher precision at the same recall every time (0.9816 against 0.9762 on
chr6 34-hap small variants; 0.4516 against 0.4145 on its SVs). The mechanism is obvious in
hindsight: few reads means a small likelihood gap, so low depth already depresses GQ. There is
nothing left for a separate cut to add.

**The two-condition pile-up guard** (`DP > 1.3 ×` local median **and** `AD/DP < 0.8`). Flags a
population 71–78% false on the rich graph and lifts SV precision 0.417 → 0.506. It reaches ~0.72
recall, and GQ thresholding to the same recall gives ~0.54. It beats *no* filter, not the
ranking already in the file.

**Maximum depth as a default.** ΔF1 against no filter, cutting on `DP` over a rolling local
median:

| cut | chr6 4-hap sm/SV | chr6 34-hap sm/SV | chr20 4-hap sm/SV | chr20 34-hap sm/SV |
|---|---|---|---|---|
| min 0.5× | −0.003 / −0.011 | −0.003 / −0.014 | −0.003 / −0.004 | −0.003 / −0.009 |
| max 2× | −0.000 / −0.099 | −0.001 / −0.043 | −0.000 / −0.043 | −0.001 / −0.027 |
| max 5× | −0.000 / −0.016 | −0.000 / **+0.022** | −0.000 / −0.009 | −0.000 / **+0.016** |
| max 20× | −0.000 / +0.000 | −0.000 / +0.016 | −0.000 / +0.000 | −0.000 / +0.008 |

It passes the GQ test in exactly one configuration — 5× the local median, structural calls,
34-haplotype graph, worth about +0.025 precision. Dominated on both 4-haplotype runs at every
threshold in the grid, and on small variants in every cell of all four. Tight cuts are far worse
than they look: at 2× the local median, GQ wins by 0.09 precision on chr6 4-hap SVs. `DP` is
emitted, so this stays a downstream option rather than a default.

### Real signal, but the sign does not transfer

**Depth as a linear ranking term.** Within one dataset, local depth ratio has AUC 0.65 against
small-variant labels (*low* depth means false) and 0.37 against SV labels (*high* depth means
false). One coefficient improves one class and damages the other.

**`BL` / `best_ln` as a global term.** The mirror image: strongest single small-variant signal
(AUC 0.79–0.84 alone), worthless to harmful for SVs (0.489 on chr6 34-hap; it *lowers* held-out
AUC from 0.719 to 0.571 on chr20 34-hap).

Both need conditioning on called-allele size before they can be used. That remains open, and it
is the concrete obstacle in front of the depth-plausibility term the design doc §5.3.3 proposed.

### Weaker variants of things that do work

**`share_resid`**, share minus its length-class baseline: about +0.003 over raw `share`, not
worth a length table.

**A Hardy–Weinberg genotype prior** reallocates precision and recall without adding
discrimination (plan §9.23).

---

## Method notes, including two bugs that produced believable wrong answers

- **Every signal is scored against GQ, never in isolation.** A weak signal usually correlates
  with correctness; the only interesting question is whether it adds to a ranking already
  available.
- **FP count at matched recall, not just AUC.** With false calls at 3% of the set, precision
  moves from 0.9859 to 0.9933 and looks like nothing, while the same change halves the false
  calls. AUC rewards the whole ranking; the operating point depends only on the tail.
- **`--min-svlen 50` against aardvark labels is meaningless.** The small-variant benchmark has
  no record ≥50 bp, so restricting its labels by size left 18 true positives. SV questions go to
  truvari.
- **An equal-weight rank sum gave opposite answers on two subsets of the same data**, because it
  dilutes a strong signal with a weak one. Fitted weights or nothing.
- **A relative dataset path printed four empty tables as a clean result** rather than erroring,
  because the caller was catching the missing-file exception in order to skip datasets that
  genuinely lack one benchmark.
- **Recall as query-side TP over the truth total is wrong for truvari**, which matches several
  base records to one query record when the caller emits an SV the benchmark decomposed — 792
  base matches against 425 query records on chr6 4-hap. It made every filter look nearly twice
  as damaging as it is.
