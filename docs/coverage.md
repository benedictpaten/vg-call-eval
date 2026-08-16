# Calling across coverage and ploidy: the titration

Stage 0 of the coverage-robustness work. The point is to find out what breaks when the caller is
given less data than the ~30x diploid it was tuned on, and to decide -- from measurement rather
than from argument -- what a coverage- and ploidy-robust quality score should be.

The headline is a **negative result that redirects the design**: every simple depth rescale of GQ
makes the picture worse once both ploidies are in view, and ploidy turns out to be a much larger
source of miscalibration than depth.

## Method

Two series, one per ploidy:

| series | contig | ploidy | source | arms |
|---|---|---|---|---|
| diploid | chr20 | 2 | 30.28x | 5, 10, 15, 20, 25, 30x |
| haploid | chrX non-PAR | 1 | 14.63x | 2.5, 5, 7.5, 10, 12.5, 14.6x |

**The two series share an x-axis in reads per haplotype.** chr20 at 30.3x across two haplotypes is
15.1 per haplotype; a male chrX at 14.6x across one is 14.6. Choosing the chrX levels at half the
chr20 levels puts the series on a common axis, so a difference between them is ploidy and not
depth. That is the only reason two series are worth running.

**Subsampling is nested.** A read is kept at level c when `crc32(name)/2^32 < c/source`, so the 5x
set is a subset of the 10x set and so on. Independent draws per level would put sampling noise into
every pairwise comparison; nesting makes them paired, which turns out to matter (see the GQ scaling
correction below). Mates share a name in this GAF, so hashing the name keeps a pair together.

**A GAF-Base database per level, not `--gaf-reads` on the subsampled GAF.** An in-memory read
source answers exactly and so has no fetch window; `local_read_rate` returns 0 for a window-less
source and the depth term switches itself off rather than inventing a rate. Titrating that way
would have measured a model with `--depth-term` silently disabled -- one of the very parameters
under study.

Reproduce with `scripts/coverage/{subsample_gaf,titrate,bench_coverage,normaliser_eval}.py|sh`.

**Controls.** The full-coverage arms reuse the whole-genome reads database, so they must reproduce
published numbers. chr20 at 30x gives ALL F1 **0.9645** against the published tier-2 `readlik`
**0.9645**; chrX at 14.6x gives **0.9362** against the whole-genome run's **0.9362**. Both exact.

**Range caveat.** The source data is 30.3x, so this covers 5-30x. The 30-40x end of the intended
range is extrapolation, not measurement.

## 1. Accuracy degrades gracefully, and differently by ploidy

| chr20 | medDP | ALL F1 | SNV F1 | Indel F1 | recall | precision |
|---|---|---|---|---|---|---|
| 5x | 5 | 0.9008 | 0.9224 | 0.8225 | 0.8414 | 0.9693 |
| 10x | 10 | 0.9419 | 0.9600 | 0.8772 | 0.9106 | 0.9754 |
| 15x | 14 | 0.9546 | 0.9709 | 0.8971 | 0.9329 | 0.9774 |
| 20x | 19 | 0.9609 | 0.9755 | 0.9096 | 0.9430 | 0.9795 |
| 30x | 29 | 0.9645 | 0.9781 | 0.9177 | 0.9495 | 0.9801 |

| chrX | medDP | ALL F1 | SNV F1 | Indel F1 | recall | precision |
|---|---|---|---|---|---|---|
| 2.5x | 2 | 0.8471 | 0.8852 | 0.7299 | 0.8113 | 0.8863 |
| 5x | 5 | 0.9097 | 0.9381 | 0.8217 | 0.9081 | 0.9114 |
| 7.5x | 7 | 0.9248 | 0.9458 | 0.8591 | 0.9247 | 0.9249 |
| 10x | 9 | 0.9311 | 0.9479 | 0.8779 | 0.9300 | 0.9322 |
| 14.6x | 14 | 0.9362 | 0.9499 | 0.8928 | 0.9335 | 0.9389 |

**Low coverage costs a diploid contig recall and a haploid contig precision.** chr20's precision
barely moves across the whole range (0.9693 to 0.9801) while recall runs 0.8414 to 0.9495. chrX
loses both, and its precision falls to 0.8863. This is the same asymmetry the chrX investigation
found at full depth (see `wgs-results.md`): under ploidy 1 a balanced pileup has no genotype that
explains it, so the model must pick one allele and produces a coin-flip call. Less depth makes
balanced pileups commoner, so haploid calling converts missing evidence into false positives where
diploid calling converts it into missing calls.

## 2. The GQ scaling law, and a correction worth keeping

Compared **across arms**, median het GQ per unit depth looks strongly superlinear on chr20 -- 2.00
at 5x rising to 3.79 at 30x, nearly a doubling. That reading is wrong, and the way it is wrong is
instructive: the median at each coverage is taken over a *different population*, because which
sites get called het changes with depth.

Paired on **identical sites** -- which the nested subsampling makes possible, the same site seen
with fewer reads rather than a different site -- the ratio is nearly flat:

| chr20, paired | 5x | 10x | 15x | 20x | 25x | 30x |
|---|---|---|---|---|---|---|
| median GQI | 16 | 30 | 50 | 70 | 91 | 113 |
| GQI / DP | 3.20 | 3.00 | 3.33 | 3.50 | 3.64 | 3.77 |

So the likelihood gap is close to linear in depth, with about 18% residual drift rather than 90%.
The scaling that is not an artifact is the **256 clamp**: it censors 0.2% of diploid calls at 30x
but **23.3% of haploid calls** at full depth. Anything normalising GQ must therefore be computed
inside the caller, on the uncensored gap, not derived from the VCF.

## 3. GQ is miscalibrated, and not because the callsets differ

Observed precision at a claimed GQ, chr20:

| cov | GQ 0-5 | 5-10 | 10-20 | 20-40 | 40-80 |
|---|---|---|---|---|---|
| 5x | 0.901 | 0.968 | 0.983 | 0.993 | 0.995 |
| 15x | 0.820 | 0.877 | 0.943 | 0.986 | 0.995 |
| 30x | 0.730 | 0.775 | 0.835 | 0.909 | 0.984 |

**The same GQ is more reliable at low coverage than at high**, which inverts the intuition. It
makes sense on reflection: at 5x a low GQ means under-powered but usually right, while at 30x a
call still marginal after 29 reads means the evidence actively conflicts -- the collapsed-paralog
signature. The two encode different kinds of uncertainty under one number.

This is not a population artifact. Restricted to the 91,106 sites called in *all six* arms the
fan-out survives almost unchanged (GQ 0-5: 0.912 at 5x against 0.755 at 30x).

The same table on chrX fans out further and sits far lower: GQ 0-5 runs 0.253 at 2.5x down to
0.110 at 14.6x. A low-GQ haploid call is almost always wrong at any coverage, where the same GQ on
a diploid contig is right 73-90% of the time.

## 4. The negative result: no simple depth rescale works

Mean precision spread across arms, lower being better calibrated:

| series | raw GQ | GQ / DP | GQ / sqrt(DP) |
|---|---|---|---|
| chr20 (diploid) | 0.101 | **0.050** | 0.083 |
| chrX (haploid) | 0.150 | 0.161 | **0.071** |
| **POOLED across ploidies** | **0.348** | 0.496 | 0.423 |

`GQ/DP` halves the spread on the diploid series and is the obvious answer if that is all you look
at. It is slightly worse on the haploid series, and **substantially worse pooled** -- it removes
the depth axis, leaves the larger ploidy axis untouched, and compresses the score range so the
ploidy gap does more damage. At a matched `GQ/DP` bucket the two ploidies sit about 0.6 apart in
precision (chr20 0.75-0.95, chrX 0.14-0.25).

Validating on chr20 alone would have shipped a field that degrades what it exists to fix. **The
pooled row is the acceptance test for any candidate**; nothing that fails it is a fix, however good
it looks on one contig.

Why ploidy dominates: at ploidy 1 the runner-up genotype is a different allele outright, so every
read discriminates fully. At ploidy 2 a het's runner-up differs on a single strand, so a read
discriminates roughly half as much. The per-read gap scale is a function of ploidy, and `1/DP`
cannot see it.

## What this changes

- **Stage 1** cannot be a rescale of GQ by depth. The remaining principled candidate is to divide
  the observed gap by the gap *achievable* at that site -- what a noise-free pileup would give
  under this site's lambda, ploidy and allele set -- which is ploidy-aware by construction. It
  needs the per-read likelihood matrix, so design it offline from `vg call --dump-likelihoods`
  before writing any C++, the way the depth-implausibility discount was designed.
- **Stage 3** must gate on that field and never on raw GQ. A gate tuned at 30x diploid would
  discard, at 5x, a population of calls that is 90% correct.
- **Stage 2** gains support: the residual low-end spread that survives normalisation is the
  conflicting-evidence population, which is what the allele-balance and `ploidy_conflict` signal
  targets directly rather than by rescaling.
- **Stage 5** should re-tune `--linkage-weight` per coverage as planned; nothing here contradicts
  that, and it is untested so far.
