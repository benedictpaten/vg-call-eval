# Calling across coverage and ploidy: the titration

> **Stale for the caller as of decide-then-render (2026-08).** Every vg figure below was measured
> before genotypes were settled ahead of record construction. That change moved the whole-genome
> autosomal numbers -- ALL F1 0.9703 -> 0.9729, Indel 0.9195 -> 0.9272, SV >=50 bp 0.5488 -> 0.5596,
> with both precision and recall improving in every class -- so the figures here understate the
> current caller by roughly that much, and any *analysis* built on which calls were wrong may have
> picked a different population. Not re-run: these arms use their own reads, truth sets and graphs, and
> re-measuring them is hours of runs that were not spent. Current numbers:
> [wgs-results.md](wgs-results.md), [pangenie-comparison.md](pangenie-comparison.md).
> What changed and what is still open: `planning/decide-then-render.md`.


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
- **Stage 5** re-tuned the three model parameters across coverage. See below.

## Stage 5: are the defaults still right at other coverages?

Each of `--linkage-weight`, `--linkage-freq-prior` and `--depth-term` was swept around its shipped
default at four corners -- low and high coverage at each ploidy -- with everything else left alone.
A coordinate sweep rather than a grid, because the defaults were fitted one at a time and the
question is whether each still holds elsewhere.

**Read the gain column, not the argmax.** An optimum a few ten-thousandths better on one arm is
noise plus a different arm; changing a default on that basis would invalidate every published
number for nothing.

### `--linkage-weight` (default 2)

| arm | 0 | 1 | 2 | 4 | 8 | best | gain |
|---|---|---|---|---|---|---|---|
| chr20 5x | 0.8546 | 0.8998 | **0.9008** | 0.9002 | 0.8987 | 2 | +0.0000 |
| chr20 30x | 0.9546 | 0.9643 | 0.9645 | **0.9646** | 0.9637 | 4 | +0.0001 |
| chrX 2.5x | 0.8471 | **0.8637** | 0.8609 | 0.8595 | 0.8593 | 1 | +0.0028 |
| chrX 14.6x | 0.9362 | **0.9438** | 0.9421 | 0.9417 | 0.9418 | 1 | +0.0018 |

**The prediction that optimal weight rises as coverage falls is refuted.** It is flat at 2 on the
diploid arms and sits at 1 on the haploid ones, with gains over the default between 0.0000 and
0.0028. The default stands.

The valuable column is `0`: linkage is worth **+0.046** at 5x diploid and **+0.017** at 2.5x
haploid, so the layer earns its place at low coverage even though its weight does not need tuning.

This sweep is also what exposed the [haploid linkage bug](#a-bug-this-sweep-found).

### `--linkage-freq-prior` (default 5)

| arm | 0 | 3 | 5 | 8 | best | gain |
|---|---|---|---|---|---|---|
| chr20 5x | 0.8960 | **0.9021** | 0.9008 | 0.8981 | 3 | +0.0012 |
| chr20 30x | 0.9593 | 0.9631 | 0.9645 | **0.9649** | 8 | +0.0003 |
| chrX 2.5x | 0.8561 | **0.8613** | 0.8609 | 0.8601 | 3 | +0.0005 |
| chrX 14.6x | 0.9385 | 0.9406 | 0.9421 | **0.9445** | 8 | +0.0024 |

**This is the one real coverage trend**, and it is consistent: low coverage prefers 3, high
coverage prefers 8, on *both* ploidies. Replication across ploidy is what makes it credible rather
than four independent argmaxes.

The direction is the opposite of the intuition behind the linkage-weight prediction. A stronger
frequency prior pulls calls toward alleles the panel carries often; where the reads are weak that
overrides them and costs rare true variants, so the prior should be *weaker* at low coverage, not
stronger.

The effect is nonetheless small -- at most 0.0024, and the default of 5 is never more than that
from the best value on any arm.

### `--depth-term` (default 0.1)

| arm | 0 | 0.05 | 0.1 | 0.2 | best | gain |
|---|---|---|---|---|---|---|
| chr20 5x | 0.9007 | 0.9008 | 0.9008 | 0.9008 | 0.1 | +0.0000 |
| chr20 30x | 0.9646 | 0.9646 | 0.9645 | 0.9645 | 0 | +0.0000 |
| chrX 2.5x | **0.8614** | 0.8609 | 0.8609 | 0.8608 | 0 | +0.0005 |
| chrX 14.6x | 0.9421 | 0.9420 | 0.9421 | 0.9421 | 0 | +0.0001 |

Flat to four decimals across the whole range at every coverage. On small variants the term does
essentially nothing, which is expected -- it exists for large deletions, where absent reads are the
evidence, and those are a small share of an F1 dominated by SNVs. This sweep says nothing about
that case and should not be read as doing so.

### The decision: neither auto-scaling nor a table

The plan left open whether to auto-scale parameters from the measured coverage or to ship a
documented table. **Neither is justified.** Two of the three parameters have optima that do not
move with coverage at all, and the third moves by at most 0.0024 F1 -- less than the spread between
adjacent values at a single coverage.

Auto-scaling would make a run's behaviour depend on its own data in a way that is hard to explain
and hard to reproduce, and it would buy under 0.003. A table would ask users to look something up
for the same. The defaults stay, and what is documented instead is the shape: `--linkage-freq-prior`
3 is worth about 0.001 below 10x and 8 is worth about 0.002 above 15x, for anyone who wants it.

### A bug this sweep found

chrX returned *identical* F1 at every `--linkage-weight`, and the VCFs were byte-identical, while
chr20's moved. But `--progress` reported 8,945 genotypes changed.

`apply_linkage_change` built the genotype it expected as `"i/j"`, which a haploid record's bare
allele can never match, so its guard rejected every haploid change. The linkage layer had been
doing the work on haploid contigs and discarding all of it -- and since phasing and the mosaic are
built from the post-linkage genotypes, the mosaic described genotypes the VCF did not contain.

Fixed in vg at the time. `apply_linkage_change` itself no longer exists -- the record is now built
from the settled genotype, so there is no line to patch and no guard to get wrong; see
`planning/decide-then-render.md`. The account above is kept because the failure mode it names,
a silent no-op that the progress counter reported as work done, is not specific to that function.
**Haploid linkage is worth +0.017 F1 at 2.5x and
+0.008 at 14.6x**, none of which chrY or non-pseudoautosomal chrX was receiving.

Three points of method, all learned the hard way in this stage:

- A resume marker that skips completed work is wrong for a sweep. The fix landed mid-sweep and the
  marker kept the pre-fix chrX arms, which then scored as though they were the fixed caller and
  were identifiable only by their file timestamps. `sweep_params.sh` now reuses a result only if it
  is **newer than the vg binary** -- a sweep measures a build, so a result from an older build is
  an answer to a different question.
- The same applies to the scoring cache. `bench_coverage.py` and `sweep_report.py` reuse a
  `.renamed.vcf.gz` if one exists, which will happily re-score a deleted-and-regenerated VCF's
  stale twin.
- In zsh, one non-matching glob aborts the whole `rm`, so a cleanup command that looks like it ran
  may have done nothing. Use `find -delete` and check.
