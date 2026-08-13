# Findings

Tier-0 only. Simulated reads mapped back to the graph they were simulated from, so absolute numbers
are optimistic; these are caller-vs-caller comparisons. See [simulation.md](simulation.md).

> **Superseded for anything but method.** Tier 2 runs real HG002 reads against the GIAB draft
> benchmark on two chromosomes and two graphs, and reverses this page's headline conclusion. Go to
> [tier2-chr20-results.md](tier2-chr20-results.md) and [tier2-chr6-results.md](tier2-chr6-results.md)
> for the accuracy numbers, [tier2-chr20-graph-comparison.md](tier2-chr20-graph-comparison.md) and
> [tier2-chr6-graph-comparison.md](tier2-chr6-graph-comparison.md) for the graph comparison, and
> [tier2-quality-signals.md](tier2-quality-signals.md) for how calls are ranked. This page is kept
> because the *method* lessons below still hold — the saturation table, the five-seed rule, and the
> two retracted claims.

## Only low depth discriminates

| Config | Result |
|---|---|
| 20 kb, 20x, 150 bp | every arm F1 = 1.0000 |
| 400 kb, 20x, 150 bp | poisson recall 0.9996, readlik-support 0.9985 — saturated, no signal |
| 400 kb, **4x, 100 bp** | separates the callers; this is the useful regime |

A configuration where every arm scores near-perfectly is a statement about the configuration.

## 400 kb, 4x, 100 bp, 5 replicate seeds

GT recall, mean over seeds, with the count of seeds where read-likelihood beat Poisson:

| Type | poisson | readlik-support | readlik-support wins |
|---|---|---|---|
| ALL | **0.9096** | 0.8991 | 0/5 |
| SNV | **0.9142** | 0.8975 | 0/5 |
| Insertion | 0.8613 | **0.9210** | **5/5** |
| Deletion | **0.9183** | 0.8929 | 0/5 |

BASEPAIR (sequence-level, all types): poisson recall 0.9372 / precision 0.9725 / F1 0.9534 against
readlik-support 0.8886 / 0.8426 / 0.8610. Read-likelihood emits materially more wrong sequence, and with much
higher variance across seeds (recall sd 0.11 vs 0.06).

### What this says

**As it stands the read-likelihood caller does not beat the Poisson caller** in the only regime where
the two can be told apart. It loses on SNVs and deletions and on sequence-level precision, and it wins
on **insertions, consistently in 5/5 seeds by about 6 points of recall** — the one clear directional
strength so far.

Three things to hold alongside that:

- ~~The Poisson baseline is winning while carrying a known bug.~~ **Tested, and this was wrong.**
  `depth_err` at `snarl_caller.cpp:602` really is malformed, but it is *inert*: its only consumer
  inside `genotype_likelihood` is commented out deliberately, and it never reaches the VCF.
  Patching it gives byte-identical calls across three 400 kb replicates. The Poisson baseline is
  therefore a fair comparison exactly as shipped, and the read-likelihood deficit above is **not**
  explained away by a bug in what it is being compared against.
- **4x is structurally unfavourable to this model.** It is depth-agnostic by design and cannot use a
  coverage anomaly as evidence, whereas the Poisson model's depth prior is doing real work at low
  coverage. The discriminating regime and the model's weakest regime are the same regime, which is
  awkward for evaluation and worth designing around.
- **A single seed said the opposite.** At seed 31 alone, read-likelihood looked marginally ahead on
  overall F1. Five seeds reversed it. Do not report single-replicate differences.

## Allele enumeration is not the bottleneck (in tier 0)

Which traversal finder each arm uses matters for interpreting the above, and no longer follows from the
flags alone. A GBZ input used to change nothing by itself — only `-g`/`-z` did — but `--read-likelihood`
now enumerates from the GBZ's haplotype panel by default wherever that panel holds at least two
haplotypes, so under that caller the *plain* invocation is the one using `GBWTTraversalFinder`:

| Arm | Traversal finder |
|---|---|
| `poisson`, `readlik-support` | `FlowTraversalFinder` — Yen's k-widest paths, node/edge weights from the pack file |
| `readlik`, `readlik-nomismap` | `GBWTTraversalFinder` — haplotypes recorded in the GBZ |

`readlik-support` exists to keep the caller comparison honest. It asks for support enumeration
explicitly (`--enumerate-support`), so poisson-vs-readlik-support still holds enumeration constant and
varies only the genotyping model. Comparing `poisson` against `readlik` varies both at once: that is a
comparison of shipped defaults, not of models, and the two questions have different answers.

The fourth arm turns out to be an unintentionally strong diagnostic. In tier 0 the GBZ is built by
`autoindex` from the truth VCF, so **its haplotypes are the true haplotypes** — that arm is enumerating
alleles from the answer. And it recovers exactly nothing extra: TP, FP and FN are identical to the
flow-enumeration arm in **5/5 seeds** (e.g. 828/77/92 at seed 31). The per-record differences are almost
entirely `1/0` versus `0/1`, i.e. allele order within an unphased genotype, which is the same call.

**Conclusion: the ~10% of truth missed at 4x is not lost in candidate generation.** It is lost in the
read evidence or in the genotyping. That is where to look for the SNV deficit.

**Limit on that conclusion:** in tier 0 the graph is constructed from the truth VCF, so it already
contains every true allele and the flow finder is choosing among candidates that include the answer.
Enumeration could easily matter on real data against an imperfect pangenome. This says enumeration is
not the bottleneck *here*, not in general.

### Open questions this raises — and what tier 2 answered

1. ~~Is the SNV recall deficit (-1.7 points, 5/5 seeds) a modelling limit or a bug?~~ **Neither, in the
   end.** It did not reproduce on real data: on chr20 and chr6, at both graph richnesses, the
   read-likelihood caller leads the Poisson caller on GT F1 for every class including SNVs. The tier-0
   deficit was a property of 4x simulated coverage, which is the one regime a depth-agnostic model is
   structurally worst in — as the third bullet above already suspected.
2. **Answered: multi-allelic sites, which is exactly the untested regime named here.** The
   read-likelihood caller's margin over the Poisson caller *widens* on a 34-haplotype graph, where
   more haplotypes mean more true alleles and more wrong ones and the question becomes whether the
   genotyper can separate them read by read. See the hap32 pages.
3. **Yes, and it was partly an artefact.** The insertion advantage survives, but the apparent
   insertion *BASEPAIR* deficit that showed up alongside it turned out to be benchmark scope: the
   `smvar` truth set holds no record >=50 bp, so a correct large insertion inside its confident
   region scores FP on every base. Size-matching both callers collapses the gap from 0.139 to 0.008.

Two further method lessons from tier 2, in the same spirit as the ones above:

- **A default that is inert on one graph is not thereby harmless.** `--mismap-max` at 0.1 looked
  irrelevant on the 4-haplotype graph and was actively wrong on the 34-haplotype one.
- **A sweep that sets a default must be scored on every benchmark the project runs.** The mismapping
  floor was first tuned on small variants alone and the setting it chose cost about 0.01 of SV F1.
