# Findings

Tier-0 only. Simulated reads mapped back to the graph they were simulated from, so absolute numbers
are optimistic; these are caller-vs-caller comparisons. See [simulation.md](simulation.md).

## Only low depth discriminates

| Config | Result |
|---|---|
| 20 kb, 20x, 150 bp | every arm F1 = 1.0000 |
| 400 kb, 20x, 150 bp | poisson recall 0.9996, readlik 0.9985 — saturated, no signal |
| 400 kb, **4x, 100 bp** | separates the callers; this is the useful regime |

A configuration where every arm scores near-perfectly is a statement about the configuration.

## 400 kb, 4x, 100 bp, 5 replicate seeds

GT recall, mean over seeds, with the count of seeds where read-likelihood beat Poisson:

| Type | poisson | readlik | readlik wins |
|---|---|---|---|
| ALL | **0.9096** | 0.8991 | 0/5 |
| SNV | **0.9142** | 0.8975 | 0/5 |
| Insertion | 0.8613 | **0.9210** | **5/5** |
| Deletion | **0.9183** | 0.8929 | 0/5 |

BASEPAIR (sequence-level, all types): poisson recall 0.9372 / precision 0.9725 / F1 0.9534 against
readlik 0.8886 / 0.8426 / 0.8610. Read-likelihood emits materially more wrong sequence, and with much
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

Which traversal finder each arm uses matters for interpreting the above, and is worth stating because
a GBZ input does *not* by itself change it — only `-g`/`-z` does:

| Arm | Traversal finder |
|---|---|
| `poisson`, `readlik`, `readlik-nomismap` | `FlowTraversalFinder` — Yen's k-widest paths, node/edge weights from the pack file |
| `readlik-gbwt-nopack` | `GBWTTraversalFinder` — haplotypes recorded in the GBZ |

So three of the four arms share support-driven enumeration, which is good experimental hygiene by
accident: the poisson-vs-readlik comparison holds enumeration constant and varies only the genotyping
model.

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

### Open questions this raises

1. Is the SNV recall deficit (-1.7 points, 5/5 seeds) a modelling limit or a bug? Given the two scoring
   bugs already found by testing, a bug should not be assumed away.
2. Where *should* this model win? Its claimed advantages are per-base quality weighting and allele
   balance without a `het_bias` knob. Neither is stressed by a low-error simulation at uniform depth.
   Higher error rates and multi-allelic sites are the untested regimes.
3. Does the insertion advantage survive on real data, and at higher depth where everything saturates?
