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

- **The Poisson baseline is winning while carrying a known bug.** `depth_err` at
  `snarl_caller.cpp:602` is always 0.0 or 1.0 rather than the real standard error. The
  `poisson-depthfix` arm needs a separately patched vg build and has not been run, so the gap may
  widen, not close.
- **4x is structurally unfavourable to this model.** It is depth-agnostic by design and cannot use a
  coverage anomaly as evidence, whereas the Poisson model's depth prior is doing real work at low
  coverage. The discriminating regime and the model's weakest regime are the same regime, which is
  awkward for evaluation and worth designing around.
- **A single seed said the opposite.** At seed 31 alone, read-likelihood looked marginally ahead on
  overall F1. Five seeds reversed it. Do not report single-replicate differences.

### Open questions this raises

1. Is the SNV recall deficit (-1.7 points, 5/5 seeds) a modelling limit or a bug? Given the two scoring
   bugs already found by testing, a bug should not be assumed away.
2. Where *should* this model win? Its claimed advantages are per-base quality weighting and allele
   balance without a `het_bias` knob. Neither is stressed by a low-error simulation at uniform depth.
   Higher error rates and multi-allelic sites are the untested regimes.
3. Does the insertion advantage survive on real data, and at higher depth where everything saturates?
