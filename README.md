# vg-call-eval

Concordance and performance testing for `vg call`, kept **outside the vg tree** so that none of its
dependencies (aardvark, truvari, a Python analysis stack) land in vg's build.

It exists to answer one question that vg's own test suite cannot: **is a change to `vg call` actually
more accurate?** vg's in-tree harness has exactly one truth-based concordance assertion, and it covers
the `-v` re-genotyping path — the default de novo path has never been measured against truth.

Implements stages 3b and 4 of the read-likelihood design. Stage 4b — fitting `read_weight` —
was **cancelled rather than completed**: the parameter provably cannot change a genotype, so the
quantity it was to be fitted against cannot move. See [docs/tier2-parameters.md](docs/tier2-parameters.md).

## Where the numbers are

| page | what it covers |
|---|---|
| [docs/tier2-chr20-results.md](docs/tier2-chr20-results.md), [docs/tier2-chr6-results.md](docs/tier2-chr6-results.md) | the full five-arm accuracy tables per chromosome, small variants and SVs |
| [docs/tier2-chr20-hap32.md](docs/tier2-chr20-hap32.md), [docs/tier2-chr6-hap32.md](docs/tier2-chr6-hap32.md) | 4-haplotype against 34-haplotype graph, the same reads remapped |
| [docs/tier2-quality-signals.md](docs/tier2-quality-signals.md) | how calls are *ranked*: `AD`, `BL`, `GQI`, the explained-share discount in `GQ`, the size-gated depth discount behind it, and the filters that turned out not to help |
| [docs/tier2-parameters.md](docs/tier2-parameters.md) | the caller's tuned parameters, re-swept twice: why `--mismap-max` moved to 0.7, why `--read-weight` was removed, and why `--mismap-min` stays at 0.02 even though the depth term dissolved the trade it was balancing |
| [docs/tier2-depth-term.md](docs/tier2-depth-term.md) | the depth term, predicted offline and then built: `--depth-term` puts the read model ahead of both Poisson arms on all four datasets, and a read counts toward depth as `1 − e_r` rather than as one read |
| [docs/tier2-sv-errors.md](docs/tier2-sv-errors.md) | what the SV errors *are*, per record: why the read model trailed on structural variants (heterozygous deletions, and nothing else), why the 34-haplotype graph's SV penalty is mostly not the caller — at matched sensitivity it is 0.021 on chr6 and zero on chr20 — and how much of "false positive" is the metric rather than the caller |
| [docs/findings.md](docs/findings.md), [docs/results.md](docs/results.md) | tier 0, superseded for accuracy but kept for its method lessons |
| [docs/simulation.md](docs/simulation.md) | how tier 0 works and what it cannot tell you |

**Design and planning documents live in [planning/](planning/)** — the caller's design, the
harness plan and its full investigation log, the characterization of `vg call` as it was, and two
outbound drafts. They sit here rather than in the vg tree so that the reasoning is next to the
evaluation it came from. Start at [planning/README.md](planning/README.md).

When a number in `planning/` disagrees with one in `docs/`, **`docs/` is right**: those are
regenerated from run artefacts, the planning documents are transcribed by hand.

## The one thing to read before quoting a number

**Tier-0 numbers are optimistic and are not absolute performance.** Reads are simulated from the graph
and mapped back to that same graph, so mapping is unrealistically easy. Tier 0 exists to compare
callers *to each other*.

A worked example of why this matters: at 20 kb and 20x, **every caller scores F1 = 1.0000**. The task
is simply too easy to discriminate. At 4x with 100 bp reads the same harness separates them. If a
configuration gives everything a perfect score, that is a statement about the configuration, not the
callers.

**Tier-2 numbers are benchmark-relative.** The GIAB truth set is a *draft*, with known errors in
homozygous regions, homopolymers and tandem repeats. More sharply: its small-variant benchmark holds
**no record at all over 50 bp**, so a correct large insertion inside its confident region scores as a
false positive on every base. Anything size-restricted above 50 bp goes to truvari against the
structural benchmark; aardvark's own `Sv*` categories are scored against the small-variant truth and
should not be read as an SV result.

## Sanity controls are not optional

A comparison harness that is subtly wrong is worse than none, because it produces confident numbers.
Two controls gate everything, and they run as tests:

```bash
PYTHONPATH=src python3 -m pytest tests/ -q
```

* **Positive** — identical inputs must score 1.0.
* **Negative** — deliberately dropping and corrupting calls must be detected. If this passes at 1.0,
  the harness cannot see errors and every number it has produced is meaningless.

## Install

Needs `vg`, `bcftools`, `samtools`/`tabix`/`bgzip`, and `aardvark`. See [docs/install.md](docs/install.md)
— note aardvark ships **only an x86_64 Linux binary**, so on macOS/ARM it must be built from source
(it is pure Rust and builds cleanly).

## Run

```bash
PYTHONPATH=src python3 -m vgcalleval.cli run \
    --out /tmp/eval --vg /path/to/vg --ref-length 60000 --depth 4 --read-length 100
```

Add `--vg-depthfix /path/to/patched/vg` to include the `poisson-depthfix` arm. The harness takes a
**binary path per arm**, so two vg builds can be compared in a single matrix — which is how the
`depth_err` bug's effect on the baseline gets quantified.

### Arms

| Arm | Purpose |
|---|---|
| `poisson` | the current default, as shipped |
| `poisson-depthfix` | with the `depth_err` one-liner patched. **Verified byte-identical output** - that bug is inert, since its only consumer in the likelihood is commented out. Retained as a control: if this arm ever diverges from `poisson`, someone has re-enabled the depth term. |
| `readlik` | the read-level likelihood caller |
| `readlik-nomismap` | `--no-mismap-term`, to measure what the mismapping term contributes |
| `readlik-gbwt-nopack` | `-z` haplotype enumeration with no pack file |

## Tier 2

Real HG002 reads against the GIAB draft benchmark on CHM13v2.0, run on a 32 GB laptop. Two
chromosomes (chr20, chr6) × two graphs (4-haplotype, 34-haplotype) × five arms, scored against both
the small-variant benchmark (aardvark) and the structural one (truvari).

```bash
# one contig, one graph: subgraph, node list, reference FASTA, truth slices, reads, pack
scripts/tier2/prep_contig.sh chr6 data/…HG002.gbz work/graph.gbz.db work/reads.gaf.db work/tier2-chr6
# arms, then the SV cross-check, then the size-matched control
python3 scripts/tier2/run_arms.py    --contig chr6 --graph … --out work/tier2-chr6/results
python3 scripts/tier2/truvari_sv.py  --contig chr6 --work  work/tier2-chr6 --label chr6-4hap
python3 scripts/tier2/size_matched.py             --results work/tier2-chr6/results …
# regenerate the pages
python3 scripts/tier2/report.py         --contig chr6
python3 scripts/tier2/compare_graphs.py --contig chr6
```

### SV error forensics

Runs on the truvari output the commands above already produced, so it needs no re-calling:

```bash
python3 scripts/tier2/patch_truvari_pysam.py          # once: truvari refine is broken under pysam 0.24
python3 scripts/tier2/sv_error_atlas.py               # per-record FP/FN/TP tables -> work/sv-atlas/
python3 scripts/tier2/sv_metric_sensitivity.py --refine
python3 scripts/tier2/sv_error_report.py              # every table in docs/tier2-sv-errors.md
python3 scripts/tier2/hetdel_mechanism.py            # the heterozygous-deletion mechanism test
python3 scripts/tier2/score_vcf.py --vcf … --label … # score any experimental VCF on BOTH benchmarks
python3 scripts/tier2/depth_term_offline.py          # Stage 0 depth-term prediction
python3 scripts/tier2/depth_count_runs.py \
    --datasets chr20-4hap chr6-4hap             # depth counted raw vs as 1 - e_r
python3 scripts/tier2/depth_gq.py \
    --tag dgrid-w0.1-f0.02-c0.7                 # DR as a GQ discount, 8 cells
python3 scripts/tier2/depth_grid.py \
    --datasets chr20-4hap                       # depth weight x mismapping floor
python3 scripts/tier2/hap32_precision.py            # why 34-hap emits more false SVs
python3 scripts/tier2/fn_decomposition.py           # misses that were called, spelled smaller
python3 scripts/tier2/sv_metric_sensitivity.py \
    --arms readlik-z --refine --refdists 500        # harmonised representation
```

A caller-side change can be put through the whole five-arm matrix without editing the arm list:

```bash
READLIK_EXTRA="--depth-term 0.1" CANARY=1 JOBS=2 scripts/tier2/refresh_all.sh
```

### What a refresh costs, and where to spend less

Measured on the depth-term refresh, and it is worth knowing before optimising the wrong
thing: **95% of the time is `vg call`.** All the aardvark, truvari and size-matched scoring
for twenty arms comes to **four minutes**, so parallelising the scoring buys nothing.

| | time | share |
|---|---|---|
| `vg call`, 20 arms | 1.36 h | 95% |
| all scoring | 4 min | 4% |
| **total** | **1.43 h** | |

Within the calling, most of it is arms that cannot answer the question:

| | time |
|---|---|
| the two Poisson arms — untouched by any read-likelihood change | **36 min (43%)** |
| `readlik` + `readlik-nomismap` — diagnostics | 32 min |
| **`readlik-z`, all four datasets — the arm that decides** | **14 min** |

So there are two tiers, and most of this project's work happened at the fast one:

* **Fast (~14 min):** `readlik-z` on all four datasets, both benchmarks. Use `param_sweep.py`
  or `depth_grid.py`. Enough for a go/no-go on a caller change.
* **Full matrix (~30 min with `CANARY=1 JOBS=2`, 86 min without):** all five arms, one build,
  every page regenerated. Needed before changing a default or publishing numbers.

`CANARY=1` re-runs one cheap Poisson arm and compares it byte for byte against the cached
copy, reusing the Poisson rows only if they are identical. That is **stronger** than
re-running them, not weaker: a mismatch means something touched shared code, which blind
re-running would have absorbed silently. `JOBS=2` runs two datasets at once — peak RSS is
6–8 GB per call, so two fit in 32 GB.

Sweep runners key their cache on `(binary, dataset, flags)` rather than on the tag, so the
same configuration under a different name is free. Tags are hand-written and collide: three
separate tags in one session here were the same experiment.

Parameter sweeps, searched on chr20 and validated on chr6 so the validation set stays held out:

```bash
python3 scripts/tier2/param_sweep.py --param mismap-max --values 0.5 0.7 0.9 0.99 \
    --param2 mismap-min --values2 0.01 0.02 0.05 --datasets chr20-34hap chr20-4hap
python3 scripts/tier2/param_sweep.py --param mismap-min --values 0.01 0.02 0.05 0.10 0.20 \
    --datasets chr20-34hap chr20-4hap
```

Every point is scored on **both** benchmarks plus the heterozygous SV class breakdown and the
genotype mix, and the surface is printed whole. There is deliberately no single objective: a
setting that buys SV F1 with small-variant genotype F1 is a judgement about what the caller is
for, and picking one number to maximise is how the mismapping cap ended up at 0.1 originally.

The mechanism test needs likelihood matrices and a traversal enumeration, both restricted to
large snarls so they are cheap:

```bash
vg call work/tier2-chr6/chr6_0_chr6.gbz -p CHM13#0#chr6 --read-likelihood -z -c 800 \
    --dump-likelihoods work/sv-atlas/chr6-large.dump.tsv --gaf-base … --gbz-base … > /dev/null
vg call work/tier2-chr6/chr6_0_chr6.gbz -p CHM13#0#chr6 -z -T -k work/tier2-chr6/chr6.pack -c 800 \
    > work/sv-atlas/chr6-trav.gaf
```

`prep_contig.sh` extracts the reference FASTA from *each* graph and stops if two graphs for the same
contig disagree, so a cross-graph comparison can never be a coordinate mismatch dressed up as an
accuracy difference.

## Status

Working: tier-0 simulation, tier-2 real data on two chromosomes and two graphs, the caller matrix,
aardvark and truvari comparison, size-matched controls, sanity controls, per-arm timing, the
quality-signal analysis behind the `GQ` change, the per-record SV error atlas, the depth-term grid
search, and the representation analysis behind the SV numbers.

Not yet built: tier 1 (vg's HGSVC fixture), and a `best_ln` filter — blocked on the sign reversal
documented in [docs/tier2-quality-signals.md](docs/tier2-quality-signals.md). Stage 4b's `read_weight`
calibration fit is not pending but **cancelled**: the parameter provably cannot change a genotype, and
has been removed from vg ([docs/tier2-parameters.md](docs/tier2-parameters.md)).

Built since, and no longer open: the **depth-plausibility term**. The read-likelihood model computed
P(reads | genotype) conditioned on the reads it was given and never asked whether that many reads should be
there, which is why collapsed-repeat pile-ups survived it and why the Poisson caller led on heterozygous
deletions above 1 kb. The length-weighted mixture fixed the *relative* weight between a genotype's haplotypes
([docs/tier2-sv-errors.md](docs/tier2-sv-errors.md)) and is the default; `--depth-term` supplies the absolute
half and puts the read model ahead of both Poisson arms on all four datasets. It is **on by default** at
`--depth-term 0.1` — see [docs/tier2-depth-term.md](docs/tier2-depth-term.md).

Still open: **`DR` in the quality field**. The term detects collapsed repeats emphatically and still cannot
outvote the read evidence at them. `DR` is emitted whether or not the term is armed, and counting each read
as `1 − e_r` rather than as one read raises its power to rank false positives above true ones from 0.51–0.55
to 0.62–0.64 on all four datasets — which is the signal a depth-implausibility discount would be built on,
as a sibling of the explained-share discount already in `GQ`.
