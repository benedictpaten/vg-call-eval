# vg-call-eval

Concordance and performance testing for `vg call`, kept **outside the vg tree** so that none of its
dependencies (aardvark, truvari, a Python analysis stack) land in vg's build.

It exists to answer one question that vg's own test suite cannot: **is a change to `vg call` actually
more accurate?** vg's in-tree harness has exactly one truth-based concordance assertion, and it covers
the `-v` re-genotyping path — the default de novo path has never been measured against truth.

Implements stages 3b, 4 and 4b of the read-likelihood design.

## Where the numbers are

| page | what it covers |
|---|---|
| [docs/tier2-chr20-results.md](docs/tier2-chr20-results.md), [docs/tier2-chr6-results.md](docs/tier2-chr6-results.md) | the full five-arm accuracy tables per chromosome, small variants and SVs |
| [docs/tier2-chr20-hap32.md](docs/tier2-chr20-hap32.md), [docs/tier2-chr6-hap32.md](docs/tier2-chr6-hap32.md) | 4-haplotype against 34-haplotype graph, the same reads remapped |
| [docs/tier2-quality-signals.md](docs/tier2-quality-signals.md) | how calls are *ranked*: `AD`, `BL`, `GQI`, the explained-share discount in `GQ`, and the filters that turned out not to help |
| [docs/tier2-parameters.md](docs/tier2-parameters.md) | the caller's tuned parameters re-swept after the mixture change: why `--mismap-max` moved to 0.7, why `--mismap-min` stays at 0.02, and why `--read-weight` was removed |
| [docs/tier2-sv-errors.md](docs/tier2-sv-errors.md) | what the SV errors *are*, per record: why the read model trails on structural variants (heterozygous deletions, and nothing else), what the 34-haplotype precision loss is made of, and how much of "false positive" is the metric rather than the caller |
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
```

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
quality-signal analysis behind the `GQ` change, and the per-record SV error atlas.

Not yet built: tier 1 (vg's HGSVC fixture), the `read_weight` calibration fit (stage 4b), and a
size-conditional depth or `best_ln` term — blocked on the sign reversal documented in
[docs/tier2-quality-signals.md](docs/tier2-quality-signals.md).

Open and specific: a **depth-plausibility term**. The read-likelihood model computes P(reads | genotype)
conditioned on the reads it is given and never asks whether that many reads should be there, which is why
collapsed-repeat pile-ups survive it and why the Poisson caller still leads on heterozygous deletions above
1 kb (0.79-0.84 against 0.44). The length-weighted mixture fixed the *relative* weight between a genotype's
haplotypes ([docs/tier2-sv-errors.md](docs/tier2-sv-errors.md)) and is now the default; absolute depth is the
remaining half. It is blocked on the same sign reversal as above — depth discriminates in opposite directions
for small variants and SVs — so it needs conditioning on called-allele size.
