# How tier-0 tests are simulated

This is the part of the harness that most needs explaining, because everything downstream — the
concordance numbers, the GQ calibration — inherits whatever assumptions are baked in here.

## The problem simulation solves

To measure a caller you need to know the right answer. Real benchmark sets give you that only
approximately, only in "confident" regions, and only for one sample. Simulation gives it exactly, and
gives one thing no real dataset can: **the true origin of every individual read**, which is what the
MAPQ calibration in `mapq.py` needs.

The cost is realism. Reads simulated from a graph are easier to map back to that graph than real reads
are, so absolute precision and recall from tier 0 will be optimistic. **Tier 0 is for comparing callers
to each other and for calibration, never for absolute numbers.** That is a hard rule, and the generated
reports repeat it in their header so a number cannot be quoted out of context.

## The design in one line

Generate a diploid sample with known phased genotypes, simulate reads from its two haplotypes, **map
those reads back with a real mapper**, and call. The truth VCF is the genotypes we generated.

## Step by step

### 1. Generate a reference and a phased truth VCF

`simulate.py` writes a random reference contig and a phased VCF for a single diploid sample, from a
seed. Everything is parameterised so a config can target a regime:

| Parameter | Controls |
|---|---|
| `ref_length` | contig size |
| `snp_rate`, `indel_rate`, `sv_rate` | variant density per class |
| `indel_size_range`, `sv_size_range` | size distribution within a class |
| `het_fraction` | proportion of heterozygous sites |
| `seed` | reproducibility |

Generating the variants ourselves rather than reusing a fixture is deliberate. It is the only way to get
**enough variants per size class to bin a calibration curve**, and the only way to dial SV density up to
probe the regime this caller is supposed to be good at. vg's own `test/small/x.fa` (1001 bp, 75
variants) is available as an extra-small smoke fixture, but 75 variants cannot support a reliability
curve.

Two constraints the generator enforces, both of which produce silently wrong truth if violated:

- **Variants may not overlap or abut.** Two variants close enough to interact make the "truth genotype"
  ambiguous — the pair could be represented several ways, and the comparison engine's haplotype replay
  may legitimately resolve it differently from how we wrote it. A minimum inter-variant spacing is
  enforced (default 20 bp, and at least the size of any indel involved).
- **The reference base must match the contig.** Trivial, and trivially easy to get wrong when
  generating indels; asserted rather than assumed.

### 2. Build the graph

```
vg construct -r ref.fa -v truth.vcf.gz -a     # alt paths, for -v regenotyping arms
vg autoindex -r ref.fa -v truth.vcf.gz -w giraffe -p idx
```

`autoindex` produces the GBZ, which carries the sample's two phased haplotypes as paths. Those paths
are what we simulate from, and they are also what `-g`/`-z` enumeration will later use as alleles — so
tier 0 exercises the pack-free GBWT path honestly.

### 3. Simulate reads from each haplotype separately

```
vg sim -x idx.giraffe.gbz -P <hap0 path> -n N/2 -l L -a -s SEED   > hap0.gam
vg sim -x idx.giraffe.gbz -P <hap1 path> -n N/2 -l L -a -s SEED+1 > hap1.gam
```

Simulating each haplotype separately, rather than sampling the graph as a whole, is what makes the read
set a genuine diploid mixture in known proportion. `-a` records each read's true alignment, which is
kept for two purposes below.

Read length, count (via target depth) and error rate are config parameters. Depth is specified rather
than read count, and converted: `N = depth * ref_length / read_length`.

### 4. Map the reads back — do not call on the simulated alignments

This step is the one most easily skipped and most important to keep.

```
vg view -X sim.gam > reads.fq          # simulated reads, as FASTQ
vg giraffe -Z idx.giraffe.gbz -f reads.fq > mapped.gam
```

Calling on `vg sim -a`'s output directly would mean genotyping perfect alignments with perfect MAPQ.
Every mapping error would be assumed away — including exactly the errors the mismapping term in the
model exists to handle, and the strand handling that a previous bug got wrong. The harness therefore
throws away the simulated alignments for calling purposes and re-maps the reads with a real mapper.

The simulated alignments are kept, though, and used for two things nothing else can provide:

- **`mapq.py`**: join mapped reads to their true origin by read name, and measure the real relationship
  between reported MAPQ and actual mismapping rate. That is the measurement behind the design doc's
  proposed MAPQ→`e_r` recalibration table, and the evidence for whether the `e_max` clamp is sensible.
- **A sanity gate**: if the mapper places a wildly implausible fraction of reads wrongly, the
  configuration is broken (wrong graph, wrong reads) and the run should fail loudly rather than produce
  a bad concordance number.

### 5. Call

`vg pack` then `vg call`, once per arm of the run matrix. The five arms:

| Arm | Purpose |
|---|---|
| `poisson` | the current default, as shipped |
| `poisson-depthfix` | with the `depth_err` one-liner patched — quantifies how much that bug distorts the baseline |
| `readlik` | the new caller as shipped: panel enumeration from the GBZ, no pack file |
| `readlik-nomismap` | `--no-mismap-term`, to measure what the mismapping term contributes |
| `readlik-support` | `--enumerate-support`, which puts enumeration back on the pack so the comparison against `poisson` varies only the genotyping model |

`poisson-depthfix` needs a separately built vg, which is why the harness takes a **binary path per arm**
rather than a single global one.

### 6. Compare

The called VCF is normalised (`bcftools norm -f ref.fa -m -any`), the sample renamed to match truth, and
compared against the truth VCF over a BED covering the whole contig — for tier 0 everything is confident
by construction, which is one more reason it is a good place to develop.

## What tier 0 cannot tell you

Stated here so it is stated somewhere:

- **Absolute precision/recall.** Reads come from the graph, so mapping is too easy.
- **Anything about real error profiles.** `vg sim`'s error model is not a sequencer.
- **Anything about reference bias from real data**, since the truth and the graph agree perfectly.
- **Whether calibration transfers.** A `read_weight` fitted here is fitted to this error model and this
  mapper. It must be re-fitted at tier 2 before being believed.

Tier 1 (vg's HGSVC fixture, real reads) partially covers the second and third of these, but it contains
only 8 non-reference SVs, so it is a regression smoke test and not a measurement. Tier 2 is where
absolute numbers come from, and it is deferred until the read source can handle that scale.
