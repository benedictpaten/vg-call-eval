# Running a whole genome, and how fast it goes

`vg call --read-likelihood` on a 34-haplotype HPRC graph, 30x reads, on a 10-core laptop with
32 GB. One contig per invocation, several contigs at once, packed under a memory budget.

## Decide-then-render: read I/O restored, CPU up

The current arm settles every genotype before building its record. The cost that mattered was read
I/O, because the arm before it bought the same coherence guarantee by re-reading the contig once per
generation:

| | inline | post-linkage descent | decide-then-render |
|---|---|---|---|
| reads fetched | 609.1 M | 903.3 M (+48.8%) | **609.7 M (1.001x)** |
| CPU user | 6.26 h | 11.45 h | **9.96 h (+59%)** |
| wall clock | 3.20 h | 4.89 h | 6.03 h |
| linkage pass, serial | 1,064 s | 950 s | 1,762 s |
| peak RSS, worst contig | 4.30 GB (chr1) | — | 6.05 GB (chr1) |

**The read penalty is gone to within 0.1%**, which was the entire point of the single sweep. It costs
+59% CPU against inline, and is 13% *cheaper* on CPU than the deferred arm while also dropping that
arm's read penalty -- so it strictly dominates the arm it replaces.

Read CPU, not wall clock, for the reasons the rest of this file gives; these runs were separately
scheduled and their contention differs.

**Per-contig peak RSS rose by about 50%** -- chr1 4.30 -> 6.05 GB, chr4 3.05 -> 4.60 GB -- because
every record's render inputs are retained until the barrier settles. chr1's retention was measured
directly at 1,063 MB against a 1.23 GB projection, so the estimate that decision rested on was 14%
conservative rather than wrong. `schedule_wgs.py`'s memory predictions do NOT model this, and happen
to remain safe only because they were already conservative (chr1 predicted 7.16 GB, actual 6.05). They
should be refitted from these numbers rather than left to that coincidence.

**Whole-genome wall clock on this machine is not a measurement of the caller**, and the cleanest
demonstration of that is the single-sweep nested arm against the inline one. Summed per-contig wall
clock went 163.8 → 213.1 minutes, +30%. Summed CPU went 457.3 → 472.0 minutes, **+3.2%** — and that
second number is the cost of the change. (The review fixes then took summed CPU back down to 459.9
minutes over the same 24 contigs, −2.6%, and the packed run end-to-end to 65.1 minutes.)

The gap is entirely six contigs that were starved of cores in the later run. Thread occupancy,
`(user + sys) / real`, is the diagnostic:

| contig | inline | single sweep | wall clock |
|---|---|---|---|
| chr17 | 2.83 | **1.09** | 4.5 → 12.6 min |
| chr22 | 3.06 | **1.35** | 5.6 → 13.6 min |
| chr11 | 2.94 | **1.46** | 5.5 → 11.8 min |
| chr13 | 2.62 | **1.49** | 5.1 → 10.1 min |
| chr14 | 2.65 | **1.49** | 7.7 → 15.3 min |
| chr16 | 2.92 | **1.55** | 4.2 → 8.3 min |
| the other 18 | 2.47–3.20 | 2.24–3.68 | within ±20% |

A contig that got one core where it previously got three takes three times as long having done the
same work, and its CPU total says so. Read wall clock as a measure of what else was running.

Earlier full runs of the same scheduler took 54.3 and 60.9 minutes end to end against 144.4 minutes
summed, a 2.66x and 2.37x packing speedup, with nothing between them but load — one had `vg`'s own
test suite competing for cores for ten of its minutes. **Treat any single whole-genome wall clock as
±10% at best, and prefer CPU time when comparing two builds.**

```bash
python3 scripts/wgs/schedule_wgs.py --work work/wgs      # calls every contig
bash    scripts/wgs/assemble_wgs.sh                      # one VCF + one mosaic
python3 scripts/wgs/bench_wgs.py --work work/wgs --out docs/wgs-results.md
```

## The mosaic, and why assembly is not `cat`

The genome mosaic is **180,858 segments over 5,037,872 sites in 14.27 MB**, at 78.9 bytes per
segment.

**92.28% of segments carry a GBWT position**, down from 99.82% before nested sites entered the
mosaic, and the shortfall is one identifiable population rather than a degradation: of the 13,960
segments without one, 12,813 are wildcard rows whose haplotype is `*` — no single panel haplotype
is named, so there is no position to record — and most are one to three sites long. That is the
phase-block fragmentation that nested ploidy-1 sites cause, which is tracked as its own problem and
is not a property of the mosaic format. Fixing the linkage layer's position and record keying (so a
nested child no longer loses its phasing to a parent at the same POS) cut that population from
13,676 to 12,813 and lifted coverage from 91.87%, which is a dent in the problem rather than a
solution to it.

Concatenating the per-contig files needs `scripts/wgs/concat_mosaic.sh`, not `cat`, because two
mosaic columns are relative to the graph that produced them:

- **`hap_index`** is the haplotype's position in *that chunk's* GBWT metadata, and the chunks do not
  agree on an ordering. Appending 24 files under one `#haplotype` table relabels haplotypes silently
  — no error, no missing data, a genome-wide file describing the wrong ancestry. The script reindexes
  on the `haplotype` (`sample#phase`) name, which is portable by construction.
- **`gbwt_offset`** cannot be fixed that way at all. It is a rank among the sequences visiting a
  node, and the whole-genome GBWT has more of them, so the same offset addresses a different path
  there. The output therefore names each contig's *own* graph in a `#contig` table instead of
  claiming a single whole-genome GBZ.

`start_node`/`end_node` do survive: `vg chunk` preserves whole-genome node IDs, which is why they
are the authoritative anchors.

The structural check worth keeping is that the two strands of every diploid contig agree on their
site total, and that the strand-0 total equals the VCF's record count — 4,742,752, exactly. chrY is
the only single-strand contig; chrX carries 298 strand-1 segments, which are its pseudoautosomal
regions arriving via `--ploidy-bed`.

## Why one contig at a time

The caller buffers every emitted record and every linkage site until the chain resolves, so peak
memory scales with the contig rather than the genome. Whole-genome in one process would need tens
of gigabytes; per contig the worst case measured is **6.1 GB** (chr3).

## The memory model, and a correction worth reading

Contigs are packed under a budget rather than run at fixed concurrency, using the truth's record
count for the contig as a predictor — known before the run, and a far better predictor of peak
memory than contig length.

```
peak GB ~ 2.25 + 11.2e-6 * emitted_records
```

Fitted on all 24 contigs of a full run. **The previous coefficient was nearly double this**, and it
made the scheduler throttle itself on a fiction:

| | old model | refitted | measured |
|---|---|---|---|
| chr1 (353,741 records) | 9.6 GB | 6.2 GB | **5.7 GB** |
| chr6 (284,529) | 8.2 GB | 5.4 GB | **5.0 GB** |
| chr20 (105,251) | 4.4 GB | 3.4 GB | **3.1 GB** |
| worst residual | 4.39 GB | **0.87 GB** | |

Every contig was overestimated, by up to 4.4 GB, so the budget refused packings that would have fit
comfortably. The model had been fitted from an earlier serial run whose memory behaviour no longer
holds.

## Thread count and concurrency

The caller uses about 3.5 of 10 cores at `-t 5`, so three concurrent jobs saturate the machine.
Lower `-t` is *more* CPU-efficient per unit of work — measured on chr20, `-t` 1/2/5 gives 0.99/1.79/
3.48 CPU for 422/247/142 s, so about 70% parallel efficiency at 5 against 90% at 2 — which suggests
running more, thinner jobs.

Measured on a six-contig subset with the corrected memory model, two replicates each:

| configuration | run 1 | run 2 | mean |
|---|---|---|---|
| **`-t 5`, 3 jobs** (shipped) | 896 s | 806 s | **851 s** |
| `-t 3`, 5 jobs | 852 s | 886 s | 869 s |
| `-t 2`, 6 jobs | 1012 s | — | 1012 s |

**`-t 5` with 3 jobs stays the default, and the second replicate is why.** On one run `-t 3`/5
looked 5% faster (852 against 896) and it was tempting to ship it. Replicated, the ordering
reverses (886 against 806): run-to-run variation of about ±45 s exceeds the difference between the
two configurations, so they are indistinguishable and the apparent gain was noise. A single
timing run on this machine cannot resolve 5%.

`-t 2` with 6 jobs is genuinely worse, and losing by 19% is well outside that noise. This is the
second time that configuration has lost, and it was worth re-testing: the first test predated the
memory refit, so the obvious suspicion was that the old model had throttled that arm below its
intended concurrency. **It had not.** The thinner-jobs reasoning is simply wrong here — past about
five concurrent jobs the per-window `posix_spawn` of a `gbz-base` process and the reopening of the
22 GB read database cost more than the extra parallelism returns.

### What the memory refit actually bought

Not speed. At three concurrent jobs the machine is CPU-bound, not memory-bound, so halving the
predicted footprint changes no scheduling decision on this hardware and the genome run takes the
same hour it did.

What it buys is that the budget now means something. The old model would have serialised a run
given a smaller `--budget-gb`, and would have blocked anyone raising `--max-jobs` on a machine with
more cores, both for no reason. A predictor wrong by 70% is worth fixing even when it is not
currently binding.

## What is *not* the bottleneck

**Not I/O**, despite appearances. Every read-fetch window spawns a `gbz-base` and reopens the read
database, which looks like enough to explain a process sitting near one CPU. Measured, it is wrong:
the caller parallelises at ~70% efficiency to `-t 5`. The spawn-per-window is real and still looks
like the bottleneck in the source; it is not the one that governs wall clock.

**Not the memory budget**, now that the model is honest. At three jobs the worst case in flight is
about 18 GB against a 24 GB budget. The remaining 8 GB of the machine is not slack either: the read
database is 22 GB and the OS page cache is doing real work, so squeezing it trades one bottleneck
for a worse one.

**Not the tail.** Largest-first scheduling puts chr1 and chr2 in the first wave; the longest single
contig is 15.4 minutes against a 60.9-minute total, so the critical path is not one slow contig.

## Resume, and why it checks the binary

Every stage skips work already done, but keyed on **freshness rather than existence**: a `.done`
marker, a compressed VCF or a cached score is reused only if it is newer than the `vg` binary or
the input that produced it.

This is not fastidiousness. Keying on existence produced, in one session, a coverage sweep that
kept pre-fix results across a rebuild and scored them as though they were the fixed caller, and a
whole-genome "refresh" that concatenated the previous run's compressed VCFs and reported numbers
byte-identical to the run it was meant to replace. Neither said anything was wrong; only file
timestamps did.

`scripts/test_harness.sh` asserts the property directly across every script that caches.
