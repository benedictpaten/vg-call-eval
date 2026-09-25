# Running a whole genome, and how fast it goes

`vg call --read-likelihood` on the 32-haplotype hap32 graph (34 panel haplotypes with the CHM13 and
GRCh38 paths), 30x reads, on a 10-core laptop with 32 GB. One contig per invocation, several contigs
at once, packed under a memory budget.

## The current run

2026-09-23, on the defaults vg `0cab3fbd4` ships (the run passed `--mismap-max 0.95` to the build
before it, the one default that commit changes):

| | short reads, 24 contigs |
|---|---|
| CPU (user + sys) | **8.58 h** (6.20 + 2.39) |
| peak RSS, worst contig | **10.0 GiB** (chr14) |
| wall clock | **86.5 min**, packed about two at a time (`--budget-gb 24 --threads 5 --max-jobs 3`) |
| slowest contig | chr21, 1,214 s wall and 2,306 CPU s, last to finish |

None of these is a serial cost. Run alone, chr20 takes 172 s at a 7.5 GiB peak RSS (the tier-2
harness); in this run it took 313 s at 3.8 GiB with two other contigs alongside. CPU and RSS come from
`/usr/bin/time -l` around each `vg call`, so CPU includes the `gbz-base` processes it waits for.

```bash
python3 scripts/wgs/schedule_wgs.py --work work/<run> --budget-gb 24 --threads 5 --max-jobs 3   # calls every contig
W=work/<run> OUT=work/<run>/HG002 bash scripts/wgs/assemble_wgs.sh                             # one VCF (and the mosaic, below)
python3 scripts/wgs/bench_wgs.py --work work/<run> --out work/<run>/score/wgs-summary.md --threads 5
```

`bench_wgs.py` scores every contig and writes chr1-22+X totals to its own summary. Never point its
`--out` at docs/wgs-results.md: that page is maintained by hand. Autosome totals, and the insertion
and deletion rows, come from `small()` and `sv()` in `scripts/bench_metrics.py` over the run's
`score/per-contig.json`. The mosaic half of `assemble_wgs.sh` rejects mosaic-version 5 files, so it
writes an empty genome mosaic (see *The mosaic* below).

## Decide-then-render (2026-08-23): read I/O restored, CPU up

This arm, measured 2026-08-23, settles every genotype before building its record, which is how
records are still built. The cost that mattered was read I/O, because the arm before it bought the
same coherence guarantee by re-reading the contig once per generation:

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
conservative rather than wrong. The scheduler's memory model is fitted on the current binary; see
*The memory model* below.

**Whole-genome wall clock on this machine is not a measurement of the caller**, and the cleanest
demonstration of that is the single-sweep nested arm against the inline one (2026-08-19). Summed
per-contig wall clock went 163.8 → 213.1 minutes, +30%. Summed CPU went 457.3 → 472.0 minutes,
**+3.2%** — and that second number is the cost of the change. (The review fixes then took summed CPU back down to 459.9
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

Earlier full runs of the same scheduler (2026-08-16 and 08-17) took 54.3 and 60.9 minutes end to
end against 144.4 minutes summed, a 2.66x and 2.37x packing speedup, with nothing between them but
load — one had `vg`'s own test suite competing for cores for ten of its minutes. **Treat any single
whole-genome wall clock as ±10% at best, and prefer CPU time when comparing two builds.**

## The mosaic, and why assembly is not `cat`

No genome mosaic has been assembled since the per-contig files moved to mosaic-version 5:
`concat_mosaic.sh` accepts only versions 3 and 4 and rejects version 5, so a run's
`HG002.mosaic.tsv` is empty and the per-contig `chr*.mosaic.tsv` files are the mosaic.

The last genome mosaic, assembled 2026-08-20, was 180,858 segments over 5,037,872 sites in
14.27 MB. 92.28% of its segments carried a GBWT position, and the shortfall was one identifiable
population rather than a degradation: of the 13,960 segments without one, 12,813 were wildcard rows
whose haplotype is `*` — no single panel haplotype is named, so there is no position to record —
and most were one to three sites long. That is the phase-block fragmentation that nested ploidy-1
sites cause, which is tracked as its own problem and is not a property of the mosaic format.

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
site total, and that the strand-0 total equals the VCF's record count; on the 2026-08-17 mosaic it
did, exactly (4,742,752). chrY is the only single-strand contig; chrX's strand-1 segments are its
pseudoautosomal regions arriving via `--ploidy-bed`.

## Why one contig at a time

The caller buffers every emitted record and every linkage site until the chain resolves, so peak
memory scales with the contig rather than the genome. Whole-genome in one process would need tens
of gigabytes; per contig the worst peak RSS in the current run is **10.0 GiB** (chr14), then chr6
8.8, chrX 8.0 and chr1 and chr4 7.9 GiB, all measured packed. Run alone at the 0.7 ceiling, chr2
reached 10.7 GiB.

## The memory model, and a correction worth reading

Contigs are packed under a budget rather than run at fixed concurrency, using the contig's
small-variant truth record count as the predictor — known before the run, and a far better
predictor of peak memory than contig length.

```
peak GiB ~ 5.57 + 13.0e-6 * truth_records
```

Refitted 2026-09-22 on 20 contigs of the current binary run one at a time at the 0.7 ceiling, where its
residuals span −1.51 GiB (chr18) to +1.91 GiB (chr14). **The model before it under-predicted every
contig**, by 2.08 to 5.65 GiB (worst chr14): it was fitted while `--max-snarl-edges` still capped
large snarls, and uncapping raised peak memory. That is the dangerous direction — it would have
packed three 10 GiB contigs into a 24 GiB budget and swapped.

Against the packed 2026-09-23 run, in GiB:

| contig (truth records) | predicted | peak RSS | peak footprint |
|---|---|---|---|
| chr1 (438,017) | 11.26 | 7.93 | 13.02 |
| chr2 (379,611) | 10.50 | 5.18 | 11.97 |
| chr6 (313,919) | 9.65 | 8.78 | 11.23 |
| chr14 (235,337) | 8.63 | **10.00** | **12.14** |
| chr20 (163,602) | 7.70 | 3.77 | 6.23 |
| chr21 (146,784) | 7.48 | 7.24 | 7.54 |
| chrX (132,387) | 7.29 | 7.95 | 9.75 |

Packed, peak RSS reads low: it came in under the prediction on 22 of 24 contigs (from −5.32 GiB for
chr2 to +1.37 for chr14, mean −2.39), and chr2, 10.7 GiB alone at 0.7, peaked at 5.18. macOS's peak memory
footprint reads the other way, over the prediction on 12 of 24 (from −1.74 GiB for chr18 to +3.51
for chr14, mean +0.26). The model is fitted to single-job RSS, so these are a check on packing, not
a refit. Only chr14 and chrX exceed it on both measures.

## Thread count and concurrency

Measured 2026-08-16 on the capped binary, under the memory model of that date, and kept as a record.
`-t 5` and `--max-jobs 3` are still what the scheduler is given, but the budget, not `--max-jobs`,
now sets concurrency: see *What the budget does now*.

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

## What the budget does now

It binds. chr1 to chr4 are predicted at 10.15 to 11.26 GiB, so three large contigs do not fit in 24
and the scheduler holds the third back: the first pair alone, chr1 and chr2, is 21.8 GiB predicted.
The 2026-09-23 run had three jobs in flight for 14.1 of its 86.5 minutes, two for 61.5 and one for
10.9, an effective concurrency of 2.03. That is why it took 86.5 minutes where the 2026-09-12 run,
on the capped binary under the old model, took 61.4 at similar CPU (8.58 against 8.80 CPU-hours).

The remaining 8 GB of the machine is not slack: the read database is 22 GB and the OS page cache is
doing real work, so squeezing it trades one bottleneck for a worse one.

## What is and is not the bottleneck

**The memory budget is**, for the first hour of the 2026-09-23 run: see above.

**The tail is.** chr21 is one of the smallest contigs by truth records (146,784, predicted 7.48 GiB),
so largest-first scheduling starts it late, at 10:38 in a run that began at 09:32. Yet it has the
third-highest CPU of any contig, 2,306 s after chr2 and chr1, takes 1,214 s of wall, and finishes
last, alone for about the final 11 minutes. Its record count does not predict its cost, and the
model has no term for it.

**Not I/O**, as measured on 2026-08-16. Every read-fetch window spawns a `gbz-base` and reopens the
read database, which looks like enough to explain a process sitting near one CPU. Measured then, it
was wrong: the caller parallelised at ~70% efficiency to `-t 5`. The spawn-per-window is real and
still looks like the bottleneck in the source. The current binary changed read-query sizing
(`8acbb43a2`) and has not been profiled, so this is unverified on it.

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
