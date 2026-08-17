# Running a whole genome, and how fast it goes

`vg call --read-likelihood` on a 34-haplotype HPRC graph, 30x reads, on a 10-core laptop with
32 GB. One contig per invocation, several contigs at once, packed under a memory budget.

**Current: 24 contigs in 54.3 minutes**, against 144.4 minutes of summed per-contig time — a
2.66x speedup. An earlier full run of the same scheduler took 60.9 minutes (2.37x). Nothing between
them was a performance change, and the faster run had `vg`'s own test suite competing for cores for
ten of its minutes, so treat the difference as run-to-run variation on a laptop and not as a gain:
whole-genome wall clock here is reproducible to about ±10%, which is worth knowing before reading
any single number as a measurement.

```bash
python3 scripts/wgs/schedule_wgs.py --work work/wgs      # calls every contig
bash    scripts/wgs/assemble_wgs.sh                      # one VCF + one mosaic
python3 scripts/wgs/bench_wgs.py --work work/wgs --out docs/wgs-results.md
```

## The mosaic, and why assembly is not `cat`

The genome mosaic is **143,365 segments over 4,742,752 sites in 11.05 MB**, 3.46 MB gzipped, at 80.8
bytes per segment. 99.82% of segments carry a GBWT position; 390 are fragment splits.

Concatenating the per-contig files needs `scripts/wgs/concat_mosaic.sh`, not `cat`, because two
mosaic v2 columns are relative to the graph that produced them:

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
