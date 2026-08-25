# Where `vg call` spends its time

Profiled on chr20, HG002 against the 32-haplotype CHM13 panel, `-t 5` on a 10-core M-series
laptop with 32 GB. Sampled with macOS `sample` at 20 ms across every thread, one 4 s window
every 12 s for the whole run, so the profile covers the run rather than a guessed-at hot spot.

The command is the one `call_wgs.sh` issues:

```
vg call chr20.gbz -p CHM13#0#chr20 -s HG002 -d 2 -t 5 --progress \
    --read-likelihood --phased --mosaic-out chr20.mosaic.tsv \
    --gaf-base reads.hap32.gaf.db --gbz-base graph.hap32.gbz.db
```

## The complaint, quantified

chr20 had roughly doubled since the middle of August. Every arm on disk, same contig, same
`-t 5`:

| run | date | wall | user+sys | CPUx |
|---|---|---|---|---|
| stage2 | 08-17 | 140.6 s | 489.5 s | 3.48 |
| coh8 | 08-18 | 161.9 s | 539.7 s | 3.33 |
| d21 | 08-18 | 154.6 s | 524.8 s | 3.39 |
| wgs-single | 08-19 | 237.1 s | 626.3 s | 2.64 |
| wgs-current | 08-20 | 307.0 s | 601.3 s | **1.96** |
| wgs-dtr | 08-22 | 338.0 s | 753.5 s | 2.23 |
| wgs-atomize | 08-24 | 373.0 s | 799.7 s | 2.14 |

Total CPU rose 63%, but wall clock rose 165%. The gap is parallelism: the speedup fell from
3.4x to 2.1x out of a possible 5. Read as Amdahl, serial time went from about 53 s to about
266 s -- **the run got five times more single-threaded**, and that, not extra work, is what
was felt.

## Where the thread budget actually goes

5 threads x 328 s is 1,640 thread-seconds. Only 32% of sampled thread time was executing:

| | share of thread samples | what it is |
|---|---|---|
| `__psynch_cvwait` | 34.7% | OpenMP threads parked -- a phase running on one thread |
| `__psynch_mutexwait` | 20.9% | blocked on `LinkageCollector::mutex` |
| `__wait4` | 16.9% | blocked on a `gbz-base` child process |
| working | 32.4% | |

By phase, from the timestamped progress log:

| phase | s | threads busy |
|---|---|---|
| graph load | 1.9 | 1 |
| snarl decomposition | 46.0 | **1 of 5** |
| top-level calling | 79.9 | 5, about half of it blocked on `gbz-base` |
| nested descent | 21.0 | 5, same |
| linkage generation 0 | 7.0 | **1 of 5** |
| linkage generations 1-5 | 58.4 | **1 of 5** |
| render sweep | 95.7 | **nominally 5, effectively 1** |
| VCF write | 10.8 | 1 |

## 1. The render sweep was quadratic, under a global lock

`renderRetainedRecords` is `#pragma omp parallel for` over the retained records, and each
record calls two `LinkageCollector` accessors: `settled_traversals` for the settled pair, and
`set_allele_map` on the way out of `emit_variant`. Both -- and five more like them -- took the
collector's one mutex and then **walked the whole `entries` vector** looking for a matching
record key.

219,316 records against 219,195 entries of 128 bytes each. Quadratic, and serialised, at once.
Through this phase 79% of thread samples are `__psynch_mutexwait` and the single thread making
progress is inside those two scans.

Fixed by an intrusive chain: `first_by_key` to the head of a key's entries, `next_same_key`
along it. First-match-in-insertion-order is preserved exactly, duplicate keys included, so no
accessor changed meaning. `entries` is only ever appended to, so the indices are stable.

**95.7 s -> 1.6 s.**

## 2. Every generation re-decoded the whole chromosome and discarded it

`resolve_generation` puts every earlier generation's sites in the chain so that this
generation's sites have transition context. An earlier site is *clamped*: a delta emission at
its settled genotype, pinned to the phase already in the VCF. Both the result loop and the
phasing loop then skip it.

Only nested sites arrive after generation 0, and nested sites are held out of the diploid runs
by construction. So generations 1 through 5 each ran a forward-backward and a Viterbi over
chr20's 192,045 top-level sites and threw every answer away -- which is why generation 5, with
**seven sites**, cost 6.11 s, the same as generation 0 with 192,045.

Fixed with a per-chain guard on whether any site in it belongs to this generation.

**58.4 s -> 25.8 s**, and peak RSS 4.05 -> 3.93 GB, since the skipped passes no longer build
192,045 `Site` objects to discard.

## Result

| | before | after |
|---|---|---|
| chr20 wall | 328.3 s | **197.2 s** |
| user+sys | 723.5 s | 593.8 s |
| peak RSS | 4.40 GB | 3.93 GB |
| records | 115,392 | 115,392 |

The VCF body and the mosaic TSV are **byte-identical**. `vg test` 12,547,761 assertions in 858
cases, `t/18_vg_call.t` 309/309.

## What is left, in order

Re-profiled after the two fixes, same method. `__psynch_mutexwait` is gone from the profile
entirely, and the subprocess wait is now the single largest consumer of the thread budget:

| | before | after |
|---|---|---|
| `__wait4` (gbz-base child) | 16.9% | **33.5%** |
| `__psynch_cvwait` (threads parked) | 34.7% | 30.9% |
| `__psynch_mutexwait` | 20.9% | **absent** |
| working | 32.4% | 35.4% |


### a. Snarl decomposition, 46 s on one thread (23% of the run now)

`IntegratedSnarlFinder::find_snarls_parallel` parks four of five threads for 46 s. It is
upstream code and the same on master.

`vg call -r snarls.pb` skips it, but not for free: the caller builds the finder with
`extra_node_weight` biasing the decomposition towards the reference path's endpoints, and
`vg snarls` does not do that, so a cached file from `vg snarls` is not the same decomposition.
Worth checking whether the difference reaches the calls before adopting it -- if it does not,
this is 46 s per contig for a one-off cost, which matters most for the A/B loop where the same
contig is called over and over.

### b. `gbz-base` subprocess round-trips, about 29% of thread time in the calling phases

Every read fetch `posix_spawn`s `gbz-base query`, which opens a 6.8 GB graph database and a
21 GB read database, extracts a subgraph, writes GAF to a temp file and exits; vg then parses
the GAF back. chr20 does this 1,446 times. A single-node query -- pure startup -- measures
0.05-0.10 s, and the profile puts about 0.19 s of blocked thread time on each spawn, so
roughly half the cost is opening databases that were open a moment ago.

Three levers, and they do three different things. Only one of them reduces the spawn count:

| | spawns | cost of each | thread blocked on it |
|---|---|---|---|
| prefetch the next window | unchanged | unchanged | hidden |
| bigger `--read-window` | fewer | higher | roughly in proportion |
| persistent worker | unchanged | much lower | lower |

- **Prefetch the next window.** This does *not* reduce spawns. Today a worker thread issues its
  query, sleeps in `waitpid` for ~0.19 s, parses the GAF and only then genotypes the window.
  Because sites are visited in node-ID order the thread knows window N+1 before it starts
  consuming window N, so it can spawn that query into a second temp file and `waitpid` only
  once window N is exhausted. Same children doing the same work -- but running alongside the
  genotyping instead of alternating with it.
- **Sweep `--read-window`.** 4096 today. This is the lever that reduces the spawn count:
  doubling the window halves it, at the price of a larger cached window and a longer per-query
  scan. One-line experiment, no code.
- **A persistent worker.** One `gbz-base` per thread, fed queries over a pipe, would remove the
  per-spawn startup -- half the cost by the single-node measurement above -- but `gbz-base
  query` is one-shot, so this needs a change in https://github.com/jltsiren/gbz-base.

**Prefetch only pays if there are cores to overlap into, and there are.** Measured two ways.
Total CPU across `vg` and every `gbz-base` child, sampled once a second through the calling
phase, runs at a **median 357% of a possible 1000%** -- between three and four of ten cores,
with the rest idle. And raising the thread count, which is the crude version of the same
overlap, works: chr20 at `-t 8` is **173.6 s against 197.2 s at `-t 5`**, with the top-level
calling phase falling 81.0 s to 58.8 s while snarl decomposition sits unchanged at 48 s.

That the children add so little CPU says they are I/O-bound on the 21 GB read database rather
than compute-bound, so the ceiling here is the disk's tolerance for concurrent random reads,
not the core count. `-t 8` is evidence that ceiling has not been reached.

One caveat on where this matters. `schedule_wgs.py` already packs several contigs onto the
machine at once, so a whole-genome run fills the cores by other means and would gain less than
these figures suggest. The single-contig case -- the A/B loop, where the same chromosome is
called over and over -- is where it is worth the most.

### c. Cutting the chain at pinned sites

Generations 1-4 still decode, because each has at least one non-nested site and one such site
keeps its whole 192,045-site chain alive. A per-site pin zeroes every state but one, so it
severs the path: a free site's answer depends only on the sites between its bracketing pins.
Decoding only those sub-chains would take generations 1-4 to near zero -- about 24 s -- but it
is a change to the model's windowing, where `phasing()` chains its windows sequentially, and it
wants its own measurement rather than being folded into a cleanup.

### d. `WindowedSiteReadSource::deliver` scans the whole window per site query

3.39 billion reads considered to deliver 24.5 million, 138 rejected for every one kept. The
bounds array is compact and the loop is already the cheap version of this, so it is only about
0.3% of samples -- but sorting a window's bounds once on fetch and binary-searching per query
would remove most of it. Small, and the least risky thing on this list.

## Redundancy

**Comments are not duplicated.** Across the branch's eleven files, comment blocks of two or
more lines that appear more than once account for 49 lines out of about 4,900 -- 1%. Header
comment density is high by design (48-68% in `linkage_model.hpp`, `allele_likelihood.hpp`,
`graph_caller.hpp`) and it is rationale rather than restatement.

**One misleading message, fixed.** `Computing snarls` / `Computed snarls` bracketed the
finder's constructor rather than the decomposition, so the log said the step was done in 1.4 s
and then went silent for 46 s.

**Duplicated code, inherited.** `FlowCaller::call_snarl_internal` and
`NestedFlowCaller::call_snarl_recursive` share four copy-pasted comment blocks and the code
around them. `NestedFlowCaller` is the upstream `-A` path, about 420 lines, unreachable from
the default caller -- out of scope for this branch to touch, and the `child_support_map` bug
already filed against it lives there.

**Options.** The branch adds 31 to `vg call`'s 70. Of those:

| option | eval harness | vg TAP tests |
|---|---|---|
| `--read-min-mapq` | 0 | 0 |
| `--atomize-blocks` | 0 | 0 |
| `--no-atomize-blocks` | 0 | 0 |
| `--linkage-scale`, `--linkage-freq-prior` | 2 | 0 |
| `--depth-count-raw`, `--depth-quality` | sweeps | 0 |

`--read-min-mapq` is referenced nowhere outside `call_main.cpp` and is the one candidate for
outright removal. `--no-atomize-blocks` is the escape hatch for a default that changed on
08-23 and should have a test rather than be deleted. `--atomize-blocks`, `--nested` and
`--phased` are no longer defaults-changing, but they are not no-ops either: each sets an
`explicit` flag that turns a silent decline into an error, which is worth keeping.
