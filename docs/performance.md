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

## Where it stands after the cleanup pass

chr20, `-t 5`, byte-identical output at every step:

| | wall |
|---|---|
| before any of this | 328.3 s |
| key index + clamped-chain guard | 197.2 s |
| compact space built once, searched not mapped | **187.6 s** |
| ... with cached snarls (`-r`) | **148.5 s** |
| ... at `-t 8` rather than `-t 5`, no `-r` | **159.1 s** |

**A caution on peak RSS.** Across six runs that produced byte-identical output, `maximum resident
set size` ranged from 3.66 to 4.77 GB. That is allocator variance, and it is wider than several of
the RSS differences quoted earlier in this file and in the commit messages. Read those as noise
unless the change is structural -- the clamped-chain guard, which stopped building 192,045 `Site`
objects five times over, is; the rest are not.

## What is left, in order

Re-profiled after the first two fixes. `__psynch_mutexwait` is gone from the profile entirely and
the subprocess wait is now the largest single consumer:

| | before | after |
|---|---|---|
| `__wait4` (gbz-base child) | 16.9% | **33.5%** |
| `__psynch_cvwait` (threads parked) | 34.7% | 30.9% |
| `__psynch_mutexwait` | 20.9% | **absent** |
| working | 32.4% | 35.4% |

### a. Snarl decomposition -- DONE, 45 s

`IntegratedSnarlFinder::find_snarls_parallel` parked four of five threads for 46 s. `vg snarls`
reproduces the caller's decomposition byte-for-byte, but only with **both** flags: `-P <ref path>`,
because `snarls_main.cpp:284` and `call_main.cpp:1424` apply the same `EXTRA_WEIGHT` to the same
first and last node of each reference path; and `-T`, because the symbolic projection keys chain
symbols on child chain boundaries and omitting trivial snarls changes the nested structure -- chr20
gained 236 records and 714 lines moved without it.

chr20 187.6 s -> 148.5 s. Wired into `prep_wgs.sh`, with `call_wgs.sh` passing `-r` when the file is
there. `vg snarls` itself costs 44.6 s, so one whole-genome pass breaks even and every run after the
first is ahead.

### b. `gbz-base` subprocess round-trips -- bounded, and smaller than it looked

33.5% of the thread budget is a worker asleep in `waitpid`. The obvious fix is to prefetch the next
window, since sites are visited in node-ID order. But the ceiling on *any* overlap is what
concurrency alone can recover, and a thread sweep finds it:

| threads | wall | CPU multiple | peak RSS |
|---|---|---|---|
| 5 | 187.6 s | 3.04 | 4.19 GB |
| 8 | **159.1 s** | 4.06 | 4.98 GB |
| 10 | 157.5 s | 4.69 | 5.09 GB |

Output identical at all three. So the whole overlap is worth about **30 s of 187.6**, and it
**saturates at 8** -- past that the limit is not cores but the read database's tolerance for
concurrent random reads, plus the ~60 s of serial phases Amdahl leaves behind.

That makes prefetching a **memory** lever rather than a speed one: `-t 8` buys the 28 s for 0.8 GB
of extra per-thread window cache, where prefetching would buy the same with one extra in-flight
query per thread. Worth building where memory is the binding constraint, which under
`schedule_wgs.py` packing contigs against a 24 GB budget it is. Not worth building to make a single
contig faster, because raising `-t` already does that for free.

**Actionable now:** use `-t 8` for single-contig work. Leave the whole-genome scheduler at 5 and let
it pack contigs, which fills the machine the same way.

Independently, a third of every query builds a GFA subgraph that goes to `/dev/null`.
`planning/gbz-base-c-api-request.md` now carries a finished ask for a reads-only query. It has not
been sent.

### c. Cutting the chain at pinned sites -- sized, not built

A per-site pin zeroes every state but one, so it severs the path: a site of this generation depends
only on the sites between its bracketing pins, and the chain could be cut there and decoded in
pieces. The caller now reports how much of each decode is live:

| generation | sites decoded | live | pinned |
|---|---|---|---|
| 1 | 209,244 | 17,199 | 192,045 |
| 2 | 211,480 | 2,236 | 209,244 |
| 3 | 211,930 | 450 | 211,480 |
| 4 | 211,952 | 22 | 211,930 |

Generations 1-4 cost 24.4 s of the 187.6, and sub-chaining would take that to roughly 2 s.

Not built, and the reason is in the details rather than the idea. The phase set is
`sites.front().position`, so it has to keep coming from the whole chain rather than from a piece of
it. And the severing rests on every pin being *accepted*, where `window_phasing` declines one whose
pair cannot spell the site's constrained genotype. Both are checkable and neither is free, and the
gate available is one chromosome.

Note also that the original framing of this was wrong: it assumed only a handful of non-nested sites
arrive after generation 0. Generation 1 has 17,199 of them.

### d. `WindowedSiteReadSource::deliver` -- dropped

3.39 billion reads considered to hand over 24.5 million sounds large and is not. The profile puts
`deliver` and `touches` together at 61 samples, **1.2 thread-seconds of 319** -- about one cycle an
iteration, which is what a streaming two-compare loop over a compact bounds array should cost.
Sorting and binary-searching would save on the order of a second. Not worth the code.

### e. The read ingestion path -- unscheduled, and the largest working cost

~12% of thread time against 2.1% for the scoring it feeds. Every read fetched -- 14.2 M on chr20 --
is parsed from GAF text into a protobuf `Alignment` with a `Path` of `Mapping`s and `Edit`s, cached,
then destroyed; the allocator attribution puts 21.7% of its samples inside `gaf_to_alignment` and
17.5% inside `get_next_record_from_gaf`. A leaner record -- node-id run, offsets, sequence, quality
-- would cut the parse, the allocation and the destruction together. Listed because the profile says
the ratio is 5:1 the wrong way, not because it should be done next: it is a read-path refactor and
wants its own plan.

## What the render-pass fix was actually worth: chr6, measured three ways

The chr20 figure in this file -- 328 s to 197 s -- understates the fix, and reporting it as "131
seconds" was the wrong unit. The defect was quadratic in the number of sites the linkage layer
holds, and chr20 holds 219,195 of them where chr6 holds 381,612 on the thin panel and 489,028 on
the rich one. The contig that showed it worst was the one nobody was re-measuring.

Ten `vg call` runs per build, both chr6 tier-2 datasets, one machine, one session, same dataset
order on both sides. Pre-session is `08a34da4b`; current is `648296d56`. **Every arm produces
byte-identical output across the two builds**, so all of what follows is speed.

Three axes, because they answer different questions and only the third is confound-free:

- **wall clock** moves with work done, waiting avoided, and page-cache state alike.
- **CPU seconds** (user+sys) drops the waiting but still moves with the machine's clock state.
- **instructions retired** is immune to frequency and cache state. If it falls, the caller genuinely
  does less. If it holds while cycles fall, the win was stalls and serialisation.

`--threads 5` throughout.

**chr6, 34-haplotype graph** — control factor **0.952**, entries vector 59.7 MB

| arm | linkage | wall | CPU | instructions | IPC | attributable | output |
|---|---|---|---|---|---|---|---|
| `poisson-z` | no | 223 → **197 s** | 597 → 541 s | 6.53 → 6.39 T (-2.2%) | 3.57 → 3.83 | **1.08x** | byte-identical |
| `readlik-support` | no | 392 → **354 s** | 1,332 → 1,242 s | 8.87 → 8.87 T (+0.0%) | 3.55 → 3.72 | **1.05x** | byte-identical |
| `readlik-nolink` | no | 286 → **292 s** | 1,009 → 1,024 s | 5.87 → 5.87 T (+0.0%) | 3.69 → 3.63 | **0.93x** | byte-identical |
| `readlik-nomismap` | **yes** | 916 → **374 s** | 1,721 → 1,183 s | 9.52 → 7.73 T (-18.8%) | 2.45 → 3.62 | **2.34x** | byte-identical |
| `readlik` | **yes** | 902 → **369 s** | 1,700 → 1,186 s | 9.51 → 7.74 T (-18.6%) | 2.47 → 3.65 | **2.33x** | byte-identical |

**chr6, 4-haplotype graph** — control factor **1.050**, entries vector 46.6 MB

| arm | linkage | wall | CPU | instructions | IPC | attributable | output |
|---|---|---|---|---|---|---|---|
| `poisson-z` | no | 159 → **163 s** | 364 → 367 s | 3.97 → 3.99 T (+0.4%) | 3.53 → 3.52 | **1.02x** | **22 of 289,002 differ** |
| `readlik-support` | no | 297 → **312 s** | 965 → 1,002 s | 6.31 → 6.29 T (-0.3%) | 3.67 → 3.52 | **1.00x** | byte-identical |
| `readlik-nolink` | no | 245 → **257 s** | 856 → 880 s | 4.72 → 4.72 T (-0.1%) | 3.42 → 3.38 | **1.00x** | byte-identical |
| `readlik-nomismap` | **yes** | 565 → **272 s** | 1,213 → 933 s | 5.94 → 5.01 T (-15.7%) | 2.39 → 3.34 | **2.18x** | byte-identical |
| `readlik` | **yes** | 564 → **275 s** | 1,202 → 934 s | 5.94 → 5.01 T (-15.8%) | 2.41 → 3.33 | **2.15x** | byte-identical |

### Reading it

**The controls are the two linkage-free read-likelihood arms, not the Poisson one.**
`readlik-support` and `readlik-nolink` hold instructions retired to **+0.0%** across the builds --
identical work, three significant figures -- while `poisson-z` wobbles 2%, because the Poisson path
walks Flow traversals whose work has a small scheduling dependence. Any instruction-count movement
on the linkage arms is therefore the fix and nothing else.

**A sequential A/B cannot avoid warming the machine for whichever side runs second.** The control
arms land between 0.93x and 1.08x, so there is a roughly +-8% band around "no change", and the
attributable column divides out the wall-clock part of it. It cannot manufacture a factor of two,
and it applies to every arm equally -- which is the point of the arms that show nothing.

**What the fix did, in proportion.** Instructions fell ~19% against a 0.0% baseline: 1.8 trillion
instructions of quadratic scan, deleted rather than rescheduled. IPC rose from ~2.46 to ~3.64,
landing on the linkage-free arms' 3.6-3.8 -- the old path was memory-bound streaming a 60 MB vector
per lookup and the new one is not. And wall clock improved about 1.7x more than CPU, because the
scan held `LinkageCollector::mutex`: four of five threads were parked in `__psynch_mutexwait`, so
removing it returns cycles *and* parallelism and the two compound. The chr20 sample profile
predicted that shape; the hardware counters confirm it without reference to it.

**The noise floor, stated rather than implied.** `readlik-nolink` came in 6 s *slower* after the fix
on provably identical instruction counts. Wall clock on this benchmark is worth about +-2% at best
even with the work pinned, which is why every claim here is quoted against a control rather than raw.

### The Poisson arm is not bit-reproducible, which is why it is not the control

`poisson-z` shows 22 of 289,002 records differing across the two builds on the 4-haplotype
dataset. That is not the fix reaching it -- neither `record_key_of` call site executes without a
linkage collector, and the arm has none. Running the *same* binary twice gives **20 differing lines
of 289,002**, so the rate is the same with the build held fixed.

What moves is confined to depth-derived floats -- `QUAL` 403.316 against 403.304, `GL`, `XD` -- with
`GT`, `AD` and `GQ` identical on every affected record. No genotype changes. The shape points at
accumulation order in the binned depth model rather than anything about calling.

It matters for two reasons. It explains why `poisson-z` was the one arm whose instruction count
wobbled (-2.2% where the read-likelihood controls held at +-0.4%), so the control factors here come
from `readlik-support` and `readlik-nolink` instead. And it means **byte-identity is not a usable
regression gate on the Poisson path**, which is worth knowing since it is the gate this branch has
leaned on all along for the read-likelihood arms -- where it does hold, exactly, across every one of
the twenty runs above.

### What this changes about how the earlier chr20 figure was reported

"chr20: 328 s to 197 s" is in this file above, and as a description of the fix it was the wrong
unit. The defect was quadratic in linkage-layer site count. chr20 holds 219,195 sites; chr6 holds
381,612 and 489,028. Quoting seconds from the smallest of the three understated the fix by more than
half and hid the thing that actually predicts its size. The ratio and the scaling term are the
reportable pair.

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
