# Parked work: re-fitting the ONT preset under `--realign`, and the PR #4990 review

Written 2026-09-14, after `--realign` shipped opt-in (vg `c60503ac3`, PR #4990). Nothing here is
started. Both parts are independent of each other and of whatever is being worked on now.

---

## Part 1 — the ONT preset was fitted against a walk that no longer runs by default

> **CLOSED, 2026-09-14, negative. Nothing moves.** The hypothesis below is wrong and the table
> immediately under it is wrong: `--gap-open` *had* been swept on the anchor-constrained walk, in
> `exact-walk.md`, before `--realign` shipped. Measured again anyway under the shipped preset, and
> both open parameters confirm their shipped values — `--gap-open 1` on the floor of its legal
> range with a 0.042 cliff above it, `--mismap-min 0.05` interior. See
> [Gate 0](#gate-0--answered-2026-09-14-negative-and-the-parameter-is-now-closed) for the numbers.
> Everything from here to that section is kept as written, so the reasoning that turned out wrong
> is on the record beside the result.


`--preset ont` sets seven things. Only **one** has been re-fitted since `--realign` became part of
it:

| preset member | fitted against | re-fitted under `--realign`? |
|---|---|---|
| `--insertion-nats 0.9` | greedy | **yes** — optimum flat 0.9-1.2, 0.9 kept |
| `--gap-open 1` | greedy | **no** |
| `--gap-extend 1` | greedy | **no** |
| `--mismap-min 0.05` | greedy | **no** |
| `--realign` | -- | (the change itself) |
| `--read-phasing` | greedy | no, but it does not touch the walk |
| `--regenotype` | greedy | no, but it does not touch the walk |

The first three are the ones that interact with the walk, and `--gap-open` is the big lever: worth
**+0.067 indel F1** on its own when the preset was built.

### The hypothesis worth testing, which is sharper than "parameters may drift"

The two walks do **not** cost an identical correspondence. Greedy charges a fresh gap **open** per
inserted read visit; the optimal walk's insertion state opens once and extends. A `k`-visit
insertion run therefore differs by

    (k - 1) * (gap_open - gap_extend)

which is **exactly zero when `gap_open == gap_extend`** — and the fitted preset is `gap-open 1
--gap-extend 1`.

So the fitted value sits precisely at the point where greedy's over-charging of multi-visit
insertion runs vanishes. That is unlikely to be a coincidence: **part of what `--gap-open 1` was
fitted to correct may be greedy's handling of insertion runs, which `--realign` now fixes
properly.** If so, the optimum under `--realign` should move *up*, because the walk no longer needs
the gap parameters to compensate for it.

That is a real, falsifiable prediction and the reason this is worth doing rather than box-ticking.

### Method

Standard: sweep on chr20, **chr6 confirms a chosen value and never chooses between candidates.**

    python3 work/ont-preset/arm.py --tag <tag> --arm ont-chr20 --vg <pinned> --preset ont \
        --gap-open <v>

1. **`--gap-open` first** — the big lever and the one the hypothesis names. Sweep 1, 2, 3, 4, 6.
   Include 1 as a control: it must reproduce `bd-c20` (indel 0.86659) exactly.
2. **`--mismap-min`** next, at whatever gap-open won. Sweep 0.02, 0.05, 0.10, 0.15.
3. **Interaction check** only if both moved: the pair may not be separable.
4. Confirm the winner on chr6 ONT, and check chr20 short reads are untouched (they should be --
   none of these are short-read defaults, but `--gap-open` is, so verify the greedy path is
   unaffected by anything adopted).

### Cost

About 6 minutes per arm (≈260s call + ≈90s scoring) on chr20 ONT, run one at a time. A 5-point
gap-open sweep is ~30 min; the whole plan including chr6 confirmation, ~1.5 hours of wall clock.

### Adoption threshold

**Do not move a fitted default for less than ~0.003 indel F1.** The `--insertion-nats` re-sweep
peaked 0.00045 above the shipped value and was correctly left alone; that is the calibration for
what counts as noise here. Neighbouring grid points on a real curve move 0.0026-0.0039.

### Gate 0 — ANSWERED 2026-09-14. Negative, and the parameter is now closed.

Run **only the `--gap-open 2` arm** first, one measurement, 6 minutes. The hypothesis predicts it
beats `--gap-open 1`. If it does not, the compensation story is wrong and the rest of the sweep is
unlikely to pay; stop and record the negative.

**It does not.** vg `3519d4980`, chr20 ONT, `--preset ont` with one flag varied, control and test
on the same binary:

| arm | indel F1 | indel P | indel R | SNV F1 | ALL F1 | indel FP |
|---|---|---|---|---|---|---|
| `--gap-open 1` (shipped) | **0.86659** | 0.8558 | 0.8777 | 0.98563 | 0.95863 | 3,318 |
| `--gap-open 2` | 0.82441 | 0.8017 | 0.8484 | 0.98543 | 0.94808 | 4,778 |
| `--mismap-min 0.02` | 0.86117 | 0.8480 | 0.8748 | 0.98532 | 0.95700 | 3,534 |
| **`--mismap-min 0.05` (shipped)** | **0.86659** | 0.8558 | 0.8777 | 0.98563 | 0.95863 | 3,318 |
| `--mismap-min 0.10` | 0.85825 | 0.8454 | 0.8715 | 0.98354 | 0.95494 | 3,594 |

`--gap-open 2` loses **0.0422** indel F1 — fourteen times the adoption threshold, in the wrong
direction, with precision *and* recall both down and 44% more indel false positives. The
homopolymer cell moves the same way: HP>=5 one-base FPs go 1,429 -> 2,433, and their share of all
FPs rises 38.5% -> 46.9%. So raising the gap scale inflates exactly the class `--gap-open 1` was
fitted to suppress, and the compensation story is wrong.

**And the premise was already refuted in this repo.** `exact-walk.md` records a gap_open sweep run
on the anchor-constrained walk — which *is* `--realign` — during the walk investigation:
0.86758 at gap_open 1 against 0.82313 at 2, with the note "all three walks share the same optimum
at gap_open 1 ... the anchored walk needs no re-tuning". The table at the top of this Part, which
says `--gap-open 1` was never re-fitted under `--realign`, is wrong. This measurement reproduces
that earlier result under the shipped preset rather than discovering it.

**`--gap-open` cannot go lower.** The arm at 0 exits immediately: `--gap-open must be between 1 and
127`, and so must `--gap-extend`. The fitted value sits on the **floor of the flag's legal range**,
with the curve rising steeply above it, so there is no sweep left to run on either gap parameter —
and the floor is right, because at `gap_open 0, gap_extend 1` a one-base gap costs nothing.

**`--mismap-min 0.05` is an interior optimum** under `--realign`, losing 0.0054 at 0.02 and 0.0083
at 0.10. Both neighbours are well outside noise and both are worse.

**Part 1 is closed. Nothing moves.** Three of the preset's four walk-adjacent members are now
confirmed at their shipped values under the shipped walk; the fourth, `--gap-extend`, is at the
same legal floor as `--gap-open` and has nowhere to go. No chr6 confirmation was run, correctly:
chr6 confirms an adopted change, and nothing is being adopted.

Two things fell out of the control arm for free. It reproduces `bd-c20` **to the digit**
(0.86659 / 0.98563 / 0.95863, FP 3,318, FN 2,407), which confirms that the anchor `gqn` fix, the
phase-ordering fix and `--anchors-hom-split` are all VCF-inert as intended — a byte gate would have
been better still, but this is the gate that was already being run.

---

## Part 2 — PR #4990 review (adamnovak, comment 5671460225)

> **STATUS 2026-09-14.** 2a **done** (vg, option table owns the ownership; the old list had
> already drifted on four live flags — see below). 2c **done** as `doc/read-likelihood-architecture.md`
> in vg (`eefa809a1`), which answers both questions and finds a third thing neither question named.
> 2b was decided by the author. 2d is scoped by 2c's "order of work" and is **not started**: its
> first item, splitting `linkage_model`, is independent of the layout; the folder move and the
> `VCFOutputCaller` → `FlowCaller` demotion both churn an open PR and want sign-off first.


Three separable things. Only the first is a defect.

### 2a. CLI options are parsed twice, by hand — the actionable one

Verified. `src/subcommand/call_main.cpp`:

- a **hand-maintained list of 56 flag strings** (`read_likelihood_only`, line ~1555) against **105
  entries** in `long_options`
- a **hand-rolled reimplementation of getopt_long's resolution** (line ~1575): exact match, else
  unambiguous prefix, because whole-token matching let `--mosaic` through as `--mosaic-out`

It is as brittle as the review says, with two near-misses in a single session: `--insertion-nats`
was silently accepted and dropped because it was never added to the list (found by reading the test
file, not by any check), and `--realign`/`--no-realign` had to be remembered into it.

**Bounded fix, which does not pre-empt the layout question:** derive the refusal set from the option
table — tag each `long_options` entry with the subsystem that owns it — so a new flag cannot drift
out of sync and the hand-rolled prefix matcher disappears. Roughly a day across ~105 options, and it
changes error-message wording, so it wants sign-off before starting.

**DONE.** The table is now `CallOption {name, has_arg, val, owner}` with `OWN_CORE` /
`OWN_READ_LIKELIHOOD`, and getopt's `struct option` array is generated from it. The parse loop
records each `val` getopt returns — in argv order, without repeats — and the check is a scan of
that, so the hand-rolled prefix matcher is gone: getopt had already done the resolution. 106
options, 61 read-likelihood-owned, 45 core. The message is unchanged except that repeats are now
deduped.

**The list had drifted, and the drift was live** — demonstrated on the shipped binary before the
change, not inferred:

| flag | before | after |
|---|---|---|
| `--mosaic-out` | refused | refused |
| `--mosaic-patch-gaps` | **accepted, silently dropped** | refused |
| `--no-mosaic-patch-gaps` | **accepted, silently dropped** | refused |
| `--no-mosaic-nested` | **accepted, silently dropped** | refused |
| `--mosaic-break-unexplained` | **accepted, silently dropped** | refused |
| `--traversals` (genuinely shared) | accepted | accepted |

Ownership had been written down in three places that could disagree — the helptext's section
headings, the refusal list, and *nothing at all* in the option table, which is the one place a new
option must be touched. Now it is in the table.

Three TAP tests guard it, including that an `OWN_CORE` option is still let through: a check that
refuses everything would pass the obvious assertion and be useless.

**And it generalises, which is the half the review actually asked for** — "so we can get errors
about using flags for a subsystem you aren't using, for free". The owner enum now carries
`OWN_ANCHORS`, `OWN_MOSAIC` and `OWN_REGENOTYPE` alongside `OWN_READ_LIKELIHOOD`, and the check is
a list of gates rather than one:

| owner | enabling flag | options | were they refused before? |
|---|---|---|---|
| `OWN_ANCHORS` | `--anchors-out` | 8 | no — silently dropped |
| `OWN_MOSAIC` | `--mosaic-out` | 4 | no — silently dropped |
| `OWN_REGENOTYPE` | `--regenotype` | 7 | no — silently dropped |

So `--anchors-min-q 30` with no `--anchors-out`, or `--regeno-temper 0.2` with no `--regenotype`,
now say so instead of doing nothing. Two details that are not cosmetic:

- **The outer gate wins.** An `--anchors-*` flag with neither `--anchors-out` nor
  `--read-likelihood` is reported against `--read-likelihood`, which is the one the user has to fix
  first; otherwise they would fix `--anchors-out` and hit a second error.
- **The regenotype gate tests `regenotype && read_phasing`, not `regenotype`.** `--preset ont` arms
  regenotyping and `--no-read-phasing` disarms it again, and that resolution happens *after* this
  check — so testing the raw flag would accept a `--regeno-*` option here and leave it inert there.
  That combination is documented and must not error, so the gate has to know about it.

Four more TAP tests, 421 total.

### 2b. Layout: folders. DECIDED (author, 2026-09-14)

> "I think we want folders for sanity. The flat directory structure has gotten too big."

It has: **297 `.cpp`/`.hpp` files in a flat `src/`.** So this is not only about the read-likelihood
family, though that family is what PR #4990 adds and the natural first tenant.

### 2c. Architecture review — answers 2 and 3, and proposes the layout

The author has asked for a review rather than an opinion, so this is a task with a deliverable, not
a question to forward. It must answer:

1. **How coupled is each subsystem to `ReadLikelihoodSnarlCaller`?** Per subsystem, not in general.
2. **Should `AlleleLikelihoodCalculator` live inside `ReadLikelihoodSnarlCaller`?**
3. **What folder layout follows**, given folders are decided.

**Preliminary probe, already done — start from this, not from zero.** Header-level dependencies
between the caller and the four subsystems adamnovak names:

| | mentions `ReadLikelihoodSnarlCaller` | includes from the read-likelihood family |
|---|---|---|
| `read_phasing.{hpp,cpp}` (385 lines) | **0** | none |
| `regenotype.{hpp,cpp}` (755 lines) | **0** | none |
| `linkage_model.{hpp,cpp}` (3,525 lines) | **0** | none |
| `anchor.{hpp,cpp}` (1,070 lines) | **0** | `site_read_source.hpp` |
| `read_likelihood_caller.hpp` | -- | includes **none** of the four |

So the preliminary answer to question 2 is **"not wedded at all"** — the subsystems do not name the
caller and, `anchor` aside, include nothing from its family. The dependency runs one way and the
composition must therefore happen in the `.cpp` files or in `call_main.cpp`.

**That is a header-level probe and is not sufficient.** The review must also check:
- shared *data types* (`CallInfo`, `AlleleReadLikelihoods`, `SnarlTraversal`) — decoupling by class
  name can hide coupling by struct
- where composition actually happens, since neither side includes the other
- whether `linkage_model` at 3,525 lines is one subsystem or several wearing one name
- what else in the flat `src/` belongs in the same folder, given the layout serves 297 files

**Candidate layouts to evaluate, not to assume:** `callers/read_likelihood/` with the subsystems
inside; or `lib/call/` + `lib/call/subsystems/`. The evidence above mildly favours subsystems as
*siblings* of the caller rather than nested inside it, since they do not depend on it — but that is
exactly what the review is for.

**Do not start 2d before this lands.**

### 2d. Structural reorganisation

Follows from 2c.

---

## Not to reopen

- **Base-level WFA** in place of the node walk: >24x the walk's CPU, unfinished on chr20 ONT;
  BiWFA refuses ends-free; the penalty transform costs ~10x. See `exact-walk.md`.
- **Forcing perfect-match pairings** as an exactness-preserving bound: not exact here, because the
  free flank breaks the global-alignment argument. It IS accuracy-neutral (-0.00015) and the
  cheapest CPU of any variant tested, so it remains a legitimate *performance* option -- but it is
  a different approximation from the band, not a strict improvement.

---

## Part 3 — the residual anchor `gqn` vs VCF `GQN` disagreement (102 snarls, 0.089%)

Mostly fixed in vg `c73fdedc7`; what is left is structural and is parked here rather than
pretended away.

### Where it stands

chr20 ONT, anchor `gqn` column joined to the VCF `GQN` field on the snarl ID:

| | before (`v6`) | after (`v7`, shipped) |
|---|---|---|
| joinable snarls | 114,625 | 114,625 |
| `gqn` != `GQN` | 8,149 (7.109%) | **102 (0.089%)** |
| sign disagreements | 4,232 | **37** |
| sign flips with \|value\| >= 0.25 | 1,548 | **16** |

### Why the last 102 survive

The VCF's `GL` is built by iterating **site** genotypes and mapping each through `site_to_scored`
into `genotype_lls` (`src/read_likelihood_caller.cpp:452-466`). So its genotype space is the one
over the **emitted** alleles. The fix already restricts the anchor's `best_other` to
`{ref_trav_idx} ∪ settled genotype`, which is the emitted set *for almost every record* -- and that
is what took 8,149 down to 102.

It is not the emitted set when the **symbolic layer collapses two traversals onto one ALT**: two
distinct traversals that spell the same ALT sequence become a single site allele, so the VCF's
genotype space is *smaller* than `{ref} ∪ called`. The residual disagreements read VCF-more-positive,
which is the signature of exactly that -- a smaller competitor set gives a larger margin.

Reproducing it needs the traversal-to-ALT map, which is built inside `emit_variant`, and anchors are
collected one line **before** `emit_variant` (`graph_caller.cpp:6401`) because that is the single
place the settled genotype and the per-read evidence coexist.

### Options, none started

1. **Leave it, documented.** 0.089%, and the sign is right on all but 37 snarls. The anchor file's
   own header already says "Anything else about the site is in the VCF, joinable on the snarl
   column" -- so the honest position is that `gqn` is a convenience copy and the VCF is normative.
2. **Collect anchors after `emit_variant`.** Then the allele map exists. Needs care: `emit_variant`
   hands the `CallInfo` on to `update_vcf_info`, so the evidence's lifetime has to be checked, and
   the ordering comment at `graph_caller.cpp:6392-6400` exists for a reason.
3. **Expose the emitted allele set** from the allele-map builder, keyed by record, and have anchor
   collection consult it. Only works if the map can be built before the anchors, which is the same
   ordering problem as (2).
4. **Replicate the collapse at anchor time.** Duplicates a rendering decision in a second place;
   rejected on principle unless (2) and (3) both prove impossible.

### Gate 0 - ANSWERED 2026-09-14, and it changes the recommendation

The gate was: characterise the 102 -- multi-allelic, nested, or SV? If they are a class a consumer
would filter out anyway, option 1 is the answer and this closes.

**There are 579, not 102, and the extra 482 are a class the earlier count silently dropped.** The
join was re-run on an anchors file and a VCF from *one* `vg call` invocation (`work/ont-preset/p3.*`,
vg `3bef00e05`) -- the pair on disk before this was v7b anchors against a v6 VCF, and although that
VCF turns out to be byte-identical to the fresh one, that was luck rather than method.

| chr20 ONT, 115,107 joinable snarls | |
|---|---|
| anchor `gqn` != VCF `GQN` | **579 (0.503%)** |
| -- VCF is `.`, anchor is a number | **482 (83.2%)** |
| -- anchor is `.`, VCF is a number | **0** |
| -- both numeric, differing | 97 (16.8%) |

The 482 are the 653 records whose VCF `GQN` is `.` *and* that join to an anchor: every single one
has a number on the anchor side. It is not a scattering, it is the whole class, and it is
one-directional. Split by why the VCF blanked:

| | | |
|---|---|---|
| **(b)** GT is diploid, so the patch ran but `achievable_phred` could not be recovered | 324 | 67.2% |
| **(a)** GT is haploid or half-missing, which the patch skips outright | 158 | 32.8% |

**These are two different defects and they point in opposite directions.**

**(a) The VCF cannot express a post-linkage `GQN` for a haploid genotype, and the anchor can.**
`apply_linkage_change` computes the re-derived `GQN` only under `called.size() == 2`, clearing
`called` the moment a GT field is `.` -- a deliberate choice, with the comment "a haploid record's
GL is indexed by allele, and conflating the two orders is how a plausible wrong number gets
written". So it blanks. But these records have real margins to report: `GQI` 116, 3 and 212 on the
first three, against anchor `gqn` 0.821, 0.094 and 0.871 computed from the live likelihoods, where
ploidy is not in doubt. **Here the anchor is right and the VCF is losing information it holds.**

**(b) The anchor reports a PRE-linkage margin on a record linkage moved.** `anchor_gqn_for` returns
a NaN meaning "use the sweep's value" for two quite different situations: linkage left the call
alone, where the sweep's value is still correct, and linkage moved it but the post-linkage margin
could not be computed, where the sweep's value is the margin of the genotype linkage moved *away
from*. The second is the defect already fixed once on this branch -- "right magnitude, WRONG SIGN"
-- surviving in the fallback path. The magnitudes are small here (0.000-0.002, because recovery
fails precisely when the pre-linkage `GQN` quantized to `0.000`), so nothing reads dramatically
wrong, but the semantics do.

The remaining **97** are the class this Part was originally written about: 44.3% carry the
symbolic-collapse signature (`AT` lists more traversals than `ALT`+1), 36.1% are multi-block, 34.0%
are sign flips, median |delta| 0.207 with 41 at or above 0.25. Direction is mixed -- 39
VCF-more-positive against 58 anchor-more-positive -- so a smaller VCF competitor set is part of it
and not all of it.

**Verdict: not option 1.** (a) and (b) are each a bounded fix in one function, and neither needs the
traversal-to-ALT map that made options 2-4 expensive. The 97 stay parked.

### Implemented

**(a) `apply_linkage_change` learns the haploid GL layout.** It computed the re-derived `GQN` only
under `called.size() == 2`, and cleared `called` on any `.` field. It now drops `.` fields rather
than abandoning the record, and decides the GL layout **by the GL's own length** -- n for haploid
against n(n+1)/2 for diploid, which is unambiguous at two or more alleles -- rather than guessing
from the genotype's shape, which is what the original comment was right to refuse. Predicted from
the existing VCF before building: it supplies a value for **all 159**, every one passing the length
check with the allele index in range.

The diploid branch deliberately keeps its older, looser condition. Under `--atomize-blocks` one
snarl emits several records that **share its snarl-level GL** while each carries only its own
block's ALTs, so `gl.size()` does not match that record's allele count and a length gate there
would newly skip every multi-block record. That is a separate problem -- 35 of the 97 remaining
disagreements are multi-block -- and not this change's to fix.

**(b) `anchor_gqn_for` stops falling back past the moved check.** One NaN used to mean both
"linkage left the call alone, so the sweep's value stands" and "linkage moved it but the margin
could not be recomputed". The first is still a fallback; every failure after the moved check now
blanks. Side effect worth knowing: `--anchors-min-gqn > 0` drops NaN rows (`anchor.cpp:368`), so
those 324 snarls now fall out of a filtered anchor file instead of passing it with a stale number.

Three TAP tests: that linkage never *blanks* a haploid record's `GQN` (under the existing
`HAP_CHANGED` antecedent, so it cannot pass vacuously), that the nested fixture offers snarls to
join on at all, and that no anchor reports a `gqn` where the VCF reports none.

### Measured, chr20 ONT, vg `6a8b10e99` + the fix

| | before | after |
|---|---|---|
| joinable snarls | 115,107 | 115,107 |
| anchor `gqn` != VCF `GQN` | **579 (0.503%)** | **97 (0.084%)** |
| -- VCF `.`, anchor a number | 482 | **0** |
| -- anchor `.`, VCF a number | 0 | **0** |
| -- both numeric, differing | 97 | 97 |

**Genotypes are untouched**: `POS`, `REF`, `ALT` and `GT` are byte-identical across the two runs,
which is what a quality-field fix has to be able to say. The residual 97 are the same 97 -- same
median |delta| 0.207, same 41 at or above 0.25, same 33 sign flips -- so nothing was traded.

**And the join was seeing 8% of the problem.** Counting what each fix actually moved, rather than
what the join could see:

| | | |
|---|---|---|
| VCF records that gained a `GQN` | **243** | 142 have an anchor row, 101 do not |
| anchor rows whose stale `gqn` became `.` | **4,005** | **340** have a VCF line, **3,665** do not |

The 3,665 are nested or off-reference snarls -- they emit anchors but no VCF line, so they had
nothing to disagree with and were invisible to every count in this Part. They were carrying a
pre-linkage margin for a genotype that had been moved, exactly like the 340 that were visible.
That is the argument for fixing it at the source rather than reconciling two columns: the column
with the error is the one a consumer reads on its own.

**What it costs a consumer filtering on `gqn`.** `--anchors-min-gqn > 0` drops a NaN row
(`anchor.cpp:368`), so those 4,005 now fall out of a filtered file. Their *stale* values were low,
which is expected -- the scale recovery fails precisely where the margin rounds to `0.000`:

| their old, stale `gqn` | |
|---|---|
| median | 0.091 |
| would have passed `--anchors-min-gqn 0.10` | 1,907 (47.6%) |
| would have passed `--anchors-min-gqn 0.25` | **423 (10.6%)** |
| would have passed `--anchors-min-gqn 0.50` | 101 (2.5%) |

So at a realistic threshold the practical loss is **423 snarls out of 172,340**, and those 423 were
passing the filter on a number that described a genotype the site is not reporting. The rest would
have been filtered out anyway.

### The obvious next increment, scoped but not done: retain `achievable_gap`

Blanking those 4,005 rows is honest but it is not the best available answer. The margin *can* be
recomputed on every one of them; what cannot be recovered is the **scale** to divide it by.
`anchor_gqn_for` reconstructs that scale as `GQI / GQN` from two quantized values, and that fails
exactly when `gq_fraction` rounds to `0.000` -- which is why the failures cluster on tiny margins.

But `achievable_gap` is computed directly in `read_likelihood_caller.cpp` (~line 230) and then
thrown away; only the ratio survives on the `CallInfo`. **Retaining it as one `double` would let
both the anchor and a post-linkage VCF normalise exactly**, with no quantization round trip, and
would turn those 4,005 `.`s into real values -- which matters, because a consumer filtering on
`--anchors-min-gqn` drops a NaN row.

Not done here, and it is not free: a new `CallInfo` field has to cross the ploidy swap in
`run_deferred_descent`, which copies only the fields it knows, and that has silently dropped a
field before. Gate it with the two-path byte-identity check that caught it last time.

### Noticed while doing it, not done: the anchor file has no provenance line

The header records `#graph`, `#reads`, `#sample`, `#mismap-min`, `#sites` and `#filters` -- enough
to reproduce the *inputs*, and nothing about the binary. So a consumer holding a `v7` anchors file
cannot tell whether its `gqn` column predates this fix, and neither can we: the format is unchanged
(same columns, same sentinel), so there is correctly no version bump to distinguish them by.

A `#vg-version` header line would close that, and it is three lines of code. It is deliberately not
in this change: the TAP suite asserts the exact sorted set of header keys, so adding one is a small
format change with its own test to update, and it belongs with whatever else the header should
carry rather than being smuggled in beside a `gqn` fix.

The chr20 deliverables rebuilt against the fix are labelled `v7c` -- a **run** label. The format is
still 7 and `scripts/check_anchors.py` is unchanged.

### Do not regress

`--anchors-min-gqn` and the `.` sentinel now mean different things than in v6: `.` is NaN ("no gap
to normalise"), a negative number is a value. `scripts/check_anchors.py` refuses v6 outright and
asserts `gqn` is `.` or within [-1, 1]. Any change here must keep those two distinct -- collapsing
them is the original bug.

---

## Part 4 — evidence-driven anchor slots: split homozygous, collapse unconfident heterozygous

**The idea, in one line.** A site has two anchor slots iff its reads are *confidently partitioned*
into two haplotypes -- whatever its genotype says.

Today slot count is hard-coded by the genotype: het gets two slots, hom gets one
(`src/anchor.cpp:347-355`), and confidence only ever *filters afterwards*. That is backwards in both
directions. A homozygous site with well-phased reads carries real haplotype information and throws
it away; a heterozygous site whose reads cannot tell the alleles apart emits two slots that are a
coin flip. **58.5% of chr20 anchor sites are single-slot** (196,400 of 335,675), and every one is a
break in the haploid run an assembler is trying to build.

### The per-read score

Mostly already built. `ReadLambda::lambda` (`src/regenotype.hpp:66`) is a per-read, cross-site
log-odds, summed over the het sites a read covers and signed by the settled frame:

    site_read_log_odds(q0, p) = log( (p*q0 + (1-p)/2) / (p*(1-q0) + (1-p)/2) )

`q0` is the read's allele responsibility split at that site, `p` the probability it came from one of
the two settled haplotypes at all. It is keyed by **read alone, not (read, block)** -- deliberately,
with the comment *"the hom case needs a lookup that has no block to offer."* The structure was built
for this.

Two properties it already has: **MAPQ is inside it**, via the `(1-p)/2` escape floored by
`--mismap-min`, so a mismapped read cannot reach a confident strand -- damped rather than dropped.
And it is **calibratable**: raw lambda runs into the hundreds because reads are treated as
independent and are not, so `calibrated_log_odds` applies a fitted temper `tau ~ 0.07-0.08`.

**What is missing, and is the substance of this part.** Per-site *genotype* confidence does not enter
the sum. Verified: nothing in `regenotype.cpp` or `read_phasing.cpp` consults `gq` or `gq_fraction`.
A marginal het call and a rock-solid one contribute identically so long as the read discriminates.

There is partial self-correction -- a site wrongly called het gives `q0 ~ 0.5`, so its term is ~0 and
it contributes nothing. That is the *benign* failure. The one it does not cover is a site
confidently het and **wrongly phased**, which contributes a confident wrong-signed term. That is
switch error, and weighting each site's term by the confidence its genotype is right is the damper
for it.

### The decision, both directions

For each site, per read, a posterior over {strand 0, strand 1}:

- **hom**: from cross-site lambda alone (the site itself has no allele signal -- both haplotypes
  carry the same allele, by definition)
- **het**: from this site's allele evidence *and* cross-site lambda, with **leave-one-out** so a site
  is not used to phase itself. `read_loo` (`src/regenotype.cpp:255`) already does that subtraction;
  at a hom site there is nothing to subtract.

Then: **split iff confidently partitioned, collapse iff not.** Reads that fail the per-read gate are
**pooled, not dropped** -- see below.

### Pooling, not dropping, is the point

Every existing confidence control here is subtractive: `--anchors-min-q` drops reads,
`--phase-min-q` drops sites, `--anchors-min-gqn` drops anchors. Collapse is *graceful degradation* --
the pin survives, its reads survive, and only the haplotype claim is withdrawn.

For an assembler that is strictly better: connectivity without phase is useful, a missing pin is not.
**A wrong split is worse than no split**: a break is expected, a confident-looking chimeric haploid
run is not. So the gate is conservative and "pooled" stays a representable state.

### Cautions

- **The score saturates.** Agreement tops out at ~95-97% however large |lambda| gets, so any phred
  has a real ceiling near 13-15 -- which is about where the existing anchor score already saturates,
  `phred(--mismap-min)` = 13.01 under `--preset ont`. A naive threshold will be overconfident.
- **`tau` is fitted only when `--regenotype` is on.** Splitting must either require it or fit its own.
- **Both slots of a split hom site carry the SAME allele**, which nothing in the writer expects.
  `site_slot_weights` length-weights by `slot_allele`, and the argmax tie-break compares allele
  indices -- both need checking with identical alleles.

### Evaluation

`scripts/tier2/anchor_purity.py` is the instrument: purity, yield and calibration against reads of
known haplotype origin, simulated separately from two haplotype paths. Its calibration check is the
one that matters -- *"an uncalibrated score is worse than none, because a downstream filter would
trust it."*

**One trap in using it.** It currently *excludes* single-slot anchors from purity, because their
majority purity is ~0.5 by construction. Split them and they enter the numerator, so a naive
before/after purity comparison is not like-for-like. The comparison must either hold the site set
fixed or report the two populations separately.

### IMPLEMENTED (vg `788470620`), with the gate result

Gate 0 passed and is now a permanent per-run self-check. Held out on chr20 ONT, leave-one-out, over
3.12M reads at heterozygous sites:

| | agreement | share of reads |
|---|---|---|
| **overall** | **94.69%** | -- |
| confident (\|lambda\| >= 2) | **95.47%** | 94.6% |
| not confident | 81.20% | 5.4% |

Against a 50% chance baseline, and 95.47% is exactly where strand confidence is documented to
saturate -- the mechanism is at its ceiling, not under-tuned.

**That measurement changed the design.** The original proposal was a per-read confidence gate.
Confident reads agree 95.5% and unconfident ones 81.2%, but only 5.4% of reads are unconfident, so
gating per read buys **0.78 points of purity for a 5.4% yield loss**. The gate is therefore
SITE-level -- two confidently placed reads on each strand -- and every read at a split site is
placed rather than dropped. Pooling beats dropping, as the principle above says.

Result, `--anchors-hom-split` on chr20 ONT:

| anchor pins | before | after |
|---|---|---|
| 1 slot (no haplotype) | 196,400 (58.5%) | **45,433 (13.5%)** |
| 2 slots | 139,275 (41.5%) | **290,138 (86.5%)** |

79,962 homozygous sites split, 9,014 left collapsed. Flag off by default and byte-identical when
off; the VCF is unchanged when on; `check_anchors.py` passes.

**Still not validated against truth.** The ~95% figure is the accuracy of the INFERENCE, measured
where an answer exists. Whether a 95%-pure haplotype label helps an assembler more than the break it
replaces is unmeasured, and needs the simulated two-haplotype read set that does not exist. Until
then the flag stays off by default and the file declares `hom-split=on` with a note naming
equal-alleles-across-two-slots as the tell.

### Not done

- **Collapsing unconfident heterozygous sites** -- the reverse direction. The site-level machinery
  is now in place (`AnchorParams::phase_min`, `phase_min_side`) and the same decision applies, but it
  is unimplemented and unmeasured.
- **Weighting each site's lambda term by its genotype confidence.** The one-line change is
  `src/regenotype.cpp:120`, multiplying `site_read_log_odds` by a new `PhaseSite::confidence`. NOT at
  `regenotype.cpp:141`, which would weight the forward sum but not the leave-one-out subtraction and
  make every site partially confirm itself. The natural source is `gq_fraction`, whose -1 sentinel
  must be translated to 1.0 and which already carries the explained-share discount, so multiplying by
  `explained_share` again would double-count.

### Order of work

1. Thread a populated `LambdaTable` to anchor-emission time (inputs are already live members).
2. Weight each site's term by its genotype confidence; re-measure switch error to confirm it helps.
3. Split hom sites behind a flag, default off; pooled slot for unassignable reads.
4. Collapse unconfident het sites, same flag family.
5. Format bump, `check_anchors.py` update, TAP and unit coverage.
6. Purity/yield/calibration against the known-origin harness, with the exclusion trap handled.

### Gate 0 -- held-out, and it needs no simulation

**The purity harness cannot run.** `anchor_purity.py` needs reads whose names encode their true
haplotype (`h1_`/`h2_`); no such read set exists anywhere in the repo, nothing builds one, and the
script has never been run -- there is no recorded output of it and no baseline to compare against.
Every read source on disk is real HG002. So the obvious validation is unavailable without first
simulating from two haplotype paths, mapping and calling: a pipeline in its own right.

A better gate is available from data already in hand, and it tests the mechanism directly:

> At **heterozygous** sites, partition the reads by cross-site lambda ALONE -- leave-one-out, so the
> site's own allele evidence is excluded -- and compare against the allele-based partition the site
> actually makes. Report the agreement rate.

A het site has an allele-derived answer. Lambda-only partitioning is exactly what a hom split would
have to rely on. So this measures the accuracy of the mechanism on held-out ground truth, with no
simulation, and it is the same comparison the hom case will make blind.

If agreement is poor, the hom split cannot be trusted at hom sites either, and steps 3-6 stop. If it
is high, the split is justified and the same number calibrates the confidence threshold.

Two known confounders to control: reads whose lambda comes from a single site have nothing left
after leave-one-out (`ReadLambda::sites == 1`), and `multi_block` reads span a phase break. Both
must be excluded and counted, not silently folded in.

### What the shipped split actually did, measured on chr20 ONT

Anchor census after `--anchors-hom-split` (vg `3519d4980`, post-fix):

| diploid homozygous snarls WITH a VCF line (35,593) | | |
|---|---|---|
| every pin split into two slots | 31,103 | **87.4%** |
| one pin split, the other single-slot | 113 | 0.3% |
| no pin split -- genuinely still collapsed | 4,377 | 12.3% |

**Splitting does not require a VCF line, and nested sites are the majority of it.** Of the 79,137
snarls showing a split (two slots carrying the SAME allele -- the tell), only **39.4% have a VCF
line**; **60.6% are `reported_inline` or `no_reference`** records, which are anchored but never
written as a line. The per-read lambda comes from the het sites the reads ALSO cross, not from the
site being split, so a nested homozygous site is partitionable exactly like any other. These are the
off-reference sites the code calls "where an assembler most needs help", and the same population the
phase-ordering bug was mislabelling.

A detail that falls out right rather than by design: `apply_read_phasing` builds `phase_sites` from
`records_for_render()`, which EXCLUDES `reported_inline` and `no_reference`. Those sites therefore
have no `PhaseSite`, the leave-one-out lookup finds nothing, and nothing is subtracted -- correct,
because a site that contributed nothing to lambda has nothing to leave out.

### The remaining single-slot pins, and a correction

45,571 pins are still single-slot:

| | pins | share |
|---|---|---|
| no VCF line (nested inline / off-reference) | 17,167 | 37.7% |
| diploid heterozygous | 14,906 | 32.7% |
| diploid homozygous -- failed the split gate | 8,415 | 18.5% |
| nested haploid (`1\|.`) | 5,083 | 11.2% |

**The heterozygous row is a PIN count and was first described as if it were a site count, which
overstated it by an order of magnitude.** At site level, of 76,861 diploid het snarls with anchors:
81.8% have both slots at every pin, **16.7% have both at one pin and one at the other** -- the
partition is present, an end pin simply held too few reads for the second slot -- and only **1.5%
(1,128 sites) have no pin with both slots**, which is the only genuinely concerning group. The
homozygous population shows the opposite shape (0.3% partial), because the split is decided once per
site and applied at both pins, whereas the het partition depends on per-pin read placement.

### Two bugs found by explaining the code, not by running it

Both were live in the first shipped version and are fixed in vg `181266514`:

- `best_slot = lo > 0.0 ? 0 : 1` sent every read with **no opinion** (`lo == 0.0`) to slot 1 --
  a haplotype claim with nothing behind it, made systematically in one direction. **48,686 read
  placements** on chr20. They are now in neither slot and counted.
- the `side0`/`side1` count included reads pinned at neither end, which the placement loop drops, so
  a site could qualify on evidence that never reached the file. Cost: 415 sites, 0.5%.

The no-opinion rate is **0.38% at hom sites against 0.19% at het** -- twice as high, which is the
expected direction, since a homozygous site is likelier to sit where few het sites are near enough
to phase from.
