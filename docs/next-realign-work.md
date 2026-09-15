# Parked work: re-fitting the ONT preset under `--realign`, and the PR #4990 review

Written 2026-09-14, after `--realign` shipped opt-in (vg `c60503ac3`, PR #4990). Nothing here is
started. Both parts are independent of each other and of whatever is being worked on now.

---

## Part 1 — the ONT preset was fitted against a walk that no longer runs by default

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

### Gate 0 — is this worth starting at all?

Run **only the `--gap-open 2` arm** first, one measurement, 6 minutes. The hypothesis predicts it
beats `--gap-open 1`. If it does not, the compensation story is wrong and the rest of the sweep is
unlikely to pay; stop and record the negative.

---

## Part 2 — PR #4990 review (adamnovak, comment 5671460225)

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

### Gate 0

Before doing any of this, characterise the 102: are they multi-allelic sites, nested sites, or SVs?
If they are a class a consumer would filter out anyway, option 1 is the answer and this closes.
One join, no vg run:

```
python3 - <<'PY'   # anchors v7 vs the VCF, on the snarl column
# see the check in this session: parse GQN from the VCF FORMAT field, gqn from the A rows,
# report snarls where they differ by more than 0.0015
PY
```

### Do not regress

`--anchors-min-gqn` and the `.` sentinel now mean different things than in v6: `.` is NaN ("no gap
to normalise"), a negative number is a value. `scripts/check_anchors.py` refuses v6 outright and
asserts `gqn` is `.` or within [-1, 1]. Any change here must keep those two distinct -- collapsing
them is the original bug.

---

## Still unparked

The **homozygous anchor split** design (splitting a hom site's single slot by implied read phase, to
lengthen haploid runs) was worked out in this session but is neither planned nor implemented here.
58.5% of chr20 anchor sites are single-slot, so it is the largest outstanding item.
