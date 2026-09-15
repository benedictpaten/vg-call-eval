<!-- DRAFT, not posted. A reply to adamnovak's review comment on vgteam/vg#4990
     (comment 5671460225). Posting a PR comment is the author's call, not the agent's,
     so this sits here until Benedict sends it or edits it. -->

Thanks — that's a useful read and it changed what I did rather than what I'd say about it.

I've written the review up as `doc/read-likelihood-architecture.md` (commit `eefa809a1`) so the
numbers are in the tree rather than in a comment. Short version, and one place where your premise
and mine were both slightly off:

**"How wedded to `ReadLikelihoodSnarlCaller` is e.g. the regenotyping system?"** — the code isn't;
the data is. None of `read_phasing`, `regenotype`, `linkage_model` or `anchor` names the caller,
and `linkage_model` and `read_phasing` include nothing from vg at all — not one `#include "..."`
between them. But every `PhaseReadEvidence` regenotyping consumes is built in exactly one place,
the `AlleleReadLikelihoods` builder. So the contract is narrow and replaceable while the supply is
singular, which I'd argue is the coupling you want: it takes a struct, not a caller.

**`AlleleLikelihoodCalculator` inside `ReadLikelihoodSnarlCaller`** — it's already private to the
caller in everything but syntax: no production file outside `read_likelihood_caller.hpp` includes
it. I'd rather not nest it, for two reasons. It's a one-method abstract base with a single
implementation, i.e. an extension point, and nesting an extension point inside its only consumer
is awkward; and three unit-test files, 180 `REQUIRE`/`CHECK` between them, construct it directly to
test the DP walks and the scoring window. Documenting the intent costs a line and buys the same
thing.

**The coupling neither of us named is the one that matters.** `VCFOutputCaller` is the base class
of four callers, and **39 of its 64 data members — 61% — exist only for the read-likelihood path**
(the linkage collector and its GBWT caches, `render_phases`, `phase_sites`, `render_lambda`, the
regenotyping state, the anchor writer, the mosaic state), along with 13 of its 43 methods. Three of
the four subclasses touch none of it. So the subsystems are tidy and the *composition* isn't — and
a folder move would leave that 61% exactly where it is. I've recorded it rather than started on it;
it's bigger than this PR.

**Folders vs enclosing widgets** — both, but less of each than it sounds. Counting namespace-scope
functions: `regenotype` 8, `symbolic_allele` 5, `anchor` 4, `read_phasing` 2, `linkage_model` 1,
and zero in the four that already have a class. Of those twenty, most are genuinely pure
(`phase_link(a, b, cap)` has no state to hold); the four worth an enclosing class are the ones that
thread the same `LambdaTable`/`PhaseSite`/params triple through every call. On layout, Benedict has
asked for folders, and `src/` earns them at 297 flat files and 206,800 lines. I'd suggest **one
flat `src/caller/`** rather than `caller/` + `subsystems/`: each subdirectory costs 16–34 lines of
Makefile plumbing, and "subsystems" asserts a hierarchy the dependency graph doesn't have —
`linkage_model` and `read_phasing` don't depend on the caller.

One thing that falls out and is independent of any of the above: **`linkage_model` is two
subsystems wearing one name** — ~2,167 lines of Li-Stephens HMM and ~1,275 of site store, and
`graph_caller` names `LinkageCollector` 29 times and `LinkageModel` only for a constant. That split
is the highest-value move in the whole reorganisation and needs no layout decision.

**"you know it's parsing half the CLI options again by hand, right?"** — yes, and it had already
drifted. Fixed in this PR. The option table now carries an owner per option and getopt's
`struct option` array is *generated* from it, so a flag can't be added to one and forgotten in the
other; the hand-rolled `getopt_long` prefix matcher is gone, because the parse loop just records
which `val`s getopt hands back and getopt had already done the resolution.

What the old list had accumulated, demonstrated on the shipped binary before I touched it:
`--mosaic-out` was refused without `--read-likelihood` while `--mosaic-patch-gaps`,
`--no-mosaic-patch-gaps`, `--no-mosaic-nested` and `--mosaic-break-unexplained` were all accepted
and silently dropped.

And it does the "for free" part you asked for. The owner enum carries `OWN_ANCHORS`, `OWN_MOSAIC`
and `OWN_REGENOTYPE` as well, so the 8 `--anchors-*`, 4 `--mosaic-*` and 7 `--regeno-*` options are
refused when their own subsystem isn't on — `--anchors-min-q 30` with no `--anchors-out` now says
so instead of doing nothing. Adding a subsystem is one row in a gate table plus one enum value.

Two details that aren't cosmetic. When *both* switches are missing the outer one is named, since
that's the one you have to fix before the other can matter. And the regenotype gate tests
`regenotype && read_phasing` rather than the raw flag, because `--preset ont --no-read-phasing`
disarms regenotyping *after* this check runs — that combination is documented and must not error,
so the gate has to know about it.

Seven TAP tests, including that an `OWN_CORE` option is still let through: a check that refused
everything would pass the obvious assertion and be useless.

Also: `regenotype.hpp` now opens with a plain-English sentence before the derivation, which is
what the "computer was high" line was pointing at. The derivation itself I'd keep — it names three
properties the unit tests assert, including that the whole correction is bit-for-bit zero at
`--regeno-temper 0`.
