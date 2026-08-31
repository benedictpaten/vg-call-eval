# Scoping the `resolve_generation` rewrite

**Ask:** make the structure as simple as possible, starting from "rewrite `resolve_generation` to be
subtree-scoped so the depth-first recursion becomes available".

**Finding that reorders it:** the subtree scoping and the recursion are not where the complexity is,
and after the enabling work neither buys anything. The complexity is that
`LinkageCollector::resolve_generation` is **869 lines** doing four separable jobs. That is the thing
to fix, and it is independent of the recursion.

## What is actually big

| | lines |
|---|---|
| `LinkageCollector::resolve_generation` | **869** |
| `FlowCaller::call_snarl_internal` | 799 |
| `FlowCaller::run_deferred_descent` | 499 |
| `VCFOutputCaller::emit_block_records` | 473 |
| `LinkageModel::window_posteriors` | 217 |

`resolve_generation` breaks down as: ~100 preamble, ~75 contig-chain build, ~244 per-parent group
build, ~319 decode loop, ~90 reporting, ~40 sort and return.

## Part A -- remove the O(generations x sites) terms

Three places do work proportional to ALL recorded sites on every call, rather than to the sites that
call decodes.

| | cost on chr20 | fix |
|---|---|---|
| entry scans in `resolve_generation` | 25 ms | **done** -- `by_generation` indexes them at `record()` time |
| `pinned_phase` rebuilt from the whole accumulated phasing | 49.6 ms | make it a collector member, updated at the single site that produces a `PhaseCall` (`phasing_out->push_back`, one call site) |
| the caller's `settled` map, same shape (`graph_caller.cpp:5128`) | ~50 ms, estimated by analogy, not measured | same |
| the caller's per-generation scan of all `pending` | 213k visits, sub-ms | bucket `pending` by generation once |

**Total saving ~100 ms of a 170 s run: 0.06%.** So Part A has no performance case on its own. Its
only justification was enabling Part C, and Part C is not recommended -- see below. What survives on
its own merits is the `pinned_phase` member, because it *deletes* a rebuild loop and a per-call map
rather than merely speeding one up. ~15 lines removed. The other two are not worth doing.

## Part B -- decompose along the seam that is already there

The decode loop already consumes a work list through three parallel arrays, and both builders
already produce exactly those. That is the seam:

```cpp
struct ChainSet {
    vector<vector<size_t>> chains;             // entry indices, one list per decode unit
    vector<const vector<double>*> context;     // entering message per chain, or null
    vector<size_t> phase_set;                  // inherited phase set, or SIZE_MAX
    std::deque<vector<double>> messages;       // OWNS what `context` points at
    unordered_map<size_t, NestedPlacement> placement;
};

ChainSet build_contig_chains(...);                     // ~75, generation 0 only
ChainSet build_parent_groups(...);                     // ~244, generation > 0
size_t   decode_chains(ChainSet&, generation, ...);    // ~319, uniform over both
void     report_generation(...);                       // ~90

size_t resolve_generation(size_t generation, bool last, vector<PhaseCall>* out) {   // ~40
    ChainSet cs = generation == 0 ? build_contig_chains(...) : build_parent_groups(...);
    size_t moved = decode_chains(cs, generation, out, last);
    report_generation(generation);
    if (out != nullptr && last) sort(*out);
    return moved;
}
```

Net lines: roughly zero. It is code motion plus a struct. What it buys:

- **`generation` stops being ambient.** It is read 24 times inside the decode loop today, mostly as
  per-entry liveness (`e.generation == generation`) and clamp (`e.generation < generation`) tests.
  After the split it is one parameter of one function.
- **The message lifetime becomes structural rather than conventional.** `deltas` must outlive the
  decode because `chain_context` holds pointers into it; today that is maintained by declaring it at
  function scope and remembering why. Getting it wrong is not hypothetical -- it was a
  use-after-free earlier in this rebuild that never crashed, and whose freed-but-intact memory
  scored BETTER than the correct code on one contig. Owning the messages inside `ChainSet` makes
  that class of bug unrepresentable.
- Two builders that produce one thing can be read side by side, which is how the drift between the
  ploidy paths was found.

**Gate:** byte-identity on chr20 and chrX, VCF and mosaic, plus `vg test` and `18_vg_call.t` alone.
It is pure motion, so byte-identity should hold exactly and any difference is a real bug.

**Risk:** the decode loop closes over a good deal of state (arenas, `model`, counters). Extracting it
means either passing that state or making these member functions -- member functions, since they all
already live on `LinkageCollector`. The parameter list is the hazard: the last time this code had a
long positional parameter list of same-typed arguments, a dropped argument compiled with everything
shifted, three times. `ChainSet` plus two named parameters is the shape that avoids repeating it.

## Part C -- the subtree recursion: not recommended

After Part A the total work is already O(sites); subtree scoping changes the ORDER, not the cost. And
the order it would change to is not simpler:

- Generation 0 must remain one whole-contig pass regardless -- that is what the Li-Stephens chain
  is. So a recursion is a special case plus a walk, where today it is one loop.
- The caller's ploidy revision is naturally per-generation: settle a generation, then revise the
  children it just made reachable. Depth-first splits that into a per-subtree interleave with no
  corresponding gain.

Cost if wanted anyway: ~5 lines in the caller once A and B are done. Cost today, without A:
`resolve_generation` does O(all sites) per call, so one call per subtree root -- thousands on chr20
-- is order 20 s against 75 ms now, and against an 8 s linkage layer.

## Recommendation

Do **Part B**, and take the `pinned_phase` member with it because it deletes code rather than merely
speeding it up. Skip the rest of Part A and all of Part C, and record at the generation loop that
depth-first is available and deliberately not taken.

Roughly 400 lines moved, ~30 added, ~15 deleted, one 869-line function becoming four. One gate run
(chr20 + chrX + TAP) covers it, and byte-identity is the whole test.
