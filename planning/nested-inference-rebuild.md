# Rebuilding nested calling around inference, not VCF

The deferred-calling machinery has grown by iteration until several mechanisms decide the same
facts and none of them can be checked against the others. This is the plan to replace it with the
smallest design that produces the same answer, organised around what is being inferred rather than
around what will be printed.

Authorised by the repository owner as a rebuild: large-scale removal and replacement is in scope,
and the measure of success is code removed, not accuracy gained.

## The algorithm, as the owner states it

For each top-level chain at ploidy 1 or 2:

1. Genotype the top-level chain as today, using linkage, giving chosen haplotype traversals (one or
   two) for each snarl in the chain.
2. For each parent snarl and its haplotype traversals: (a) symbolic allele comparison to decide the
   genotypes to REPORT, which is an emission concern only; (b) compute the copy count of each nested
   chain, which may exceed 2 -- count and log those.
3. Recursively genotype each nested chain at its copy count, returning to step 1.

## Settled decisions

- **Greedy, permanently.** A parent settles before its children and is never revisited. Child
  evidence never informs a parent. Today is already greedy in this direction, so nothing is lost.
- **A nested chain's phase is expressed in its parent's frame.** Strand identity composes up the
  tree. The owner calls this the essence of the inference.
- **The output is each strand as a sequence of recombinations of the GBWT haplotypes**, including
  recombinations inside nested chains. The mosaic is the product; VCF is a projection of it.
- **Windowing and phasing behaviour does not change.** Reimplement for clarity if it comes out
  cleaner, gated on byte-identity. No parallelism inside a chain's decode: independent windows put
  the two decodings in unrelated frames, and reconciling them after the fact is a vote that gets
  noisier as the panel grows -- it would manufacture a switch at every seam, which is the artifact
  the seam pin exists to prevent. Cutting at a saturated gap is exact but `rho` reaches 1 only
  asymptotically, so there are no such gaps to cut at.
- **Copy count above 2 is counted, not represented.** The hidden state is an ordered PAIR of panel
  haplotypes; triples would be 42,875 states at m=35, a third `rho`, and a different GL layout.
- **Both conditionings are testable from one binary**: a child conditioned on its parent's hard
  called pair, versus on the parent's posterior over pairs with the copy count still hard-called.

## What the measurements already settled

**Linkage between nested chains is worth nothing.** Decoding each nested chain alone against
decoding all of a parent's children together, on three contigs:

| | chr20 | chr6 | chr17 |
|---|---|---|---|
| ALL / SNV / Indel F1 | 0 | +1.8e-6 | -1.3e-5 |
| SV F1 | -3.3e-4 | 0 | +8.4e-4 |

Every movement is one call, and the signs disagree across contigs. chr17 carries 273 multi-crossing
children where the other two carry none, so the case with the best chance of mattering does not.

**Therefore sibling chains have no transition between them, and their relative order and spacing are
unobservable.** That retires inter-chain ordering, orientation and distance entirely.

**Read retention across the recursion already exists**, in the read source rather than the caller:
descent is synchronous on the fetching thread and the cache is keyed by node-ID window, so a child's
reads hit the parent's resident window. Going back to the reads after settling was measured at
**+57.8%** (22.4M against 14.2M) because a scattered second pass re-fetches whole windows to reach a
handful of sites.

**A chain-unit read query would break that.** A top-level chain's node-ID span routinely exceeds the
4,096-node window, and a straddling query bypasses the cache and caches nothing. Scoring stays
snarl-unit inside a chain-unit loop.

**"Score only what the reads touch" buys 3-5%.** 96.3% of nested chains have discriminating reads
(94.2% with off-reference chains admitted). Worth doing as a cheap skip at descent -- it avoids the
`panel_alleles` query too -- but it cannot drive the architecture.

## Distances

- **Top-level chain**: reference positions, unchanged. The reference is a real coordinate there.
- **Within a nested chain**: the length of the settled traversal through the preceding snarl. Chain
  snarls are adjacent by construction, so this is intrinsic and needs no reference. Per-strand
  distances genuinely differ here, because the two haplotypes take different-length paths -- unlike
  the between-chain case, where 90% were exactly the reference difference.
- **Between nested chains**: none. There is no transition.
- **Parent to child**: the child's offset along the parent's settled traversal, per strand, from
  `Entry::frame_offset`. Linkage should decay across a large parent snarl, so containment is not the
  same as zero separation. To be tested against a no-transition arm and a fixed-distance control,
  which separates "decay matters" from "any transition matters".

## The plan

Steps V1-V3 are validations the plan itself needs before it can be executed.

| | step | gate |
|---|---|---|
| V1 | **DONE, and it refuted the figure.** 250 sites was said to retain 55% pair correlation, needing ~1,237 for 0.05. Measured on chr20: mean -log P(no switch) = 0.0121 per step, so **4.8% retained at 250 sites and 247 needed for 0.05**. The shipped margin is almost exactly right. (Per strand rather than per pair it is 22% retained and ~494 needed; the pair is what the window decodes, so 250 is correct on the measure that matches the model.) **Step 4 is therefore deleted.** | measured; instrument arm byte-identical |
| V2 | **DONE, settled by reading.** Descent already visits EVERY child of every snarl: a child no called allele reaches is not skipped, it sets `retain_only` and descends anyway, precisely because linkage may move the parent onto it. So collect's enumeration is already broad enough for greedy recursion. The one narrowing is the reference gate (12,486 chains on chr20), which step 8 removes. | none needed |
| V3 | Restate step 5's gate. "Zero disagreements" is unachievable on current code because the `ploidy == 2 && parent_trav >= 0` population is a known defect. The gate becomes zero AFTER that fix, making it a prerequisite rather than an aside. | plan text only |
| 0 | Instruments: pin declines, read-evidence split, copy-count histogram | published, no logic change |
| 1 | Retire the stage-14 frame instrumentation, THEN delete `frame_reversed`, `frame_end`, `frame_total`, `n_reads`. The three frame fields are not dead as claimed: `frame_end`/`frame_total` are read in the cross-parent branch of that instrumentation, into a `frame_gap` the next line discards. `frame_offset` STAYS -- it is the parent-to-child distance. | byte-identical |
| 2 | Delete inter-chain ordering and spacing | byte-identical with per-chain grouping on |
| 3 | Split `panel_alleles` out of record time | byte-identical; expect a speedup |
| ~~4~~ | ~~Expand the margin~~ -- **deleted**, V1 shows 250 is correct | -- |
| 5 | `relate()` -- copy count, carrier and strand derived from the parent's settled pair -- landed as a CHECK | zero disagreements, after the V3 prerequisite |
| 6 | Replace the barrier with the recursion; delete generations, `PendingRecord`, `respecify` | byte-identical if 5 reads zero |
| 7 | Strand composition in the parent's frame; the two conditioning arms; the three parent-to-child distance arms | accuracy, chr20 and chr6 |
| 8 | Emission split: `chain_reported_inline` stops gating descent | record set identical in CHROM/POS/ID |
| 9 | The mosaic as the product: strands as recombination sequences, including inside nested chains | new output; no gate on old behaviour |

Steps 1-3 and 6 should be byte-identical. With step 4 deleted, **only step 7 moves numbers** --
the two conditioning arms and the three parent-to-child distance arms. Everything else in the
rebuild is gated on producing the identical answer.

## Expected size

Roughly 2,600-3,200 lines removed against 900-1,200 added, plus 31 functions, ~97 struct fields, 18
named concepts and ~103 counters. The three largest blocks counted directly:
`run_deferred_descent` 715 lines, the per-strand nested pass 816, the grouping and its vetoes the
rest.

## What this is not

It is not a performance or accuracy project. Everything measured today that survived was a deletion:
inter-chain linkage is worthless, and the off-reference population costs slightly more than it gives.
The value here is the code that stops existing, and that should not be read as anything else.
