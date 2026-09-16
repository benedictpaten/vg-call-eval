# Plan: read assignment, the mapQ lever, and anchor output control

Supersedes the ordering in [phase-confidence-plan.md](phase-confidence-plan.md), whose first item
(split thresholds) is done. Written after a measurement pass that found a real defect and refuted
three of my own explanations for it, so the diagnosis step below is a gate, not a formality.

## Status: all phases done

| phase | outcome |
|---|---|
| 1 diagnose Λ = 0 | done -- 97.8% reached no phase site, 2.2% summed to nothing, 0 phase breaks |
| 2a coin | done -- 5,964 chr20 / 21,861 chr6 reads stop vanishing; see [read-assignment-results.md](read-assignment-results.md) |
| 2b chain gap | **not needed** -- phase sites 76,135 against 76,077 two-slot snarls, so the chain is not excluding het sites |
| 3 Λ == 0 site gate | done, byte-identical above q = 0 |
| 4 mapQ lever | done and **rejected**: see [mapq-for-phase-confidence.md](mapq-for-phase-confidence.md) |
| 5 output flags | dropped |
| 6 read_strand hoist | done, byte-identical |

## Done

`--split-min-q` fitted to **0.5** on chr20's phased-run length, confirmed on chr6 held out
(read-walkable run N50 2,130,682 -> 2,883,125 bp, the two arms' VCFs byte-identical), shipped in vg
`27ac3b8cd`. See [phased-run-lengths.md](phased-run-lengths.md).

## The defect this plan addresses

A read with no cross-site opinion is dropped at a **split homozygous** site (`src/anchor.cpp:492`)
but kept at het sites and at unsplit homozygous ones, because the normal responsibility loop never
consults lambda. Its trail through the file then reads present-absent-present, which downstream is
the signature of a deletion.

Measured on chr20, no-split against hom-split:

| | placements |
|---|---|
| (read, snarl) pairs lost | 46,663 |
| belonging to **5,964 reads that vanish from the file entirely** | 44,017 |
| belonging to reads still present | 2,646 |
| — **INTERIOR: present, absent, present** | **428**, on 146 reads |
| — terminal | 942 |
| — at a snarl with no VCF line | 1,002 |

The 428 sit in gaps of median 2,391 bp, p90 10,347, max 14,392, and **163 span >= 10 kb**. 88.8% of
the losses on kept reads are at snarls where BOTH slots survive, so they are genuine no-opinion
discards and not `min_reads` casualties.

## What Λ is, exactly

Per phase site the read touches (`src/graph_caller.cpp:6098`):

```
r0 = (1-e)*w0*rel(r,a0)     r1 = (1-e)*w1*rel(r,a1)
q0 = r0/(r0+r1)             p  = (r0+r1)/(r0+r1+e)
```

and the site contributes `log(a/b)` with `a = p*q0 + (1-p)/2`, `b = p*(1-q0) + (1-p)/2`
(`src/regenotype.cpp:103`). Λ is the sum over sites, tempered by τ ≈ 0.08. **Λ = 0 means no net
opinion**, and it is what `--split-min-q` thresholds.

Λ = 0 has distinct causes that the fix must tell apart:

1. **No phase site at all.** The read touches none. Nothing to know.
2. **Zero-valued contributions**, `q0 = 0.5`: the read matches both alleles equally -- the
   homopolymer case, where the two alleles differ by a length the read cannot resolve. Nothing to
   know.
3. **Exact cancellation** across sites: contradictory evidence. Nothing to know, arguably.
4. **The site was excluded from the chain.** Five guards at `src/graph_caller.cpp:6058-6084` drop a
   het site from `phase_sites`: no PhaseCall, `ploidy != 2` or `trav_first == trav_second`, no call
   info, no phase evidence, or a PhaseCall naming a traversal the matrix lacks. A read whose het
   sites are ALL excluded is absent from the lambda table entirely. **Something IS knowable and the
   chain is not using it**, and from inside the table this is indistinguishable from cause 1.
5. `multi_block`: the read spans a phase break, so its two halves label strands independently.
   Reported as 0 on both chr20 and chr6, but not structurally guaranteed.

**NOT a cause at a homozygous site: `sites == 0` after leave-one-out.** I listed it earlier and it is
wrong. Leave-one-out only subtracts when the site being asked about is itself in `phase_sites`, and
`phase_sites` excludes every record with `trav_first == trav_second` (`src/graph_caller.cpp:6065`) --
which is exactly what a homozygous site is. So nothing is subtracted and `sites` is not decremented.
Verified by three independent adversarial reads, 0/3 refuted; the codebase already states it at
`src/graph_caller.cpp:6344`: "for it the leave-one-out is automatic, because it contributed nothing
to Lambda in the first place." The value the split gate sees is the read's full tempered log-odds.

That is why a single global Λ for assignment is not a simplification: **it is already the behaviour.**
Removing the leave-one-out from this path would change no output. What it would change is legibility,
and that is worth something -- the shared accessor is what made me report a cause that cannot fire.

**Causes 1-3 are truly unassignable; a coin is right for them. Cause 4 is a chain gap; a coin would
throw away real evidence to paper over it.** That distinction is the whole design.

## Phase 1 -- diagnose Λ = 0 (BLOCKING)

Count the causes rather than infer them. Three of my inferences from output files were wrong:
tie-grade het scores (refuted: median 9.00, same as every read), low site reliability (refuted:
10.21 vs 10.13), and `min_reads` (refuted: 88.8% of losses leave both slots standing).

- A counter per guard at `graph_caller.cpp:6058-6084`, so exclusions are attributable.
- At each split homozygous site, classify every Λ = 0 read into causes 1-5 and report the histogram.
- Record, per excluded het site, whether any read there had `|q0 - 0.5| > epsilon` -- that is the
  test for "knowable but unused", and it is what separates cause 4 from cause 2.

**Gate: every Λ = 0 read at a split site attributed to exactly one cause, summing to the
`hom_split_no_opinion` counter.** No fix is written before this holds.

## Phase 2 -- assignment

**2a. Truly unassignable (causes 1-3, and 5 within a block): deterministic per-read coin.**
The 5,964 reads with no het anchor are proven safe -- measured, 0.00% of them touch a het anchor, so
no coin can contradict an assignment they do not have. Constraints:
- **Deterministic and per READ**, seeded from the read name's hash parity. Not a PRNG: a
  thread-dependent one would differ between runs and, worse, between a site's two pins.
- Applied at every site the read crosses, so the read lies wholly on one haplotype and creates no
  switch.
- `multi_block` reads get one coin per block, not one globally.

**2b. Cause 4: fix the chain, do not randomise.** Scope decided by Phase 1's histogram.

**Gates:** interior holes -> 0 by construction, confirmed not asserted; whatshap switch error against
T2T-Q100 no worse on chr20 and on chr6 held out; phased-run N50 no worse; and the hole measurement
added to `scripts/check_anchors.py` so the invariant is enforced rather than measured once.

## Phase 3 -- the Λ == 0 site gate (independent, trivial)

The splittability gate at `src/anchor.cpp:388` does **not** skip `lo == 0.0`, while the read
placement 100 lines below does. At q = 0 that pads `side1` with every uninformative read, since
`lo > 0.0` is strict. Measured: q = 0 against q = 1e-9 moves 15 sites of 89,407. Byte-identical at
any q > 0. Make the gate skip `lo == 0.0` explicitly, matching the placement.

## Phase 4 -- the mapQ lever, `--phase-mismap-min`

`allele_likelihood.cpp:333` stores the already-floored mismap and `:1514` hands that same value to
the phase path, so phasing inherits the genotyping clamp. Retain the raw probability, add a
phase-only floor, default it to the genotyping value, **gate on byte-identity at the default**, then
sweep 0.05 -> 0.01 -> 0.001 -> 0 against phased-run N50, held-out agreement, and whatshap switch
error.

The prediction to falsify is that τ absorbs it: 94.67% of chr20 ONT alignments are MAPQ 60 and 96.60%
sit at or below the 0.05 floor, so the floor rescales a near-constant, and τ is fitted per run.
Against that, agreement saturates at 95-97% however large |Λ| gets -- but that ceiling was measured
on a Λ built with the clamp, so it may be an artefact of the thing under test.

**Independent of Phase 2, and not a substitute for it.** A lower floor raises `p` and so the
MAGNITUDE of each term, which moves reads across the `--split-min-q` threshold -- but it cannot move
a term that is exactly zero, because `q0 = 0.5` gives `a = b` whatever `p` is. Causes 1, 2 and 4
survive it untouched.

## Phase 5 -- DROPPED

Both output-control flags are rejected, for the same reason as the defect this plan exists to fix:
they drop a read at some anchors and not others.

- `--anchors-max-reads N` caps per ANCHOR, so a read surviving the cap at a shallow anchor and losing
  it at a deep one reads present-absent-present -- worse than the current defect, because the holes
  would correlate with depth and cluster in high-coverage regions. The only hole-free form is
  whole-read subsampling, which belongs in read selection upstream.
- `--anchors-phased-reads` has the same trap in subtler form: Λ's leave-one-out is per site, so a
  per-(read, site) filter punches holes and only the per-read form is safe.

And the per-read provenance marker is dropped too: the reads it would mark are exactly the reads with
no phase evidence, so it is a distinction without a difference for any consumer. **No format bump.**

## Phase 6 -- build `read_strand` only where it is read (pure performance, output-identical)

`src/graph_caller.cpp:1507` builds the per-read log-odds vector at EVERY site whenever
`--anchors-hom-split` is on, but `build_site_anchors` dereferences it only inside `if (hom)`
(`src/anchor.cpp:378`, `389`, `493`). So at het, haploid and non-diploid sites -- and at every site
the function early-returns from -- the vector is built and never read.

And at a het site that IS a phase site the leave-one-out branch fires, which rebuilds the site's whole
own-map once per read: **Θ(n²)**, against Θ(n) at a hom site where the lookup misses. Worse, the
self-check at `src/graph_caller.cpp:1548` then recomputes the same values for the same reads, so a het
phase site pays that quadratic twice, once for a vector nobody reads.

Gate the construction on `hom && !haploid && genotype.size() == 2`, mirroring `src/anchor.cpp:378`
exactly. Three adversarial reads agree nothing changes: both dereferences are inside the hom branch,
and `hom_split`, `hom_unsplit` and `hom_split_no_opinion` all live there too.

**Constraint:** `src/anchor.cpp:379` gates on `read_strand->size() == evidence.reads.size()` and
indexes by the evidence index, so hom sites must still get a FULL-LENGTH vector. A filtered one makes
the gate decline silently and every hom site reports `hom_unsplit`.

**Gate: VCF and anchor file byte-identical, plus wall-clock.**
