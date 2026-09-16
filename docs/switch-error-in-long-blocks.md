# Do the long phase blocks contain switches?

Yes -- about 41 real ones on chr20 -- and genotype confidence does not remove them. Measured with
whatshap against the T2T-Q100 phased truth, vg `ec0a6ecee`, chr20 ONT.

## Switches and flips are different errors and must be read apart

`whatshap compare` reports `all_switches`, which counts a FLIP as two switches, and
`all_switchflips` as `S/F`, the decomposition into true switches and flips. A flip is one site on
the wrong haplotype with the phase immediately resuming; a switch is the phase changing and staying
changed. Only the second shortens the usable haplotype.

| stratum | assessed pairs | raw switch rate | **true switches** | flips |
|---|---|---|---|---|
| all het sites | 58,803 | 0.3826% | **55** | 85 |
| GQ >= 20 | 55,353 | 0.2186% | **45** | 38 |
| GQ >= 40 | 50,411 | 0.1765% | **41** | 24 |
| GQ >= 60 | 42,541 | 0.1669% | **41** | 15 |
| GQ >= 80 | 29,677 | 0.1719% | **33** | 9 |
| GQN >= 0.3 | 49,511 | 0.1434% | **35** | 18 |

Check: 55 + 2x85 = 225 = the raw switch count, and the identity holds on every row.

**Flips are a genotype-quality problem and switches are not.** Flips fall nine-fold with confidence,
85 -> 9. True switches fall by a third at most, 55 -> 33, and are flat across GQ 40 to 60. So
filtering on genotype confidence cleans up isolated noise and leaves the real phase changes where
they are.

## What that means for range

chr20 is one 66.2 Mb block, and roughly 41 true switches sit inside it: **one every ~1.6 Mb**.

That is shorter than the 2,798,104 bp read-walkable phased-run N50 measured from the anchors, and the
gap is the point. **An anchor run measures where the chain is UNBROKEN, not where it is CORRECT.** It
cannot see a switch, because a switch leaves every anchor in place and only changes which haplotype
they name. The two numbers are both true and they answer different questions; for "how far can a
haplotype actually be followed", 1.6 Mb is the honest figure.

## Dropping unreliable het sites: measured and REFUTED

The proposal was that an unreliable heterozygous site admitted to the chain contributes a noisy term
to every crossing read's log-odds, so raising `--phase-min-q` might buy cleaner cross-site phase.
`--phase-min-q` gates a site on its own mean read score and is fitted at 8.5.

| `--phase-min-q` | reliable het sites | sites split | agreement | run N50 | **true switches** | flips |
|---|---|---|---|---|---|---|
| **8.5 (default)** | 60,203 | 85,678 | 94.7102% | 2,798,104 | **55** | 85 |
| 10 | 3,016 | 84,790 | 94.8604% | **3,273,415** | **106** | 171 |
| 12 | 0 | 70,592 | 79.4231% | 804,382 | **1,572** | 296 |
| 15 | 0 | 70,592 | 79.4231% | 804,382 | 1,572 | 296 |

**Every metric short of switch error says 10 is better, and switch error says it is twice as bad.**
At 10 the run N50 rises 17%, the per-read held-out agreement rises, and the split count barely moves
-- while the true switches double, 55 -> 106. The same doubling survives the confidence filter, 41 ->
84 at GQ >= 40. The longer runs are longer and wronger.

Note the cliff: 8.5 -> 10 cuts reliable sites from 60,203 to 3,016, because the reliability
distribution is packed just above the default (chr20 median 8.98, ceiling phred(--mismap-min) =
13.01). At 12 nothing is reliable, the chain falls back to the panel, and switch error goes to 3.7%.

**So the unreliable hets are load-bearing.** They are individually noisy and collectively the thing
that keeps the chain dense enough that no single link carries a long stretch on its own. `8.5` stays.

## The methodological point

This is the second time in this line of work that a range metric pointed the wrong way. Range is a
property of where the chain is unbroken; correctness is a property of the labels it carries. **A
phasing change must be gated on switch error against a truth, not on run length**, and run length
alone would have adopted `--phase-min-q 10` here.
