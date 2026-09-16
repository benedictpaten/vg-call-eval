# Results: placing reads that cannot be phased

The defect, the fix, and what is left. vg `ec0a6ecee`. chr20 fits, chr6 is held out.

## The defect

A read with no cross-site opinion was dropped at a split homozygous site, but kept at het sites and
at unsplit homozygous ones, because neither of those consults the phase. Its trail through the anchor
file therefore read present, absent, present -- and a read that disappears between two anchors and
comes back is, downstream, the shape of a deletion.

## Why Λ = 0, measured rather than assumed

New counters attribute every zero. chr20:

| cause | count | assignable? |
|---|---|---|
| reached no phase site at all | 69,446 (97.8%) | yes -- nothing to know |
| reached one but the contributions summed to nothing | 1,586 (2.2%) | yes -- an allele pair the read cannot resolve |
| spanned a phase break | 0 | **no** -- it has evidence, on both sides, and they are not comparable |
| no calibrated table | 0 | n/a |

Corroborated: the run reports **76,135 phase sites** against **76,077 two-slot anchor snarls**, which
match to within 58 in the direction of phase sites with no anchors. So "reached no phase site" means
the read touches no heterozygous site, NOT that the chain silently dropped sites that had something
to say. That distinction is the one that decides whether a coin is honest.

Three earlier explanations for these reads were tested and refuted: tie-grade allele scores (their
het-anchor scores are median 9.00, identical to every read), low site reliability (10.21 against
10.13), and `min_reads` (88.8% of the losses are at snarls where both slots survive).

## The fix

Assign rather than drop, by the read name's hash parity. Sound because the site is HOMOZYGOUS: both
slots spell the same allele, so the read's sequence fits either haplotype and only the label is
arbitrary. Deterministic and **per read**, so the same read lands on the same haplotype at every site
and at both pins of each -- it lies wholly on one haplotype and creates no switch. A per-site or
per-pin random choice would create exactly the switches the discard was trying to avoid.

Reads spanning a phase break are still dropped, and are marked apart with a NaN so the writer can
tell them from honest zeros.

## Measured

| | chr20 discard | chr20 coin | chr6 discard | chr6 coin |
|---|---|---|---|---|
| reads present in both files | 72,345 | **79,385** | 207,494 | **229,355** |
| (read, snarl) placements lost | 4,630 | **587** | 7,163 | **1,070** |
| **INTERIOR, deletion-shaped** | 461 | **251** | 644 | **353** |
| terminal | 1,516 | 16 | 3,736 | 18 |
| **holes >= 10 kb** | 163 | **8** | 97 | **0** |
| median gap | 2,150 bp | 701 bp | 1,754 bp | 736 bp |
| max gap | 14,392 bp | 14,033 bp | 17,762 bp | 9,191 bp |

The headline is the first row: **5,964 reads on chr20 and 21,861 on chr6 stopped disappearing from
the output entirely.** That loss was never measured before and does not show up as a hole, because a
read that is never present cannot look deleted.

**Genotype-neutral on both contigs**: each coin arm's VCF is byte-identical to its no-split control.

## What is left, and why it is a different problem

251 interior holes on chr20 and 353 on chr6 remain. They are not the phase discard. Classified:

| | chr20 |
|---|---|
| a split slot fell below `--anchors-reads` | 291 (49.6%) |
| both slots survive | 290 (49.4%) |
| snarl absent from the hom-split file | 6 (1.0%) |

Splitting a homozygous site halves its reads between two slots, and a slot below `--anchors-reads`
(default 2) is dropped whole, taking its reads. That is a pre-existing per-slot filter interacting
with splitting, not a phase decision -- and the coin made it *better*, because coined reads now fill
slots that used to fall below the threshold. Total placement loss went from 88,591 R rows to 974,
0.0076%.

Fixing the rump properly means declining to split a site whose slots would come out thin, which needs
the per-slot counts before the split is committed. Not done here. The severity is bounded: no hole
over 10 kb survives on chr6 and eight do on chr20.

## Tests

`t/18_vg_call.t` 437 pass, including five new: the coin is deterministic across runs; no read is
placed in both slots of one snarl; `--phase-mismap-min` at the inherited value changes nothing, is
accepted at 0, and is refused above `--mismap-max`. Unit `[anchor]` 148 assertions, including one
rewritten to assert that the quiet read gets two rows, one per pin, on the SAME slot -- a count alone
would be satisfied by the per-pin randomisation that is the failure worth testing for.
