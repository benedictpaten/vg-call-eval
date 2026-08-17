# The SVs vg writes nothing for are mostly hidden by snarl scope, not misjudged

[sv-fn-mechanism.md](sv-fn-mechanism.md) isolates 1,230 autosomal truth SVs that PanGenie called, vg
missed, and vg emitted no record within 100 bp of. Because PanGenie found them the panel demonstrably
carries the alleles, so this was the population where the read model looked most clearly at fault.

It mostly is not the read model. **Every one of these events is reachable in the graph; the default
top-level-only snarl scope is what hides them.**

Measured on chr20 (48 loci) and chr6 (101), which together hold 149 of the 1,230 — 12%. Reproduce with:

```bash
python3 scripts/wgs/nocall_anatomy.py --loci work/nocall/loci.tsv \
    --allvcf work/nocall/chr20.all.vcf --contig chr20 --out docs/sv-nocall-chr20.md
```

## Why the trail goes cold in a normal run

`vg call` emits a record only where it calls non-reference, so a site confidently genotyped hom-ref
and a site the caller never examined produce identical output: nothing. `-a/--genotype-snarls` emits
reference calls too, which separates them, and `-A/--all-snarls` additionally descends into nested
snarls.

## What the caller actually did

Classified within truvari's own 500 bp refdist. An allele counts as *comparable* when it changes
length in the truth event's direction by 0.5x to 2x its size — both halves matter, since matching on
the largest change at a snarl would let a bubble offering a 5 kb allele stand in for a 64 bp
deletion, and matching without direction would let an insertion cover a deletion.

Top-level (`-a`), on both contigs, and the same chr20 loci re-examined with nesting on (`-a -A`):

| what the caller did | chr20 | chr6 | both | chr20 nested |
|---|---|---|---|---|
| no snarl within 500 bp | 18 (37.5%) | 52 (51.5%) | **70 (47.0%)** | **0 (0.0%)** |
| snarl there, no comparable allele offered | 17 (35.4%) | 33 (32.7%) | **50 (33.6%)** | 8 (16.7%) |
| comparable allele offered, called non-reference | 9 (18.8%) | 10 (9.9%) | 19 (12.8%) | 33 (68.8%) |
| comparable allele offered, called hom-ref | 4 (8.3%) | 6 (5.9%) | **10 (6.7%)** | 7 (14.6%) |

**80.6% of the population is upstream of the likelihood** — no snarl, or a snarl that never offered a
comparable allele — and the two contigs agree closely on that split. **6.7% is the model weighing
evidence and getting it wrong.**

**All 18 of chr20's "no snarl" cases are nested bubbles.** With nesting enabled the category empties
completely and 40 of 48 loci have a comparable allele on offer. Nothing about the graph or the panel
was missing; the default traversal scope never presented these sites to the model.

## The ten the model did get wrong

Small but unambiguous, and the only bucket that is a scoring failure. Quality fields from the
chr20 four:

| quantity | median | min | max |
|---|---|---|---|
| DP (reads at the site) | 41 | 7 | 58 |
| alleles offered | 2 | 2 | 4 |
| share of AD on a non-reference allele | 0.433 | 0 | 0.758 |
| GQ of the hom-ref call | 2 | 2 | 3 |
| GQN | 0.027 | 0.021 | 0.27 |
| DR | 0.371 | 0.086 | 0.484 |

Ample depth, 43% of reads supporting the alternate allele — a textbook heterozygote — called hom-ref
at GQ 2. `DR` near 0.37 says the observed depth is well under what the called genotype predicts, so
the caller's own diagnostics already flag these as implausible: a `--min-confidence` gate would mark
them, not rescue them. Ten such loci across the two contigs, 6.7% of the population.

## Does nested calling recover the recall?

No — and confirming a diagnosis is not the same as having a fix. Both nested configurations were
scored end to end through the identical `bench_wgs.py` path: same truth, same truvari invocation,
same `--sizemin/--sizefilt`, same aardvark call.

| chr20 | SV TP | SV FP | SV FN | SV recall | SV precision | **SV F1** | small-variant F1 |
|---|---|---|---|---|---|---|---|
| default (top-level) | 375 | 367 | 390 | 0.4902 | 0.4986 | **0.4944** | **0.9646** |
| nested (`-a -A`) | 392 | 459 | 373 | **0.5124** | 0.4562 | 0.4827 | 0.9671 |
| `--top-down` | 369 | 430 | 396 | 0.4824 | 0.4543 | 0.4679 | 0.9470 |

**Neither is a win, and `--top-down` is worse than doing nothing.**

`-a -A` buys the recall the diagnosis predicts — 17 more true calls — and pays 92 false ones for it,
because calling every snarl independently emits a nested event *alongside* its parent. Two records
describe one event and truvari can match only one, so the duplicate scores as a false positive.

`--top-down` propagates genotypes from parent to child instead, which was the obvious candidate for
getting the recall without the duplication. It loses on every axis, including **recall below the
default** (0.4824 against 0.4902) and small-variant F1 down 0.0176. Constraining a child by its
parent's genotype evidently costs more calls than the extra scope adds.

So the fix is not a flag. Converting this diagnosis into recall needs nested calls emitted *in place
of* their parents where they explain the event, which is a change to how nested snarls are
represented on output rather than a change to scoring or to a traversal limit.

## What this changes about the earlier conclusion

[pangenie-comparison.md](pangenie-comparison.md) attributes the residual SV recall gap to the read
model, on the grounds that the panel carried alleles vg did not call. The panel did, and vg did not —
but the model was never shown them. **80.6% of this population is snarl scope and allele enumeration;
6.7% is the likelihood weighing evidence and getting it wrong.**

That is a different kind of problem from the mixture-weight defect in
[tier2-sv-errors.md](tier2-sv-errors.md), and it will not be fixed by anything in the scoring model.

## Caveats

- **149 of 1,230 loci (12%)**, from two contigs that agree closely on the split. The remaining 20
  contigs have not been re-called, and the accuracy comparison is chr20 only.
- The `-A` numbers count what vg *called*, at loci that were false negatives in the default run.
  The truvari row above is the only part that has been scored end to end.
- Truth SVs cluster heavily in this set — 31 of the 48 have another truth SV within 500 bp — and
  truvari matches one-to-one, so some of these can never all be matched at once regardless of scope.
