# The ONT MAPQ floor is 5, and the switch-error case for it was wrong

> **Current state (2026-09-23).** `--preset ont` still sets `--read-min-mapq 5`, on the
> reliable-het case below, which was measured on vg `05acbfc6a` (2026-09-17). ONT is insensitive to
> both MAPQ levers: on the current preset, dropping the floor, at a cap of 0.7 or 0.95, moves no F1
> on chr20 or chr6 by more than 2e-5 and no truvari SV call (`work/wgs-run5`). Short reads keep
> `--read-min-mapq 0`, with `--mismap-max` now 0.95. Over chr1-22+X, with both arms at the
> then-default cap of 0.7, MQ5 raises short-read SV F1 by +0.0067 but adds 3,290 net SNV errors,
> while the 0.95 cap gains SV F1 +0.0021 with small variants unchanged within noise. Against the
> current default, that MQ5 arm still leads on SV F1 by 0.0047 and trails on SNV F1 by 0.0005, both
> significant (`work/wgs-run3/compare_mq5_vs_mm095.out`; [tier2-parameters.md](tier2-parameters.md)).

Supersedes the switch-error argument in [mapq-filter-and-the-clamp.md](mapq-filter-and-the-clamp.md)
and the `--read-min-mapq 10` row in
[anchor-defaults-and-chr20-delivery.md](anchor-defaults-and-chr20-delivery.md). Shipped as vg
`05acbfc6a` on PR #4990.

## What was wrong

`--read-min-mapq 10` was set on chr20 going 55 -> 45 true switches with F1 flat. Three problems,
found by re-measuring rather than by re-reading:

1. **There is no MAPQ 10 arm.** The scan ran 0/5/11/20/30/60 and the headline 45 is the MAPQ 20 arm.
   `mq_base.vcf` is md5-identical to `sq11.vcf`, so the shipped setting was the MAPQ 11 arm, which
   reads 49 raw and is the *worst* arm under a pinned site set.
2. **It is one locus.** Under a pairwise pin the arms read 35 vs 20 inside chr20:26-27 Mb and
   188 vs 185 over the whole rest of the contig. The discordant positions there are not independent
   -- 21 of 29 fall in three clusters inside 227 kb, some 3 bp apart -- so collapsing at any merge
   distance >= 1 kb takes McNemar from p=0.006 to p=0.42.
3. **It does not replicate.** chr6 raw is 459 all-switches against 458, and the "63 -> 60" that was
   quoted is `63S/198F -> 60S/199F`: switch/flip reclassification at constant total error.

A fourth thing was wrong on my side of it. I "corrected" the above by intersecting all six arms'
correctly-genotyped site lists and re-scoring on that -- and that test is biased toward the null by
construction. 260 of the 294 sites it deletes lie in chr20:26-30 Mb, the only region where the arms
differ at all, and 58.8% of the exclusions are driven by MAPQ 60, an arm not in an MAPQ 0 vs 20
comparison. **Pin pairwise, to the two arms actually being compared.** The right way to check for the
denominator artefact in the first place is simply to compare the denominators: on ONT they move
0.006% while the switch count moves 18%, so there was never an artefact to correct.

## What 5 is set on

Reliable het count, measured with off-reference nesting on, which is what anchor output uses:

| MAPQ 0 -> 5 | chr20 | chr6 hold-out |
|---|---|---|
| reliable hets | 60,695 -> 61,863 | 163,916 -> 164,438 |
| chain breaks | 1,269 -> 1,432 | 1,503 -> 1,534 |
| ALL F1 | 0.95825 -> 0.95826 | 0.96520 -> 0.96486 |
| TP / FP | +0 / -1 | -89 / +101 |
| true switches | 68 -> 56 | 66 -> 69 |

The reliable-het rise holds on both contigs. Calling accuracy is exact on chr20 and moves 0.03%
relative on chr6 -- reproducible, since `vg call` is byte-identical on a repeat run, but far too
small to decide anything either way. **The contigs disagree on switch error**, so no switch claim is
made and none should be added later.

5 rather than 11: the chain-break cost per reliable het recovered roughly doubles above 5.

## Off-reference nesting has a switch cost, stated separately

Isolated for the first time here -- same binary, MAPQ 0 on both sides, the only difference being
`--anchors-out`:

| | chr20 | chr6 |
|---|---|---|
| true switches | 55 -> **68** | 63 -> **66** |
| ALL F1 | -0.00022 | +0.00002 |

Nesting is not optional -- off-reference variants have to be in the anchors file, and no other
output carries them -- so this is a fact to record, not a trade to weigh. It does mean the earlier
"the branch takes chr20 55 -> 46" was misleading: nesting raises 55 -> 68 and the MAPQ filter pulled
it back, rather than the branch reducing switches on its own.

## Method note, for next time

- Compare denominators before building a fixed-set test; if they match, the raw counts are fine.
- If a fixed set is needed, pin **pairwise**.
- Localise before believing an aggregate: one megabase produced this entire effect.
- Collapse non-independent positions before any paired test.
