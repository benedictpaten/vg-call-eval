# What happened to the SVs vg wrote nothing for (chr6)

101 loci, re-called with `-a/--genotype-snarls` so reference calls are emitted too.
Match window 500 bp (truvari's refdist). An allele counts as comparable when it
changes length in the same direction as the truth event and by 0.5x to 2x
its size.

| what the caller did | n | share |
|---|---|---|
| no snarl | 52 | 51.5% |
| no comparable allele offered | 33 | 32.7% |
| offered and called, truvari rejected it | 10 | 9.9% |
| offered and called hom-ref | 6 | 5.9% |

## By truth variant size

| size | no snarl | no comparable allele offered | offered and called, truvari rejected it | offered and called hom-ref |
|---|---|---|---|---|
| 50-100 | 21 | 16 | 4 | 0 |
| 100-300 | 16 | 14 | 3 | 2 |
| 300-700 | 9 | 2 | 0 | 2 |
| 700+ | 6 | 1 | 3 | 2 |

## The hom-ref calls: was there evidence to read?

| quantity | median | min | max |
|---|---|---|---|
| DP (reads at the site) | 36 | 26 | 362 |
| alleles offered | 1 | 1 | 4 |
| share of AD on a non-ref allele | 0.411 | 0.0352 | 0.754 |
| GQ of the hom-ref call | 2 | 1 | 23 |
| GQN | 0.029 | 0.004 | 0.277 |
| DR | 0.501 | 0.289 | 0.662 |

A hom-ref call with ample DP and a real share of AD on the alternate allele is the
model mis-weighing evidence it had. One with DP near zero is a site the reads never
covered, which no scoring change reaches.
