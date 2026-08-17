# What happened to the SVs vg wrote nothing for (chr20)

48 loci, re-called with `-a/--genotype-snarls` so reference calls are emitted too.
Match window 500 bp (truvari's refdist). An allele counts as comparable when it
changes length in the same direction as the truth event and by 0.5x to 2x
its size.

| what the caller did | n | share |
|---|---|---|
| no snarl | 18 | 37.5% |
| no comparable allele offered | 17 | 35.4% |
| offered and called, truvari rejected it | 9 | 18.8% |
| offered and called hom-ref | 4 | 8.3% |

## By truth variant size

| size | no snarl | no comparable allele offered | offered and called, truvari rejected it | offered and called hom-ref |
|---|---|---|---|---|
| 50-100 | 8 | 12 | 6 | 1 |
| 100-300 | 4 | 3 | 2 | 3 |
| 300-700 | 2 | 1 | 1 | 0 |
| 700+ | 4 | 1 | 0 | 0 |

## The hom-ref calls: was there evidence to read?

| quantity | median | min | max |
|---|---|---|---|
| DP (reads at the site) | 41 | 7 | 58 |
| alleles offered | 2 | 2 | 4 |
| share of AD on a non-ref allele | 0.433 | 0 | 0.758 |
| GQ of the hom-ref call | 2 | 2 | 3 |
| GQN | 0.027 | 0.021 | 0.27 |
| DR | 0.371 | 0.086 | 0.484 |

A hom-ref call with ample DP and a real share of AD on the alternate allele is the
model mis-weighing evidence it had. One with DP near zero is a site the reads never
covered, which no scoring change reaches.
