# Why vg misses the SVs PanGenie finds: no call, or an unmatched call?

Over the 1,896 autosomal truth SVs that vg misses and PanGenie finds, asking what vg
actually wrote within 100 bp of each.

| what vg emitted at the locus | n | share |
|---|---|---|
| record of comparable size present | 883 | 46.6% |
| record present but far too small | 743 | 39.2% |
| no record within 100 bp | 259 | 13.7% |
| only a same-length substitution present | 11 | 0.6% |

`comparable size` means a record whose largest allele-length change is at least half
the truth variant's size, so the event was written but truvari declined the match --
a representation failure, and one that is scored twice, as a false negative here and as
an unmatched false positive in the same run.

## By truth variant size

| size | no record | comparable record | same-length substitution only | too small |
|---|---|---|---|---|
| 50-100 | 74 | 494 | 1 | 244 |
| 100-300 | 100 | 306 | 4 | 328 |
| 300-700 | 42 | 55 | 4 | 99 |
| 700+ | 43 | 28 | 2 | 72 |

## Loci where the only nearby call is a same-length substitution

| contig | truth pos | truth type | truth size | vg record pos | vg allele len |
|---|---|---|---|---|---|
| chr1 | 116,445,238 | DEL | 128 | 116,445,302 | 52 |
| chr2 | 860,334 | INS | 383 | 860,314 | 98 |
| chr6 | 156,230,969 | DEL | 75 | 156,230,899 | 73 |
| chr9 | 88,463,150 | INS | 532 | 88,463,173 | 134 |
| chr11 | 121,210,820 | INS | 385 | 121,210,778 | 71 |
| chr14 | 30,138,035 | INS | 108 | 30,137,952 | 50 |
