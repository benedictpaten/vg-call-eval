# Why vg misses the SVs PanGenie finds: no call, or an unmatched call?

Over the 2,630 autosomal truth SVs that vg misses and PanGenie finds, asking what vg
actually wrote within 100 bp of each.

| what vg emitted at the locus | n | share |
|---|---|---|
| no record within 100 bp | 1,230 | 46.8% |
| record of comparable size present | 889 | 33.8% |
| record present but far too small | 488 | 18.6% |
| only a same-length substitution present | 23 | 0.9% |

`comparable size` means a record whose largest allele-length change is at least half
the truth variant's size, so the event was written but truvari declined the match --
a representation failure, and one that is scored twice, as a false negative here and as
an unmatched false positive in the same run.

## By truth variant size

| size | no record | comparable record | same-length substitution only | too small |
|---|---|---|---|---|
| 50-100 | 482 | 501 | 4 | 172 |
| 100-300 | 465 | 302 | 9 | 211 |
| 300-700 | 177 | 56 | 9 | 60 |
| 700+ | 106 | 30 | 1 | 45 |

## Loci where the only nearby call is a same-length substitution

| contig | truth pos | truth type | truth size | vg record pos | vg allele len |
|---|---|---|---|---|---|
| chr1 | 718,349 | DEL | 58 | 718,260 | 77 |
| chr2 | 1,870,164 | INS | 81 | 1,870,071 | 1,428 |
| chr5 | 44,416,030 | DEL | 52 | 44,416,058 | 53 |
| chr6 | 156,230,969 | DEL | 75 | 156,230,899 | 73 |
| chr8 | 62,875 | DEL | 133 | 62,860 | 67 |
| chr8 | 62,942 | DEL | 133 | 62,860 | 67 |
