# Stage 0: would symbolic alleles collapse the records that are losing variants? (chr20)

Snarl inventory from the `-a -A` run: 222,797 snarls. Production records projected:
105,251.

A record *collapses* when every ALT's symbolic allele equals the reference traversal's, so the
whole record is reference at this level and its differences belong to nested sites. A record is
*non-leaf* when any allele's traversal crosses a child snarl at all -- the precondition for
symbolising anything.

| | n | share of its row group |
|---|---|---|
| all emitted records | 105,251 | |
| ... non-leaf | 1,569 | 1.5% |
| ... would collapse to reference | 606 | 0.6% |
| **large records (REF >=50 bp)** | 885 | |
| ... non-leaf | 677 | 76.5% |
| ... would collapse | 295 | 33.3% |
| **same-length alleles >=50 bp (the substitution FPs)** | 171 | |
| ... non-leaf | 169 | 98.8% |
| ... would collapse | 109 | 63.7% |
| **records swallowing a missed SNV** | 62 | |
| ... non-leaf | 60 | 96.8% |
| ... would collapse | 48 | 77.4% |

## Missed SNVs recovered by collapsing

| | n |
|---|---|
| swallowed SNVs total | 1,041 |
| swallowed SNVs in a collapsing record | 1,015 |
| swallowed SNVs in a non-leaf record | 1,037 |
|   of those PanGenie called | 819 |

## Examples of collapsing records

| POS | REF bp | alleles | missed SNVs inside |
|---|---|---|---|
| 3,845,329 | 3,695 | 3 | 5 |
| 3,862,587 | 2,688 | 3 | 2 |
| 3,872,014 | 2,830 | 3 | 2 |
| 5,784,980 | 1,897 | 3 | 14 |
| 8,534,844 | 4,447 | 2 | 12 |
