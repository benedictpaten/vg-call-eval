# Stage 0: would symbolic alleles collapse the records that are losing variants? (chr6)

Snarl inventory from the `-a -A` run: 491,182 snarls. Production records projected:
284,529.

A record *collapses* when every ALT's symbolic allele equals the reference traversal's, so the
whole record is reference at this level and its differences belong to nested sites. A record is
*non-leaf* when any allele's traversal crosses a child snarl at all -- the precondition for
symbolising anything.

| | n | share of its row group |
|---|---|---|
| all emitted records | 284,529 | |
| ... non-leaf | 2,884 | 1.0% |
| ... would collapse to reference | 990 | 0.3% |
| **large records (REF >=50 bp)** | 1,381 | |
| ... non-leaf | 1,039 | 75.2% |
| ... would collapse | 440 | 31.9% |
| **same-length alleles >=50 bp (the substitution FPs)** | 230 | |
| ... non-leaf | 228 | 99.1% |
| ... would collapse | 147 | 63.9% |
| **records swallowing a missed SNV** | 156 | |
| ... non-leaf | 153 | 98.1% |
| ... would collapse | 123 | 78.8% |

## Missed SNVs recovered by collapsing

| | n |
|---|---|
| swallowed SNVs total | 3,304 |
| swallowed SNVs in a collapsing record | 3,133 |
| swallowed SNVs in a non-leaf record | 3,301 |
|   of those PanGenie called | 2,942 |

## Examples of collapsing records

| POS | REF bp | alleles | missed SNVs inside |
|---|---|---|---|
| 93,877 | 10,506 | 3 | 16 |
| 4,987,353 | 1,464 | 2 | 1 |
| 5,476,367 | 2,382 | 2 | 3 |
| 7,675,124 | 28,153 | 3 | 51 |
| 8,208,065 | 1,451 | 2 | 1 |
