# How much of vg's SNV recall loss is reachable at all

Every SNV aardvark marks FN for vg, over 22 contig(s): 142,707 variants.

The panel VCF is a deconstruction of the same graph vg calls on, so a truth variant absent
from it is not carried by the graph and no caller change could produce it.

| | n | share |
|---|---|---|
| not in the panel (graph cannot express it) | 56,900 | **39.9%** |
| in the panel, swallowed by a large vg allele | 55,222 | **38.7%** |
| in the panel, its own site, still missed | 30,585 | **21.4%** |

## Split by what PanGenie made of the same variant

| | PanGenie TP | PanGenie FN | no verdict |
|---|---|---|---|
| not in the panel (graph cannot express it) | 458 | 56,442 | 0 |
| in the panel, swallowed by a large vg allele | 47,885 | 7,337 | 0 |
| in the panel, its own site, still missed | 19,023 | 11,562 | 0 |

A variant absent from the panel that PanGenie nonetheless called TP would be a
contradiction, so that cell is a check on the method rather than a result.
