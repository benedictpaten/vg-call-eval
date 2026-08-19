# Why truvari refused a comparable nearby call

Restricted to vg-only false negatives with a vg record of comparable size within 100 bp
(678 of them). Proximity was not the obstacle -- truvari's refdist default is 500 bp.

| reason | n | share | of which another truth SV within 500 bp |
|---|---|---|---|
| call was compared and rejected on similarity | 486 | 71.7% | 331 (68%) |
| call was matched to a different truth variant | 192 | 28.3% | 192 (100%) |

A call matched to a *different* truth variant is not a caller error at all. Truvari matches
one-to-one, so where truth SVs cluster, one call can satisfy only one of them and the
remainder are false negatives no matter how good the call is. That share of the recall gap
is not recoverable by changing the model, and should be excluded before quoting how much
is.

A call compared and rejected on similarity is a genuine disagreement about the variant's
sequence or size, and is the population where a better model or a better traversal could
change the outcome.
