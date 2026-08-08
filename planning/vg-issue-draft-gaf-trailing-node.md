# Draft vg issue: `vg convert -G` drops a trailing path node with zero reference length

**Status: DRAFT, not filed.** Would go to <https://github.com/vgteam/vg/issues>. Searched first —
nothing existing covers it; #4166 is adjacent (a different GAF conversion limitation) but distinct.

---

**Title:** `vg convert -G` drops a trailing path node whose mapping consumes no reference, so GAM→GAF→GAM
is not path-stable

Round-tripping an alignment through GAF loses a node from the end of the path when the read's final
mapping consumes **zero reference bases** — a pure trailing insertion. The alignment's reference span is
unchanged, but the recorded path is shorter, so any tool that asks "which nodes does this read touch"
gets a different answer depending on which format the reads arrived in.

### Reproducer

Uses only data already in the repo (`test/small/x.fa`, `test/small/x.vcf.gz`):

```bash
vg autoindex -r x.fa -v x.vcf.gz -w giraffe -p idx -t 4
vg sim -x idx.giraffe.gbz -n 2000 -l 100 -s 42 -e 0.01 -i 0.002 -a > sim.gam
vg view -a -X sim.gam > sim.fq
vg giraffe -Z idx.giraffe.gbz -f sim.fq -t 4 > mapped.gam
vg convert -G mapped.gam idx.giraffe.gbz | grep -v '^@' > mapped.gaf
vg convert -F mapped.gaf idx.giraffe.gbz > roundtrip.gam
# then compare the node list of each mapping path before and after
```

5 of 1872 mapped reads come back with a shorter path. For example:

```
before: ... 141 139 138 137 135
after : ... 141 139 138 137
```

### What is happening

The dropped node's mapping is a pure insertion — `from_length` 0, `to_length` > 0:

```
GAM : ... node 133 [(14,14,"")]   node 134 [(0,6,"TTATCT")]
GAF path column: >...>132>133      <- 134 absent
```

It is the **writer**, not the reader: the `path` column emitted by `vg convert -G` already lacks the
node, so `vg convert -F` has nothing to rebuild it from.

The rule appears to be exact. Of 1872 mapped reads, 5 have a final mapping consuming no reference bases,
and those 5 are precisely the 5 whose path changes. It never affects a leading mapping and never a
mapping mid-path — only the last one.

### Why it matters

Mostly it does not, which is why it has probably gone unnoticed: the aligned reference interval, CIGAR,
and score are all unaffected. It matters for anything deriving per-node read membership from alignment
paths. We hit it building read-level genotyping for `vg call`
([#4990](https://github.com/vgteam/vg/pull/4990)): a read that overlapped a site *only* through the
dropped node stops being evidence for it, so `DP` falls by one at a handful of sites and the genotype
likelihoods shift accordingly. Same reads, same graph, same caller — different numbers depending on
whether they were supplied as GAM or as GAF.

### Not part of this report

Round-tripping also re-partitions insertions that sit exactly on a node boundary (42 further reads in the
run above): GAM attaches such an insertion to the *following* mapping, GAF to the *preceding* one. That
one is **not** a defect — GAF's difference string has no way to express which side a zero-reference
insertion sits on, and it makes no difference to the alignment. Mentioning it only to keep it out of
scope.

### Open question for you

Two places this could reasonably be fixed, and we do not have a view on which you would prefer:

1. `vg convert -G` keeps the node in the path even though it contributes no reference bases; or
2. `vg giraffe` does not emit a trailing mapping that consumes no reference in the first place — it is
   arguably describing a path step the read does not actually align to.

Either way there is currently no test covering GAM→GAF→GAM path stability, which is probably the more
useful thing to add.
