# Mosaic v5: each strand as one contiguous walk

The mosaic's job is to state the sample's genome as a recombination of panel haplotypes. v4 states it
as a list of **sites**; a consumer wanting the **walk** cannot get there from the file. This is the
plan to fix that, with the format change kept as small as the measurements allow.

## What the measurements found

Everything below is chr20, HG002, 34-haplotype panel, and every step held the VCF byte-identical.

**v4 segments do not chain.** Only 236 of 7,875 consecutive pairs met at a shared node -- 3.0%. A
segment spans its first site's start node to its *last site's end node*, so everything between one
segment's last site and the next segment's first site is named by nothing: median 1,041 bp, p95
8.8 kb, **17.2 Mb in total, 27% of the contig**.

**`.` is a VCF notion that had leaked into a graph file, and it was most of the problem.** A `.` row
said "this strand carries no sequence at this record". In graph terms that strand is not empty
there -- it traverses the parent's other allele, which bypasses the child snarl, possibly across a
deletion edge -- so the site is simply not on its walk. Worse, emitting the row *cut the other
strand's run in three*: 419 rows, **351 of them (84%) flanked by the same haplotype on both sides**,
interrupting a walk where nothing had happened.

Segmenting each strand over the sites **it actually traverses** removed all 419 rows, merged 351
pairs of runs, took segments from 7,882 to 7,112, and made **821 of the 917 then-remaining gaps stop
existing**. They were never gaps in a walk; they were manufactured by the segmentation.

**The extension ladder, complete:**

| | boundaries | share |
|---|---|---|
| extend right -- carry this segment's haplotype to the next segment's first node | 6,902 | 97.1% |
| extend left -- carry the *next* segment's haplotype back to this one's last node | 171 | 2.4% |
| reference patch | 37 | 0.5% |
| break | **0** | -- |

Contiguity went **3.0% -> 99.5%**, and the reference closes all 37 of the remainder, so on this graph
a thread is contiguous end to end with no breaks at all.

**Extend-left must not be nested inside extend-right.** It asks whether the *next* segment's
haplotype reaches back, so it does not need this segment to be walkable. It was inside the
`p != invalid_edge()` guard and therefore never tried for the 57 boundaries whose left segment has no
position; hoisting it recovered 50 of them. The distinction only surfaced from asking what
"no position to extend from" meant as against "neither haplotype spans the gap".

**Clipping is the ordinary case, not an anomaly.** Only **2 of the 34** panel haplotypes are
contiguous across chr20; `GRCh38#0` is stored in 9 subpaths. So a haplotype that cannot be carried
across a gap is expected, and the reference -- one of the two contiguous paths -- is what covers it.

**Where the 37 are.** 31 distinct loci, 6 of them appearing on both strands: **28 loci
pericentromeric (25.8-32.5 Mb), 2 subtelomeric, 1 elsewhere**. 17 of the 37 span more than 10x
chr20's average node density and the worst two are 9,912x and 31,748x -- 329,080 nodes over 1,071
reference bases. The reference closes them, and through satellite tangle it is a poor proxy for the
sample. The file should say what it filled and how big it was; the judgement belongs to the reader.

**Two things remain untested rather than working.** No segment on chr20 or chrX carries a reverse
orientation -- 0 of 15,008 -- so inversions are unexercised. And row order is *reference* order,
which is not walk order for a strand traversing an inversion, so segmentation may be cutting runs in
the wrong frame there. Both need a fixture, not an argument.

## The format

Three semantic changes, no new columns and no new record types. Version 5.

```
#mosaic-version 5
#graph      work/wgs-tt/chr20/chr20.gbz
#sample     HG002
#reference  CHM13#0#chr20
#decoding   constrained-viterbi
#patch      reference            <- or 'none'; what was done with gaps
#haplotype  0   CHM13#0
#H  contig strand ref_start ref_end start_node end_node hap_index haplotype sites gbwt_node gbwt_offset nested_sites max_depth

H   chr20  0  22      533     114818871  114819170  4    recombination#19  5   229637742  19  0  0
H   chr20  0  945     2268    114819170  114820704  9    recombination#30  4   229638340  9   1  1
H   chr20  1  26433521 26434564 117486629 117486715 ref  CHM13#0           .   234973258  3   .  .
```

1. **`end_node` is where the next segment begins**, wherever this segment's haplotype could be
   carried that far -- so consecutive segments of one strand share a node and the strand is one walk.
   This is the change that makes the file describe a path. It is also why the version must move: a
   v4 consumer reading a v5 file would under-read every segment.

2. **`.` never appears as a haplotype.** A strand that does not traverse a site has no row for it.
   `*` is unaffected and still means "the strand traverses this, and the panel cannot name a
   haplotype for it" -- a real stretch of walk with no label.

3. **`hap_index` may be the literal `ref`**, marking a stretch filled with the reference because
   neither flanking haplotype could be carried across it. Such a row carries a real
   `gbwt_node`/`gbwt_offset` (so it is walkable like any other) and `.` for `sites`, `nested_sites`
   and `max_depth`, because it covers no called site. Its span is `start_node`..`end_node` and
   `ref_start`..`ref_end`, so a consumer can weigh a 22 bp fill against 1 kb of satellite without a
   new column.

**Which haplotype covers a gap is arbitrary and the file says so.** No called site lies in the
stretch between two segments, so nothing distinguishes the left haplotype from the right one there
and a recombination anywhere inside is equally consistent. Extending right is a convention, not an
inference: the crossover is *bracketed* by the gap, not located within it. A `#note` states this,
because a consumer who reads the boundary as the crossover point will over-trust it.

**The invariant the format now carries:** for consecutive segments of one strand,
`end_node[i] == start_node[i+1]`, with no exceptions once patching is on. That is one pass to check
and it is what "translatable into a path object" reduces to.

## Implementation

Most of it is built and unpushed. What remains is small.

**Done, measured, uncommitted:** per-strand site lists (change 2); extend right and extend left with
earned checks and the position moving with an extended start; the head-clipped bail-out fix (+80
walkable sites); the reference index found by name rather than assumed to be 0.

**To write:**

| | |
|---|---|
| the patch row itself | ~15 lines in `emit_row`: where the gap was not closed and patching is on and the reference spans it, write a second row with `hap_index = ref`. `out` is sequential, so no restructuring. |
| `--mosaic-patch-gaps` / `--no-mosaic-patch-gaps` | **on by default**. Off leaves the gap and reports it, which is what a consumer who would rather see a break than an assertion wants. |
| version 5 + the three `#note`s and `#patch` | header only |
| delete the dead `.` path | the `StrandKind::Empty` branch in `emit_row` is now unreachable -- 0 rows emitted on chr20 -- along with `empty_segments` and the header note describing `.` |

Breaks are **not** implemented. There is no population for them on these graphs (0 of 7,110), and the
case that needs them is off-reference/gref, where there is no reference to patch with. Designed but
unbuilt is the right state: PanSN subpaths (`HG002#0#chr20[offset]`) with a `T` row carrying a break
reason, matching how the graph already stores its own clipped haplotypes.

## Tests

**The acceptance test is the round trip.** `test/mosaic_to_path.py` expands each thread -- slicing
each segment out of its haplotype's walk and concatenating -- and checks the result is an exact walk
in the graph, every step a real edge. Every other property is a way for this to fail, so it subsumes
them. Written, and passing on the small fixture; it needs `vg paths -A` and `vg view -g`, so no new
subcommand.

**Structural invariants, one pass each, in `18_vg_call.t`:**

- `end_node[i] == start_node[i+1]` for every consecutive pair of a strand -- **the** invariant, and
  now assertable because the extension rule makes it true
- no row has a named haplotype and a missing position unless it is a patch
- every non-patch segment covers at least one site
- the reported counts (extended right, extended left, patched, gap) equal what is in the file

**Two of the four current assertions are withdrawn**, both because they encoded v4 accidents:

- *both strands account for the same number of sites* -- true in v4 only because `.` rows padded the
  non-traversing strand to match. Excluding them the strands genuinely differ (113,209 against
  113,034 on chr20), and under v5 there are no `.` rows to pad with.
- *segments are in reference order within a strand* -- passes because the phasing is sorted into
  reference order before segmentation, so it is near-tautological; and it would be **wrong** for a
  thread traversing an inversion, where walk order is not reference order.

**The fixture is the remaining gap, and it now has three reasons to exist.** A hand-built graph with
a clipped haplotype and an inversion, small enough for CI, gives regression tests for: the
head-clipped bail-out; the reference patch path (which no CI fixture reaches, since the small graph
has no gaps); and walk-order-against-reference-order under inversion, which nothing currently
exercises -- 0 reverse-orientation segments in 15,008 across chr20 and chrX.

## Gate

chr20 and chrX, VCF byte-identical throughout -- this is a mosaic-only change. Then contiguity at
**100%** on both (from 3.0%), the round trip passing on the fixture, `vg test`, and `18_vg_call.t`
alone.
