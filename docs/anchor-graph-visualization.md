# Seeing the anchor graph

`scripts/anchor_graph_gfa.py` condenses an anchor file into a GFA that Bandage can open.

```
python3 scripts/anchor_graph_gfa.py \
    --anchors HG002.chr20.ont.anchors.homsplit.tsv \
    --vcf     HG002.chr20.ont.vcf \
    --region  chr20:26000000-31000000 \
    --gfa out.gfa --csv out.csv
```

Then in Bandage: *File > Load graph* on the GFA, *File > Load CSV* on the companion, and colour by
`Class`. 16 s for all of chr20.

## The model

**One pin is one bidirected node.** A pin is one `A` row -- a (boundary node, snarl, slot) triple --
and a read passes through it entering one side and leaving the other, with the row's `strand` saying
which way round. The pin rather than the site-slot is the unit because not every site-slot keeps
both its pins: on chr20, 311,820 have two and 31,190 have one.

Links come from reads. Sorting one read's pins by `offset` gives its path; each consecutive pair is
a link, weighted by how many reads carry it (`RC:i:` on the `L` line, so Bandage can draw thickness).

## What you are looking at

A phased site has two slots, so it contributes **two parallel strands**, one per haplotype. A site
that failed to split has one slot that every read goes through, so the two strands **pinch into a
single node**. Those pinches are the run breaks from
[anchor-defaults-and-chr20-delivery.md](anchor-defaults-and-chr20-delivery.md) -- `Class=collapsed`
in the CSV.

The ladder is real and visible in the raw links. Over 50 kb at chr20:5.00 Mb the slot-0 chain runs
26 -> 23 -> 28 -> 20 -> 20 -> 31 -> 31 reads and the slot-1 chain 16 -> 13 -> 14 -> 8 -> 14 -> 13
beside it.

## Why it does not condense to two strands

Because the read evidence genuinely branches, and that is the interesting part rather than a defect
of the script. In the same window, a link carrying 6 reads crosses from the slot-0 chain to the
slot-1 chain against 20 staying put, and one node splits its outgoing reads 5 against 7. Those are
sites where the read partition is muddy -- the same quantity `--split-min-q` gates on -- so pruning
them away would delete the thing worth looking at.

`--min-edge-frac` (default 0.15) is the "all reads continue or terminate" rule made robust: a link
is dropped only when it is a minority continuation at **both** of the sides it joins, so it adapts
to local depth instead of guessing an absolute count, and a thin-but-only continuation survives.
Raise it to simplify the picture, lower it to see every branch. `--min-edge-reads` is a hard floor
and `--min-unitig-pins` drops fragments.

## chr20, hom-split anchors, at the defaults

| | segments | links | file |
|---|---|---|---|
| whole contig | 94,541 | 170,350 | 8.5 MB -- Bandage will open it, slowly |
| whole contig, `--min-unitig-pins 4` | 30,935 | 17,479 | 1.7 MB -- the practical whole-contig view |
| `chr20:26-31 Mb`, `--min-unitig-pins 2` | 3,882 | 5,670 | 293 KB -- the pericentromeric band |

654,830 pins condense 6.9x. Largest unitig 1,161 pins spanning 171.6 kb; median span 1.53 kb.
Class balance is even, 45,231 hap0 against 45,479 hap1, with 3,206 collapsed and 625 haploid.

The pericentromeric band is visibly worse, which is consistent with everything else about it: median
unitig span 0.28 kb against 1.53 kb contig-wide, and haploid nodes are 7.0% of the band against 0.7%
contig-wide -- a 10x enrichment.

## Annotations

`LN:i:` is the reference span in bases of the chain a segment represents, falling back to the pin
count for off-reference chains with no coordinate. `DP:f:` is mean read depth over its pins. The CSV
carries `Name,Class,Pins,SpanBp,Depth,Locus`; `Locus` reads `off-ref` for a chain no reference path
positions, which is a quick way to find the sites that exist only in the anchors.
