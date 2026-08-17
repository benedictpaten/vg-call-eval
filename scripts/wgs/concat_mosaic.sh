#!/usr/bin/env bash
# Concatenate per-contig mosaic v2 files into one genome-wide mosaic.
#
# This is not `cat`, and mosaic v2 is why. Two columns are meaningful only relative to the graph
# that produced them:
#
#   * `hap_index` is the haplotype's position in *that chunk's* GBWT metadata, and the chunks do not
#     agree on an ordering -- chr20 puts recombination#10 at 13, another contig will not. Emitting
#     24 files' worth of rows under one #haplotype table would silently relabel haplotypes. We
#     reindex on the `haplotype` (sample#phase) column, which is portable by construction; that is
#     what it is for.
#   * `gbwt_offset` is a rank among the sequences visiting that node in the chunk's GBWT. The
#     whole-genome GBWT has more sequences at the same node, so the same offset addresses a
#     different path there. Not portable, and not made portable by renaming.
#
# So the output names each contig's *own* graph in a #contig table rather than claiming a single
# whole-genome GBZ, and says which columns resolve against it. What does survive concatenation is
# `start_node`/`end_node`: `vg chunk` preserves whole-genome node IDs -- chr20's segments start
# around node 114.8M, not renumbered from 1 -- so the node anchors address the same nodes in the
# full graph. That is why they, and not the reference coordinates, are the authoritative anchors.
#
# Usage: concat_mosaic.sh SAMPLE OUT.tsv CONTIG:FILE [CONTIG:FILE ...]
set -euo pipefail

[ $# -ge 3 ] || { echo "usage: $0 SAMPLE OUT.tsv CONTIG:FILE [CONTIG:FILE ...]" >&2; exit 2; }
SAMPLE=$1; OUT=$2; shift 2

for spec in "$@"; do
    F=${spec#*:}
    [ -s "$F" ] || { echo "missing or empty mosaic: $F" >&2; exit 1; }
done

# Two passes over each file: collect the per-contig graph/reference and the union of the haplotype
# panels, then re-emit the rows with hap_index remapped onto that union.
printf '%s\n' "$@" | tr ':' ' ' | awk -v sample="$SAMPLE" '
BEGIN { g = 0; n = 0 }   # awk would otherwise index gname[""] on the first haplotype
{ contig[++n] = $1; file[$1] = $2 }
END {
    for (i = 1; i <= n; i++) {
        C = contig[i]
        while ((getline line < file[C]) > 0) {
            if (line !~ /^#/) continue
            split(line, f, "\t")
            if (f[1] == "#mosaic-version" && f[2] != "2") {
                print "expected mosaic-version 2 in " file[C] ", got " f[2] > "/dev/stderr"; exit 1
            }
            if (f[1] == "#graph")     graph[C] = f[2]
            if (f[1] == "#reference") ref[C]   = f[2]
            if (f[1] == "#haplotype" && !(f[3] in gidx)) { gidx[f[3]] = g; gname[g] = f[3]; g++ }
        }
        close(file[C])
        if (!(C in graph)) { print "no #graph line in " file[C] > "/dev/stderr"; exit 1 }
        if (!(C in ref))   { print "no #reference line in " file[C] > "/dev/stderr"; exit 1 }
    }
    print "#mosaic-version\t2"
    print "#sample\t" sample
    print "#decoding\tconstrained-viterbi"
    print "#note\tgenome-wide file: assembled from per-contig runs, so there is no single graph. gbwt_node/gbwt_offset resolve against the #contig graph for that row s own contig and only that one, because an offset is a rank among the sequences at a node and the whole-genome GBWT has more of them."
    print "#note\tstart_node/end_node are whole-genome node IDs and stay valid against any graph containing them; ref_start/ref_end are advisory, in that contig s #contig reference coordinate system."
    print "#note\tsegments are maximal runs on one panel haplotype; walk the haplotype from start_node to end_node to reconstruct it. * means the panel does not explain that strand there. Haploid contigs carry strand 0 only."
    print "#note\ta segment never spans a GBWT fragment boundary, so one position walks the whole of it; a haplotype in several fragments yields several segments."
    print "#note\thap_index is reindexed onto the union panel below and is NOT any contig s own index. haplotype (sample#phase) is the portable identifier."
    print "#C\tcontig\tgraph\treference"
    for (i = 1; i <= n; i++) print "#contig\t" contig[i] "\t" graph[contig[i]] "\t" ref[contig[i]]
    for (j = 0; j < g; j++) print "#haplotype\t" j "\t" gname[j]
    print "#H\tcontig\tstrand\tref_start\tref_end\tstart_node\tend_node\thap_index\thaplotype\tsites\tgbwt_node\tgbwt_offset"
    for (i = 1; i <= n; i++) {
        C = contig[i]
        while ((getline line < file[C]) > 0) {
            if (line ~ /^#/) continue
            # Data rows carry a leading "H" marker, so their fields line up with the #H header:
            # field 8 is hap_index and field 9 the haplotype name, in both.
            nf = split(line, f, "\t")
            if (nf != 12 || f[1] != "H") {
                print "expected an H row of 12 columns in " file[C] ", got " nf " starting " f[1] > "/dev/stderr"
                exit 1
            }
            if (f[9] != "*") {
                if (!(f[9] in gidx)) {
                    print "haplotype " f[9] " on " C " is absent from that contig s own panel" > "/dev/stderr"
                    exit 1
                }
                f[8] = gidx[f[9]]
            }
            out = f[1]
            for (k = 2; k <= 12; k++) out = out "\t" f[k]
            print out
        }
        close(file[C])
    }
}' > "$OUT"
