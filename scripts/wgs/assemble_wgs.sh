#!/usr/bin/env bash
# Assemble the per-contig calls into one genome-wide VCF and one mosaic.
#
# The mosaic can legitimately be concatenated under a single header naming the *whole-genome* GBZ,
# even though every segment was produced from a per-contig chunk. `vg chunk` preserves whole-genome
# node IDs -- chr20's segments start around node 114.8M, not renumbered from 1 -- so the node
# anchors address the same nodes in the full graph. Each chunk's own header names only its chunk,
# which would be misleading in a genome-wide file, so it is replaced rather than kept.
set -euo pipefail
cd "$(dirname "$0")/../.."

W=${W:-work/wgs}
OUT=${OUT:-work/wgs/HG002}
GBZ=${GBZ:-data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz}
SAMPLE=${SAMPLE:-HG002}
CONTIGS=${CONTIGS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"}

step() { echo "[$(date +%H:%M:%S)] $*"; }

# --- VCF -------------------------------------------------------------------
list=()
for C in $CONTIGS; do
    V="$W/$C/$C.vcf"
    [ -s "$V" ] || { echo "missing calls for $C" >&2; exit 1; }
    if [ ! -s "$V.gz.tbi" ]; then
        bgzip -f -c "$V" > "$V.gz"
        tabix -f -p vcf "$V.gz"
    fi
    list+=("$V.gz")
done

step "concat ${#list[@]} contigs"
# -a rather than plain concat: each contig's header declares only its own contig, so the files are
# not a simple ordered split of one header.
bcftools concat -a -Oz -o "$OUT.vcf.gz" "${list[@]}"
bcftools index -f -t "$OUT.vcf.gz"
echo "  $(bcftools view -H "$OUT.vcf.gz" | wc -l | tr -d ' ') records -> $OUT.vcf.gz"

# --- mosaic ----------------------------------------------------------------
step "concat mosaics"
{
    printf '#mosaic-version\t1\n'
    printf '#graph\t%s\n' "$GBZ"
    printf '#sample\t%s\n' "$SAMPLE"
    printf '#decoding\tconstrained-viterbi\n'
    printf '#note\tsegments are maximal runs on one panel haplotype; walk the haplotype from start_node to end_node to reconstruct it. * means the panel does not explain that strand there. Haploid contigs carry strand 0 only.\n'
    printf '#H\tcontig\tstrand\tref_start\tref_end\tstart_node\tend_node\thap_index\thaplotype\tsites\n'
    for C in $CONTIGS; do
        M="$W/$C/$C.mosaic.tsv"
        [ -s "$M" ] || { echo "missing mosaic for $C" >&2; exit 1; }
        grep -v "^#" "$M"
    done
} > "$OUT.mosaic.tsv"
echo "  $(grep -vc '^#' "$OUT.mosaic.tsv") segments -> $OUT.mosaic.tsv ($(du -h "$OUT.mosaic.tsv" | cut -f1))"

# Structural check: every strand's segments must tile that contig's sites exactly, which is the
# property that makes the file a description of a genome rather than a list of observations.
awk -F'\t' '!/^#/ {n[$2"/"$3] += $10} END {for (k in n) print k, n[k]}' "$OUT.mosaic.tsv" \
    | sort > "$OUT.mosaic.tiling.txt"
echo "  per contig/strand site totals in $OUT.mosaic.tiling.txt"
echo "ASSEMBLE_DONE"
