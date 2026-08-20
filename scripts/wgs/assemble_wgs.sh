#!/usr/bin/env bash
# Assemble the per-contig calls into one genome-wide VCF and one mosaic.
#
# Concatenating mosaics is not just `cat`. Two of the columns are meaningful
# only relative to the graph that produced them:
#
#   * `hap_index` is the haplotype's position in *that chunk's* GBWT metadata, and the chunks do not
#     agree on an ordering. Emitting 24 files' worth of rows under one #haplotype table would
#     silently relabel haplotypes. We reindex on the `haplotype` (sample#phase) column instead,
#     which is portable by construction -- that is what it is for.
#   * `gbwt_offset` is a rank among the sequences visiting that node in the chunk's GBWT. The
#     whole-genome GBWT has more sequences at the same node, so the same offset addresses a
#     different path there. The offsets are not portable and cannot be made portable by renaming.
#
# So the genome-wide file names each contig's *own* graph in a #contig table rather than claiming a
# single whole-genome GBZ, and says which columns resolve against it. What does survive
# concatenation is `start_node`/`end_node`: `vg chunk` preserves whole-genome node IDs -- chr20's
# segments start around node 114.8M, not renumbered from 1 -- so the node anchors address the same
# nodes in the full graph, and they are the authoritative anchors precisely because of that.
set -euo pipefail
cd "$(dirname "$0")/../.."

W=${W:-work/wgs}
# Derived from W by default, never fixed at work/wgs. Setting W alone used to redirect the inputs
# while the output stayed on work/wgs -- so assembling any other arm silently overwrote the
# --no-nested baseline's whole-genome VCF and mosaic with the other arm's data, in place, with a
# success message naming the path it had just clobbered. Override OUT explicitly to split them.
OUT=${OUT:-$W/HG002}
GBZ=${GBZ:-data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz}
SAMPLE=${SAMPLE:-HG002}
CONTIGS=${CONTIGS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"}

step() { echo "[$(date +%H:%M:%S)] $*"; }

# --- VCF -------------------------------------------------------------------
list=()
for C in $CONTIGS; do
    V="$W/$C/$C.vcf"
    [ -s "$V" ] || { echo "missing calls for $C" >&2; exit 1; }
    # Freshness, not existence. Keying on "is there a .tbi" silently concatenated the previous
    # run's compressed VCFs after a full recall: every .vcf was new, every .vcf.gz was a week old,
    # and the refreshed whole-genome numbers came out byte-identical to the run they were meant to
    # replace. Nothing in the output said so.
    if [ ! -s "$V.gz.tbi" ] || [ "$V" -nt "$V.gz" ]; then
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
specs=()
for C in $CONTIGS; do specs+=("$C:$W/$C/$C.mosaic.tsv"); done
bash scripts/wgs/concat_mosaic.sh "$SAMPLE" "$OUT.mosaic.tsv" "${specs[@]}"
echo "  $(grep -vc '^#' "$OUT.mosaic.tsv") segments -> $OUT.mosaic.tsv ($(du -h "$OUT.mosaic.tsv" | cut -f1))"

# Structural check: every strand's segments must tile that contig's sites exactly, which is the
# property that makes the file a description of a genome rather than a list of observations.
awk -F'\t' '!/^#/ {n[$2"/"$3] += $10} END {for (k in n) print k, n[k]}' "$OUT.mosaic.tsv" \
    | sort > "$OUT.mosaic.tiling.txt"
echo "  per contig/strand site totals in $OUT.mosaic.tiling.txt"
echo "ASSEMBLE_DONE"
