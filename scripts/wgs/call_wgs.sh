#!/usr/bin/env bash
# Whole-genome calling, one contig at a time, with the right ploidy for a male sample.
#
# One contig per invocation is what keeps this on a laptop. The caller buffers every emitted record
# and every linkage site until the chain is resolved, so peak memory scales with the contig, not
# with the genome: chr6 peaked at 3.4 GB, and chr1 is 1.5x that.
#
# Ploidy. HG002 is male, so chrY is haploid throughout and chrX is haploid *except* in the
# pseudoautosomal regions, where X and Y recombine and the sample carries two copies. vg's ploidy
# is per contig, so PAR cannot be expressed in one run -- chrX is therefore called twice and
# spliced. The boundaries below are not looked up: they are read off the T2T-Q100 truth itself,
# which is haploid at 120,704 chrX records and diploid at 11,683, in exactly two blocks.
#
# The seam this creates at each PAR boundary is an artefact of the run, not biology: linkage and
# the mosaic restart there. That is worth knowing when reading chrX's phasing.
set -euo pipefail
cd "$(dirname "$0")/../.."

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
export PATH="$(dirname "$VG"):$PATH"
W=${W:-work/wgs}
READS_DB=${READS_DB:-work/reads.hap32.gaf.db}
GRAPH_DB=${GRAPH_DB:-work/graph.hap32.gbz.db}
REF_SAMPLE=${REF_SAMPLE:-CHM13}
SAMPLE=${SAMPLE:-HG002}
THREADS=${THREADS:-5}
CONTIGS=${CONTIGS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"}

# PAR1 and PAR2 on CHM13 chrX, as the truth set draws them.
PAR1_END=2394370
PAR2_START=153926003

if ! command -v gbz-base >/dev/null; then
    GBZ=$(find /private/tmp/claude-501 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
    [ -n "$GBZ" ] || { echo "gbz-base not found" >&2; exit 1; }
    export PATH="$(dirname "$GBZ"):$PATH"
fi
step() { echo "[$(date +%H:%M:%S)] $*"; }

call_one() {   # contig ploidy outprefix
    local C=$1 PLOIDY=$2 OUT=$3
    local D="$W/$C"
    /usr/bin/time -l "$VG" call "$D/$C.gbz" \
        -p "${REF_SAMPLE}#0#${C}" -s "$SAMPLE" -d "$PLOIDY" -t "$THREADS" --progress \
        --read-likelihood --phased --mosaic-out "$OUT.mosaic.tsv" \
        --gaf-base "$READS_DB" --gbz-base "$GRAPH_DB" \
        > "$OUT.vcf" 2> "$OUT.log"
    local secs rss
    secs=$(grep -E "^ *[0-9.]+ real" "$OUT.log" | awk '{print $1}')
    rss=$(grep "maximum resident set size" "$OUT.log" | awk '{printf "%.1f", $1/2^30}')
    echo "    $(grep -vc '^#' "$OUT.vcf") records, ${secs}s, ${rss} GB peak"
}

for C in $CONTIGS; do
    D="$W/$C"
    [ -s "$D/$C.gbz" ] || { echo "missing subgraph for $C -- run prep_wgs.sh" >&2; exit 1; }

    if [ -s "$D/$C.done" ]; then
        echo "[$(date +%H:%M:%S)] $C: already called, skipping"
        continue
    fi

    case $C in
        chrY)
            step "$C: haploid"
            call_one "$C" 1 "$D/$C"
            ;;
        chrX)
            # Two passes over the same subgraph, spliced on the PAR boundaries. Both emit the same
            # sites, so the splice is a region filter rather than a merge.
            step "$C: haploid pass (non-PAR)"
            call_one "$C" 1 "$D/$C.hap"
            step "$C: diploid pass (PAR)"
            call_one "$C" 2 "$D/$C.dip"

            step "$C: splice PAR"
            for f in "$D/$C.hap" "$D/$C.dip"; do
                bgzip -f -c "$f.vcf" > "$f.vcf.gz"
                tabix -f -p vcf "$f.vcf.gz"
            done
            bcftools view -r "${C}:1-${PAR1_END},${C}:${PAR2_START}-" \
                -Oz -o "$D/$C.par.vcf.gz" "$D/$C.dip.vcf.gz"
            bcftools view -t "^${C}:1-${PAR1_END},^${C}:${PAR2_START}-" \
                -Oz -o "$D/$C.nonpar.vcf.gz" "$D/$C.hap.vcf.gz"
            bcftools index -f -t "$D/$C.par.vcf.gz"
            bcftools index -f -t "$D/$C.nonpar.vcf.gz"
            bcftools concat -a -Oz -o "$D/$C.vcf.gz" "$D/$C.par.vcf.gz" "$D/$C.nonpar.vcf.gz"
            bcftools index -f -t "$D/$C.vcf.gz"
            gzcat "$D/$C.vcf.gz" > "$D/$C.vcf"
            # Mosaic: PAR from the diploid pass, the rest from the haploid one.
            {
                grep "^#" "$D/$C.hap.mosaic.tsv"
                awk -F'\t' -v e=$PAR1_END -v s=$PAR2_START \
                    '!/^#/ && ($4 <= e || $4 >= s)' "$D/$C.dip.mosaic.tsv"
                awk -F'\t' -v e=$PAR1_END -v s=$PAR2_START \
                    '!/^#/ && $4 > e && $4 < s' "$D/$C.hap.mosaic.tsv"
            } > "$D/$C.mosaic.tsv"
            echo "    spliced: $(grep -vc '^#' "$D/$C.vcf") records, $(grep -vc '^#' "$D/$C.mosaic.tsv") mosaic segments"
            ;;
        *)
            step "$C: diploid"
            call_one "$C" 2 "$D/$C"
            ;;
    esac
    touch "$D/$C.done"
done
echo "CALL_DONE"
