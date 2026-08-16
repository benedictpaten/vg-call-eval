#!/usr/bin/env bash
# Whole-genome calling, one contig at a time, with the right ploidy for a male sample.
#
# One contig per invocation is what keeps this on a laptop. The caller buffers every emitted record
# and every linkage site until the chain is resolved, so peak memory scales with the contig, not
# with the genome: chr6 peaked at 3.4 GB, and chr1 is 1.5x that.
#
# Ploidy. HG002 is male, so chrY is haploid throughout and chrX is haploid *except* in the
# pseudoautosomal regions, where X and Y recombine and the sample carries two copies. That is now
# one run with --ploidy-bed chrX.par.bed, which takes ploidy per region.
#
# It used to be two runs spliced on the PAR boundaries, because vg's ploidy was per contig. The
# splice was also wrong: `bcftools view -t "^chrX:153926003-"` does not exclude an open-ended range
# the way `-r` includes one, so 190 haploid records leaked into PAR2 and the concatenated VCF
# carried 190 duplicated positions at contradicting ploidies. Nothing published was affected --
# the T2T-Q100 confident regions end at 153,910,814, before PAR2 begins, so those records were
# never scored -- but it is exactly the kind of error a splice invites and a ploidy BED cannot make.
#
# Linkage and the mosaic still break at each ploidy boundary. That is a property of the boundary
# rather than of the splice: there is no haplotype correspondence to carry across a ploidy change.
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

call_one_bed() {   # contig ploidy outprefix [ploidy-bed]
    local C=$1 PLOIDY=$2 OUT=$3 BED=${4:-}
    local D="$W/$C"
    local extra=()
    [ -n "$BED" ] && extra=(--ploidy-bed "$BED")
    # ${extra[@]+"${extra[@]}"}, not "${extra[@]}". Under `set -u`, bash 3.2 -- which is what
    # /bin/bash still is on macOS -- treats an empty array expansion as an unbound variable and
    # aborts. Every diploid contig died this way on the first whole-genome rerun: only chrX passes
    # a BED, so only chrX had exercised this function since it gained the argument.
    /usr/bin/time -l "$VG" call "$D/$C.gbz" \
        -p "${REF_SAMPLE}#0#${C}" -s "$SAMPLE" -d "$PLOIDY" -t "$THREADS" --progress \
        --read-likelihood --phased --mosaic-out "$OUT.mosaic.tsv" ${extra[@]+"${extra[@]}"} \
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

    # Resume, but only past work done by *this* binary. A marker that just says "called"
    # will happily carry a result from before a fix across a rebuild, and it did: the
    # coverage sweep kept pre-fix chrX arms that then scored as though they were the fixed
    # caller, identifiable only by their file timestamps.
    if [ -f "$D/$C.done" ]; then
        if [ "$D/$C.done" -nt "$VG" ]; then
            echo "[$(date +%H:%M:%S)] $C: already called by this binary, skipping"
            continue
        fi
        echo "[$(date +%H:%M:%S)] $C: called by an older binary, recalling"
        rm -f "$D/$C.done"
    fi

    case $C in
        chrY)
            step "$C: haploid"
            call_one_bed "$C" 1 "$D/$C"
            ;;
        chrX)
            # One pass, ploidy from the BED. Both PAR blocks come out diploid and the interior
            # haploid, with the boundaries exact -- verified against the old two-pass output:
            # every one of the 97,068 non-PAR sites matches it genotype for genotype.
            step "$C: haploid with diploid PAR (--ploidy-bed)"
            call_one_bed "$C" 1 "$D/$C" "$(dirname "$0")/chrX.par.bed"
            ;;
        *)
            step "$C: diploid"
            call_one_bed "$C" 2 "$D/$C"
            ;;
    esac
    grep -vc "^#" "$D/$C.vcf" > "$D/$C.done"
done
echo "CALL_DONE"
