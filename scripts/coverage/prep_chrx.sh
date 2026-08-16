#!/usr/bin/env bash
# Extract chrX reads so the titration has a haploid axis alongside chr20's diploid one.
#
# Why chrX and not simply chr20 called at -d 1: chr20 is genuinely diploid in HG002, so calling it
# haploid would be scored against a truth that disagrees everywhere, and the result would measure
# the mismatch rather than the model. chrX outside the pseudoautosomal regions is really haploid in
# this sample and the T2T-Q100 truth is haploid there too (120,704 haploid records against 11,683
# diploid ones, in exactly two blocks).
#
# **The two axes line up on reads per haplotype, which is the unit that matters.** chr20 carries
# 30.3x across two haplotypes, so 15 reads per haplotype; a male chrX carries about half the
# genome-wide depth across one, so also about 15 per haplotype at full coverage. Choosing chrX
# levels at half the chr20 levels therefore puts the two series on a common x-axis, and any
# remaining difference between them is ploidy rather than depth. That is the whole point of having
# two series.
#
# Only the non-PAR part is used. The PAR is diploid and would mix the two ploidies inside one
# arm -- the same reason the whole-genome run calls chrX twice and splices.
set -euo pipefail
cd "$(dirname "$0")/../.."

W=${W:-work/coverage/chrX}
SRC=${SRC:-work/wgs/chrX}
GRAPH_DB=${GRAPH_DB:-work/graph.hap32.gbz.db}
READS_DB=${READS_DB:-work/reads.hap32.gaf.db}
VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
export PATH="$(dirname "$VG"):$PATH"
if ! command -v gbz-base >/dev/null; then
    d=$(find /private/tmp/claude-501 -maxdepth 8 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
    [ -n "$d" ] || { echo "gbz-base not found" >&2; exit 1; }
    export PATH="$(dirname "$d"):$PATH"
fi
step() { echo "[$(date +%H:%M:%S)] $*"; }

mkdir -p "$W"
[ -s "$SRC/chrX.gbz" ] || { echo "missing $SRC/chrX.gbz -- run the wgs prep" >&2; exit 1; }

# 1. Node ID list for the read extraction. gbz-base errors on a node absent from the graph and
#    this ID space is sparse, so the list has to be explicit (as prep_contig.sh found).
if [ ! -s "$W/chrX_all_nodes.txt" ]; then
    step "node ID list"
    vg convert -f "$SRC/chrX.gbz" | awk '$1=="S"{print $2}' | sort -n > "$W/chrX_all_nodes.txt"
fi
echo "  nodes: $(wc -l < "$W/chrX_all_nodes.txt")"

# 2. Reads, out of the whole-genome GAF-Base rather than by streaming the 27 GB gzip.
if [ ! -s "$W/chrX.reads.gaf" ]; then
    step "extract reads"
    python3 scripts/tier2/extract_reads_from_db.py \
        --nodes "$W/chrX_all_nodes.txt" \
        --gaf-base "$READS_DB" --gbz-base "$GRAPH_DB" \
        --out "$W/chrX.reads.gaf" \
        --tmp "/tmp/gafbase_extract_coverage_chrX.gaf"
fi
echo "  reads: $(wc -l < "$W/chrX.reads.gaf")"
step "PREP_CHRX_DONE"
