#!/usr/bin/env bash
# SUPERSEDED, AND NO LONGER RUNNABLE. This compared whole-traversal against unique-content
# length weights. Unique content won, became the default, and the flags this script passes have
# since been removed from vg: --length-weighted-mixture went when the weight became the default,
# and --length-weight-whole-traversal went once the comparison was settled. Kept as the record of
# how the measurement was run, not as something to run. The result is in
# docs/tier2-sv-errors.md and in the set_unique_lengths comment in vg.
# Call the read-likelihood arm under each mixture weight, timed, on one graph.
#
# Three arms per dataset:
#   default   flat 1/|G|                              (the shipped model)
#   whole     w_h from whole traversal length         (the first weighted version)
#   unique    w_h from sequence unique to each allele (the sharpened version)
#
# The baseline is re-run rather than reusing the logged 242 s, because that number
# came from a different build on a differently loaded machine and a runtime claim
# built on it would not mean anything.
set -euo pipefail

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
THREADS=${THREADS:-5}

GBZ_BASE_BIN=${GBZ_BASE_BIN:-$(command -v gbz-base || true)}
if [ -z "$GBZ_BASE_BIN" ]; then
    GBZ_BASE_BIN=$(find /private/tmp/claude-501 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
fi
[ -n "$GBZ_BASE_BIN" ] || { echo "gbz-base not found" >&2; exit 1; }

run() {
    local ds=$1 contig=$2 gbzdb=$3 gafdb=$4 label=$5
    shift 5
    local w="$REPO/work/$ds"
    local out="$w/results/readlik-z-${label}.vcf.gz"
    if [ -s "${out}.tbi" ] && [ "${FORCE:-0}" != "1" ]; then
        echo "skip $ds $label (exists)"; return
    fi
    echo "=== $ds $label ==="
    /usr/bin/time -l "$VG" call "$w/${contig}_0_${contig}.gbz" \
        -p "CHM13#0#${contig}" -t "$THREADS" --read-likelihood -z \
        --gaf-base "$REPO/work/$gafdb" --gbz-base "$REPO/work/$gbzdb" \
        --gaf-base-binary "$GBZ_BASE_BIN" "$@" \
        2> "$w/results/readlik-z-${label}.time.log" \
        | bgzip -c > "$out"
    tabix -f -p vcf "$out"
    local secs rss
    secs=$(awk '/real/{print $1}' "$w/results/readlik-z-${label}.time.log" | tail -1)
    rss=$(awk '/maximum resident set size/{printf "%.1f", $1/1073741824}' "$w/results/readlik-z-${label}.time.log")
    echo "  $ds $label: $(bcftools index -n "$out") variants, ${secs}s, ${rss} GB peak"
}

for spec in "${DATASETS:-tier2-chr6:chr6:graph.gbz.db:reads.gaf.db}"; do
    IFS=: read -r ds contig gbzdb gafdb <<< "$spec"
    run "$ds" "$contig" "$gbzdb" "$gafdb" baseline
    run "$ds" "$contig" "$gbzdb" "$gafdb" whole  --length-weighted-mixture --length-weight-whole-traversal
    run "$ds" "$contig" "$gbzdb" "$gafdb" unique --length-weighted-mixture
done
echo MARKER_MIXTURE_RUNS_DONE
