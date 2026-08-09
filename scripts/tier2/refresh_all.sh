#!/usr/bin/env bash
# Regenerate every tier-2 number from one vg build: all five arms on all four
# datasets, both benchmarks, the size-matched control, then the pages.
#
# Why all five arms and not just the three that changed. The Poisson arms are
# untouched by a read-likelihood change, so re-running them looks like waste -- but
# the accuracy pages put all five in one table, and a table whose rows come from
# different builds is the failure this harness exists to prevent. It has already
# happened once here. One build, one pass, or the numbers are not comparable.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
THREADS=${THREADS:-5}

# run_arms.py shells out to gbz-base by name and does not take a binary path, so it
# has to be on PATH. A non-interactive shell does not inherit the login PATH, and the
# failure is a clean 40 s error on the read arms only -- a half-finished matrix that
# looks like a caller problem.
if ! command -v gbz-base >/dev/null; then
    GBZ=$(find /private/tmp/claude-501 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
    [ -n "$GBZ" ] || { echo "gbz-base not found" >&2; exit 1; }
    export PATH="$(dirname "$GBZ"):$PATH"
fi
echo "vg:        $VG"
echo "gbz-base:  $(command -v gbz-base)"

# label : work subdir : contig : gbz-base db : reads db
DATASETS=${DATASETS:-"chr20-4hap:tier2-chr20:chr20:graph.gbz.db:reads.gaf.db
chr20-34hap:tier2-chr20-hap32:chr20:graph.hap32.gbz.db:reads.hap32.gaf.db
chr6-4hap:tier2-chr6:chr6:graph.gbz.db:reads.gaf.db
chr6-34hap:tier2-chr6-hap32:chr6:graph.hap32.gbz.db:reads.hap32.gaf.db"}

while IFS=: read -r label sub contig gbzdb gafdb; do
    [ -n "$label" ] || continue
    W="$REPO/work/$sub"
    echo "### $label arms"
    python3 "$REPO/scripts/tier2/run_arms.py" \
        --vg "$VG" --contig "$contig" --threads "$THREADS" \
        --graph "$W/${contig}_0_${contig}.gbz" --pack "$W/${contig}.pack" \
        --gaf-base "$REPO/work/$gafdb" --gbz-base "$REPO/work/$gbzdb" \
        --reference "$W/${contig}.fa" \
        --truth-vcf "$W/truth.${contig}.smvar.vcf.gz" \
        --truth-bed "$W/truth.${contig}.smvar.bed" \
        --out "$W/results"

    echo "### $label truvari"
    python3 "$REPO/scripts/tier2/truvari_sv.py" --contig "$contig" --work "$W" \
        --label "$label" --arms poisson poisson-z readlik readlik-nomismap readlik-z

    echo "### $label size-matched"
    python3 "$REPO/scripts/tier2/size_matched.py" --results "$W/results" \
        --truth-vcf "$W/truth.${contig}.smvar.vcf.gz" \
        --truth-bed "$W/truth.${contig}.smvar.bed" \
        --reference "$W/${contig}.fa" --threads "$THREADS"
done <<< "$DATASETS"

echo "### pages"
for contig in chr20 chr6; do
    python3 "$REPO/scripts/tier2/report.py"         --contig "$contig"
    python3 "$REPO/scripts/tier2/compare_graphs.py" --contig "$contig"
done
echo MARKER_REFRESH_ALL_DONE
