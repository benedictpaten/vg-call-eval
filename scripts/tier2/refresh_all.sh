#!/usr/bin/env bash
# Regenerate every tier-2 number from one vg build: all five arms on all four
# datasets, both benchmarks, the size-matched control, then the pages.
#
# Why all five arms and not just the three that changed. The Poisson arms are
# untouched by a read-likelihood change, so re-running them looks like waste -- but
# the accuracy pages put all five in one table, and a table whose rows come from
# different builds is the failure this harness exists to prevent. It has already
# happened once here. One build, one pass, or the numbers are not comparable.
#
# CANARY=1 gets the same guarantee for a fifth of the cost. The two Poisson arms are 43%
# of the run time -- 36 minutes of 86 -- so instead of re-running them on faith, run one
# cheap Poisson arm and compare it byte for byte against the cached copy. Identical means
# the cached rows carry over; different means something touched shared code, which blind
# re-running would have absorbed without telling you. Verification, not assumption.
#
# JOBS=2 runs that many datasets at once. Peak RSS is 6-8 GB per vg call, so two fit in
# 32 GB with room to spare; the arms within a dataset stay serial because they share the
# results directory.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
THREADS=${THREADS:-5}
# Extra flags for the three read-likelihood arms only, so a caller-side change can be
# measured across the whole matrix before it becomes a default. Empty by default.
READLIK_EXTRA=${READLIK_EXTRA:-}
# Flags that only work with -z, so they go to readlik-z alone. --linkage-weight is refused by
# vg call without -z, so putting it in READLIK_EXTRA kills the two support-enumeration arms.
READLIK_Z_EXTRA=${READLIK_Z_EXTRA:-}

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

ARM_FILTER=()
if [ "${CANARY:-0}" = "1" ]; then
    if python3 "$REPO/scripts/tier2/canary.py" --vg "$VG" --threads "$THREADS"; then
        ARM_FILTER=(--only readlik readlik-nomismap readlik-z)
        echo "canary passed: reusing cached Poisson arms"
    else
        echo "canary failed or unavailable: running every arm"
    fi
fi
JOBS=${JOBS:-1}

# label : work subdir : contig : gbz-base db : reads db
DATASETS=${DATASETS:-"chr20-4hap:tier2-chr20:chr20:graph.gbz.db:reads.gaf.db
chr20-34hap:tier2-chr20-hap32:chr20:graph.hap32.gbz.db:reads.hap32.gaf.db
chr6-4hap:tier2-chr6:chr6:graph.gbz.db:reads.gaf.db
chr6-34hap:tier2-chr6-hap32:chr6:graph.hap32.gbz.db:reads.hap32.gaf.db"}

run_dataset() {
    local label=$1 sub=$2 contig=$3 gbzdb=$4 gafdb=$5
    W="$REPO/work/$sub"
    echo "### $label arms"
    python3 "$REPO/scripts/tier2/run_arms.py" \
        --vg "$VG" --contig "$contig" --threads "$THREADS" \
        --graph "$W/${contig}_0_${contig}.gbz" --pack "$W/${contig}.pack" \
        --gaf-base "$REPO/work/$gafdb" --gbz-base "$REPO/work/$gbzdb" \
        --reference "$W/${contig}.fa" \
        --truth-vcf "$W/truth.${contig}.smvar.vcf.gz" \
        --truth-bed "$W/truth.${contig}.smvar.bed" \
        --out "$W/results" ${READLIK_EXTRA:+--readlik-extra "$READLIK_EXTRA"} \
        ${READLIK_Z_EXTRA:+--readlik-z-extra "$READLIK_Z_EXTRA"} \
        ${ARM_FILTER[@]+"${ARM_FILTER[@]}"}

    echo "### $label truvari"
    python3 "$REPO/scripts/tier2/truvari_sv.py" --contig "$contig" --work "$W" \
        --label "$label" --arms poisson poisson-z readlik readlik-nomismap readlik-z

    echo "### $label size-matched"
    python3 "$REPO/scripts/tier2/size_matched.py" --results "$W/results" \
        --truth-vcf "$W/truth.${contig}.smvar.vcf.gz" \
        --truth-bed "$W/truth.${contig}.smvar.bed" \
        --reference "$W/${contig}.fa" --threads "$THREADS"
}

# PIDs of the backgrounded datasets, waited on individually at the end. A failure inside a
# background job is otherwise invisible: `set -e` does not cross the boundary, so each exit
# status has to be collected explicitly.
#
# Not `wait -n`. That needs bash 4.3, and macOS ships 3.2 -- where it fails as an invalid
# option, which `|| exit 1` turns into an exit *after* the pool is already full, leaving
# orphaned `vg call` children running against a dead parent. Polling `jobs -rp` is portable, and
# waiting on every PID rather than whichever one `wait -n` reaped is strictly better failure
# detection anyway.
PIDS=""
while IFS=: read -r label sub contig gbzdb gafdb; do
    [ -n "$label" ] || continue
    if [ "$JOBS" -gt 1 ]; then
        run_dataset "$label" "$sub" "$contig" "$gbzdb" "$gafdb" \
            > "$REPO/work/refresh-$label.log" 2>&1 &
        PIDS="$PIDS $!"
        while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
    else
        run_dataset "$label" "$sub" "$contig" "$gbzdb" "$gafdb"
    fi
done <<< "$DATASETS"
if [ "$JOBS" -gt 1 ]; then
    FAILED=0
    for pid in $PIDS; do wait "$pid" || FAILED=1; done
    for f in "$REPO"/work/refresh-*.log; do echo "--- $f"; tail -3 "$f"; done
    [ "$FAILED" -eq 0 ] || { echo "a dataset failed; see work/refresh-*.log" >&2; exit 1; }
fi

echo "### pages"
for contig in chr20 chr6; do
    python3 "$REPO/scripts/tier2/report.py"         --contig "$contig"
    python3 "$REPO/scripts/tier2/compare_graphs.py" --contig "$contig"
done
echo MARKER_REFRESH_ALL_DONE
