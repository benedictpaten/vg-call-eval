#!/usr/bin/env bash
# One-off whole-genome read database for tier 2 (plan §9.7).
#
# sort -> construct through a FIFO, so the ~199 GB uncompressed sorted GAF is
# never written to disk. Verified to work: gaf-base construct reads sequentially.
set -o pipefail
cd "$(dirname "$0")/.."

GBZ=data/hprc-v2.1-mc-chm13-eval.HG002.gbz
GAF=data/hprc-v2.1-mc-chm13-eval.HG002.gaf.gz
OUT=work/reads.gaf.db
FIFO=work/sorted.fifo
GB=/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-vg-planning/65260683-9f63-4117-b57b-52dccfd3b871/scratchpad/gbztools/bin

rm -f "$FIFO" "$OUT"
mkfifo "$FIFO"

echo "[$(date +%H:%M:%S)] starting gaf-base sort -> construct"
"$GB/gaf-base" sort "$GAF" -o "$FIFO" -p > work/sort.log 2>&1 &
SORT_PID=$!

"$GB/gaf-base" construct "$FIFO" -r "$GBZ" -o "$OUT" --overwrite > work/construct.log 2>&1
CONSTRUCT_RC=$?

wait $SORT_PID; SORT_RC=$?
rm -f "$FIFO"

echo "[$(date +%H:%M:%S)] sort rc=$SORT_RC construct rc=$CONSTRUCT_RC"
if [ $CONSTRUCT_RC -ne 0 ] || [ $SORT_RC -ne 0 ]; then
    echo "FAILED"; tail -5 work/sort.log work/construct.log; exit 1
fi
ls -la "$OUT" | awk '{printf "reads.gaf.db: %.2f GB\n", $5/1073741824}'
tail -4 work/construct.log
echo "DONE"
