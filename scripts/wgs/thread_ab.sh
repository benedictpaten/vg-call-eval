#!/usr/bin/env bash
# Does lowering -t cost time, and does it save memory?
#
# The claim under test: the caller is I/O-bound -- every read-fetch window spawns a gbz-base and
# reopens a 22 GB SQLite database -- so it uses about one CPU whatever -t says, and the threads
# mainly buy per-thread read caches. If that holds, -t can be lowered nearly for free and the
# memory freed spends better on running more contigs at once.
#
# It is a hypothesis, not a measurement, which is why this exists. The failure mode if it is wrong
# is a scheduler that packs more contigs while each runs several times slower, which would look
# like a win on paper and lose on the clock.
#
# One contig, one build, machine otherwise idle, replicated and interleaved so drift spreads across
# the arms rather than landing on whichever ran last.
set -euo pipefail
cd "$(dirname "$0")/../.."

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
export PATH="$(dirname "$VG"):$PATH"
W=${W:-work/wgs}
C=${C:-chr20}
READS_DB=${READS_DB:-work/reads.hap32.gaf.db}
GRAPH_DB=${GRAPH_DB:-work/graph.hap32.gbz.db}
REPS=${REPS:-2}
THREAD_SET=${THREAD_SET:-"1 2 5"}
OUT=${OUT:-/tmp/thread_ab.tsv}

if ! command -v gbz-base >/dev/null; then
    GBZ=$(find /private/tmp/claude-501 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
    [ -n "$GBZ" ] || { echo "gbz-base not found" >&2; exit 1; }
    export PATH="$(dirname "$GBZ"):$PATH"
fi
[ -s "$W/$C/$C.gbz" ] || { echo "no subgraph for $C" >&2; exit 1; }

: > "$OUT"
for rep in $(seq 1 "$REPS"); do
    for t in $THREAD_SET; do
        /usr/bin/time -l "$VG" call "$W/$C/$C.gbz" \
            -p "CHM13#0#${C}" -s HG002 -d 2 -t "$t" \
            --read-likelihood --phased \
            --gaf-base "$READS_DB" --gbz-base "$GRAPH_DB" \
            > /tmp/thread_ab.vcf 2> /tmp/thread_ab.err
        secs=$(grep -E "^ *[0-9.]+ real" /tmp/thread_ab.err | awk '{print $1}')
        rss=$(grep "maximum resident set size" /tmp/thread_ab.err | awk '{printf "%.2f", $1/2^30}')
        cpu_u=$(grep -E "^ *[0-9.]+ real" /tmp/thread_ab.err | awk '{print $3}')
        cpu_s=$(grep -E "^ *[0-9.]+ real" /tmp/thread_ab.err | awk '{print $5}')
        # CPU seconds over wall seconds is the number that settles the argument: at 1.0 the
        # threads are doing nothing, whatever -t was set to.
        eff=$(awk -v u="$cpu_u" -v s="$cpu_s" -v w="$secs" 'BEGIN{printf "%.2f", (u+s)/w}')
        printf "%s\t%s\t%s\t%s\t%s\n" "$rep" "$t" "$secs" "$rss" "$eff" >> "$OUT"
        echo "rep$rep -t $t: ${secs}s  ${rss} GB  ${eff} CPU"
    done
done
echo "THREAD_AB_DONE  -> $OUT"
awk -F'\t' '{n[$2]++; s[$2]+=$3; m[$2]+=$4; c[$2]+=$5}
     END {printf "\n%-4s %10s %10s %10s\n", "-t", "mean s", "mean GB", "mean CPU";
          for (k in n) printf "%-4s %10.1f %10.2f %10.2f\n", k, s[k]/n[k], m[k]/n[k], c[k]/n[k]}' \
    "$OUT" | sort -n
