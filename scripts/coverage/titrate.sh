#!/usr/bin/env bash
# Build a GAF-Base database per coverage level, then call at each.
#
# Why a database per level and not `--gaf-reads` on the subsampled GAF. An in-memory read source
# answers site queries exactly and so has no fetch window; `local_read_rate` returns 0 for a source
# with no window, and the depth term switches itself off rather than inventing a rate (see
# allele_likelihood.hpp). A titration run that way would measure a model with `--depth-term`
# silently disabled -- which is precisely one of the parameters whose coverage behaviour we are
# here to measure. So: databases, at the cost of a build per level.
#
# The full-coverage arm reuses the whole-genome database rather than rebuilding one. It is the same
# reads, it is already built, and using it means the 30x arm is byte-for-byte the configuration
# that produced the existing tier-2 numbers -- which is what makes "full coverage must reproduce
# tier-2 exactly" a real check rather than a re-measurement of something slightly different.
#
# Reads databases are built against the *whole-genome* GBZ, and calls use the whole-genome
# gbz-base. Node IDs are global (vg chunk preserves them), so a per-contig graph database would buy
# nothing and would add a way for the two to disagree.
set -euo pipefail
cd "$(dirname "$0")/../.."

CONTIG=${CONTIG:-chr20}
W=${W:-work/coverage/$CONTIG}
GBZ=${GBZ:-data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz}
GRAPH_DB=${GRAPH_DB:-work/graph.hap32.gbz.db}
FULL_READS_DB=${FULL_READS_DB:-work/reads.hap32.gaf.db}
SUB=${SUB:-work/tier2-chr20-hap32/chr20_0_chr20.gbz}
REF_SAMPLE=${REF_SAMPLE:-CHM13}
SAMPLE=${SAMPLE:-HG002}
PLOIDY=${PLOIDY:-2}
THREADS=${THREADS:-5}
LEVELS=${LEVELS:-"5 10 15 20 25"}
FULL_TAG=${FULL_TAG:-30}

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
export PATH="$(dirname "$VG"):$PATH"

# Two binaries are needed and they are not interchangeable: `gbz-base` *queries* a database (vg
# call shells out to it per read window) and the gaf sorter/builder *creates* one. Scratchpads get
# cleaned, and gaf-base had already vanished once while gbz-base survived -- so each is located
# independently rather than assuming one implies the other.
#
# The crate ships under two binary layouts depending on how it was built. `cargo install gbz-base`
# produces separate `gafsort`/`gaf2db` executables; the repo build produces one `gaf-base` with
# `sort`/`construct` subcommands. Both are supported here because this machine has one of each.
#
# **Pin the version.** The whole-genome database was built by 0.5.1 and crates.io now serves 0.6.0.
# A 0.6.0 builder could write a database that the 0.5.1 gbz-base doing the querying cannot read,
# and nothing would say so until the numbers came out wrong. Rebuild the builder with:
#     cargo install gbz-base --version 0.5.1 --root <dir>
find_tool() {   # name -> directory containing it, or empty
    local n=$1 p
    command -v "$n" >/dev/null && { dirname "$(command -v "$n")"; return; }
    p=$(find /private/tmp/claude-501 -maxdepth 8 -name "$n" -type f -perm +111 2>/dev/null | head -1)
    [ -n "$p" ] && dirname "$p"
}
d=$(find_tool gbz-base)
[ -n "$d" ] || { echo "gbz-base not found -- vg call needs it to query the database" >&2; exit 1; }
export PATH="$d:$PATH"
echo "querying with $(command -v gbz-base) ($(gbz-base --version 2>&1 | awk '{print $2}'))"

# Builder: prefer the split binaries, fall back to the subcommand form.
if d=$(find_tool gaf2db) && [ -n "$d" ]; then
    export PATH="$d:$PATH"
    GAF_SORT="gafsort"; GAF_BUILD="gaf2db"
elif d=$(find_tool gaf-base) && [ -n "$d" ]; then
    export PATH="$d:$PATH"
    GAF_SORT="gaf-base sort"; GAF_BUILD="gaf-base construct"
else
    echo "no GAF-Base builder found (gaf2db or gaf-base)." >&2
    echo "  cargo install gbz-base --version 0.5.1 --root <dir>   # 0.5.1, not latest" >&2
    exit 1
fi
echo "building with $GAF_BUILD"
step() { echo "[$(date +%H:%M:%S)] $*"; }

# --- 1. A database per subsampled level -------------------------------------
for C in $LEVELS; do
    DB="$W/reads.${C}x.gaf.db"
    GAF="$W/$CONTIG.${C}x.gaf"
    [ -s "$GAF" ] || { echo "missing $GAF -- run subsample_gaf.py first" >&2; exit 1; }
    if [ -s "$DB" ]; then
        step "$C x: database exists, skipping"
        continue
    fi
    step "$C x: sort -> construct"
    # Through a FIFO so the sorted GAF is never written to disk, as the whole-genome build does.
    FIFO="$W/sorted.${C}x.fifo"
    rm -f "$FIFO"; mkfifo "$FIFO"
    $GAF_SORT "$GAF" -o "$FIFO" -p > "$W/sort.${C}x.log" 2>&1 &
    SORT_PID=$!
    $GAF_BUILD "$FIFO" -r "$GBZ" -o "$DB" --overwrite > "$W/construct.${C}x.log" 2>&1
    RC=$?
    wait $SORT_PID; SRC=$?
    rm -f "$FIFO"
    if [ $RC -ne 0 ] || [ $SRC -ne 0 ]; then
        echo "FAILED building $DB (construct=$RC sort=$SRC)" >&2
        tail -5 "$W/sort.${C}x.log" "$W/construct.${C}x.log" >&2
        exit 1
    fi
    ls -la "$DB" | awk '{printf "    db: %.2f GB\n", $5/1073741824}'
done

# --- 2. Call at each level, plus full ----------------------------------------
call_one() {   # tag readsdb
    local TAG=$1 DB=$2
    local OUT="$W/call.${TAG}x"
    if [ -s "$OUT.vcf.gz" ]; then
        step "$TAG x: called already, skipping"
        return
    fi
    step "$TAG x: vg call"
    /usr/bin/time -l "$VG" call "$SUB" \
        -p "${REF_SAMPLE}#0#${CONTIG}" -s "$SAMPLE" -d "$PLOIDY" -t "$THREADS" --progress \
        --read-likelihood \
        --gaf-base "$DB" --gbz-base "$GRAPH_DB" \
        > "$OUT.vcf" 2> "$OUT.log"
    bgzip -f "$OUT.vcf" && tabix -f -p vcf "$OUT.vcf.gz"
    local secs rss
    secs=$(grep -E "^ *[0-9.]+ real" "$OUT.log" | awk '{print $1}')
    rss=$(grep "maximum resident set size" "$OUT.log" | awk '{printf "%.1f", $1/2^30}')
    echo "    $(bcftools index -n "$OUT.vcf.gz") records, ${secs}s, ${rss} GB peak"
}

for C in $LEVELS; do
    call_one "$C" "$W/reads.${C}x.gaf.db"
done
# DO_FULL=0 to work on the subsampled arms alone; the full arm reuses the whole-genome database
# and is the control that ties this series back to the published tier-2 number.
if [ "${DO_FULL:-1}" = "1" ]; then
    call_one "$FULL_TAG" "$FULL_READS_DB"
fi

step "TITRATE_DONE"
