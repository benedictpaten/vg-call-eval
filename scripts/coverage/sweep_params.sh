#!/usr/bin/env bash
# Sweep one parameter across the four titration corners, so a default can be chosen from
# measurement rather than from the 30x diploid case it was originally fitted on.
#
# The corners are low and high coverage at each ploidy: chr20 at 5x and 30x, chrX non-PAR at 2.5x
# and 14.6x. Those four span both axes at their extremes, which is where a parameter that is only
# right in the middle will show it. The intermediate levels exist and can be added, but they cost
# wall clock and tell you about the interior of a curve whose ends are the question.
#
# One parameter at a time, around the shipped defaults, rather than a grid. A full grid over
# linkage weight, frequency prior and depth term at four arms is several hundred calls; the
# defaults were fitted one at a time and the question here is whether each still holds at other
# coverages, which a coordinate sweep answers.
set -uo pipefail
cd "$(dirname "$0")/../.."

PARAM=${PARAM:?e.g. --linkage-weight}
VALUES=${VALUES:?e.g. "0 1 2 4 8"}
TAG=${TAG:?short name for the output directory, e.g. lw}
W=${W:-work/coverage/sweep/$TAG}

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
export PATH="$(dirname "$VG"):$PATH"
if ! command -v gbz-base >/dev/null; then
    d=$(find /private/tmp/claude-501 -maxdepth 8 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
    [ -n "$d" ] || { echo "gbz-base not found" >&2; exit 1; }
    export PATH="$(dirname "$d"):$PATH"
fi
mkdir -p "$W"

# arm : ploidy : subgraph : reads db
ARMS=${ARMS:-"chr20.5:2:work/tier2-chr20-hap32/chr20_0_chr20.gbz:work/coverage/chr20/reads.5x.gaf.db
chr20.30:2:work/tier2-chr20-hap32/chr20_0_chr20.gbz:work/reads.hap32.gaf.db
chrX.2.5:1:work/wgs/chrX/chrX.gbz:work/coverage/chrX/reads.2.5x.gaf.db
chrX.14.6:1:work/wgs/chrX/chrX.gbz:work/reads.hap32.gaf.db"}

while IFS= read -r spec; do
    [ -z "$spec" ] && continue
    IFS=: read -r ARM PLOIDY SUB DB <<< "$spec"
    CONTIG=${ARM%%.*}
    for V in $VALUES; do
        OUT="$W/${ARM}.${TAG}${V}"
        # Skip only what is newer than the binary that would produce it. A sweep is a way of
        # measuring a build, so a result from an older build is not a cached answer to this
        # question -- it is a different one. This is not hypothetical: the haploid linkage fix
        # landed mid-sweep, the resume marker happily kept the pre-fix chrX arms, and they scored
        # as though they were the fixed caller. They were identifiable only by their timestamps.
        if [ -s "$OUT.vcf.gz" ]; then
            if [ "$OUT.vcf.gz" -nt "$VG" ]; then
                echo "[$(date +%H:%M:%S)] $ARM $PARAM $V: done, skipping"
                continue
            fi
            echo "[$(date +%H:%M:%S)] $ARM $PARAM $V: older than $VG, recomputing"
            rm -f "$OUT.vcf.gz" "$OUT.vcf.gz.tbi" "$OUT.renamed.vcf.gz" "$OUT.renamed.vcf.gz.tbi"
            rm -rf "$W/aardvark.${ARM}.${V}"
        fi
        echo "[$(date +%H:%M:%S)] $ARM $PARAM $V"
        "$VG" call "$SUB" -p "CHM13#0#${CONTIG}" -s HG002 -d "$PLOIDY" -t 5 \
            --read-likelihood "$PARAM" "$V" \
            --gaf-base "$DB" --gbz-base work/graph.hap32.gbz.db \
            > "$OUT.vcf" 2> "$OUT.log" \
          && bgzip -f "$OUT.vcf" && tabix -f -p vcf "$OUT.vcf.gz" \
          && echo "    $(bcftools index -n "$OUT.vcf.gz") records" \
          || { echo "    FAILED"; tail -3 "$OUT.log"; }
    done
done <<< "$ARMS"
echo SWEEP_DONE
