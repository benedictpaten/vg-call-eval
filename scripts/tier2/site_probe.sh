#!/usr/bin/env bash
# Re-call with `-a/--genotype-snarls` so that *reference* calls are emitted too.
#
# !! THIS DOES NOT ANSWER THE QUESTION IT WAS BUILT FOR. Kept as a recorded negative
# !! result; sv_error_atlas.py no longer reads its output.
# !!
# !! `-a` changes the snarl decomposition, not merely which genotypes are printed. On
# !! chr6 4-hap, `poisson-z` calls 47 of 48 large heterozygous deletions in its normal
# !! run and its own `-a` probe carries a matching allele at only 26 of them -- the loss
# !! falls precisely on the large nested sites the probe exists to interrogate. Bulk
# !! agreement hides it completely: 286,557 non-reference records either way and 93.9%
# !! position agreement. Believing it would have turned a genotyper defect into a
# !! graph-content story.
# !!
# !! To ask "was this allele enumerated?", try `-T/--traversals`, which reports the
# !! candidate traversals without genotyping at all.
#
# Why this exists. `vg call` emits no 0/0 records at all -- every record in a normal
# arm VCF is a non-reference call. So for a truth SV the caller missed, the VCF cannot
# distinguish "the site was never offered" (a graph/enumeration limitation) from "the
# site was genotyped and reference won" (a genotyper limitation). That is exactly the
# split the SV error investigation turns on, and `-a` is the only way to get it.
#
# These runs are a *probe*, not an arm: the non-reference calls they contain should be
# identical to the corresponding arm VCF, and sv_error_atlas.py asserts that before
# using them. If it ever fails, the probe is describing a different experiment.
set -euo pipefail

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
THREADS=${THREADS:-6}

# `--gaf-base` shells out to `gbz-base`, which is not installed system-wide here. A
# non-interactive shell does not get the login PATH, so this fails 40 s in with a clear
# error -- but only for the read-likelihood arms, leaving the Poisson ones to succeed
# and the run to look half-finished rather than misconfigured.
GBZ_BASE_BIN=${GBZ_BASE_BIN:-$(command -v gbz-base || true)}
if [ -z "$GBZ_BASE_BIN" ]; then
    GBZ_BASE_BIN=$(find /private/tmp/claude-501 -name gbz-base -type f -perm +111 2>/dev/null | head -1)
fi
[ -n "$GBZ_BASE_BIN" ] || { echo "gbz-base not found; set GBZ_BASE_BIN" >&2; exit 1; }
echo "using gbz-base: $GBZ_BASE_BIN"

run() {
    local ds=$1 contig=$2 gbzdb=$3 gafdb=$4 arm=$5
    shift 5
    local w="$REPO/work/$ds"
    local out="$w/results/${arm}.allsnarls.vcf"
    if [ -s "${out}.gz" ]; then echo "skip $ds $arm (exists)"; return; fi
    echo "=== $ds $arm ==="
    /usr/bin/time -l "$VG" call "$w/${contig}_0_${contig}.gbz" \
        -p "CHM13#0#${contig}" -t "$THREADS" --progress -a "$@" \
        > "$out" 2> "$w/results/${arm}.allsnarls.log" || {
            echo "FAILED $ds $arm"; tail -5 "$w/results/${arm}.allsnarls.log"; return; }
    bgzip -f "$out"
    tabix -f -p vcf "${out}.gz"
    echo "  $(bcftools index -n "${out}.gz") records"
}

# readlik-z: haplotype enumeration, no pack. poisson-z: haplotype enumeration, needs pack.
for spec in \
    "tier2-chr20        chr20 graph.gbz.db        reads.gaf.db" \
    "tier2-chr20-hap32  chr20 graph.hap32.gbz.db  reads.hap32.gaf.db" \
    "tier2-chr6         chr6  graph.gbz.db        reads.gaf.db" \
    "tier2-chr6-hap32   chr6  graph.hap32.gbz.db  reads.hap32.gaf.db" \
; do
    set -- $spec
    ds=$1 contig=$2 gbzdb=$3 gafdb=$4
    run "$ds" "$contig" "$gbzdb" "$gafdb" readlik-z \
        --read-likelihood -z --gaf-base "$REPO/work/$gafdb" --gbz-base "$REPO/work/$gbzdb" \
        --gaf-base-binary "$GBZ_BASE_BIN"
    run "$ds" "$contig" "$gbzdb" "$gafdb" poisson-z \
        -z -k "$REPO/work/$ds/${contig}.pack"
done
echo MARKER_SITE_PROBE_DONE
