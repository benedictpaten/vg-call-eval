#!/usr/bin/env bash
# H3 (plan §9.19): does the mismapping floor want to be higher on a richer graph?
#
# §9.15 set --mismap-min to 0.01 by measuring on the 4-haplotype graph. If the floor
# is really standing in for "a read fits an allele it did not come from", then a graph
# with 34 haplotypes -- where the mismapping term already suppresses ~10x more calls --
# should want a larger value. If the optimum moves, the floor should scale with graph
# complexity rather than be a constant, and §9.15's single number is incomplete.
#
# The 4-haplotype sweep is re-run at the same points, because the two must be compared
# on the same vg build: the read path changed substantially in 44fd008, and while that
# was verified byte-identical, "verified" is worth less than "measured side by side".
set -o pipefail
cd "$(dirname "$0")/../.."
V=/Users/benedictpaten/CLionProjects/vg/bin/vg
export PATH="/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/gbztools/bin:$PATH"

run_one() {  # <tag> <workdir> <gaf-base> <gbz-base> <floor>
    local tag=$1 W=$2 GAFDB=$3 GBZDB=$4 E=$5
    local out=$W/results
    # Cache the *call* only. Skipping the comparison too meant one bad aardvark flag
    # silently produced a sweep with no metrics in it, and a re-run repaired nothing.
    if [ -s "$out/$tag.vcf.gz" ]; then
        echo "[$(date +%H:%M:%S)] $tag: call cached"
    else
    echo "[$(date +%H:%M:%S)] $tag: --mismap-min $E"
    $V call "$W/chr20_0_chr20.gbz" -p "CHM13#0#chr20" -z --read-likelihood \
       --gaf-base "$GAFDB" --gbz-base "$GBZDB" \
       --mismap-min "$E" -t 6 > "$out/$tag.vcf" 2> "$out/$tag.log" || return 1
    # The harness renames the reference contig so the VCF matches the truth's naming.
    python3 - "$out/$tag.vcf" <<'PY'
import sys
p=sys.argv[1]
src=open(p).read().replace("ID=CHM13#0#chr20","ID=chr20")
out=[]
for line in src.split("\n"):
    if line.startswith("#") or not line: out.append(line); continue
    f=line.split("\t")
    if f[0]=="CHM13#0#chr20": f[0]="chr20"
    out.append("\t".join(f))
open(p,"w").write("\n".join(out))
PY
    bgzip -f -c "$out/$tag.vcf" > "$out/$tag.vcf.gz" && tabix -f -p vcf "$out/$tag.vcf.gz"
    rm -f "$out/$tag.vcf"
    fi
    aardvark compare -r "$W/chr20.fa" -t "$W/truth.chr20.smvar.vcf.gz" \
       -q "$out/$tag.vcf.gz" -b "$W/truth.chr20.smvar.bed" -o "$out/aardvark-$tag" \
       --truth-sample HG002 --query-sample SAMPLE --compare-label "$tag" \
       --min-variant-gap 50 --max-branch-factor 50 --enable-record-basepair-metrics \
       --threads 6 > "$out/aardvark-$tag.log" 2>&1
    echo "[$(date +%H:%M:%S)] $tag done"
}

for E in 0.005 0.01 0.02 0.05 0.10 0.20; do
    tag="sweep-$E"
    run_one "$tag" work/tier2-chr20-hap32 work/reads.hap32.gaf.db work/graph.hap32.gbz.db "$E"
done
echo "MARKER_HAP32_SWEEP_DONE"

for E in 0.005 0.01 0.02 0.05 0.10 0.20; do
    tag="sweep-$E"
    run_one "$tag" work/tier2-chr20 work/reads.gaf.db work/graph.gbz.db "$E"
done
echo "MARKER_SWEEP_ALL_DONE"
