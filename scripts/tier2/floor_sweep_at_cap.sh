#!/usr/bin/env bash
# Re-sweep --mismap-min at the *new* --mismap-max default of 0.5, on both graphs.
#
# The earlier floor sweep was run with the cap at its old 0.1, which is no longer the
# operating point, and the two knobs are not independent: at cap 0.5 floor 0.10 was
# already worse than 0.05, and floor 0.20 was only expressible at all once the cap rose
# (vg enforces min <= max). Carrying "0.05 is optimal" across a cap change would be
# assuming exactly the kind of transfer that made 0.1 the wrong cap in the first place.
#
# Points already computed during the cap investigation are cached and skipped.
set -o pipefail
cd "$(dirname "$0")/../.."
V=/Users/benedictpaten/CLionProjects/vg/bin/vg
export PATH="/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/gbztools/bin:$PATH"

run_one() {  # <workdir> <gaf-base> <gbz-base> <floor>
    local W=$1 GAFDB=$2 GBZDB=$3 FLOOR=$4
    local out=$W/results tag="cap0.5-fl$FLOOR"
    if [ -s "$out/$tag.vcf.gz" ]; then
        echo "[$(date +%H:%M:%S)] $W $tag: call cached"
    else
        echo "[$(date +%H:%M:%S)] $W $tag"
        $V call "$W/chr20_0_chr20.gbz" -p "CHM13#0#chr20" -z --read-likelihood \
           --gaf-base "$GAFDB" --gbz-base "$GBZDB" \
           --mismap-max 0.5 --mismap-min "$FLOOR" -t 6 \
           > "$out/$tag.vcf" 2> "$out/$tag.log" || { echo "  FAILED: $(tail -1 "$out/$tag.log")"; return 1; }
        python3 - "$out/$tag.vcf" <<'PY'
import sys
p=sys.argv[1]
src=open(p).read().replace("ID=CHM13#0#chr20","ID=chr20")
o=[]
for line in src.split("\n"):
    if line.startswith("#") or not line: o.append(line); continue
    f=line.split("\t")
    if f[0]=="CHM13#0#chr20": f[0]="chr20"
    o.append("\t".join(f))
open(p,"w").write("\n".join(o))
PY
        bgzip -f -c "$out/$tag.vcf" > "$out/$tag.vcf.gz" && tabix -f -p vcf "$out/$tag.vcf.gz"
        rm -f "$out/$tag.vcf"
    fi
    [ -f "$out/aardvark-$tag/summary.tsv" ] && { echo "  compare cached"; return 0; }
    aardvark compare -r "$W/chr20.fa" -t "$W/truth.chr20.smvar.vcf.gz" \
       -q "$out/$tag.vcf.gz" -b "$W/truth.chr20.smvar.bed" -o "$out/aardvark-$tag" \
       --truth-sample HG002 --query-sample SAMPLE --compare-label "$tag" \
       --min-variant-gap 50 --max-branch-factor 50 --enable-record-basepair-metrics \
       --threads 6 > "$out/aardvark-$tag.log" 2>&1
    echo "[$(date +%H:%M:%S)] $W $tag done"
}

for FLOOR in 0.01 0.02 0.03 0.05 0.10 0.20; do
    run_one work/tier2-chr20-hap32 work/reads.hap32.gaf.db work/graph.hap32.gbz.db "$FLOOR"
done
for FLOOR in 0.01 0.02 0.03 0.05 0.10 0.20; do
    run_one work/tier2-chr20 work/reads.gaf.db work/graph.gbz.db "$FLOOR"
done
echo MARKER_FLOOR_AT_CAP_DONE
