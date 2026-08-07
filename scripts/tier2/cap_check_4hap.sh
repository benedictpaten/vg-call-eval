#!/usr/bin/env bash
# H4 follow-on (plan §9.19): sweep --mismap-max, the cap, on the 32-haplotype graph.
#
# Why the cap and not the floor. At the SNV sites that go wrong, **29.2% of reads have
# MAPQ < 10** against 1.85% at the sites that go right. Phred 10 is p=0.1, which is the
# default cap -- so at those sites almost a third of the evidence is being told it is at
# most 10% likely to be mismapped, when the mapper is saying it could easily be somewhere
# else entirely. A MAPQ 0 read gets e_r = 0.1 rather than 1.0.
#
# §9.14 swept this knob on the 4-haplotype graph and found it inert, because only 6.3% of
# reads there sat at MAPQ <= 9. That conclusion does not transfer: the richer graph moves
# far more reads into the range where the cap is the thing deciding, not the MAPQ.
#
# The floor is swept alongside because vg enforces min <= max, so raising the cap also
# unlocks floors that were previously rejected -- 0.20 failed in the floor sweep for
# exactly that reason.
set -o pipefail
cd "$(dirname "$0")/../.."
V=/Users/benedictpaten/CLionProjects/vg/bin/vg
export PATH="/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-2026-07-vg-call-refactor/e061d8cf-7996-49a4-986b-3686374ef250/scratchpad/gbztools/bin:$PATH"
W=work/tier2-chr20
out=$W/results

run_one() {  # <cap> <floor>
    local CAP=$1 FLOOR=$2
    local tag="cap$CAP-fl$FLOOR"
    if [ -s "$out/$tag.vcf.gz" ]; then
        echo "[$(date +%H:%M:%S)] $tag: call cached"
    else
        echo "[$(date +%H:%M:%S)] $tag: --mismap-max $CAP --mismap-min $FLOOR"
        $V call "$W/chr20_0_chr20.gbz" -p "CHM13#0#chr20" -z --read-likelihood \
           --gaf-base work/reads.gaf.db --gbz-base work/graph.gbz.db \
           --mismap-max "$CAP" --mismap-min "$FLOOR" -t 6 \
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
    aardvark compare -r "$W/chr20.fa" -t "$W/truth.chr20.smvar.vcf.gz" \
       -q "$out/$tag.vcf.gz" -b "$W/truth.chr20.smvar.bed" -o "$out/aardvark-$tag" \
       --truth-sample HG002 --query-sample SAMPLE --compare-label "$tag" \
       --min-variant-gap 50 --max-branch-factor 50 --enable-record-basepair-metrics \
       --threads 6 > "$out/aardvark-$tag.log" 2>&1
    echo "[$(date +%H:%M:%S)] $tag done"
}

# Cap alone, holding the floor at the value the floor sweep preferred.
for CAP in 0.5; do run_one "$CAP" 0.05; done
# Floors the old cap forbade, now that it is raised.
for CAP in 0.5; do run_one "$CAP" 0.01; done
echo "MARKER_CAP4HAP_DONE"
