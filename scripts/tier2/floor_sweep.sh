#!/usr/bin/env bash
set -o pipefail
cd "$(dirname "$0")"
V=/Users/benedictpaten/CLionProjects/vg/bin/vg
export PATH="/private/tmp/claude-501/-Users-benedictpaten-My-Drive-papers-and-projects-2026-vg-planning/65260683-9f63-4117-b57b-52dccfd3b871/scratchpad/gbztools/bin:$PATH"
for E in 0.01 0.05; do
  tag="fl$E"
  echo "[$(date +%H:%M:%S)] calling --mismap-min $E"
  $V call chr20_0_chr20.gbz -p CHM13#0#chr20 -z --read-likelihood \
     --gaf-base ../reads.gaf.db --gbz-base ../graph.gbz.db --read-window 1024 \
     --mismap-min $E -t 6 > results/$tag.vcf 2> results/$tag.log
  bgzip -f -c results/$tag.vcf > results/$tag.vcf.gz && tabix -f -p vcf results/$tag.vcf.gz
  aardvark compare -r chr20.fa -t truth.chr20.smvar.vcf.gz -q results/$tag.vcf.gz \
     -b truth.chr20.smvar.bed -o results/aardvark-$tag --truth-sample HG002 \
     --query-sample SAMPLE --compare-label $tag --min-variant-gap 50 \
     --max-branch-factor 50 --threads 6 > results/aardvark-$tag.log 2>&1
  echo "[$(date +%H:%M:%S)] $tag done"
done
echo ALLDONE
