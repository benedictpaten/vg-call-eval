#!/usr/bin/env bash
# Run the five arms on chr20 of the 32-haplotype graph, then both benchmarks.
#
# Everything is the same as the 4-haplotype run except the graph and the alignments,
# which is the point: the truth slices, the confident-region BEDs and the reference
# FASTA are shared (and prep_hap32_chr20.sh checks the FASTA really is identical).
#
# The two runs are not a clean single-variable comparison and should not be described
# as one. Reads mapped to the 4-haplotype graph cannot be used with this one -- node
# IDs differ -- so the 32-haplotype arm uses its own alignments. Graph and alignment
# move together. That is what anyone adopting a richer graph would actually do, but it
# means "the graph is worth X" is not a claim this design can support on its own.
set -euo pipefail
cd "$(dirname "$0")/../.."

NEW=work/tier2-chr20-hap32
OUT=$NEW/results

for f in "$NEW/chr20_0_chr20.gbz" "$NEW/chr20.pack" "$NEW/chr20.fa" \
         work/reads.hap32.gaf.db work/graph.hap32.gbz.db; do
    [ -s "$f" ] || { echo "missing: $f -- run prep_hap32_chr20.sh first"; exit 1; }
done
mkdir -p "$OUT"

echo "=== small variants: five arms ==="
python3 scripts/tier2/run_arms.py \
    --graph "$NEW/chr20_0_chr20.gbz" --pack "$NEW/chr20.pack" \
    --gaf-base work/reads.hap32.gaf.db --gbz-base work/graph.hap32.gbz.db \
    --reference "$NEW/chr20.fa" \
    --truth-vcf "$NEW/truth.chr20.smvar.vcf.gz" --truth-bed "$NEW/truth.chr20.smvar.bed" \
    --out "$OUT"

echo
echo "=== structural variants ==="
python3 scripts/tier2/compare_sv.py \
    --results "$OUT" \
    --truth-vcf "$NEW/truth.chr20.stvar.vcf.gz" --truth-bed "$NEW/truth.chr20.stvar.bed" \
    --reference "$NEW/chr20.fa"

echo
echo "=== size-matched insertion check ==="
python3 scripts/tier2/size_matched.py \
    --results "$OUT" \
    --truth-vcf "$NEW/truth.chr20.smvar.vcf.gz" --truth-bed "$NEW/truth.chr20.smvar.bed" \
    --reference "$NEW/chr20.fa"

echo "DONE"
