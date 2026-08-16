#!/usr/bin/env bash
# Lay an external call set out so bench_wgs.py scores it by exactly the same code path as ours.
#
# The point of the exercise is a comparison, so the scoring must not differ in any respect: same
# truth VCFs, same confident-region BEDs, same reference FASTA, same aardvark and truvari
# invocations, same per-contig-then-sum aggregation, same chrY exclusion. Rather than write a
# second scorer and hope it matches, this builds a work directory of the shape bench_wgs.py
# already expects -- per-contig calls beside symlinks to the truth files the vg run used -- so the
# existing script runs over it unchanged.
#
# Two adjustments the input needs, both of which would otherwise make the comparison wrong:
#
#   * **Hom-ref records are dropped.** PanGenie genotypes a fixed panel and emits a record for
#     every site in it, including 0/0. vg call emits only non-reference calls (absent -a). Leaving
#     them in would hand the comparison a large population of records that are not calls at all.
#   * **Contig headers are added.** The file has none, so it cannot be indexed or sliced.
#
# What is deliberately *not* adjusted: the allele representation. The input is biallelic-split
# where ours is multiallelic. aardvark compares by local haplotype rather than by record, so this
# is the difference it is built to absorb; normalising it by hand would be the more dangerous move.
set -euo pipefail
cd "$(dirname "$0")/../.."

SRC=${SRC:?path to the external VCF}
W=${W:-work/pangenie}
REF_WORK=${REF_WORK:-work/wgs}
SAMPLE=${SAMPLE:-HG002}
CONTIGS=${CONTIGS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"}

step() { echo "[$(date +%H:%M:%S)] $*"; }
mkdir -p "$W"

# 1. Contig header lines, taken from the same FASTAs the truth is scored against, so the
#    coordinate system is asserted rather than assumed.
if [ ! -s "$W/contigs.txt" ]; then
    step "building contig header from $REF_WORK/*/*.fa.fai"
    : > "$W/contigs.txt"
    for C in $CONTIGS; do
        FAI="$REF_WORK/$C/$C.fa.fai"
        [ -s "$FAI" ] || { echo "missing $FAI" >&2; exit 1; }
        awk -v c="$C" '$1==c {printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "$FAI" >> "$W/contigs.txt"
    done
fi
echo "  $(wc -l < "$W/contigs.txt") contig lines"

# 2. Non-reference records only, with the header repaired, sorted and indexed.
if [ ! -s "$W/calls.vcf.gz" ]; then
    step "filtering to non-reference calls and indexing"
    # grep, not `head -n -1`: BSD head on macOS rejects a negative count, and the failure is a
    # usage error rather than anything that looks like a VCF problem.
    bcftools view -h "$SRC" 2>/dev/null | grep -v "^#CHROM" > "$W/hdr.txt"
    cat "$W/contigs.txt" >> "$W/hdr.txt"
    bcftools view -h "$SRC" 2>/dev/null | grep "^#CHROM" >> "$W/hdr.txt"
    # GT="alt" keeps anything with a non-reference allele and drops 0/0 and ./. .
    {
        cat "$W/hdr.txt"
        bcftools view -H -i 'GT="alt"' "$SRC" 2>/dev/null
    } | bcftools sort -Oz -o "$W/calls.vcf.gz" -
    bcftools index -f -t "$W/calls.vcf.gz"
fi
echo "  kept $(bcftools index -n "$W/calls.vcf.gz") non-reference records"

# 3. Per-contig layout: the calls sliced, the truth and reference symlinked from the vg run so
#    both call sets are scored against byte-identical inputs.
for C in $CONTIGS; do
    mkdir -p "$W/$C"
    if [ ! -s "$W/$C/$C.vcf.gz" ]; then
        bcftools view -r "$C" -Oz -o "$W/$C/$C.vcf.gz" "$W/calls.vcf.gz"
        bcftools index -f -t "$W/$C/$C.vcf.gz"
    fi
    for f in "truth.$C.smvar.vcf.gz" "truth.$C.smvar.vcf.gz.tbi" "truth.$C.smvar.bed" \
             "truth.$C.stvar.vcf.gz" "truth.$C.stvar.vcf.gz.tbi" "truth.$C.stvar.bed" \
             "$C.fa" "$C.fa.fai"; do
        [ -e "$REF_WORK/$C/$f" ] && ln -sf "$PWD/$REF_WORK/$C/$f" "$W/$C/$f"
    done
    printf "  %-6s %8s records\n" "$C" "$(bcftools index -n "$W/$C/$C.vcf.gz")"
done
step "PREP_EXTERNAL_DONE"
