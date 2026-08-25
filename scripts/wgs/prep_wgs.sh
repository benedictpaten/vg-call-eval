#!/usr/bin/env bash
# Whole-genome prep: one GBZ subgraph per contig, plus the reference FASTA and truth slices.
#
# Deliberately *not* prep_contig.sh. That script also extracts a per-contig GAF and builds a pack,
# which cost 12 GB and 3 GB for chr6 alone -- 270 GB extrapolated to a genome, against 300 GB free.
# Both exist only for support enumeration, and panel enumeration is now the default under
# --read-likelihood, so neither is needed here. What remains is ~200 MB a contig.
#
# The reads stay where they are: the whole-genome GAF-Base and GBZ-Base databases are queried
# directly by node ID, and `vg chunk` preserves whole-genome node IDs, so a subgraph's IDs address
# the same reads the full graph would. That is also what keeps the emitted mosaic readable against
# the whole-genome GBZ rather than only against the chunk it came from.
set -euo pipefail
cd "$(dirname "$0")/../.."

VG=${VG:-$HOME/CLionProjects/vg/bin/vg}
export PATH="$(dirname "$VG"):$PATH"
GBZ=${GBZ:-data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz}
W=${W:-work/wgs}
REF_SAMPLE=${REF_SAMPLE:-CHM13}
TRUTH=data/truth
THREADS=${THREADS:-6}
CONTIGS=${CONTIGS:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"}

mkdir -p "$W"
for tool in vg bcftools samtools; do
    command -v $tool >/dev/null || { echo "not on PATH: $tool"; exit 1; }
done
step() { echo "[$(date +%H:%M:%S)] $*"; }

for CONTIG in $CONTIGS; do
    D="$W/$CONTIG"
    mkdir -p "$D"
    SUB="$D/${CONTIG}.gbz"

    if [ ! -s "$SUB" ]; then
        step "$CONTIG: chunk"
        vg chunk -x "$GBZ" --gbz --contig "$CONTIG" -b "$D/chunk" -t "$THREADS"
        mv "$(ls "$D/chunk_0_"*.gbz | head -1)" "$SUB"
        rm -f "$D/chunk_"*.gbz
    fi

    # Snarls, cached. `vg call` decomposes the contig itself when it is not given any, and that
    # decomposition is single-threaded: 46 s of a 197 s chr20 run, with four of five threads
    # parked. Doing it once here and passing -r takes chr20 to 148.5 s.
    #
    # The two flags are both load-bearing, and the output is byte-identical only with both:
    #
    #   -T  include trivial snarls. `vg call`'s own SnarlManager has them and the symbolic
    #       projection keys chain symbols on child chain boundaries, so omitting them changes the
    #       nested structure. Measured: without -T, chr20 gained 236 records and 714 lines moved.
    #   -P  upweight the reference path's tips. snarls_main.cpp and call_main.cpp apply the same
    #       EXTRA_WEIGHT to the same first and last node of each reference path, so this
    #       reproduces the caller's decomposition rather than merely a decomposition.
    #
    # Regenerated when the graph is newer, on the same reasoning as call_wgs.sh's .done marker: a
    # stale snarl file would silently describe a different graph.
    if [ ! -s "$D/$CONTIG.snarls.pb" ] || [ "$SUB" -nt "$D/$CONTIG.snarls.pb" ]; then
        step "$CONTIG: snarls"
        vg snarls -T -P "${REF_SAMPLE}#0#${CONTIG}" -t "$THREADS" "$SUB" > "$D/$CONTIG.snarls.pb"
    fi

    if [ ! -s "$D/$CONTIG.fa.fai" ]; then
        step "$CONTIG: reference FASTA"
        vg paths -x "$SUB" -F -Q "${REF_SAMPLE}#0#${CONTIG}" \
            | sed "1s/.*/>$CONTIG/" > "$D/$CONTIG.fa"
        samtools faidx "$D/$CONTIG.fa"
    fi

    for kind in smvar stvar; do
        if [ ! -s "$D/truth.$CONTIG.$kind.vcf.gz.tbi" ]; then
            step "$CONTIG: truth $kind"
            bcftools view -r "$CONTIG" -Oz \
                -o "$D/truth.$CONTIG.$kind.vcf.gz" \
                "$TRUTH/CHM13v2.0_HG2-T2TQ100-V1.1_${kind}.vcf.gz"
            bcftools index -t -f "$D/truth.$CONTIG.$kind.vcf.gz"
            awk -v c="$CONTIG" '$1==c' \
                "$TRUTH/CHM13v2.0_HG2-T2TQ100-V1.1_${kind}.benchmark.bed" \
                > "$D/truth.$CONTIG.$kind.bed"
        fi
    done
    echo "  $CONTIG: $(du -h "$SUB" | cut -f1) subgraph, $(du -h "$D/$CONTIG.snarls.pb" | cut -f1) snarls, $(bcftools view -H "$D/truth.$CONTIG.smvar.vcf.gz" | wc -l | tr -d ' ') smvar truth records"
done
echo "PREP_DONE  total $(du -sh "$W" | cut -f1)"
