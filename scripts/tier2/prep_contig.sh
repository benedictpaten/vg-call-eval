#!/usr/bin/env bash
# Build one contig's tier-2 artefacts for one graph. Generalises prep_hap32_chr20.sh,
# which was chr20- and hap32-specific.
#
#   prep_contig.sh <contig> <graph.gbz> <graph.gbz.db> <reads.gaf.db> <workdir>
#
# The reference FASTA is extracted from *this* graph rather than reused, then compared
# against any sibling work directory for the same contig. Two graphs that disagreed about
# the reference sequence would be scored against different truth coordinates, and every
# accuracy number would be quietly wrong -- so a mismatch stops the run rather than
# warning.
set -euo pipefail
cd "$(dirname "$0")/../.."

CONTIG=${1:?contig, e.g. chr6}
GBZ=${2:?graph .gbz}
GRAPH_DB=${3:?gbz-base .db}
READS_DB=${4:?gaf-base .db}
W=${5:?work directory}
REF_SAMPLE=${REF_SAMPLE:-CHM13}
REF_PATH="${REF_SAMPLE}#0#${CONTIG}"
TRUTH=data/truth

mkdir -p "$W"
for tool in vg gbz-base bcftools samtools; do
    command -v $tool >/dev/null || { echo "not on PATH: $tool"; exit 1; }
done
step() { echo "[$(date +%H:%M:%S)] $*"; }

# 1. Contig subgraph, as GBZ so the GBWT haplotypes -z needs survive the cut.
if [ ! -s "$W/${CONTIG}_0_${CONTIG}.gbz" ]; then
    step "$W: vg chunk $CONTIG"
    vg chunk -x "$GBZ" --gbz --contig "$CONTIG" -b "$W/$CONTIG" -t 6
fi
SUB=$(ls "$W/${CONTIG}_0_"*.gbz | head -1)
echo "  subgraph: $SUB ($(du -h "$SUB" | cut -f1))"

# 2. Node ID list. gbz-base errors on a node that is not in the graph and this ID space
#    is sparse, so the read extraction needs the list explicitly.
if [ ! -s "$W/${CONTIG}_all_nodes.txt" ]; then
    step "$W: node ID list"
    vg convert -f "$SUB" | awk '$1=="S"{print $2}' | sort -n > "$W/${CONTIG}_all_nodes.txt"
fi
echo "  nodes: $(wc -l < "$W/${CONTIG}_all_nodes.txt")"

# 3. Reference FASTA from this graph.
if [ ! -s "$W/$CONTIG.fa" ]; then
    step "$W: reference FASTA"
    vg paths -x "$SUB" -F -Q "$REF_PATH" | sed "1s/.*/>$CONTIG/" > "$W/$CONTIG.fa"
    samtools faidx "$W/$CONTIG.fa"
fi

# 4. Truth slices. Reference-based, so identical for every graph -- but sliced per
#    work directory so a run is self-contained.
for kind in smvar stvar; do
    if [ ! -s "$W/truth.$CONTIG.$kind.vcf.gz" ]; then
        step "$W: truth $kind"
        bcftools view -r "$CONTIG" -Oz \
            -o "$W/truth.$CONTIG.$kind.vcf.gz" \
            "$TRUTH/CHM13v2.0_HG2-T2TQ100-V1.1_$kind.vcf.gz"
        bcftools index -t -f "$W/truth.$CONTIG.$kind.vcf.gz"
        awk -v c="$CONTIG" '$1==c' "$TRUTH/CHM13v2.0_HG2-T2TQ100-V1.1_$kind.benchmark.bed" \
            > "$W/truth.$CONTIG.$kind.bed"
    fi
done
echo "  truth: $(bcftools view -H "$W/truth.$CONTIG.smvar.vcf.gz" | wc -l) smvar records"

# 5. Reads for this contig, out of the whole-genome GAF-Base.
if [ ! -s "$W/$CONTIG.reads.gaf" ]; then
    step "$W: extract reads"
    python3 scripts/tier2/extract_reads_from_db.py \
        --nodes "$W/${CONTIG}_all_nodes.txt" \
        --gaf-base "$READS_DB" --gbz-base "$GRAPH_DB" \
        --out "$W/$CONTIG.reads.gaf" \
        --tmp "/tmp/gafbase_extract_$(basename "$W").gaf"
fi

# 6. Pack, needed by four of the five arms for allele enumeration.
if [ ! -s "$W/$CONTIG.pack" ]; then
    step "$W: vg pack"
    vg pack -x "$SUB" -a "$W/$CONTIG.reads.gaf" -o "$W/$CONTIG.pack" -t 6
fi

# 7. Cross-check the reference against any sibling directory for the same contig.
for other in work/tier2-$CONTIG work/tier2-$CONTIG-hap32; do
    [ "$other" = "$W" ] && continue
    [ -s "$other/$CONTIG.fa" ] || continue
    if cmp -s "$W/$CONTIG.fa" "$other/$CONTIG.fa"; then
        echo "  reference FASTA identical to $other -- comparisons are like for like"
    else
        echo "  reference FASTA DIFFERS from $other. Stopping: the two runs would be"
        echo "  scored against different sequence and every accuracy number would be suspect."
        exit 1
    fi
done
step "$W: DONE"
