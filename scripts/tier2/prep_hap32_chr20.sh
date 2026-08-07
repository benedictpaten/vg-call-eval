#!/usr/bin/env bash
# Build the chr20 artefacts for the 32-haplotype graph, mirroring plan §9.7 steps 1-6.
#
# Only the graph and its alignments change: the truth slices and the confident-region
# BEDs are reference-based and are reused as they are. The reference FASTA is *checked*
# rather than reused -- it is extracted from this graph and diffed against the one taken
# from the 4-haplotype graph, because if the two CHM13 paths ever disagreed, aardvark
# would be comparing against a sequence the caller never saw and every number would be
# quietly wrong.
#
# Run after work/build_hap32_dbs.sh finishes. Steps 1-3 need only the GBZ, but they peak
# around 10 GB, so they are kept out of the database build's way rather than overlapped.
set -euo pipefail
cd "$(dirname "$0")/../.."

GBZ=data/hprc-v2.1-mc-chm13-eval.HG002.hap32.gbz
GRAPH_DB=work/graph.hap32.gbz.db
READS_DB=work/reads.hap32.gaf.db
OLD=work/tier2-chr20
NEW=work/tier2-chr20-hap32
REF_PATH="CHM13#0#chr20"

mkdir -p "$NEW"
command -v vg >/dev/null || { echo "vg not on PATH"; exit 1; }
command -v gbz-base >/dev/null || { echo "gbz-base not on PATH"; exit 1; }

step() { echo "[$(date +%H:%M:%S)] $*"; }

# 1. chr20 subgraph, as GBZ so the GBWT haplotypes -z needs survive the cut.
if [ ! -s "$NEW/chr20_0_chr20.gbz" ]; then
    step "vg chunk: chr20 subgraph"
    vg chunk -x "$GBZ" --gbz --contig chr20 -b "$NEW/chr20" -t 6
fi
ls -la "$NEW"/chr20*.gbz

# 2. Node ID list, for the read extraction below. gbz-base errors on a node that is not
#    in the graph, and this ID space is sparse, so the list has to be explicit.
if [ ! -s "$NEW/chr20_all_nodes.txt" ]; then
    step "node ID list"
    vg convert -f "$NEW/chr20_0_chr20.gbz" | awk '$1=="S"{print $2}' | sort -n \
        > "$NEW/chr20_all_nodes.txt"
fi
echo "nodes: $(wc -l < "$NEW/chr20_all_nodes.txt")"

# 3. Reference FASTA from *this* graph, then checked against the 4-haplotype one.
if [ ! -s "$NEW/chr20.fa" ]; then
    step "reference FASTA"
    vg paths -x "$NEW/chr20_0_chr20.gbz" -F -Q "$REF_PATH" \
        | sed "1s/.*/>chr20/" > "$NEW/chr20.fa"
    samtools faidx "$NEW/chr20.fa"
fi
if cmp -s "$NEW/chr20.fa" "$OLD/chr20.fa"; then
    echo "reference FASTA: identical to the 4-haplotype graph's -- comparisons are like for like"
else
    echo "reference FASTA: DIFFERS from the 4-haplotype graph's. Stopping: the two runs would"
    echo "  be scored against different sequence, and every accuracy number would be suspect."
    exit 1
fi

# 4. chr20 reads out of the whole-genome GAF-Base.
if [ ! -s "$NEW/chr20.reads.gaf" ]; then
    step "extract chr20 reads from the GAF-Base"
    python3 scripts/tier2/extract_reads_from_db.py \
        --nodes "$NEW/chr20_all_nodes.txt" \
        --gaf-base "$READS_DB" --gbz-base "$GRAPH_DB" \
        --out "$NEW/chr20.reads.gaf" \
        --tmp /tmp/gafbase_extract_hap32.gaf
fi

# 5. Pack, needed by four of the five arms for allele enumeration.
if [ ! -s "$NEW/chr20.pack" ]; then
    step "vg pack"
    vg pack -x "$NEW/chr20_0_chr20.gbz" -a "$NEW/chr20.reads.gaf" -o "$NEW/chr20.pack" -t 6
fi
ls -la "$NEW/chr20.pack"

# 6. Truth slices are reference-based, so the 4-haplotype run's are correct here too.
for f in truth.chr20.smvar.vcf.gz truth.chr20.smvar.vcf.gz.tbi truth.chr20.smvar.bed \
         truth.chr20.stvar.vcf.gz truth.chr20.stvar.vcf.gz.tbi truth.chr20.stvar.bed; do
    [ -e "$NEW/$f" ] || ln -s "../tier2-chr20/$f" "$NEW/$f"
done

step "DONE"
