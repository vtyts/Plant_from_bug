#!/usr/bin/env bash
# Main pipeline to extract plant barcode sequences from insect Illumina data.
set -euo pipefail

PLANT_BARCODE_DIR=${1:-"plant_barcodes"}
INSECT_FASTQ_DIR=${2:-"plant_genes_Nov25/data"}
OUTPUT_DIR=${3:-"results"}

THREADS=${THREADS:-32}
EVALUE=${EVALUE:-1e-3}
FASTQ_SUFFIX=${FASTQ_SUFFIX:-"_R1_R2.fastq.gz"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

required_bins=(python3 blastn makeblastdb)
for bin in "${required_bins[@]}"; do
    if ! command -v "$bin" >/dev/null 2>&1; then
        echo "Missing dependency: $bin"
        exit 1
    fi
done

if [[ ! -d "$PLANT_BARCODE_DIR" ]]; then
    echo "Plant barcode directory not found: $PLANT_BARCODE_DIR"
    exit 1
fi

if [[ ! -d "$INSECT_FASTQ_DIR" ]]; then
    echo "Insect FASTQ directory not found: $INSECT_FASTQ_DIR"
    exit 1
fi

FASTAS_DIR="$OUTPUT_DIR/fastas"
BLASTDB_DIR="$OUTPUT_DIR/blastdbs"
BLAST_DIR="$OUTPUT_DIR/blast"
UNIQUE_DIR="$OUTPUT_DIR/unique"

mkdir -p "$FASTAS_DIR" "$BLASTDB_DIR" "$BLAST_DIR/matK" "$BLAST_DIR/rbcL" "$UNIQUE_DIR"

echo "==> Converting FASTQ libraries to FASTA"
python3 "$SCRIPT_DIR/prepare_fastas.py" \
    --input "$INSECT_FASTQ_DIR" \
    --output "$FASTAS_DIR" \
    --suffix "$FASTQ_SUFFIX"

combine_barcodes() {
    local gene=$1
    local pattern=$2
    local target="$OUTPUT_DIR/${gene}_barcodes.fasta"

    shopt -s nullglob
    local files=("$PLANT_BARCODE_DIR"/$pattern)
    shopt -u nullglob

    if [[ ${#files[@]} -eq 0 ]]; then
        echo "No barcode files matching '$pattern' in $PLANT_BARCODE_DIR"
        exit 1
    fi

    cat "${files[@]}" >"$target"
    echo "==> Combined ${gene} barcodes -> $target"
}

combine_barcodes "matK" "matK_*.fasta"
combine_barcodes "rbcL" "rbcL*.fasta"

run_blast() {
    local gene=$1
    local query="$OUTPUT_DIR/${gene}_barcodes.fasta"
    local outdir="$BLAST_DIR/$gene"
    shopt -s nullglob
    local fasta_files=("$FASTAS_DIR"/*.fasta)
    shopt -u nullglob
    if [[ ${#fasta_files[@]} -eq 0 ]]; then
        echo "No FASTA files produced in $FASTAS_DIR"
        exit 1
    fi

    for fasta in "${fasta_files[@]}"; do
        local sample
        sample=$(basename "$fasta" .fasta)
        local db_prefix="$BLASTDB_DIR/${sample}"
        makeblastdb -in "$fasta" -dbtype nucl -out "$db_prefix" >/dev/null
        blastn \
            -query "$query" \
            -db "$db_prefix" \
            -evalue "$EVALUE" \
            -num_threads "$THREADS" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
            -out "$outdir/${sample}.tsv"
        echo "==> ${gene} vs ${sample} complete"
    done

    cat "$outdir"/*.tsv >"$BLAST_DIR/${gene}_all.tsv"
    echo "==> Aggregated ${gene} hits -> $BLAST_DIR/${gene}_all.tsv"
}

run_blast "matK"
run_blast "rbcL"

echo "==> Deriving unique hits"
python3 "$SCRIPT_DIR/collect_unique_hits.py" \
    --blast-tsv "$BLAST_DIR/matK_all.tsv" \
    --unique-fasta "$UNIQUE_DIR/matK_unique_hits.fasta" \
    --unique-table "$UNIQUE_DIR/matK_unique_hits.tsv"

python3 "$SCRIPT_DIR/collect_unique_hits.py" \
    --blast-tsv "$BLAST_DIR/rbcL_all.tsv" \
    --unique-fasta "$UNIQUE_DIR/rbcL_unique_hits.fasta" \
    --unique-table "$UNIQUE_DIR/rbcL_unique_hits.tsv"

cat <<EOF
Pipeline complete.
- Unique matK hits: $UNIQUE_DIR/matK_unique_hits.fasta
- Unique rbcL hits: $UNIQUE_DIR/rbcL_unique_hits.fasta

To identify plant taxa via GenBank nt, run for each file:
  export NCBI_EMAIL="you@example.com"
  ./scripts/blast_nt_hits.sh results/unique/matK_unique_hits.fasta results/nt/matK_vs_nt
EOF
