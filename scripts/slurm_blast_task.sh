#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <manifest.tsv> <gene> <query_fasta> <blast_dir> <blastdb_dir>"
    echo "Tip: run via sbatch --array ... slurm_blast_task.sh manifest gene query outdir blastdb"
    exit 1
fi

MANIFEST=$1
GENE=$2
QUERY_FASTA=$3
BLAST_DIR=$4
BLASTDB_DIR=$5

TASK_ID=${SLURM_ARRAY_TASK_ID:-${TASK_INDEX:-}}
if [[ -z "$TASK_ID" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set (did you forget --array?)."
    echo "For local testing, set TASK_INDEX to the zero-based line number."
    exit 1
fi

LINE_NUMBER=$((TASK_ID + 1))
if [[ ! -f "$MANIFEST" ]]; then
    echo "Manifest not found: $MANIFEST"
    exit 1
fi

ENTRY=$(sed -n "${LINE_NUMBER}p" "$MANIFEST" || true)
if [[ -z "$ENTRY" ]]; then
    echo "No manifest entry for task $TASK_ID (line $LINE_NUMBER); exiting."
    exit 0
fi

IFS=$'\t' read -r SAMPLE FASTA_PATH <<<"$ENTRY"
if [[ -z "$SAMPLE" || -z "$FASTA_PATH" ]]; then
    echo "Malformed manifest entry on line $LINE_NUMBER: $ENTRY"
    exit 1
fi

if [[ ! -f "$FASTA_PATH" ]]; then
    echo "FASTA not found for sample $SAMPLE: $FASTA_PATH"
    exit 1
fi

mkdir -p "$BLAST_DIR" "$BLASTDB_DIR"

DB_PREFIX="$BLASTDB_DIR/${SAMPLE}"
makeblastdb -in "$FASTA_PATH" -dbtype nucl -out "$DB_PREFIX" >/dev/null

BLAST_THREADS=${BLAST_THREADS:-${THREADS:-32}}
EVALUE=${EVALUE:-1e-3}

blastn \
    -query "$QUERY_FASTA" \
    -db "$DB_PREFIX" \
    -evalue "$EVALUE" \
    -num_threads "$BLAST_THREADS" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
    -out "$BLAST_DIR/${SAMPLE}.tsv"

echo "[ok] ${GENE} vs ${SAMPLE}"
