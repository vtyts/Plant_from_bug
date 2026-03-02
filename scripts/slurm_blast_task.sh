#!/usr/bin/env bash
set -euo pipefail
module load ncbi-blast/2.17.0+

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
SKIP_DB=$6

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

echo "[start] ${GENE} vs ${SAMPLE}"

mkdir -p "$BLAST_DIR"
DB_PREFIX="$BLASTDB_DIR/${SAMPLE}"

if (( $SKIP_DB == 0 )); then
    mkdir -p "$BLASTDB_DIR"
    makeblastdb -in "$FASTA_PATH" -dbtype nucl -out "$DB_PREFIX" #>/dev/null
fi

sleep 10

echo "[db] makeblastdb finished"

BLAST_THREADS=${BLAST_THREADS:-3}
EVALUE=${EVALUE:-1e-3}
MAX_RETRIES=${MAX_RETRIES:-2}
RETRY_DELAY=${RETRY_DELAY:-300}

OUT_FILE="$BLAST_DIR/${SAMPLE}.tsv"
FAIL_FILE="${BLAST_DIR}/FAILED_BLASTS_${SAMPLE}.tab"

ATTEMPT=1
SUCCESS=0
 
while [ $ATTEMPT -le $MAX_RETRIES ]

do

    echo "[blastn] BLAST attempt ${ATTEMPT} for ${GENE} vs ${SAMPLE}"

    rm -f "${OUT_FILE}"

    set +e
    blastn \
        -query "$QUERY_FASTA" \
        -db "$DB_PREFIX" \
        -evalue "$EVALUE" \
        -num_threads "$BLAST_THREADS" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \
        -out "$OUT_FILE"
    status=$?
    set -e

    if [[ $status -eq 0 && -f "$OUT_FILE" ]]; then
        SUCCESS=1
        break
    fi

    echo "[warn] BLAST failed for ${SAMPLE} (exit code: $status)"

    if [ $ATTEMPT -le $MAX_RETRIES ]; then
        echo "Waiting ${RETRY_DELAY}s before retry..."
        sleep $RETRY_DELAY
    fi

    ATTEMPT=$(($ATTEMPT + 1))
done

if [ $SUCCESS -eq 0 ] ; then
    echo "[error] BLAST failed after ${MAX_RETRIES} attempts: ${GENE} vs ${SAMPLE}"
    echo -e "${GENE}\t${SAMPLE}" > "$FAIL_FILE"
    exit 0
fi

echo "[ok] ${GENE} vs ${SAMPLE}"
