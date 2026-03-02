#!/usr/bin/env bash
set -euo pipefail
module load ncbi-blast/2.17.0+

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <unique_hits_fasta> <output_prefix>"
    exit 1
fi

UNIQUE_FASTA=$1
OUT_PREFIX=$2
OUT_FILE="${OUT_PREFIX}.tsv"

EVALUE=${EVALUE:-1e-3}
MAX_TARGET_SEQS=${MAX_TARGET_SEQS:-25}
TAX_FILTER=${TAX_FILTER:-"Viridiplantae[ORGN]"}
MAX_RETRIES=${MAX_RETRIES:-2}
RETRY_DELAY=${RETRY_DELAY:-300}

if [[ ! -f "$UNIQUE_FASTA" ]]; then
    echo "Unique hits FASTA not found: $UNIQUE_FASTA"
    exit 1
fi

if [[ -z "${NCBI_EMAIL:-}" ]]; then
    echo "Please set NCBI_EMAIL before running remote BLAST."
    exit 1
fi

API_FLAG=()
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    API_FLAG=(-api_key "$NCBI_API_KEY")
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

ATTEMPT=1
SUCCESS=0

while [ $ATTEMPT -le $MAX_RETRIES ]

do

    echo "==> BLAST attempt $ATTEMPT for $UNIQUE_FASTA"

    rm -f "$OUT_FILE"

    set +e
    blastn \
        -query "$UNIQUE_FASTA" \
        -db nt \
        -remote \
        -out "$OUT_FILE" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms stitle" \
        -max_target_seqs "$MAX_TARGET_SEQS" \
        -evalue "$EVALUE" \
        -entrez_query "$TAX_FILTER" \
        "${API_FLAG[@]}"
    status=$?
    set -e 

    if [[ $status -eq 0 && -s "$OUT_FILE" ]]; then
        SUCCESS=1
        break
    fi

    echo "[warn] BLAST failed (exit code: $status) or output empty."

    if [ $ATTEMPT -le $MAX_RETRIES ]; then
        echo "Waiting ${RETRY_DELAY}s before retry..."
        sleep $RETRY_DELAY
    fi
    ATTEMPT=$(($ATTEMPT + 1))
    
done

if [ $SUCCESS -eq 0 ] ; then
    echo "[error] BLAST failed after ${MAX_RETRIES} attempts: $UNIQUE_FASTA"
    exit 0
fi

echo "[ok] Remote BLAST complete: $OUT_FILE"
