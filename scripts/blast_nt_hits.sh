#!/usr/bin/env bash
set -euo pipefail
module load ncbi-blast/2.16.0+

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <unique_hits_fasta> <output_prefix>"
    echo "Example: $0 results/unique/matK_unique_hits.fasta results/nt/matK_vs_nt"
    exit 1
fi

UNIQUE_FASTA=$1
OUT_PREFIX=$2
EVALUE=${EVALUE:-1e-3}
MAX_TARGET_SEQS=${MAX_TARGET_SEQS:-25}
TAX_FILTER=${TAX_FILTER:-"Viridiplantae[ORGN]"}

if [[ ! -f "$UNIQUE_FASTA" ]]; then
    echo "Unique hits FASTA not found: $UNIQUE_FASTA"
    exit 1
fi

if [[ -z "${NCBI_EMAIL:-}" ]]; then
    echo "Please set the NCBI_EMAIL environment variable before running remote BLAST."
    echo "export NCBI_EMAIL='you@example.com'"
    exit 1
fi

API_FLAG=()
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    API_FLAG=(-api_key "$NCBI_API_KEY")
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

blastn \
    -query "$UNIQUE_FASTA" \
    -db nt \
    -remote \
    -out "${OUT_PREFIX}.tsv" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms stitle" \
    -max_target_seqs "$MAX_TARGET_SEQS" \
    -evalue "$EVALUE" \
    -entrez_query "$TAX_FILTER" \
    "${API_FLAG[@]}"

echo "[ok] Remote BLAST complete: ${OUT_PREFIX}.tsv"
