#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <manifest.tsv> <prepare_fastas.py>"
    exit 1
fi

MANIFEST=$1
CONVERTER=$2

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

IFS=$'\t' read -r SAMPLE FASTQ_PATH FASTA_PATH <<<"$ENTRY"
if [[ -z "$SAMPLE" || -z "$FASTQ_PATH" || -z "$FASTA_PATH" ]]; then
    echo "Malformed manifest entry on line $LINE_NUMBER: $ENTRY"
    exit 1
fi

if [[ ! -f "$FASTQ_PATH" ]]; then
    echo "FASTQ not found for sample $SAMPLE: $FASTQ_PATH"
    exit 1
fi

python3 "$CONVERTER" \
    --fastq-file "$FASTQ_PATH" \
    --output-fasta "$FASTA_PATH" \
    --force

echo "[ok] FASTQ -> FASTA complete for ${SAMPLE}"
