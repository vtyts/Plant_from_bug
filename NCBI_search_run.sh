#!/bin/bash
CODE=$1
MARK=$2
DATE=$3
echo "This script will run blast_nt_hits.sh for batches of samples from results/unique for a marker of choice"
echo "The library code (00-00-) is $1"
echo "The plant marker is $2"
echo "The analysis folder is $3"
cd "./$DATE/results/unique/by_sample/$MARK"
for file in $(ls $CODE*.fasta)
    do
    bash ../../../../../scripts/blast_nt_hits.sh "${file%%_*}_${MARK}_unique_hits.fasta" "../../nt/${file%%_*}_${MARK}_vs_nt"
    done
    