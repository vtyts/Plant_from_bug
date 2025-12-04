#!/usr/bin/env bash
# Main pipeline to extract plant barcode sequences from insect Illumina data.
set -euo pipefail
module load ncbi-blast/2.16.0+

PLANT_BARCODE_DIR=${1:-"plant_barcodes"}
DATASET_ARG=${2:-"plant_genes_Nov25"}
DATASET_ARG=${DATASET_ARG%/}

THREADS=${THREADS:-32}
EVALUE=${EVALUE:-1e-3}
FASTQ_SUFFIX=${FASTQ_SUFFIX:-"_R1_R2.fastq.gz"}
GENOME_BATCH_SIZE=${GENOME_BATCH_SIZE:-10}
SLURM_WAIT_POLL=${SLURM_WAIT_POLL:-30}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"

abs_path() {
    local target=$1
    if [[ "$target" = /* ]]; then
        printf "%s\n" "$target"
    else
        printf "%s/%s\n" "$(pwd)" "${target#./}"
    fi
}

declare -a SBATCH_EXTRA_ARGS=()
if [[ -n "${SLURM_PARTITION:-}" ]]; then
    SBATCH_EXTRA_ARGS+=("--partition" "$SLURM_PARTITION")
fi
if [[ -n "${SLURM_ACCOUNT:-}" ]]; then
    SBATCH_EXTRA_ARGS+=("--account" "$SLURM_ACCOUNT")
fi
if [[ -n "${SLURM_TIME:-}" ]]; then
    SBATCH_EXTRA_ARGS+=("--time" "$SLURM_TIME")
fi
if [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
    SBATCH_EXTRA_ARGS+=("--mem-per-cpu" "$SLURM_MEM_PER_CPU")
fi
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
    SBATCH_EXTRA_ARGS+=("--cpus-per-task" "$SLURM_CPUS_PER_TASK")
fi
if [[ -n "${SLURM_QOS:-}" ]]; then
    SBATCH_EXTRA_ARGS+=("--qos" "$SLURM_QOS")
fi
if [[ -n "${SLURM_SBATCH_OPTS:-}" ]]; then
    # shellcheck disable=SC2206
    extra_opts=($SLURM_SBATCH_OPTS)
    SBATCH_EXTRA_ARGS+=("${extra_opts[@]}")
fi

if [[ -d "$DATASET_ARG" && -d "$DATASET_ARG/data" ]]; then
    DATASET_ROOT="$DATASET_ARG"
    INSECT_FASTQ_DIR="$DATASET_ARG/data"
elif [[ -d "$DATASET_ARG" ]]; then
    INSECT_FASTQ_DIR="$DATASET_ARG"
    if [[ "$DATASET_ARG" == */data ]]; then
        DATASET_ROOT="${DATASET_ARG%/data}"
    else
        DATASET_ROOT="$DATASET_ARG"
    fi
else
    echo "Insect dataset directory not found: $DATASET_ARG"
    exit 1
fi

if [[ -n "${3:-}" ]]; then
    if [[ "$3" = /* ]]; then
        OUTPUT_DIR="$3"
    else
        OUTPUT_DIR="$DATASET_ROOT/$3"
    fi
else
    OUTPUT_DIR="$DATASET_ROOT/results"
fi
OUTPUT_DIR=$(abs_path "$OUTPUT_DIR")

required_bins=(python3 blastn makeblastdb sbatch squeue)
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

echo "==> Dataset root: $DATASET_ROOT"
echo "==> FASTQ source: $INSECT_FASTQ_DIR"

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

shopt -s nullglob
FASTAS=("$FASTAS_DIR"/*.fasta)
shopt -u nullglob
if [[ ${#FASTAS[@]} -eq 0 ]]; then
    echo "No FASTA files produced in $FASTAS_DIR"
    exit 1
fi

MANIFEST_DIR="$OUTPUT_DIR/manifests"
MANIFEST_FILE="$MANIFEST_DIR/insect_fastas.tsv"
mkdir -p "$MANIFEST_DIR"
: >"$MANIFEST_FILE"
for fasta in "${FASTAS[@]}"; do
    sample=$(basename "$fasta" .fasta)
    printf "%s\t%s\n" "$sample" "$fasta" >>"$MANIFEST_FILE"
done
MANIFEST_FILE=$(abs_path "$MANIFEST_FILE")
TOTAL_SAMPLES=${#FASTAS[@]}
echo "==> Prepared manifest for $TOTAL_SAMPLES FASTA libraries"

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

wait_for_job() {
    local job_id=$1
    if [[ -z "$job_id" ]]; then
        echo "Missing Slurm job ID to wait on."
        exit 1
    fi

    while true; do
        local remaining
        remaining=$(squeue --noheader -j "$job_id" 2>/dev/null || true)
        if [[ -z "$remaining" ]]; then
            break
        fi
        sleep "$SLURM_WAIT_POLL"
    done

    if command -v sacct >/dev/null 2>&1; then
        local state
        state=$(sacct -j "${job_id}" --format=State --noheader 2>/dev/null | awk 'NR==1 {gsub(/[[:space:]]/, "", $1); print $1}')
        case "$state" in
        ""|PENDING|RUNNING|COMPLETED|COMPLETING|COMPLETED*)
            ;;
        CANCELLED*|FAILED|TIMEOUT|NODE_FAIL|PREEMPTED)
            echo "Slurm job ${job_id} finished with state ${state}"
            exit 1
            ;;
        *)
            ;;
        esac
    fi
}

submit_blast_array() {
    local gene=$1
    local query_abs
    query_abs=$(abs_path "$OUTPUT_DIR/${gene}_barcodes.fasta")
    local outdir_abs
    outdir_abs=$(abs_path "$BLAST_DIR/$gene")
    local blastdb_abs
    blastdb_abs=$(abs_path "$BLASTDB_DIR")
    mkdir -p "$outdir_abs"

    local total_tasks
    total_tasks=$(wc -l <"$MANIFEST_FILE" | tr -d ' ')
    if (( total_tasks == 0 )); then
        echo "No entries found in manifest $MANIFEST_FILE"
        exit 1
    fi

    local array_spec="0-$((total_tasks - 1))"
    local concurrency_note="$GENOME_BATCH_SIZE"
    if (( GENOME_BATCH_SIZE > 0 )); then
        array_spec="${array_spec}%${GENOME_BATCH_SIZE}"
    else
        concurrency_note="unlimited"
    fi

    local sbatch_cmd=(sbatch --parsable --export=ALL --array="$array_spec" --job-name "blast_${gene}" --output "$outdir_abs/slurm-%A_%a.out")
    if [[ ${#SBATCH_EXTRA_ARGS[@]} -gt 0 ]]; then
        sbatch_cmd+=("${SBATCH_EXTRA_ARGS[@]}")
    fi
    sbatch_cmd+=("$SCRIPT_DIR/slurm_blast_task.sh" "$MANIFEST_FILE" "$gene" "$query_abs" "$outdir_abs" "$blastdb_abs")

    echo "==> Launching Slurm array for ${gene} (${total_tasks} samples; max ${concurrency_note} concurrent)"
    local job_id
    job_id=$("${sbatch_cmd[@]}") || {
        echo "Failed to submit Slurm job for ${gene}"
        exit 1
    }
    echo "==> Slurm job ${job_id} submitted for ${gene}"
    wait_for_job "$job_id"
    echo "==> Slurm job ${job_id} completed for ${gene}"

    shopt -s nullglob
    local tsv_files=("$outdir_abs"/*.tsv)
    shopt -u nullglob
    if [[ ${#tsv_files[@]} -eq 0 ]]; then
        echo "No BLAST TSV files created in $outdir_abs for ${gene}"
        exit 1
    fi

    cat "${tsv_files[@]}" >"$BLAST_DIR/${gene}_all.tsv"
    echo "==> Aggregated ${gene} hits -> $BLAST_DIR/${gene}_all.tsv"
}

submit_blast_array "matK"
submit_blast_array "rbcL"

echo "==> Deriving unique hits"
python3 "$SCRIPT_DIR/collect_unique_hits.py" \
    --blast-tsv "$BLAST_DIR/matK_all.tsv" \
    --unique-fasta "$UNIQUE_DIR/matK_unique_hits.fasta" \
    --unique-table "$UNIQUE_DIR/matK_unique_hits.tsv"

python3 "$SCRIPT_DIR/collect_unique_hits.py" \
    --blast-tsv "$BLAST_DIR/rbcL_all.tsv" \
    --unique-fasta "$UNIQUE_DIR/rbcL_unique_hits.fasta" \
    --unique-table "$UNIQUE_DIR/rbcL_unique_hits.tsv"

echo "==> Cleaning up intermediate FASTA files and BLAST databases"
rm -rf "$FASTAS_DIR" "$BLASTDB_DIR"

cat <<EOF
Pipeline complete.
- Unique matK hits: $UNIQUE_DIR/matK_unique_hits.fasta
- Unique rbcL hits: $UNIQUE_DIR/rbcL_unique_hits.fasta

To identify plant taxa via GenBank nt, run for each file:
  export NCBI_EMAIL="you@example.com"
  ./scripts/blast_nt_hits.sh results/unique/matK_unique_hits.fasta results/nt/matK_vs_nt
EOF
