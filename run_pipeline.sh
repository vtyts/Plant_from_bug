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
FASTQ_BATCH_SIZE=${FASTQ_BATCH_SIZE:-10}
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
MANIFEST_DIR="$OUTPUT_DIR/manifests"

mkdir -p "$FASTAS_DIR" "$BLASTDB_DIR" "$BLAST_DIR/matK" "$BLAST_DIR/rbcL" "$UNIQUE_DIR" "$MANIFEST_DIR"

FASTQ_MANIFEST="$MANIFEST_DIR/fastq_inputs.tsv"
: >"$FASTQ_MANIFEST"

shopt -s nullglob
FASTQ_FILES=("$INSECT_FASTQ_DIR"/*"$FASTQ_SUFFIX")
shopt -u nullglob
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "No FASTQ files ending with '$FASTQ_SUFFIX' found in $INSECT_FASTQ_DIR"
    exit 1
fi

for fastq in "${FASTQ_FILES[@]}"; do
    sample=$(basename "$fastq")
    sample=${sample%"$FASTQ_SUFFIX"}
    fasta_path="$FASTAS_DIR/${sample}.fasta"
    printf "%s\t%s\t%s\n" "$sample" "$(abs_path "$fastq")" "$fasta_path" >>"$FASTQ_MANIFEST"
done
FASTQ_MANIFEST=$(abs_path "$FASTQ_MANIFEST")
TOTAL_FASTQ=${#FASTQ_FILES[@]}
echo "==> Prepared FASTQ manifest for $TOTAL_FASTQ libraries"

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

submit_fastq_conversion() {
    local manifest=$1
    local total_tasks=$2
    if (( total_tasks == 0 )); then
        echo "No FASTQ entries to convert."
        return
    fi

    local array_spec="0-$((total_tasks - 1))"
    local concurrency_note="$FASTQ_BATCH_SIZE"
    if (( FASTQ_BATCH_SIZE > 0 )); then
        array_spec="${array_spec}%${FASTQ_BATCH_SIZE}"
    else
        concurrency_note="unlimited"
    fi

    local log_dir
    log_dir=$(abs_path "$FASTAS_DIR/slurm_fastq")
    mkdir -p "$log_dir"

    local sbatch_cmd=(sbatch --parsable --export=ALL --array="$array_spec" --job-name "fastq2fa" --output "$log_dir/slurm-%A_%a.out")
    if [[ ${#SBATCH_EXTRA_ARGS[@]} -gt 0 ]]; then
        sbatch_cmd+=("${SBATCH_EXTRA_ARGS[@]}")
    fi
    sbatch_cmd+=("$SCRIPT_DIR/slurm_fastq_to_fasta.sh" "$manifest" "$SCRIPT_DIR/prepare_fastas.py")

    echo "==> Launching Slurm array for FASTQ -> FASTA (${total_tasks} libraries; max ${concurrency_note} concurrent)"
    local job_id
    job_id=$("${sbatch_cmd[@]}") || {
        echo "Failed to submit Slurm job for FASTQ conversion"
        exit 1
    }
    echo "==> Slurm job ${job_id} submitted for FASTQ conversion"
    wait_for_job "$job_id"
    echo "==> FASTQ conversion job ${job_id} completed"
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

}

derive_per_sample_unique_hits() {
    local gene=$1
    local gene_blast_dir="$BLAST_DIR/$gene"
    local sample_out_dir="$UNIQUE_DIR/by_sample/$gene"

    mkdir -p "$sample_out_dir"

    shopt -s nullglob
    local sample_files=("$gene_blast_dir"/*.tsv)
    shopt -u nullglob

    if [[ ${#sample_files[@]} -eq 0 ]]; then
        echo "No per-sample BLAST TSV files found for ${gene} in $gene_blast_dir"
        exit 1
    fi

    local processed=0
    for tsv in "${sample_files[@]}"; do
        local sample
        sample=$(basename "$tsv")
        sample=${sample%.tsv}
        if [[ ! -s "$tsv" ]]; then
            echo "[skip] ${gene} BLAST output for ${sample} has no hits; skipping."
            continue
        fi
        python3 "$SCRIPT_DIR/collect_unique_hits.py" \
            --blast-tsv "$tsv" \
            --unique-fasta "$sample_out_dir/${sample}_${gene}_unique_hits.fasta" \
            --unique-table "$sample_out_dir/${sample}_${gene}_unique_hits.tsv"
        ((processed++))
    done
    echo "==> Deduplicated per-sample ${gene} hits (${processed} samples) -> $sample_out_dir"
}

submit_fastq_conversion "$FASTQ_MANIFEST" "$TOTAL_FASTQ"

MANIFEST_FILE="$MANIFEST_DIR/insect_fastas.tsv"
: >"$MANIFEST_FILE"
while IFS=$'\t' read -r sample _ fasta_path; do
    [[ -z "${sample:-}" ]] && continue
    if [[ ! -f "$fasta_path" ]]; then
        echo "Missing FASTA output for $sample at $fasta_path"
        exit 1
    fi
    printf "%s\t%s\n" "$sample" "$fasta_path" >>"$MANIFEST_FILE"
done <"$FASTQ_MANIFEST"
MANIFEST_FILE=$(abs_path "$MANIFEST_FILE")
TOTAL_SAMPLES=$(wc -l <"$MANIFEST_FILE" | tr -d ' ')
if (( TOTAL_SAMPLES == 0 )); then
    echo "No FASTA entries found in $MANIFEST_FILE"
    exit 1
fi
echo "==> Prepared manifest for $TOTAL_SAMPLES FASTA libraries"

submit_blast_array "matK"
submit_blast_array "rbcL"

echo "==> Deriving per-sample unique hits"
derive_per_sample_unique_hits "matK"
derive_per_sample_unique_hits "rbcL"

echo "==> Cleaning up intermediate FASTA files and BLAST databases"
rm -rf "$FASTAS_DIR" "$BLASTDB_DIR"

cat <<EOF
Pipeline complete.
- Per-sample matK hits: $UNIQUE_DIR/by_sample/matK
- Per-sample rbcL hits: $UNIQUE_DIR/by_sample/rbcL

To identify plant taxa via GenBank nt, run for each file:
  export NCBI_EMAIL="you@example.com"
  ./scripts/blast_nt_hits.sh results/unique/by_sample/matK/<sample>_matK_unique_hits.fasta results/nt/matK_vs_nt
EOF
