# Plant_from_bug

Pipeline for extracting plant barcodes (matK and rbcL) from insect Illumina
metagenomes, collapsing unique hits, and identifying candidate host plants via
GenBank `nt`.

## Repository Layout

- `plant_barcodes/`: provide reference barcode FASTA files.
  - matK files must match `matK_*.fasta`
  - rbcL files must match `rbcL*.fasta`
- `plant_genes_Nov25/data/`: place insect metagenomes (`*_R1_R2.fastq.gz`)
- `run_pipeline.sh`: orchestrates the workflow end-to-end
- `scripts/prepare_fastas.py`: converts compressed FASTQ libraries to FASTA
- `scripts/collect_unique_hits.py`: collapses BLAST output to unique hits
- `scripts/blast_nt_hits.sh`: helper to BLAST unique hits against GenBank `nt`
- `scripts/slurm_blast_task.sh`: Slurm array worker for per-sample BLAST jobs
- `scripts/slurm_fastq_to_fasta.sh`: Slurm array worker for FASTQ to FASTA conversion

## Dependencies

- Python 3.8+
- NCBI BLAST+ (`blastn`, `makeblastdb`)
- GNU coreutils (sed, cat, etc.; available on most Linux distros)
- Slurm client utilities (`sbatch`, `squeue`; `sacct` is optional but used when present)

## Running the Pipeline

```bash
# optional: customize threads/e-value
module load ncbi-blast/2.16.0+       # if not already in your shell startup
export THREADS=32
export EVALUE=1e-3
# optional: Slurm/resource tuning
export FASTQ_BATCH_SIZE=10          # max concurrent FASTQ->FASTA jobs
export GENOME_BATCH_SIZE=10        # max concurrent insect genomes
export BLAST_THREADS=3             # threads per BLAST task
export SLURM_MEM_PER_CPU=10G              # whatever per-core memory you need
# (you can mix in other options, e.g. SLURM_PARTITION, SLURM_TIME, etc.)
# export SLURM_PARTITION=general   # uncomment to target a partition

# Analyses stay inside each dataset directory (e.g. plant_genes_Nov25/results)
bash run_pipeline.sh plant_barcodes plant_genes_Nov25

# For a different batch (e.g. December genomes)
bash run_pipeline.sh plant_barcodes plant_genes_Dec25
```

Notes:

- Pass the dataset root (`plant_genes_Nov25`). The script looks for `data/`
  inside that directory; if you instead provide the `data/` path directly, it is
  detected automatically.
- Outputs land in `<dataset>/results` by default. Provide a third argument (e.g.
  `analysis_run2`) to place them in `<dataset>/analysis_run2`. Absolute paths
  are honored as-is.
- `run_pipeline.sh` submits Slurm job arrays for both the FASTQ conversion and BLAST stages; run it from a
  login/submit node with access to your shared filesystem.

What happens:

1. `*_R1_R2.fastq.gz` libraries are streamed into `results/fastas/*.fasta`
   through a Slurm array (10 concurrent conversions by default)
2. All matK and rbcL barcode references are concatenated per gene
3. Each FASTA becomes its own `makeblastdb` target
4. Slurm job arrays process the FASTA manifest ~10 samples at a time (tunable)
   and run `blastn` (default e-value `1e-3`) for matK/rbcL queries
5. BLAST results are combined per gene and deduplicated into unique hits

Key outputs (relative to each dataset directory):

- `results/blast/`: raw BLAST tables per sample and combined `*_all.tsv`
- `results/unique/matK_unique_hits.fasta` / `rbcL_unique_hits.fasta`: deduplicated hit sequences
- `results/unique/*.tsv`: metadata for each unique hit (pident, bitscore, etc.)
- Intermediates (`results/fastas`, `results/blastdbs`) are deleted automatically
  at the end of each run to save space—rerunning the pipeline regenerates them.
- Slurm stdout/stderr for each BLAST task lands in
  `results/blast/<gene>/slurm-<job>_<task>.out`.
- Slurm stdout/stderr for each FASTQ conversion task lands in
  `results/fastas/slurm_fastq/slurm-<job>_<task>.out`.
- Manifests in `results/manifests/` track FASTQ inputs and FASTA outputs used by the arrays.

### Slurm scheduling details

- `FASTQ_BATCH_SIZE` (default `10`) controls how many FASTQ conversion jobs run
  simultaneously. Set to `0` or a negative value to allow unlimited concurrency.
- `GENOME_BATCH_SIZE` (default `10`) caps how many insect genomes run
  simultaneously. Set to `0` or a negative value to allow unlimited concurrency.
- `BLAST_THREADS` (or `THREADS`) controls the per-task `blastn -num_threads`.
  Ensure `GENOME_BATCH_SIZE * BLAST_THREADS` fits within your allocation.
- Optional environment variables forwarded to `sbatch`:
  `SLURM_PARTITION`, `SLURM_TIME`, `SLURM_MEM_PER_CPU`,
  `SLURM_CPUS_PER_TASK`, `SLURM_QOS`, plus any extra flags via
  `SLURM_SBATCH_OPTS` (space-separated string).
- `SLURM_WAIT_POLL` (seconds, default `30`) sets how often the pipeline polls
  `squeue` to wait for arrays to finish.
- Set `TASK_INDEX` to a zero-based value if you need to run
  `scripts/slurm_blast_task.sh` manually for debugging outside of Slurm.

## Identifying Plant Taxa via GenBank nt

Remote BLAST requires an email (and benefits from an NCBI API key).

```bash
export NCBI_EMAIL="you@example.com"
# optional
#export NCBI_API_KEY="XXXX"

./scripts/blast_nt_hits.sh results/unique/matK_unique_hits.fasta results/nt/matK_vs_nt
./scripts/blast_nt_hits.sh results/unique/rbcL_unique_hits.fasta results/nt/rbcL_vs_nt
```

Outputs are tab-delimited tables enriched with taxonomy columns (`staxids`,
`sscinames`, `sskingdoms`, `stitle`). Adjust `MAX_TARGET_SEQS` or `TAX_FILTER`
(defaults to `Viridiplantae[ORGN]`) through environment variables if needed.

## Customization Tips

- Change `FASTQ_SUFFIX` if your libraries follow a different naming scheme.
- Override `THREADS` to match the number of CPU cores available.
- Use `scripts/prepare_fastas.py --force ...` to re-generate FASTA files.
- Tweak `EVALUE` from the environment when invoking `run_pipeline.sh` or
  `scripts/blast_nt_hits.sh`.
- `GENOME_BATCH_SIZE`, `BLAST_THREADS`, and the `SLURM_*` environment variables
  let you shape how aggressively the BLAST stage submits work to your cluster.
