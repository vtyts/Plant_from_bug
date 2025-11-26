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

## Dependencies

- Python 3.8+
- NCBI BLAST+ (`blastn`, `makeblastdb`)
- GNU coreutils (sed, cat, etc.; available on most Linux distros)

## Running the Pipeline

```bash
# optional: customize threads/e-value
export THREADS=32
export EVALUE=1e-3

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

What happens:

1. `*_R1_R2.fastq.gz` libraries are streamed into `results/fastas/*.fasta`
2. All matK and rbcL barcode references are concatenated per gene
3. Each FASTA becomes its own `makeblastdb` target
4. `blastn` (e-value 1e-3, multithreaded) runs matK/rbcL queries against every sample
5. BLAST results are combined per gene and deduplicated into unique hits

Key outputs (relative to each dataset directory):

- `results/blast/`: raw BLAST tables per sample and combined `*_all.tsv`
- `results/unique/matK_unique_hits.fasta` / `rbcL_unique_hits.fasta`: deduplicated hit sequences
- `results/unique/*.tsv`: metadata for each unique hit (pident, bitscore, etc.)
- Intermediates (`results/fastas`, `results/blastdbs`) are deleted automatically
  at the end of each run to save space—rerunning the pipeline regenerates them.

## Identifying Plant Taxa via GenBank nt

Remote BLAST requires an email (and benefits from an NCBI API key).

```bash
export NCBI_EMAIL="you@example.com"
# optional
export NCBI_API_KEY="XXXX"

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
