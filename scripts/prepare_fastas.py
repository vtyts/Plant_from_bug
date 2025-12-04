#!/usr/bin/env python3
"""
Stream-convert compressed FASTQ libraries into FASTA files.

Each Illumina library is expected to be stored as "<sample>_R1_R2.fastq.gz".
The script reads the file in chunks to avoid loading everything into memory
and writes a FASTA file with the matching sample name.

You can either convert an entire directory (default behavior) or a single
FASTQ file by providing the `--fastq-file` / `--output-fasta` pair.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterator


def iter_fastq_records(handle: gzip.GzipFile) -> Iterator[tuple[str, str]]:
    """Yield (header_without_at, sequence) tuples from a FASTQ stream."""
    while True:
        header = handle.readline()
        if not header:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()

        if not (seq and plus and qual):
            raise ValueError("FASTQ file ended unexpectedly; check input integrity.")

        header_str = header.decode().strip()
        if not header_str.startswith("@"):
            raise ValueError(f"Malformed FASTQ header: {header_str}")

        seq_str = seq.decode().strip()
        yield header_str[1:].split()[0], seq_str


def convert_file(src: Path, dest: Path, force: bool) -> None:
    """Convert a gzip-compressed FASTQ file to FASTA."""
    if dest.exists() and not force:
        print(f"[skip] {dest} already exists (use --force to overwrite)")
        return

    dest.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fh_in, dest.open("w") as fh_out:
        record_count = 0
        for header, seq in iter_fastq_records(fh_in):
            fh_out.write(f">{header}\n{seq}\n")
            record_count += 1
        print(f"[ok] {src.name} -> {dest.name} ({record_count} records)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert *_R1_R2.fastq.gz files to FASTA for downstream BLAST."
    )
    parser.add_argument(
        "--input",
        type=Path,
        help="Directory containing *_R1_R2.fastq.gz files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Destination directory for FASTA files.",
    )
    parser.add_argument(
        "--suffix",
        default="_R1_R2.fastq.gz",
        help="File name suffix used to detect FASTQ inputs (default: %(default)s).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing FASTA files instead of skipping them.",
    )
    parser.add_argument(
        "--fastq-file",
        type=Path,
        help="Convert this single FASTQ.gz file instead of scanning a directory.",
    )
    parser.add_argument(
        "--output-fasta",
        type=Path,
        help="Destination FASTA path when using --fastq-file.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    single_mode = args.fastq_file is not None or args.output_fasta is not None
    if single_mode:
        if args.fastq_file is None or args.output_fasta is None:
            raise SystemExit("--fastq-file and --output-fasta must be provided together.")
        if not args.fastq_file.exists():
            raise SystemExit(f"FASTQ file not found: {args.fastq_file}")
        convert_file(args.fastq_file, args.output_fasta, args.force)
        return

    if args.input is None or args.output is None:
        raise SystemExit(
            "Provide --input/--output for directory conversion or "
            "--fastq-file/--output-fasta for single-file conversion."
        )

    if not args.input.is_dir():
        raise SystemExit(f"Input directory not found: {args.input}")

    fastq_files = sorted(args.input.glob(f"*{args.suffix}"))
    if not fastq_files:
        raise SystemExit(
            f"No files ending with '{args.suffix}' found in {args.input}. "
            "Check the directory path or suffix."
        )

    for fastq_gz in fastq_files:
        sample_name = fastq_gz.name[: -len(args.suffix)]
        fasta_path = args.output / f"{sample_name}.fasta"
        convert_file(fastq_gz, fasta_path, args.force)


if __name__ == "__main__":
    main()
