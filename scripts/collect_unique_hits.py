#!/usr/bin/env python3
"""
Collapse BLAST tabular output into unique hit sets and emit FASTA + summary TSV.

Expected BLAST format:
6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from textwrap import wrap
from typing import Dict, List

BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "sseq",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate unique BLAST hits and companion FASTA/TSV outputs."
    )
    parser.add_argument(
        "--blast-tsv",
        required=True,
        type=Path,
        help="Combined BLAST tabular file (outfmt 6 with sseq column).",
    )
    parser.add_argument(
        "--unique-fasta",
        required=True,
        type=Path,
        help="Output FASTA path for unique sequences.",
    )
    parser.add_argument(
        "--unique-table",
        required=True,
        type=Path,
        help="Output TSV path describing each unique hit.",
    )
    parser.add_argument(
        "--allow-empty",
        action="store_true",
        help="Permit BLAST inputs with zero hits and emit empty outputs instead of exiting.",
    )
    return parser.parse_args()


def load_hits(blast_tsv: Path) -> List[Dict[str, str]]:
    hits: List[Dict[str, str]] = []
    with blast_tsv.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line in reader:
            if not line:
                continue
            if len(line) != len(BLAST_COLUMNS):
                raise ValueError(
                    f"Expected {len(BLAST_COLUMNS)} columns but found {len(line)} "
                    f"in line: {line}"
                )
            hits.append(dict(zip(BLAST_COLUMNS, line)))
    return hits


def select_unique_hits(hits: List[Dict[str, str]]) -> Dict[str, Dict[str, str]]:
    """Return a mapping of sequence -> best hit metadata (highest bitscore)."""
    unique: Dict[str, Dict[str, str]] = {}
    for hit in hits:
        seq = hit["sseq"]
        bitscore = float(hit["bitscore"])
        current = unique.get(seq)
        if current is None or bitscore > float(current["bitscore"]):
            unique[seq] = hit
    return unique


def write_fasta(records: Dict[str, Dict[str, str]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w") as fh:
        for idx, (sequence, hit) in enumerate(records.items(), start=1):
            header = (
                f"hit{idx}|query:{hit['qseqid']}|subject:{hit['sseqid']}|"
                f"pident:{hit['pident']}|bitscore:{hit['bitscore']}"
            )
            fh.write(f">{header}\n")
            for chunk in wrap(sequence, width=80):
                fh.write(f"{chunk}\n")


def write_table(records: Dict[str, Dict[str, str]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["unique_id", *BLAST_COLUMNS])
        for idx, (sequence, hit) in enumerate(records.items(), start=1):
            unique_id = f"hit{idx}"
            writer.writerow([unique_id, *[hit[col] for col in BLAST_COLUMNS[:-1]], sequence])


def main() -> None:
    args = parse_args()
    if not args.blast_tsv.exists():
        raise SystemExit(f"BLAST file not found: {args.blast_tsv}")

    hits = load_hits(args.blast_tsv)
    if not hits:
        if args.allow_empty:
            write_fasta({}, args.unique_fasta)
            write_table({}, args.unique_table)
            print(
                f"[warn] No hits found in {args.blast_tsv}. "
                f"Emitted empty outputs ({args.unique_fasta}, {args.unique_table})"
            )
            return
        raise SystemExit(f"No hits found in {args.blast_tsv}.")

    unique_hits = select_unique_hits(hits)
    write_fasta(unique_hits, args.unique_fasta)
    write_table(unique_hits, args.unique_table)
    print(
        f"[ok] Unique hits: {len(unique_hits)} "
        f"(FASTA: {args.unique_fasta}, TSV: {args.unique_table})"
    )


if __name__ == "__main__":
    main()
