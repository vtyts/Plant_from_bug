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
from typing import Dict, List, Sequence

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


def should_replace(existing: Dict[str, str], candidate: Dict[str, str]) -> bool:
    """Return True if candidate should replace existing for the same subject ID."""
    exist_len = int(existing["length"])
    cand_len = int(candidate["length"])
    if cand_len > exist_len:
        return True
    if cand_len < exist_len:
        return False

    # Tie-break on bitscore (higher is better), then e-value (lower is better)
    exist_bits = float(existing["bitscore"])
    cand_bits = float(candidate["bitscore"])
    if cand_bits > exist_bits:
        return True
    if cand_bits < exist_bits:
        return False

    exist_eval = float(existing["evalue"])
    cand_eval = float(candidate["evalue"])
    if cand_eval < exist_eval:
        return True
    if cand_eval > exist_eval:
        return False

    # Preserve the most recent candidate when all metrics are equal.
    return True


def select_unique_hits(hits: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Return the best hit per subject (sseqid), preferring the longest alignment."""
    unique: Dict[str, Dict[str, str]] = {}
    for hit in hits:
        subject = hit["sseqid"]
        current = unique.get(subject)
        if current is None or should_replace(current, hit):
            unique[subject] = hit
    return list(unique.values())


def write_fasta(records: Sequence[Dict[str, str]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w") as fh:
        for idx, hit in enumerate(records, start=1):
            header = (
                f"hit{idx}|query:{hit['qseqid']}|subject:{hit['sseqid']}|"
                f"pident:{hit['pident']}|bitscore:{hit['bitscore']}"
            )
            fh.write(f">{header}\n")
            for chunk in wrap(hit["sseq"], width=80):
                fh.write(f"{chunk}\n")


def write_table(records: Sequence[Dict[str, str]], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["unique_id", *BLAST_COLUMNS])
        for idx, hit in enumerate(records, start=1):
            unique_id = f"hit{idx}"
            writer.writerow([unique_id, *[hit[col] for col in BLAST_COLUMNS[:-1]], hit["sseq"]])


def main() -> None:
    args = parse_args()
    if not args.blast_tsv.exists():
        raise SystemExit(f"BLAST file not found: {args.blast_tsv}")

    hits = load_hits(args.blast_tsv)
    if not hits:
        print(f"[skip] No hits found in {args.blast_tsv}; skipping unique outputs.")
        return

    unique_hits = select_unique_hits(hits)
    if not unique_hits:
        print(f"[skip] No unique hits found in {args.blast_tsv}; skipping unique outputs.")
        return

    write_fasta(unique_hits, args.unique_fasta)
    write_table(unique_hits, args.unique_table)
    print(
        f"[ok] Unique hits: {len(unique_hits)} "
        f"(FASTA: {args.unique_fasta}, TSV: {args.unique_table})"
    )


if __name__ == "__main__":
    main()
