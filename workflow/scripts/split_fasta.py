#!/usr/bin/env python3
"""
Split FASTA file into batches while preserving sequence order and IDs.

Output: batch FASTA files + batch_manifest.txt.
"""

import argparse
import os
import sys
from pathlib import Path

from Bio import SeqIO


def _write_manifest(out_dir: Path, entries: list[tuple[str, list[str]]]) -> None:
    with open(out_dir / "batch_manifest.txt", "w") as f:
        f.write("# batch_id\tseq_count\tseq_ids\n")
        for batch_id, seq_ids in entries:
            f.write(f"{batch_id}\t{len(seq_ids)}\t{','.join(seq_ids)}\n")


def split_fasta(
    input_fasta: Path,
    output_dir: Path,
    batch_size: int,
) -> tuple:
    """Return (n_batches, n_sequences) after writing batch FASTAs."""
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_idx = 0
    current_batch: list = []
    current_ids: list[str] = []
    batch_entries: list[tuple[str, list[str]]] = []
    total = 0

    for record in SeqIO.parse(str(input_fasta), "fasta"):
        total += 1

        current_batch.append(record)
        current_ids.append(record.id)

        if len(current_batch) == batch_size:
            bid = f"{batch_idx:05d}"
            SeqIO.write(current_batch, str(output_dir / f"batch_{bid}.fa"), "fasta")
            batch_entries.append((bid, current_ids))
            batch_idx += 1
            current_batch, current_ids = [], []

    if current_batch:
        bid = f"{batch_idx:05d}"
        SeqIO.write(current_batch, str(output_dir / f"batch_{bid}.fa"), "fasta")
        batch_entries.append((bid, current_ids))
        batch_idx += 1

    _write_manifest(output_dir, batch_entries)

    return batch_idx, total


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_fasta", help="Input FASTA file")
    p.add_argument("output_dir", help="Directory for batch FASTA files")
    p.add_argument("batch_size", type=int, help="Sequences per batch")
    args = p.parse_args()

    if not os.path.isfile(args.input_fasta):
        print(f"Error: input file not found: {args.input_fasta}", file=sys.stderr)
        sys.exit(1)
    if args.batch_size <= 0:
        print("Error: batch_size must be positive", file=sys.stderr)
        sys.exit(1)

    n_batches, total = split_fasta(
        Path(args.input_fasta),
        Path(args.output_dir),
        args.batch_size,
    )
    print(f"Split {total} sequences into {n_batches} batches")


if __name__ == "__main__":
    main()
