#!/usr/bin/env python3
"""
Split FASTA file into batches while preserving sequence order and IDs.
Sequences exceeding --max-seq-len are routed to --oversized-dir as
single-sequence FASTAs and excluded from the normal batches.

Output: batch FASTA files + batch_manifest.txt (+ oversized FASTAs and their
manifest when --max-seq-len is given).
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

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
    max_seq_len: Optional[int] = None,
    oversized_dir: Optional[Path] = None,
) -> tuple:
    """Split FASTA, routing oversized sequences to oversized_dir.

    Returns (n_batches, n_sequences, n_oversized).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_idx = 0
    oversized_idx = 0
    current_batch: list = []
    current_ids: list[str] = []
    batch_entries: list[tuple[str, list[str]]] = []
    oversized_entries: list[tuple[str, list[str]]] = []
    skip_rows: list[tuple[str, int]] = []
    total = 0

    for record in SeqIO.parse(str(input_fasta), "fasta"):
        total += 1

        if max_seq_len is not None and len(record.seq) > max_seq_len:
            oversized_dir.mkdir(parents=True, exist_ok=True)
            ob_id = f"oversized_{oversized_idx:05d}"
            SeqIO.write([record], str(oversized_dir / f"batch_{ob_id}.fa"), "fasta")
            oversized_entries.append((ob_id, [record.id]))
            skip_rows.append((record.id, len(record.seq)))
            oversized_idx += 1
            continue

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

    if oversized_dir is not None:
        oversized_dir.mkdir(parents=True, exist_ok=True)
        _write_manifest(oversized_dir, oversized_entries)
        skip_tsv = oversized_dir.parent / "scanfold_oversized_skip.tsv"
        with open(skip_tsv, "w") as f:
            f.write("seq_id\tseq_len\tsource_batch\n")
            for seq_id, seq_len in skip_rows:
                f.write(f"{seq_id}\t{seq_len}\tinitial_split\n")

    return batch_idx, total, oversized_idx


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_fasta", help="Input FASTA file")
    p.add_argument("output_dir", help="Directory for batch FASTA files")
    p.add_argument("batch_size", type=int, help="Sequences per batch")
    p.add_argument(
        "--max-seq-len",
        type=int,
        default=None,
        metavar="N",
        help="Route sequences longer than N nt to --oversized-dir",
    )
    p.add_argument(
        "--oversized-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Directory for single-sequence oversized FASTAs (required with --max-seq-len)",
    )
    args = p.parse_args()

    if (args.max_seq_len is None) != (args.oversized_dir is None):
        p.error("--max-seq-len and --oversized-dir must be used together")
    if not os.path.isfile(args.input_fasta):
        print(f"Error: input file not found: {args.input_fasta}", file=sys.stderr)
        sys.exit(1)
    if args.batch_size <= 0:
        print("Error: batch_size must be positive", file=sys.stderr)
        sys.exit(1)

    n_batches, total, n_oversized = split_fasta(
        Path(args.input_fasta),
        Path(args.output_dir),
        args.batch_size,
        args.max_seq_len,
        args.oversized_dir,
    )
    print(f"Split {total} sequences into {n_batches} batches ({n_oversized} oversized)")


if __name__ == "__main__":
    main()
