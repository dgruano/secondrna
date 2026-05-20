#!/usr/bin/env python3
"""
Split FASTA file into batches while preserving sequence order and IDs.
Output: batch FASTA files + manifest mapping sequence IDs to batches.
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO


def split_fasta(input_fasta, output_dir, batch_size):
    """
    Split FASTA into batches, create manifest with sequence IDs.

    Args:
        input_fasta: Path to input FASTA file
        output_dir: Directory to write batch files
        batch_size: Number of sequences per batch

    Returns:
        List of batch IDs created, total sequence count
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_id = 0
    batch_count = 0
    current_batch_fasta = []
    current_batch_ids = []
    total_sequences = 0
    manifest_lines = []

    # Read FASTA sequentially and group into batches
    for record in SeqIO.parse(input_fasta, "fasta"):
        total_sequences += 1
        current_batch_fasta.append(record)
        current_batch_ids.append(record.id)
        batch_count += 1

        # Write batch when full
        if batch_count == batch_size:
            batch_file = output_dir / f"batch_{batch_id:05d}.fa"
            SeqIO.write(current_batch_fasta, str(batch_file), "fasta")

            # Record in manifest: batch_id, seq_count, comma-separated seq IDs
            manifest_line = f"{batch_id:05d}\t{len(current_batch_ids)}\t{','.join(current_batch_ids)}"
            manifest_lines.append(manifest_line)

            batch_id += 1
            current_batch_fasta = []
            current_batch_ids = []
            batch_count = 0

    # Write final partial batch if any sequences remain
    if current_batch_fasta:
        batch_file = output_dir / f"batch_{batch_id:05d}.fa"
        SeqIO.write(current_batch_fasta, str(batch_file), "fasta")
        manifest_line = f"{batch_id:05d}\t{len(current_batch_ids)}\t{','.join(current_batch_ids)}"
        manifest_lines.append(manifest_line)
        batch_id += 1

    # Write manifest
    manifest_file = output_dir / "batch_manifest.txt"
    with open(manifest_file, "w") as f:
        f.write("# batch_id\tseq_count\tseq_ids\n")
        for line in manifest_lines:
            f.write(line + "\n")

    return batch_id, total_sequences


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input_fasta> <output_dir> <batch_size>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    batch_size = int(sys.argv[3])

    if not os.path.isfile(input_fasta):
        print(f"Error: Input file not found: {input_fasta}", file=sys.stderr)
        sys.exit(1)

    if batch_size <= 0:
        print(f"Error: batch_size must be positive, got {batch_size}", file=sys.stderr)
        sys.exit(1)

    try:
        num_batches, total_seqs = split_fasta(input_fasta, output_dir, batch_size)
        print(f"Successfully split {total_seqs} sequences into {num_batches} batches")
        print(f"Batch size: {batch_size}")
        print(f"Output: {output_dir}")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
