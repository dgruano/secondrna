#!/usr/bin/env python3
"""
Merge per-batch rG4detector CSV results using streaming to minimize memory.

Processes batches in original FASTA order from manifest without loading all data into memory.
Supports both prediction (1 row/seq) and detection (2 rows/seq) formats.
"""

import csv
import os
import sys
from pathlib import Path


def read_manifest(manifest_file):
    """Read batch manifest and return batch_id -> (batch_dir, seq_ids)."""
    batches = {}
    with open(manifest_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                batch_id = parts[0]
                batch_dir = parts[1]
                seq_ids = [s.strip() for s in parts[2].split(",") if s.strip()]
                batches[batch_id] = (batch_dir, seq_ids)
    return batches


def read_nonempty_rows(csv_file):
    """Read CSV file and filter out empty rows."""
    with open(csv_file, newline="") as f:
        rows = []
        for row in csv.reader(f):
            if row and any(len(cell) > 0 and cell.strip() for cell in row):
                rows.append(row)
        return rows


def is_detection_format(rows):
    """Check if rows represent detection format (2-line blocks per sequence)."""
    return len(rows) >= 2 and len(rows[1]) > 0 and rows[1][0].strip() == ""


def merge_detection(manifest_file, csv_files, output_file):
    """Merge detection format (2 rows per sequence, original order)."""
    batches = read_manifest(manifest_file)
    batch_map = {Path(f).parent.name: f for f in csv_files}

    seq_count = 0
    with open(output_file, "w", newline="") as out:
        writer = csv.writer(out)

        for batch_id, (_, seq_ids) in batches.items():
            csv_path = batch_map.get(batch_id)
            if not csv_path:
                raise ValueError(f"Batch {batch_id} not found in provided CSV files")

            rows = read_nonempty_rows(csv_path)
            if len(rows) % 2 != 0:
                raise ValueError(f"Malformed detection CSV (odd number of rows): {csv_path}")

            record_count = len(rows) // 2
            if record_count != len(seq_ids):
                raise ValueError(
                    f"Detection record count mismatch in {csv_path}: "
                    f"expected {len(seq_ids)}, got {record_count}"
                )

            for i in range(record_count):
                header_row = rows[i * 2]
                score_row = rows[i * 2 + 1]

                seq_id = header_row[0]
                if seq_id != seq_ids[i]:
                    raise ValueError(
                        f"Detection order mismatch in {csv_path}: "
                        f"expected {seq_ids[i]}, got {seq_id}"
                    )

                writer.writerow(header_row)
                writer.writerow(score_row)
                seq_count += 1

    print(f"Merged {seq_count} sequences to {output_file}")


def merge_prediction(manifest_file, csv_files, output_file):
    """Merge prediction format (1 row per sequence, original order)."""
    batches = read_manifest(manifest_file)
    batch_map = {Path(f).parent.name: f for f in csv_files}

    seq_count = 0
    writer = None
    with open(output_file, "w", newline="") as out:
        for batch_id, (_, seq_ids) in batches.items():
            csv_path = batch_map.get(batch_id)
            if not csv_path:
                raise ValueError(f"Batch {batch_id} not found in provided CSV files")

            rows = read_nonempty_rows(csv_path)
            if not rows:
                continue

            header = rows[0]
            data_rows = rows[1:]

            if len(data_rows) != len(seq_ids):
                raise ValueError(
                    f"Prediction row count mismatch in {csv_path}: "
                    f"expected {len(seq_ids)}, got {len(data_rows)}"
                )

            if writer is None:
                writer = csv.DictWriter(out, fieldnames=header)
                writer.writeheader()

            for idx, row in enumerate(data_rows):
                seq_id = row[0]
                if seq_id != seq_ids[idx]:
                    raise ValueError(
                        f"Prediction order mismatch in {csv_path}: "
                        f"expected {seq_ids[idx]}, got {seq_id}"
                    )

                record = dict(zip(header, row))
                writer.writerow(record)
                seq_count += 1

    print(f"Merged {seq_count} sequences to {output_file}")


def main():
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <manifest> <output> <csv_file1> [csv_file2 ...]")
        sys.exit(1)

    manifest_file = sys.argv[1]
    output_file = sys.argv[2]
    csv_files = sys.argv[3:]

    if not os.path.isfile(manifest_file):
        print(f"Error: Manifest file not found: {manifest_file}", file=sys.stderr)
        sys.exit(1)

    for csv_file in csv_files:
        if not os.path.isfile(csv_file):
            print(f"Error: CSV file not found: {csv_file}", file=sys.stderr)
            sys.exit(1)

    try:
        # Detect format from first batch CSV
        first_rows = read_nonempty_rows(csv_files[0])
        is_detection = is_detection_format(first_rows)

        if is_detection:
            merge_detection(manifest_file, csv_files, output_file)
        else:
            merge_prediction(manifest_file, csv_files, output_file)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
