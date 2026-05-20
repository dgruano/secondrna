#!/usr/bin/env python3
"""
Merge per-batch rG4detector CSV results into a single output file.

Prediction mode writes a normal row-oriented CSV, while detection mode writes
2-line blocks per sequence. This script preserves the original FASTA order for
both formats.
"""

import csv
import os
import sys
from pathlib import Path


def read_manifest(manifest_file):
    """Read batch manifest and return batch_id -> ordered sequence IDs."""
    batch_ids = {}
    with open(manifest_file, "r") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                batch_ids[parts[0]] = parts[2].split(",") if parts[2] else []
    return batch_ids


def read_nonempty_rows(csv_file):
    with open(csv_file, newline="") as handle:
        return [row for row in csv.reader(handle) if row and any(cell.strip() for cell in row)]


def infer_batch_id(csv_file):
    return Path(csv_file).parent.name


def is_detection_format(rows):
    return len(rows) >= 2 and len(rows[1]) > 0 and rows[1][0].strip() == ""


def read_prediction_batch(csv_file, expected_ids):
    rows = read_nonempty_rows(csv_file)
    if not rows:
        return {}

    header = rows[0]
    data_rows = rows[1:]
    if expected_ids and len(data_rows) != len(expected_ids):
        raise ValueError(
            f"Prediction row count mismatch in {csv_file}: expected {len(expected_ids)}, got {len(data_rows)}"
        )

    records = {}
    for index, row in enumerate(data_rows):
        seq_id = row[0]
        if expected_ids and seq_id != expected_ids[index]:
            raise ValueError(
                f"Prediction order mismatch in {csv_file}: expected {expected_ids[index]}, got {seq_id}"
            )
        records[seq_id] = dict(zip(header, row))
    return records


def read_detection_batch(csv_file, expected_ids):
    rows = read_nonempty_rows(csv_file)
    if len(rows) % 2 != 0:
        raise ValueError(f"Malformed detection CSV (odd number of rows): {csv_file}")

    record_count = len(rows) // 2
    if expected_ids and record_count != len(expected_ids):
        raise ValueError(
            f"Detection record count mismatch in {csv_file}: expected {len(expected_ids)}, got {record_count}"
        )

    records = {}
    for index in range(record_count):
        header_row = rows[index * 2]
        score_row = rows[index * 2 + 1]
        seq_id = header_row[0]
        if expected_ids and seq_id != expected_ids[index]:
            raise ValueError(
                f"Detection order mismatch in {csv_file}: expected {expected_ids[index]}, got {seq_id}"
            )
        records[seq_id] = (header_row, score_row)
    return records


def read_batch_csvs(csv_files, manifest_batches):
    """Read all per-batch CSVs and return a format-aware record map."""
    data_dict = {}
    output_format = None

    for csv_file in csv_files:
        try:
            batch_id = infer_batch_id(csv_file)
            expected_ids = manifest_batches.get(batch_id)
            if expected_ids is None:
                raise ValueError(f"Batch ID {batch_id} is missing from the manifest")

            rows = read_nonempty_rows(csv_file)
            current_format = "detection" if is_detection_format(rows) else "prediction"

            if output_format is None:
                output_format = current_format
            elif output_format != current_format:
                raise ValueError(
                    f"Mixed rG4detector output formats are not supported: {output_format} and {current_format}"
                )

            if current_format == "detection":
                batch_records = read_detection_batch(csv_file, expected_ids)
            else:
                batch_records = read_prediction_batch(csv_file, expected_ids)

            overlap = set(data_dict).intersection(batch_records)
            if overlap:
                raise ValueError(
                    f"Duplicate sequence IDs encountered in batch results: {sorted(overlap)[:10]}"
                )

            data_dict.update(batch_records)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}", file=sys.stderr)
            raise

    return data_dict, output_format


def merge_results(manifest_file, csv_files, output_file):
    """Merge per-batch CSVs in original sequence order."""
    manifest_batches = read_manifest(manifest_file)
    original_seq_ids = []
    for batch_ids in manifest_batches.values():
        original_seq_ids.extend(batch_ids)
    print(f"Original sequence count: {len(original_seq_ids)}")

    data_dict, output_format = read_batch_csvs(csv_files, manifest_batches)
    print(f"Sequences in batch CSVs: {len(data_dict)}")

    missing = [sid for sid in original_seq_ids if sid not in data_dict]
    if missing:
        print(f"Error: {len(missing)} sequences missing from batch results", file=sys.stderr)
        print(f"Missing IDs: {missing[:10]}...", file=sys.stderr)
        sys.exit(1)

    extra = [sid for sid in data_dict if sid not in original_seq_ids]
    if extra:
        print(f"Warning: {len(extra)} unexpected sequences in batch results", file=sys.stderr)

    if output_format == "detection":
        with open(output_file, "w", newline="") as handle:
            writer = csv.writer(handle)
            for seq_id in original_seq_ids:
                header_row, score_row = data_dict[seq_id]
                writer.writerow(header_row)
                writer.writerow(score_row)
        print(f"Merged {len(original_seq_ids)} sequences to {output_file}")
        return

    output_rows = [data_dict[seq_id] for seq_id in original_seq_ids]
    if output_rows:
        header = list(output_rows[0].keys())
        with open(output_file, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header)
            writer.writeheader()
            writer.writerows(output_rows)
    print(f"Merged {len(output_rows)} sequences to {output_file}")


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
        merge_results(manifest_file, csv_files, output_file)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
