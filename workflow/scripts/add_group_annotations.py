#!/usr/bin/env python3
"""
Add group annotations to rG4detector summary results.

Merges a summary CSV (e.g., from rg4detector_summarize.py) with external
group/class annotations (e.g., protein-coding vs lncRNA classifications).

Usage:
    python add_group_annotations.py \
        --summary rg4_summary.csv \
        --annotations group_labels.tsv \
        --output rg4_summary_annotated.csv \
        [--id-column transcript_id] \
        [--annotation-columns group,class,biotype] \
        [--separator \\t]

Example annotation file format (TSV):
    seq_ID      group       biotype
    TX001       protein     protein_coding
    TX002       lncRNA      lncRNA
    ...

The script will add all annotation columns to the summary while preserving
all original summary columns.
"""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add group annotations to rG4detector summary results."
    )
    parser.add_argument(
        "--summary",
        required=True,
        help="Summary CSV file (output from rg4detector_summarize.py)",
    )
    parser.add_argument(
        "--annotations",
        required=True,
        help="Annotation file (TSV, CSV, or similar) with group/class labels",
    )
    parser.add_argument(
        "--output", required=True, help="Output CSV file with annotations added"
    )
    parser.add_argument(
        "--id-column",
        default="transcript_id",
        help="Name of ID column in summary file (default: transcript_id)",
    )
    parser.add_argument(
        "--annotation-id-column",
        default=None,
        help="Name of ID column in annotation file (default: auto-detect first column)",
    )
    parser.add_argument(
        "--separator", default="\t", help="Delimiter for annotation file (default: tab)"
    )
    parser.add_argument(
        "--annotation-columns",
        default=None,
        help="Comma-separated list of annotation columns to add (default: all columns except ID)",
    )
    parser.add_argument(
        "--keep-unmatched",
        action="store_true",
        help="Keep summary rows without annotations (default: only keep matched rows)",
    )
    parser.add_argument(
        "--on-duplicate",
        default="first",
        choices=["first", "last", "error"],
        help="How to handle duplicate IDs in annotation file (default: first)",
    )
    return parser.parse_args()


def load_summary(path, id_column):
    """Load summary CSV and set ID as index."""
    df = pd.read_csv(path)

    if id_column not in df.columns:
        raise ValueError(
            f"ID column '{id_column}' not found in summary file. "
            f"Available columns: {list(df.columns)}"
        )

    df = df.set_index(id_column)
    return df


def load_annotations(path, id_column, separator, on_duplicate):
    """Load annotation file and set ID as index."""
    # Try to infer the separator if not explicitly provided
    if separator == "\t":
        df = pd.read_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep=separator)

    # Auto-detect ID column if not provided
    if id_column is None:
        id_column = df.columns[0]
        print(f"[add_annotations] Using first column '{id_column}' as ID column")

    if id_column not in df.columns:
        raise ValueError(
            f"ID column '{id_column}' not found in annotation file. "
            f"Available columns: {list(df.columns)}"
        )

    # Handle duplicates
    if df[id_column].duplicated().any():
        n_dups = df[id_column].duplicated().sum()
        print(f"[add_annotations] Found {n_dups} duplicate IDs in annotation file")

        if on_duplicate == "first":
            df = df.drop_duplicates(subset=[id_column], keep="first")
        elif on_duplicate == "last":
            df = df.drop_duplicates(subset=[id_column], keep="last")
        else:  # error
            raise ValueError(
                f"Duplicate IDs found in annotation file and on_duplicate='error'"
            )

    df = df.set_index(id_column)
    return df


def select_annotation_columns(annotation_df, requested_columns):
    """Select specific columns from annotation dataframe."""
    if requested_columns is None:
        # Use all columns except the index (ID column)
        return annotation_df

    cols = [c.strip() for c in requested_columns.split(",")]
    missing = set(cols) - set(annotation_df.columns)

    if missing:
        raise ValueError(
            f"Requested columns {missing} not found in annotation file. "
            f"Available: {list(annotation_df.columns)}"
        )

    return annotation_df[cols]


def merge_summary_with_annotations(summary_df, annotation_df, keep_unmatched):
    """Merge summary and annotations on index."""
    if keep_unmatched:
        merged = summary_df.join(annotation_df, how="left")
        n_matched = annotation_df.index.isin(summary_df.index).sum()
        n_unmatched = len(summary_df) - n_matched
        print(
            f"[add_annotations] Matched {n_matched:,} rows, "
            f"{n_unmatched:,} rows without annotations"
        )
    else:
        merged = summary_df.join(annotation_df, how="inner")
        n_matched = len(merged)
        n_unmatched = len(summary_df) - n_matched
        print(
            f"[add_annotations] Matched {n_matched:,} rows, "
            f"dropped {n_unmatched:,} rows without annotations"
        )

    return merged


def main():
    args = parse_args()

    print(f"[add_annotations] Loading summary from {args.summary} ...")
    summary_df = load_summary(args.summary, args.id_column)
    print(f"[add_annotations] Loaded {len(summary_df):,} summary rows")

    print(f"[add_annotations] Loading annotations from {args.annotations} ...")
    annotation_df = load_annotations(
        args.annotations, args.annotation_id_column, args.separator, args.on_duplicate
    )
    print(f"[add_annotations] Loaded {len(annotation_df):,} annotation rows")

    # Select specific annotation columns if requested
    annotation_df = select_annotation_columns(annotation_df, args.annotation_columns)
    print(f"[add_annotations] Using annotation columns: {list(annotation_df.columns)}")

    # Merge
    print(f"[add_annotations] Merging summary with annotations ...")
    merged_df = merge_summary_with_annotations(
        summary_df, annotation_df, args.keep_unmatched
    )

    # Write output
    print(f"[add_annotations] Writing output to {args.output} ...")
    merged_df.reset_index().to_csv(args.output, index=False)

    print(f"[add_annotations] Done. {len(merged_df):,} rows written")


if __name__ == "__main__":
    main()
