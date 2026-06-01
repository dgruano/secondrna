import csv
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths

# ---------------------------------------------------------------------------
# Low-level CSV helpers
# ---------------------------------------------------------------------------


def read_nonempty_rows(path):
    """Read CSV file and return only non-empty rows."""
    with open(path, newline="", encoding="utf-8") as f:
        return [
            row
            for row in csv.reader(f)
            if row and any(len(cell) > 0 and cell.strip() for cell in row)
        ]


def is_detection_format(rows):
    """Return True if rows represent block (detection) format — 2 rows per sequence."""
    return len(rows) >= 2 and len(rows[1]) > 0 and rows[1][0].strip() == ""


def peek_format(path):
    """Read first two non-empty lines of a CSV and return True if block (detection) format."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        lines = []
        for row in reader:
            if row and any(cell.strip() for cell in row):
                lines.append(row)
                if len(lines) == 2:
                    break
    return is_detection_format(lines)


# ---------------------------------------------------------------------------
# Detection CSV streamers
# ---------------------------------------------------------------------------


def stream_detection(path):
    """
    Stream block-format detection CSV one transcript at a time.

    File format (block / detection format):
        Two consecutive non-empty rows per transcript:
          Row 1 (header row):  transcript_id, nt_1, nt_2, ..., nt_N
          Row 2 (score row):   <empty>, score_1, score_2, ..., score_N
        Blank rows between transcripts are ignored.
        Positions are assigned as 1-based integers (1 … N).

        Example:
            ENST00000001,A,C,G,U,A
            ,0.12,0.87,0.94,0.03,0.21
            ENST00000002,G,G,A,C
            ,0.55,0.61,0.08,0.44

    NOTE: This is the output format of the rG4detector tool.

    Yields (transcript_id, scores_array, positions_array) without loading
    the full file into memory. Peak memory: O(max_transcript_length).
    """
    with open(path, newline="\n", encoding="utf-8") as handle:
        header_line = None
        for line in handle:
            line = line.strip()
            if not line:
                continue

            if header_line is None:
                header_line = line
                continue

            header_fields = header_line.split(",")

            transcript_id = header_fields[0].strip()

            scores = np.fromstring(line.strip(","), sep=",", dtype=np.float32)

            yield transcript_id, scores
            header_line = None


# ---------------------------------------------------------------------------
# Annotation loader
# ---------------------------------------------------------------------------


def load_annotation(path):
    """
    Load optional BED-like file: transcript_id, start, end, feature
    (tab-separated, no header). Returns None if path is None.
    """
    if path is None:
        return None
    df = pd.read_csv(
        path, sep="\t", header=None, names=["transcript_id", "start", "end", "feature"]
    )
    return df


# ---------------------------------------------------------------------------
# Processing and summarization
# ---------------------------------------------------------------------------
def count_peaks_from_file(file_path: str, th: float = 1.56) -> pd.DataFrame:
    """
    Finds and counts peaks above a given threshold for each transcript in the detection file.
    Returns a DataFrame with columns "transcript_id" and "peak_count".

    Args:
        file_path: Path to the detection CSV file.
        th: Threshold for peak detection (default: 1.56).

    Returns:
        pd.DataFrame: A DataFrame with columns "transcript_id" and "peak_count".
    """
    tx_ids = []
    peak_counts = []

    for transcript_id, scores in stream_detection(file_path):
        scores_array = np.array(scores)
        peaks, _ = find_peaks(scores_array, height=th)
        peak_count = len(peaks)
        tx_ids.append(transcript_id)
        peak_counts.append(peak_count)

    df = pd.DataFrame({"transcript_id": tx_ids, "peak_count": peak_counts})
    return df


def peak_file_to_df(detection_file: str, threshold: float) -> pd.DataFrame:
    """
    Reads a detection CSV file, finds peaks above the given threshold,
    and returns a DataFrame with one row per peak, including columns:
        - transcript_id
        - peak_position (0-based)
        - peak_strand
        - peak_score

    Args:
        - detection_file (str): Path to the detection CSV file.
        - threshold (float): Threshold for peak detection.

    Returns:
        pd.DataFrame: A DataFrame with columns "transcript_id", "peak_position", "peak_strand", "peak_score".
    """
    records = []
    for transcript_id, scores in stream_detection(detection_file):
        scores_array = np.array(scores)
        peaks, properties = find_peaks(scores_array, height=threshold)
        df = pd.DataFrame(
            {
                "transcript_id": transcript_id,
                "peak_position": peaks,
                "peak_score": properties["peak_heights"],
            }
        )
        records.append(df)
    df = pd.concat(records, ignore_index=True)
    # Peak strand is + if score >0, else -
    df["peak_strand"] = np.where(df["peak_score"] > 0, "+", "-")
    df["peak_score"] = df["peak_score"].abs()
    return df
