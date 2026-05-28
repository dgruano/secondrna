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
            row for row in csv.reader(f)
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


def stream_detection_blocks(path):
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
    seq_num = 0
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header_row = None
        for row in reader:
            if not row or not any(cell.strip() for cell in row):
                continue
            if header_row is None:
                header_row = row
            else:
                score_row = row
                seq_num += 1

                transcript_id = header_row[0].strip()
                if not transcript_id:
                    transcript_id = f"sequence_{seq_num}"

                scores = np.array(
                    [float(v) for v in (c.strip() for c in score_row[1:]) if v],
                    dtype=np.float32,
                )
                nucleotides = [c.strip() for c in header_row[1:] if c.strip()]
                length = min(len(nucleotides), len(scores)) if nucleotides else len(scores)
                positions = np.arange(1, length + 1, dtype=np.int32)
                scores = scores[:length]

                yield transcript_id, scores, positions
                header_row = None



def stream_detection_table(path, chunk_size=500_000):
    """
    Stream table-format detection CSV one transcript at a time using chunked reads.

    File format (table format):
        Standard CSV with a header row followed by one row per nucleotide position.
        Required columns (case-insensitive, stripped; first match wins):
          ID column:    sequence_id | transcript_id | name | id | seq_id
          Position col: position | pos | nucleotide | nt
          Score col:    score | prediction | rg4_score | value
        Rows for the same transcript must be grouped (or will be reassembled
        across chunk boundaries). Position values are used as-is (int32).

        Example:
            transcript_id,position,score
            ENST00000001,1,0.12
            ENST00000001,2,0.87
            ENST00000002,1,0.55

    Yields (transcript_id, scores_array, positions_array).
    Peak memory: O(chunk_size + max_transcript_length).
    """
    id_candidates    = ["sequence_id", "transcript_id", "name", "id", "seq_id"]
    pos_candidates   = ["position", "pos", "nucleotide", "nt"]
    score_candidates = ["score", "prediction", "rg4_score", "value"]

    def resolve(candidates, cols):
        for c in candidates:
            if c in cols:
                return c
        raise ValueError(f"Cannot find column — expected one of {candidates}, got {list(cols)}")

    header_df = pd.read_csv(path, nrows=0)
    header_df.columns = header_df.columns.str.strip().str.lower()
    id_col    = resolve(id_candidates,    header_df.columns)
    pos_col   = resolve(pos_candidates,   header_df.columns)
    score_col = resolve(score_candidates, header_df.columns)

    pending_id     = None
    pending_scores = []
    pending_pos    = []

    for chunk in pd.read_csv(path, chunksize=chunk_size,
                              usecols=[id_col, pos_col, score_col],
                              dtype={score_col: np.float32, pos_col: np.int32}):
        chunk.columns = chunk.columns.str.strip().str.lower()
        chunk = chunk.rename(columns={id_col: "transcript_id",
                                      pos_col: "position",
                                      score_col: "score"})

        for tx_id, group in chunk.groupby("transcript_id", sort=False):
            group = group.sort_values("position")
            scores = group["score"].values
            positions = group["position"].values

            if tx_id == pending_id:
                pending_scores.append(scores)
                pending_pos.append(positions)
            else:
                if pending_id is not None:
                    yield pending_id, np.concatenate(pending_scores), np.concatenate(pending_pos)

                pending_id     = tx_id
                pending_scores = [scores]
                pending_pos    = [positions]

    if pending_id is not None:
        yield pending_id, np.concatenate(pending_scores), np.concatenate(pending_pos)


def stream_detection(path):
    """
    Yield (transcript_id, scores, positions) one transcript at a time.
    Automatically selects block or table streaming path.
    """
    if peek_format(path):
        yield from stream_detection_blocks(path)
    else:
        yield from stream_detection_table(path)


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
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["transcript_id", "start", "end", "feature"])
    return df
