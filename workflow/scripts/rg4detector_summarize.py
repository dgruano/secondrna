#!/usr/bin/env python3
"""
rG4detector Detection Mode — Per-Transcript Summarization
Computes peak-based, distribution, and density metrics from
rG4detector_detection.csv output.

Usage:
    python rg4detector_summarize.py \
        --input rG4detector_detection.csv \
        --output rg4_summary.csv \
        [--threshold_low 1.0] \
        [--threshold_high 2.0] \
        [--min_peak_prominence 0.3] \
        [--min_peak_width 20] \
        [--annotation annotation.bed]  # optional, for regional breakdown
"""

import argparse
import csv
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths
import warnings
warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Summarize rG4detector detection output per transcript.")
    parser.add_argument("--input",  required=True, help="rG4detector_detection.csv path")
    parser.add_argument("--output", required=True, help="Output summary CSV path")
    parser.add_argument("--threshold_low",  type=float, default=1.0,
                        help="Score threshold for candidate rG4 peaks (default: 1.0)")
    parser.add_argument("--threshold_high", type=float, default=2.0,
                        help="Score threshold for high-confidence rG4 peaks (default: 2.0)")
    parser.add_argument("--min_peak_prominence", type=float, default=0.3,
                        help="Minimum peak prominence for peak calling (default: 0.3)")
    parser.add_argument("--min_peak_width", type=int, default=20,
                        help="Minimum peak width in nucleotides (default: 20)")
    parser.add_argument("--annotation", default=None,
                        help="Optional BED file with feature annotations (transcript_id, start, end, feature_name)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Core metrics per transcript
# ---------------------------------------------------------------------------
def summarize_transcript(transcript_id, scores, positions,
                         threshold_low, threshold_high,
                         min_prominence, min_width,
                         annotation_df=None):
    scores = np.array(scores, dtype=float)
    positions = np.array(positions, dtype=int)
    length = len(scores)

    # --- Global distribution metrics ---
    mean_score   = scores.mean()
    max_score    = scores.max()
    std_score    = scores.std()
    p90          = np.percentile(scores, 90)
    p95          = np.percentile(scores, 95)
    p99          = np.percentile(scores, 99)

    # --- Coverage metrics ---
    frac_above_low  = (scores > threshold_low).mean()
    frac_above_high = (scores > threshold_high).mean()

    # --- Peak calling (all peaks above baseline 0.5) ---
    peaks_idx, props = find_peaks(
        scores,
        height=threshold_low,
        prominence=min_prominence,
        width=min_width
    )

    n_peaks_low  = len(peaks_idx)
    n_peaks_high = int((scores[peaks_idx] > threshold_high).sum()) if n_peaks_low > 0 else 0

    # --- Per-peak metrics (top peak only + aggregates) ---
    if n_peaks_low > 0:
        peak_heights     = scores[peaks_idx]
        peak_prominences = props["prominences"]

        # Width at half-maximum (relative to base width)
        widths_half, _, _, _ = peak_widths(scores, peaks_idx, rel_height=0.5)

        top_peak_idx         = peaks_idx[np.argmax(peak_heights)]
        top_peak_score       = scores[top_peak_idx]
        top_peak_position    = int(positions[top_peak_idx])
        top_peak_frac        = top_peak_position / (positions[-1] if positions[-1] > 0 else 1)
        top_peak_width       = float(widths_half[np.argmax(peak_heights)])
        top_peak_prominence  = float(peak_prominences[np.argmax(peak_heights)])

        mean_peak_score      = float(peak_heights.mean())
        mean_peak_prominence = float(peak_prominences.mean())
        mean_peak_width      = float(widths_half.mean())

        # Peak area: sum of scores within each peak region
        peak_areas = []
        for i, (lo, hi) in enumerate(zip(
            props["left_bases"].astype(int),
            props["right_bases"].astype(int)
        )):
            peak_areas.append(scores[lo:hi+1].sum())
        top_peak_area  = float(peak_areas[np.argmax(peak_heights)])
        mean_peak_area = float(np.mean(peak_areas))

        # Inter-peak distance
        if n_peaks_low > 1:
            inter_peak_distances = np.diff(positions[peaks_idx])
            mean_inter_peak_dist = float(inter_peak_distances.mean())
            min_inter_peak_dist  = float(inter_peak_distances.min())
        else:
            mean_inter_peak_dist = np.nan
            min_inter_peak_dist  = np.nan

    else:
        top_peak_score = top_peak_position = top_peak_frac = np.nan
        top_peak_width = top_peak_prominence = top_peak_area = np.nan
        mean_peak_score = mean_peak_prominence = mean_peak_width = mean_peak_area = np.nan
        mean_inter_peak_dist = min_inter_peak_dist = np.nan

    # --- Density (peaks per kb) ---
    transcript_len_nt = int(positions[-1] - positions[0]) if length > 1 else length
    rg4_density_per_kb = (n_peaks_low / (transcript_len_nt / 1000)) if transcript_len_nt > 0 else np.nan

    # --- Optional regional breakdown (requires annotation BED) ---
    regional_counts = {}
    if annotation_df is not None and n_peaks_low > 0:
        tx_ann = annotation_df[annotation_df["transcript_id"] == transcript_id]
        for _, row in tx_ann.iterrows():
            feat = row["feature"]
            in_region = (positions[peaks_idx] >= row["start"]) & (positions[peaks_idx] < row["end"])
            regional_counts[f"n_peaks_{feat}"] = int(in_region.sum())

    result = {
        "transcript_id":          transcript_id,
        "transcript_length_nt":   transcript_len_nt,
        # Distribution
        "mean_score":             round(mean_score, 4),
        "max_score":              round(max_score, 4),
        "std_score":              round(std_score, 4),
        "p90_score":              round(p90, 4),
        "p95_score":              round(p95, 4),
        "p99_score":              round(p99, 4),
        # Coverage
        "frac_above_low":         round(frac_above_low, 4),
        "frac_above_high":        round(frac_above_high, 4),
        # Peak counts
        "n_peaks_above_low":      n_peaks_low,
        "n_peaks_above_high":     n_peaks_high,
        "rg4_density_per_kb":     round(rg4_density_per_kb, 4) if not np.isnan(rg4_density_per_kb) else np.nan,
        # Top peak
        "top_peak_score":         round(top_peak_score, 4) if not np.isnan(top_peak_score) else np.nan,
        "top_peak_position_nt":   top_peak_position,
        "top_peak_position_frac": round(top_peak_frac, 4) if not np.isnan(top_peak_frac) else np.nan,
        "top_peak_width_nt":      round(top_peak_width, 1) if not np.isnan(top_peak_width) else np.nan,
        "top_peak_prominence":    round(top_peak_prominence, 4) if not np.isnan(top_peak_prominence) else np.nan,
        "top_peak_area":          round(top_peak_area, 2) if not np.isnan(top_peak_area) else np.nan,
        # Aggregate peak stats
        "mean_peak_score":        round(mean_peak_score, 4) if not np.isnan(mean_peak_score) else np.nan,
        "mean_peak_prominence":   round(mean_peak_prominence, 4) if not np.isnan(mean_peak_prominence) else np.nan,
        "mean_peak_width_nt":     round(mean_peak_width, 1) if not np.isnan(mean_peak_width) else np.nan,
        "mean_peak_area":         round(mean_peak_area, 2) if not np.isnan(mean_peak_area) else np.nan,
        "mean_inter_peak_dist_nt":round(mean_inter_peak_dist, 1) if not np.isnan(mean_inter_peak_dist) else np.nan,
        "min_inter_peak_dist_nt": round(min_inter_peak_dist, 1) if not np.isnan(min_inter_peak_dist) else np.nan,
    }
    result.update(regional_counts)
    return result


# ---------------------------------------------------------------------------
# I/O helpers — streaming, no full-file load
# ---------------------------------------------------------------------------
def _peek_format(path):
    """Read only first 2 non-empty lines to detect block vs. table format."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        lines = []
        for row in reader:
            if row and any(cell.strip() for cell in row):
                lines.append(row)
                if len(lines) == 2:
                    break
    return len(lines) >= 2 and len(lines[1]) > 0 and lines[1][0].strip() == ""


def _stream_detection_blocks(path):
    """
    Stream block-format detection CSV one transcript at a time.

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


def _stream_detection_table(path, chunk_size=500_000):
    """
    Stream table-format detection CSV one transcript at a time using chunked reads.

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

    # Read one row to resolve column names
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
                # Transcript spans chunk boundary — accumulate
                pending_scores.append(scores)
                pending_pos.append(positions)
            else:
                if pending_id is not None:
                    all_scores = np.concatenate(pending_scores)
                    all_pos    = np.concatenate(pending_pos)
                    yield pending_id, all_scores, all_pos

                pending_id     = tx_id
                pending_scores = [scores]
                pending_pos    = [positions]

    # Flush last transcript
    if pending_id is not None:
        yield pending_id, np.concatenate(pending_scores), np.concatenate(pending_pos)


def stream_detection(path):
    """
    Yield (transcript_id, scores, positions) one transcript at a time.
    Automatically selects block or table streaming path.
    """
    if _peek_format(path):
        yield from _stream_detection_blocks(path)
    else:
        yield from _stream_detection_table(path)


def load_annotation(path):
    """
    Optional BED-like file: transcript_id, start, end, feature
    (tab-separated, no header)
    """
    if path is None:
        return None
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["transcript_id", "start", "end", "feature"])
    return df


# ---------------------------------------------------------------------------
# Main — streaming, writes output incrementally
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    ann = load_annotation(args.annotation)

    print(f"[rg4_summarize] Streaming {args.input} ...")

    writer = None
    n_processed = 0

    with open(args.output, "w", newline="") as out_fh:
        for tx_id, scores, positions in stream_detection(args.input):
            rec = summarize_transcript(
                transcript_id  = tx_id,
                scores         = scores,
                positions      = positions,
                threshold_low  = args.threshold_low,
                threshold_high = args.threshold_high,
                min_prominence = args.min_peak_prominence,
                min_width      = args.min_peak_width,
                annotation_df  = ann,
            )

            if writer is None:
                writer = csv.DictWriter(out_fh, fieldnames=list(rec.keys()))
                writer.writeheader()

            writer.writerow(rec)
            n_processed += 1

            if n_processed % 10_000 == 0:
                print(f"[rg4_summarize] Processed {n_processed:,} transcripts ...")

    print(f"[rg4_summarize] Done. {n_processed:,} transcripts written to {args.output}")


if __name__ == "__main__":
    main()
