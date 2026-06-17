#!/usr/bin/env python3
"""
rG4 Peak Counter — Per-Transcript Peak Counts and Summary Statistics

Streams detection.csv, calls peaks above a threshold, and writes:
  1. A per-transcript counts CSV  (--output-counts)
  2. An aggregate statistics TSV  (--output-stats)

Optionally accepts a binary class-label file (--labels) to run
association tests between peak presence/count and class membership:
  - Chi-square test + Cramer's V + Odds Ratio  (has_peak vs real)
  - Mann-Whitney U + VDA (A12)                 (peak_count by real)

Usage:
    python rg4_count_peaks.py \\
        --input detection.csv \\
        --output-counts peak_counts.csv \\
        --output-stats  peak_stats.tsv \\
        [--threshold 1.56] \\
        [--labels binary_class_table.tsv] \\
        [--log rg4_count_peaks.log]

Snakemake script: block usage (auto-detected via snakemake global):
    rule rg4_count_peaks:
        input:
            detection = "results/{dataset}/rg4detector/detection.csv",
            labels    = "data/{dataset}/labels.tsv"   # optional
        output:
            counts = "results/{dataset}/rg4detector/peak_counts.csv",
            stats  = "results/{dataset}/rg4detector/peak_stats.tsv"
        log:
            "logs/{dataset}/rg4_count_peaks.log"
        benchmark:
            "benchmarks/{dataset}/rg4_count_peaks.tsv"
        params:
            threshold = 1.56
        script:
            "scripts/rg4_count_peaks.py"
"""

import argparse
import csv
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.rg4_parsing import stream_detection

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------


def setup_logging(log_file=None, level=logging.INFO):
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=level,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=handlers,
    )
    return logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Count rG4 peaks per transcript and compute summary statistics.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", required=True, help="rG4detector detection.csv path"
    )
    parser.add_argument(
        "--output-counts",
        required=True,
        dest="output_counts",
        help="Output CSV: per-transcript peak counts",
    )
    parser.add_argument(
        "--output-stats",
        required=True,
        dest="output_stats",
        help="Output TSV: aggregate and class-association statistics",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.56,
        help="Peak height threshold (default: 1.56)",
    )
    parser.add_argument(
        "--labels",
        default=None,
        help="Optional TSV with seq_ID and real columns for class-association tests",
    )
    parser.add_argument("--log", default=None, help="Log file path")
    parser.add_argument("--verbose", "-v", action="store_true")
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def count_peaks(detection_path, threshold, logger):
    """
    Stream detection.csv and return per-transcript counts DataFrame
    with additional per-transcript statistics.

    """
    records = []
    n = 0
    for transcript_id, scores in stream_detection(detection_path):
        peaks, _ = find_peaks(scores, height=threshold)
        records.append(
            {
                "transcript_id": transcript_id,
                "transcript_length": len(scores),
                "rg4_peak_count": len(peaks),
                "rg4_peak_density": (
                    len(peaks) / len(scores) if len(scores) > 0 else 0.0
                ),
                "rg4_peak_rel_position": (
                    (sum(p / len(scores) for p in peaks) / len(peaks))
                    if len(peaks) > 0
                    else -1.0
                ),
            }
        )
        n += 1
        if n % 10_000 == 0:
            logger.info(f"Processed {n:,} transcripts ...")
    logger.info(f"Finished: {n:,} transcripts total")
    df = pd.DataFrame(records)
    df["has_peak"] = df["rg4_peak_count"] > 0
    return df


def per_bin_stats(peaks, n_bins=10):
    """
    Per-transcript do:
        - Cut transcript in 10 equal-length bins
        - Compute peak count/density per bin
        - Return data
    Args:
        peaks: list of peak positions (0-based)
        n_bins: number of bins to divide transcript into
    Returns:
        bin_counts: list of counts per bin
        bin_densities: list of densities per bin

    NOTE: Preliminary function, not used for now.
    """
    if len(peaks) == 0:
        return [0] * n_bins, [0.0] * n_bins
    max_pos = max(peaks)
    bin_edges = np.linspace(0, max_pos + 1, n_bins + 1)
    bin_counts, _ = np.histogram(peaks, bins=bin_edges)
    bin_lengths = np.diff(bin_edges)
    bin_densities = bin_counts / bin_lengths
    return bin_counts.tolist(), bin_densities.tolist()


def global_stats(df):
    """Compute aggregate statistics from counts DataFrame."""
    rows = [
        ("n_transcripts", len(df)),
        ("n_with_peaks", int(df["has_peak"].sum())),
        ("frac_with_peaks", round(df["has_peak"].mean(), 6)),
        ("mean_peak_count", round(df["rg4_peak_count"].mean(), 4)),
        ("median_peak_count", round(df["rg4_peak_count"].median(), 4)),
        ("max_peak_count", int(df["rg4_peak_count"].max())),
        ("mean_transcript_length", round(df["transcript_length"].mean(), 1)),
    ]
    return rows


def class_association_stats(df, labels_path, logger):
    """Load labels and run chi-square + Mann-Whitney tests."""
    from scipy.stats import chi2_contingency, mannwhitneyu

    labels = pd.read_csv(labels_path, sep="\t")
    if "seq_ID" not in labels.columns or "real" not in labels.columns:
        logger.error("--labels file must have columns: seq_ID, real")
        sys.exit(1)

    # Strip version suffix from ENST IDs to allow merge
    merged = df.copy()
    merged.index = merged["transcript_id"].str.split("|").str[0]
    merged = merged.merge(
        labels[["seq_ID", "real"]].set_index("seq_ID"),
        left_index=True,
        right_index=True,
        how="inner",
    )
    logger.info(f"Merged {len(merged):,} transcripts with class labels")

    merged["real"] = merged["real"].astype(bool)

    # Chi-square: has_peak vs real
    contingency = pd.crosstab(merged["has_peak"], merged["real"]).reindex(
        index=[False, True], columns=[False, True], fill_value=0
    )
    chi2, p_chi2, dof, _ = chi2_contingency(contingency)
    n = contingency.to_numpy().sum()
    cramers_v = np.sqrt(
        chi2 / (n * min(contingency.shape[0] - 1, contingency.shape[1] - 1))
    )

    a = contingency.loc[True, True]
    b = contingency.loc[True, False]
    c = contingency.loc[False, True]
    d = contingency.loc[False, False]
    if min(a, b, c, d) == 0:
        odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        or_note = "_haldane_corrected"
    else:
        odds_ratio = (a * d) / (b * c)
        or_note = ""

    # Mann-Whitney U: peak_count by real
    peak_true = merged.loc[merged["real"], "rg4_peak_count"]
    peak_false = merged.loc[~merged["real"], "rg4_peak_count"]
    u_stat, p_mwu = mannwhitneyu(peak_true, peak_false, alternative="two-sided")
    vda = u_stat / (len(peak_true) * len(peak_false))

    rows = [
        ("n_merged", len(merged)),
        ("n_labels_real_true", int(merged["real"].sum())),
        ("n_labels_real_false", int((~merged["real"]).sum())),
        # Contingency table
        ("has_peak_T_real_T", int(a)),
        ("has_peak_T_real_F", int(b)),
        ("has_peak_F_real_T", int(c)),
        ("has_peak_F_real_F", int(d)),
        # Chi-square
        ("chi2", round(chi2, 4)),
        ("chi2_dof", dof),
        ("chi2_pvalue", f"{p_chi2:.3e}"),
        ("cramers_v", round(cramers_v, 6)),
        (f"odds_ratio{or_note}", round(odds_ratio, 4)),
        # Mann-Whitney
        ("mannwhitneyu_U", round(u_stat, 1)),
        ("mannwhitneyu_pvalue", f"{p_mwu:.3e}"),
        ("vda_a12", round(vda, 6)),
    ]
    return rows


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        args = argparse.Namespace()
        args.input = smk.input.detection
        args.output_counts = smk.output.counts
        args.output_stats = smk.output.stats
        args.threshold = getattr(smk.params, "threshold", 1.56)
        labels_val = getattr(smk.input, "labels", None)
        args.labels = labels_val if labels_val else None
        args.log = smk.log[0] if smk.log else None
        args.verbose = False
    else:
        args = parse_args()

    logger = setup_logging(args.log, logging.DEBUG if args.verbose else logging.INFO)

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)

    Path(args.output_counts).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_stats).parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Streaming {input_path} (threshold={args.threshold}) ...")
    df = count_peaks(str(input_path), args.threshold, logger)

    # Write per-transcript counts
    df.to_csv(args.output_counts, index=False)
    logger.info(f"Peak counts written to {args.output_counts}")

    # Collect statistics rows
    stat_rows = global_stats(df)
    if args.labels:
        logger.info(f"Loading class labels from {args.labels} ...")
        stat_rows += class_association_stats(df, args.labels, logger)

    with open(args.output_stats, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["statistic", "value"])
        writer.writerows(stat_rows)
    logger.info(f"Statistics written to {args.output_stats}")


if __name__ == "__main__":
    main()
