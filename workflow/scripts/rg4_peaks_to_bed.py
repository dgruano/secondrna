#!/usr/bin/env python3
"""
rG4 Peaks to BED — Export Per-Transcript Peaks as BED6

Streams detection.csv, calls peaks above a threshold, and writes a
BED6 file (chrom, start, end, name, score, strand) where each row is
one peak:
  - chrom  : transcript_id
  - start  : 0-based peak position
  - end    : start + 1
  - name   : transcript_id
  - score  : absolute peak score
  - strand : + (score > 0) or - (score <= 0)

Usage:
    python rg4_peaks_to_bed.py \\
        --input detection.csv \\
        --output peaks.bed \\
        [--threshold 1.56] \\
        [--log rg4_peaks_to_bed.log]

Snakemake script block usage (auto-detected via snakemake global):
    rule rg4_peaks_to_bed:
        input:
            detection = "results/{dataset}/detection.csv"
        output:
            bed = "results/{dataset}/peaks.bed"
        log:
            "logs/rg4_peaks_to_bed_{dataset}.log"
        params:
            threshold = 1.56
        script:
            "scripts/rg4_peaks_to_bed.py"
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent.parent))
from utils.rg4_parsing import peak_file_to_df

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
        description="Export rG4 peaks to BED6 format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", required=True, help="rG4detector detection.csv path"
    )
    parser.add_argument("--output", "-o", required=True, help="Output BED6 file path")
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.56,
        help="Peak height threshold (default: 1.56)",
    )
    parser.add_argument("--log", default=None, help="Log file path")
    parser.add_argument("--verbose", "-v", action="store_true")
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def peaks_to_bed6(detection_path, threshold, logger):
    """Stream detection.csv and return peaks as a BED6 DataFrame."""
    logger.info(f"Calling peaks (threshold={threshold}) ...")
    df = peak_file_to_df(detection_path, threshold)
    logger.info(
        f"Found {len(df):,} peaks across {df['transcript_id'].nunique():,} transcripts"
    )

    # TODO: Map transcripts to their original chromosomal coordinates if needed (requires additional metadata)
    df["end"] = df["peak_position"] + 1

    # Reorder to BED6: chrom, start, end, name, score, strand
    df = df[
        [
            "transcript_id",
            "peak_position",
            "end",
            "transcript_id",
            "peak_score",
            "peak_strand",
        ]
    ].copy()
    df.columns = ["chrom", "start", "end", "name", "score", "strand"]
    return df


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        args = argparse.Namespace()
        args.input = smk.input.detection
        args.output = smk.output.bed
        args.threshold = getattr(smk.params, "threshold", 1.56)
        args.log = smk.log[0] if smk.log else None
        args.verbose = False
    else:
        args = parse_args()

    logger = setup_logging(args.log, logging.DEBUG if args.verbose else logging.INFO)

    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Streaming {input_path} ...")
    df = peaks_to_bed6(str(input_path), args.threshold, logger)

    df.to_csv(args.output, sep="\t", index=False, header=False)
    logger.info(f"BED file written to {args.output} ({len(df):,} peaks)")


if __name__ == "__main__":
    main()
