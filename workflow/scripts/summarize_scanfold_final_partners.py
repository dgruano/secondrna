#!/usr/bin/env python3
"""
Compute per-transcript FinalPartners statistics from one or more directories.

Each directory is auto-detected as either a tarball collection (contains *.tar.gz)
or a flat/recursive tree of *.ScanFold.FinalPartners.txt files.  Results are
concatenated into a single TSV with a 'source_dir' column identifying the origin.

Usage (CLI):
    python summarize_scanfold_final_partners.py \\
        --dirs /mnt/cbib/LNClassifier/RNA_ScanFold2.0 \\
                results/gencode.v47.repeat.simple/scanfold2_gpu \\
        --output results/final_partners_stats.tsv

Snakemake usage:
    script: "../scripts/summarize_scanfold_final_partners.py"
    (expects smk.input as a list of directories, smk.output[0], smk.log[0],
     and optionally smk.params.n_bins)
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.utils.scanfold_lib import (
    summarize_final_partners_dir,
    summarize_final_partners_tarballs,
)


def setup_logging(log_file=None):
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        handlers=handlers,
    )
    return logging.getLogger(__name__)


def _has_tarballs(directory: Path) -> bool:
    """Return True if the directory contains any *.tar.gz files."""
    return any(directory.glob("*.tar.gz"))


def process_directory(
    directory: Path, n_bins: int, threads: int, logger: logging.Logger
) -> pd.DataFrame:
    """Summarize one directory, auto-detecting tarballs vs. flat files."""
    directory = Path(directory)
    if not directory.is_dir():
        raise NotADirectoryError(f"Not a directory: {directory}")

    if _has_tarballs(directory):
        logger.info(f"[tarball] {directory}")
        df = summarize_final_partners_tarballs(
            directory, n_bins=n_bins, threads=threads, logger=logger
        )
        # Drop 'tarball' provenance column — source_dir covers it at this level
        df = df.drop(columns="tarball", errors="ignore")
    else:
        logger.info(f"[files]   {directory}")
        df = summarize_final_partners_dir(directory, n_bins=n_bins, logger=logger)

    logger.info(f"          {len(df):,} transcripts")
    df.insert(0, "source_dir", str(directory))
    return df


def run(
    dirs: list[str],
    output: str,
    n_bins: int = 10,
    threads: int = 1,
    log_file: str | None = None,
) -> None:
    logger = setup_logging(log_file)

    frames = []
    for d in dirs:
        frames.append(process_directory(Path(d), n_bins, threads, logger))

    result = pd.concat(frames)
    logger.info(f"Total: {len(result):,} transcripts from {len(dirs)} director(y/ies)")

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.reset_index().to_csv(out_path, sep="\t", index=False)
    logger.info(f"Saved {out_path}  ({out_path.stat().st_size / 1e6:.1f} MB)")


def parse_args(args=None):
    p = argparse.ArgumentParser(
        description="Compute per-transcript FinalPartners stats from one or more directories.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--dirs",
        nargs="+",
        required=True,
        metavar="DIR",
        help="One or more directories (tarball collections or FinalPartners file trees).",
    )
    p.add_argument("--output", required=True, help="Output TSV path.")
    p.add_argument(
        "--n-bins",
        type=int,
        default=10,
        dest="n_bins",
        help="Number of positional bins (default: 10).",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Parallel workers for tarball processing.",
    )
    p.add_argument("--log", help="Log file path.")
    return p.parse_args(args)


def main():
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        args = argparse.Namespace(
            dirs=list(smk.input),
            output=smk.output[0],
            n_bins=getattr(smk.params, "n_bins", 10),
            threads=smk.threads,
            log=smk.log[0] if smk.log else None,
        )
    else:
        args = parse_args()

    run(
        args.dirs,
        args.output,
        n_bins=args.n_bins,
        threads=args.threads,
        log_file=args.log,
    )


if __name__ == "__main__":
    main()
