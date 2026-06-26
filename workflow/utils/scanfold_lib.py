"""Utilities for loading and exploring ScanFold2 output files."""

from __future__ import annotations

import io
import re
import tarfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import IO, Iterable

import numpy as np
import pandas as pd

SCANFOLD_COLS = ["Start", "End", "Temperature", "NativeMFE", "Z-score", "p-value", "ED"]


def _read_scanfold(src: str | Path | IO[bytes], usecols: list[str]) -> pd.DataFrame:
    """Core reader shared by path-based and file-object-based loaders."""
    return pd.read_csv(src, sep="\t", index_col=0, usecols=list(usecols))


def load_scanfold(path: str | Path, usecols: list[str] | None = None) -> pd.DataFrame:
    """Load a ScanFold2 tab-separated output file.

    Parameters
    ----------
    path:
        Path to the .csv (tab-separated) file produced by ScanFold2.
    usecols:
        Columns to keep beyond the index. Defaults to SCANFOLD_COLS
        (Start, End, Temperature, NativeMFE, Z-score, p-value, ED).

    Returns
    -------
    DataFrame with the requested columns; row index matches ScanFold window index.
    """
    if usecols is None:
        usecols = SCANFOLD_COLS
    return _read_scanfold(path, usecols)


def load_scanfold_dir(
    directory: str | Path,
    pattern: str = "*.csv",
    usecols: list[str] | None = None,
) -> dict[str, pd.DataFrame]:
    """Load all ScanFold2 files in a directory.

    Returns a dict mapping transcript filename stem to its DataFrame.
    """
    directory = Path(directory)
    return {
        p.stem: load_scanfold(p, usecols=usecols)
        for p in sorted(directory.glob(pattern))
    }


def concat_scanfold(
    directory: str | Path,
    pattern: str = "*.csv",
    usecols: list[str] | None = None,
) -> pd.DataFrame:
    """Concatenate all ScanFold2 files in a directory into a single DataFrame.

    Adds a 'transcript_id' column from the filename stem.
    """
    frames = []
    for path in sorted(Path(directory).glob(pattern)):
        df = load_scanfold(path, usecols=usecols)
        df.insert(0, "transcript_id", path.stem)
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No files matching '{pattern}' in {directory}")
    return pd.concat(frames, ignore_index=True)


def concat_scanfold_tarballs(
    directory: str | Path,
    tarball_pattern: str = "*.tar.gz",
    csv_suffix: str = ".csv",
    usecols: list[str] | None = None,
) -> pd.DataFrame:
    """Read ScanFold2 CSV files from all tar.gz archives in a directory.

    Streams each CSV directly from the archive without extracting to disk.
    Adds 'tarball' (archive filename stem) and 'transcript_id' (CSV stem) columns.

    Parameters
    ----------
    directory:
        Directory to scan for .tar.gz files.
    tarball_pattern:
        Glob pattern for archive files.
    csv_suffix:
        Extension used to identify ScanFold CSV members inside each archive.
    usecols:
        Columns to keep. Defaults to SCANFOLD_COLS.

    Returns
    -------
    Concatenated DataFrame across all archives and all CSVs within them.
    """
    if usecols is None:
        usecols = SCANFOLD_COLS

    frames = []
    archives = sorted(Path(directory).glob(tarball_pattern))
    if not archives:
        raise FileNotFoundError(f"No files matching '{tarball_pattern}' in {directory}")

    for archive in archives:
        tarball_stem = archive.name.removesuffix(".tar.gz")
        with tarfile.open(archive, "r:gz") as tf:
            members = [
                m for m in tf.getmembers() if m.name.endswith(csv_suffix) and m.isfile()
            ]
            for member in members:
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                transcript_id = Path(member.name).stem
                df = _read_scanfold(fobj, usecols)
                df.insert(0, "tarball", tarball_stem)
                df.insert(1, "transcript_id", transcript_id)
                frames.append(df)

    if not frames:
        raise ValueError(
            f"No '{csv_suffix}' files found inside tarballs in {directory}"
        )
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Per-transcript analysis
# ---------------------------------------------------------------------------


def _count_peaks(mask: np.ndarray) -> int:
    """Count contiguous True runs in a boolean array."""
    if not mask.any():
        return 0
    # A new peak starts wherever we transition False→True
    padded = np.concatenate([[False], mask])
    return int(np.count_nonzero(np.diff(padded) == 1))


def compute_scanfold_stats(df: pd.DataFrame, n_bins: int = 10) -> pd.Series:
    """Per-transcript ScanFold2 statistics (vectorized).

    Parameters
    ----------
    df:
        Single-transcript slice from ``concat_scanfold`` / ``groupby``.
        Required columns: ``NativeMFE``, ``Z-score``, ``ED``, ``End``.
        Rows must be in ascending Start order (preserved by concat_scanfold).
    n_bins:
        Number of equal-width sequence bins for positional counts (default 10).

    Returns
    -------
    ``pd.Series`` of scalar statistics — one row per transcript when used
    inside ``summarize_transcripts``.

    Notes
    -----
    "Nucleotides with Z < threshold" counts window *positions* (one window
    per nucleotide in a step-1 sliding window).  "Peaks" counts contiguous
    runs of consecutive windows below the threshold.
    """
    z = df["Z-score"].to_numpy()
    mfe = df["NativeMFE"].to_numpy()
    ed = df["ED"].to_numpy()
    end = df["End"].to_numpy()

    transcript_len = int(end.max())
    n = len(z)

    mask1 = z < -1
    mask2 = z < -2

    # Assign each window to one of n_bins equal positional buckets
    bin_idx = np.floor(np.arange(n) * n_bins / n).astype(int).clip(0, n_bins - 1)

    stats: dict = {
        # Z-score
        "zscore_min": float(z.min()),
        "zscore_mean": float(z.mean()),
        "zscore_std": float(z.std(ddof=1)),
        # MFE
        "mfe_min": float(mfe.min()),
        "mfe_mean": float(mfe.mean()),
        "mfe_std": float(mfe.std(ddof=1)),
        # ED
        "ed_min": float(ed.min()),
        "ed_mean": float(ed.mean()),
        "ed_std": float(ed.std(ddof=1)),
        # Nucleotide counts
        "n_nt_z_lt_minus1": int(mask1.sum()),
        "n_nt_z_lt_minus2": int(mask2.sum()),
        # Peak counts
        "n_peaks_z_lt_minus1": _count_peaks(mask1),
        "n_peaks_z_lt_minus2": _count_peaks(mask2),
        # Relative position [0, 1] of global minimum (anchored to window End)
        "rel_pos_zscore_min": float(end[z.argmin()] / transcript_len),
        "rel_pos_mfe_min": float(end[mfe.argmin()] / transcript_len),
        "rel_pos_ed_min": float(end[ed.argmin()] / transcript_len),
    }

    # Per-bin counts
    for b in range(n_bins):
        sel = bin_idx == b
        stats[f"n_nt_z_lt_minus1_bin{b + 1}"] = int(mask1[sel].sum())
        stats[f"n_nt_z_lt_minus2_bin{b + 1}"] = int(mask2[sel].sum())
        stats[f"n_peaks_z_lt_minus1_bin{b + 1}"] = _count_peaks(mask1[sel])
        stats[f"n_peaks_z_lt_minus2_bin{b + 1}"] = _count_peaks(mask2[sel])

    return pd.Series(stats)


def summarize_transcripts(all_df: pd.DataFrame, n_bins: int = 10) -> pd.DataFrame:
    """Compute per-transcript statistics for a concatenated ScanFold2 DataFrame.

    Parameters
    ----------
    all_df:
        Output of ``concat_scanfold`` or ``concat_scanfold_tarballs``.
        Must have a ``transcript_id`` column; rows within each transcript
        must be in ascending Start order (guaranteed by both loaders).
        For large datasets, cast ``transcript_id`` to ``category`` dtype
        before calling to reduce memory usage.
    n_bins:
        Passed through to ``compute_scanfold_stats``.

    Returns
    -------
    DataFrame indexed by ``transcript_id``, one row per transcript,
    columns matching those returned by ``compute_scanfold_stats``.
    """
    import gc

    records: list[pd.Series] = []
    index: list[str] = []
    # Explicit iteration instead of groupby.apply() to avoid pandas' internal
    # bookkeeping overhead, which can double peak memory on large DataFrames.
    for transcript_id, group in all_df.groupby("transcript_id", sort=False):
        records.append(compute_scanfold_stats(group, n_bins=n_bins))
        index.append(transcript_id)

    gc.collect()
    return pd.DataFrame(records, index=pd.Index(index, name="transcript_id"))


# ---------------------------------------------------------------------------
# FinalPartners.txt — per-nucleotide consensus pairing output
# ---------------------------------------------------------------------------

# Suffix shared by all FinalPartners files (used for discovery and ID stripping)
_FP_SUFFIX = ".ScanFold.FinalPartners.txt"

# Regex matching data rows: "nt-42:" or "nt-42*:"
_FP_ROW_RE = re.compile(r"^nt-(\d+)(\*?):")


def _parse_final_partners_lines(lines: Iterable[str]) -> pd.DataFrame:
    """Parse an iterable of text lines from a FinalPartners.txt file.

    Accepts any iterable of strings — an open file handle, a list, or a
    ``io.TextIOWrapper`` wrapping a tarball member.

    Returns
    -------
    DataFrame with columns:
        pos          int   nucleotide position (1-based)
        bp_j         int   partner position; equals pos when unpaired
        avgMFE       float average MFE of windows supporting this pairing
        avgZ         float average Z-score of windows supporting this pairing
        avgED        float average Ensemble Diversity of those windows
        competition  bool  True when the row label carries a ``*`` flag
        paired       bool  True when pos != bp_j
    """
    rows = []
    for line in lines:
        m = _FP_ROW_RE.match(line)
        if m is None:
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 6:
            continue
        rows.append(
            (
                int(m.group(1)),  # pos
                int(parts[2]),  # bp_j
                float(parts[3]),  # avgMFE
                float(parts[4]),  # avgZ
                float(parts[5]),  # avgED
                m.group(2) == "*",  # competition flag
            )
        )
    df = pd.DataFrame(
        rows, columns=["pos", "bp_j", "avgMFE", "avgZ", "avgED", "competition"]
    )
    df["paired"] = df["pos"] != df["bp_j"]
    return df


def load_final_partners(path: str | Path) -> pd.DataFrame:
    """Load a single ScanFold2 FinalPartners.txt file.

    Parameters
    ----------
    path:
        Path to ``*.ScanFold.FinalPartners.txt``.

    Returns
    -------
    DataFrame — see ``_parse_final_partners_lines`` for column descriptions.
    """
    with open(path) as fh:
        return _parse_final_partners_lines(fh)


def compute_final_partners_stats(df: pd.DataFrame, n_bins: int = 10) -> pd.Series:
    """Compute per-transcript statistics from a FinalPartners DataFrame (vectorized).

    Parameters
    ----------
    df:
        Output of ``load_final_partners`` or ``_parse_final_partners_lines``.
        Required columns: pos, bp_j, avgMFE, avgZ, avgED, competition, paired.
    n_bins:
        Number of equal-width positional bins for per-bin counts (default 10).

    Returns
    -------
    ``pd.Series`` of scalar statistics:

    Summary statistics (min / mean / median / max):
        avgZ_{min,mean,median,max}, avgMFE_{...}, avgED_{...}

    Sequence-level counts and coverage (fraction of total nucleotides):
        length                total nucleotides
        n_paired              count of paired bases (bp_i != bp_j)
        frac_paired           n_paired / length
        n_competition         count of competition-flagged bases
        frac_competition      n_competition / length
        n_z_lt_minus1         count with avgZ < -1
        frac_z_lt_minus1      n_z_lt_minus1 / length
        n_z_lt_minus2         count with avgZ < -2
        frac_z_lt_minus2      n_z_lt_minus2 / length

    Relative position of global extremes [0 = 5' end, 1 = 3' end]:
        rel_pos_avgZ_min, rel_pos_avgMFE_min, rel_pos_avgED_min

    Per-bin fractions (suffix _bin1 … _bin{n_bins}):
        These could be either:
        - 1) fraction of nucleotides in that bin that meet the criterion
        - 2) fraction of nucleotides that meet the criterion in the total sequence that fall in that bin
        The second option answers better the question "where in the transcript do these features occur?",
        and is therefore closer to our question on feature distribution.
        We call them prop_*_bin{n} for proportion, to avoid confusion with the sequence-level fractions above.
        prop_paired_bin*, prop_competition_bin*, prop_z_lt_minus1_bin*, prop_z_lt_minus2_bin*
    """
    z = df["avgZ"].to_numpy()
    mfe = df["avgMFE"].to_numpy()
    ed = df["avgED"].to_numpy()
    paired = df["paired"].to_numpy()
    competition = df["competition"].to_numpy()

    n = len(df)
    if n == 0:
        return pd.Series({"length": 0}, dtype=float)

    mask_z1 = z < -1
    mask_z2 = z < -2

    # Relative position of each nucleotide: pos is 1-based, scale to [1/n … 1]
    rel_pos = df["pos"].to_numpy() / n

    stats: dict = {
        # Sequence length
        "length": n,
        # Z-score summary
        "avgZ_min": float(z.min()),
        "avgZ_mean": float(z.mean()),
        "avgZ_median": float(np.median(z)),
        "avgZ_max": float(z.max()),
        # MFE summary
        "avgMFE_min": float(mfe.min()),
        "avgMFE_mean": float(mfe.mean()),
        "avgMFE_median": float(np.median(mfe)),
        "avgMFE_max": float(mfe.max()),
        # ED summary
        "avgED_min": float(ed.min()),
        "avgED_mean": float(ed.mean()),
        "avgED_median": float(np.median(ed)),
        "avgED_max": float(ed.max()),
        # Paired bases
        "n_paired": int(paired.sum()),
        "frac_paired": float(paired.mean()),
        # Competition-flagged bases
        "n_competition": int(competition.sum()),
        "frac_competition": float(competition.mean()),
        # Bases below Z-score thresholds
        "n_z_lt_minus1": int(mask_z1.sum()),
        "frac_z_lt_minus1": float(mask_z1.mean()),
        "n_z_lt_minus2": int(mask_z2.sum()),
        "frac_z_lt_minus2": float(mask_z2.mean()),
        # Relative position [0, 1] of the global minimum for each metric
        "rel_pos_avgZ_min": float(rel_pos[z.argmin()]),
        "rel_pos_avgMFE_min": float(rel_pos[mfe.argmin()]),
        "rel_pos_avgED_min": float(rel_pos[ed.argmin()]),
    }

    # Per-bin counts: divide the sequence into n_bins equal positional slices
    bin_idx = np.floor(np.arange(n) * n_bins / n).astype(int).clip(0, n_bins - 1)
    for b in range(n_bins):
        sel = bin_idx == b
        stats[f"prop_paired_bin{b + 1}"] = (
            int(paired[sel].sum()) / stats["n_paired"] if stats["n_paired"] > 0 else 0
        )
        stats[f"prop_competition_bin{b + 1}"] = (
            int(competition[sel].sum()) / stats["n_competition"]
            if stats["n_competition"] > 0
            else 0
        )
        stats[f"prop_z_lt_minus1_bin{b + 1}"] = (
            int(mask_z1[sel].sum()) / stats["n_z_lt_minus1"]
            if stats["n_z_lt_minus1"] > 0
            else 0
        )
        stats[f"prop_z_lt_minus2_bin{b + 1}"] = (
            int(mask_z2[sel].sum()) / stats["n_z_lt_minus2"]
            if stats["n_z_lt_minus2"] > 0
            else 0
        )

    return pd.Series(stats)


def summarize_final_partners_dir(
    directory: str | Path,
    pattern: str = f"**/*{_FP_SUFFIX}",
    n_bins: int = 10,
    logger=None,
) -> pd.DataFrame:
    """Compute per-transcript FinalPartners stats for all files under a directory.

    Files are processed one at a time (streamed) so memory usage is bounded
    by the size of a single transcript, not the entire dataset.

    The default glob pattern (``**/*.ScanFold.FinalPartners.txt``) recurses into
    sub-directories, which covers both flat output directories and batch layouts
    such as ``results/<sample>/scanfold2_gpu/batch_*/``.

    Parameters
    ----------
    directory:
        Root directory to search.
    pattern:
        Glob pattern relative to *directory*. Override for non-recursive search
        (e.g. ``"*.ScanFold.FinalPartners.txt"``).
    n_bins:
        Passed to ``compute_final_partners_stats``.

    Returns
    -------
    DataFrame indexed by ``transcript_id`` (filename stem), one row per transcript.
    """
    directory = Path(directory)
    files = sorted(directory.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching '{pattern}' under {directory}")

    records: list[pd.Series] = []
    index: list[str] = []
    for path in files:
        transcript_id = path.name.removesuffix(_FP_SUFFIX)
        df = load_final_partners(path)
        if df.empty:
            if logger:
                logger.warning(f"[skip] empty FinalPartners: {path}")
            continue
        records.append(compute_final_partners_stats(df, n_bins=n_bins))
        index.append(transcript_id)

    return pd.DataFrame(records, index=pd.Index(index, name="transcript_id"))


def _process_archive(args: tuple) -> list[tuple[str, str, pd.Series]]:
    """Worker: extract stats from one .tar.gz; returns list of (transcript_id, tarball_stem, stats)."""
    archive_path, n_bins = args
    archive = Path(archive_path)
    tarball_stem = archive.name.removesuffix(".tar.gz")
    results = []
    with tarfile.open(archive, "r:gz") as tf:
        members = [
            m for m in tf.getmembers() if m.name.endswith(_FP_SUFFIX) and m.isfile()
        ]
        for member in members:
            fobj = tf.extractfile(member)
            if fobj is None:
                continue
            lines = io.TextIOWrapper(fobj, encoding="utf-8")
            df = _parse_final_partners_lines(lines)
            if df.empty:
                continue
            transcript_id = Path(member.name).name.removesuffix(_FP_SUFFIX)
            results.append(
                (
                    transcript_id,
                    tarball_stem,
                    compute_final_partners_stats(df, n_bins=n_bins),
                )
            )
    return results


def summarize_final_partners_tarballs(
    directory: str | Path,
    tarball_pattern: str = "*.tar.gz",
    n_bins: int = 10,
    threads: int = 1,
    logger=None,
) -> pd.DataFrame:
    """Compute per-transcript FinalPartners stats from all tarballs in a directory.

    Each FinalPartners file is streamed directly from its archive without
    extracting to disk, keeping memory proportional to one transcript at a time.

    Parameters
    ----------
    directory:
        Directory containing ``.tar.gz`` archives.
    tarball_pattern:
        Glob pattern for archive files.
    n_bins:
        Passed to ``compute_final_partners_stats``.

    Returns
    -------
    DataFrame indexed by ``transcript_id``, one row per transcript.
    Includes a ``tarball`` column with the archive filename stem (without
    ``.tar.gz``), useful when the same transcript_id appears in multiple archives.
    """
    archives = sorted(Path(directory).glob(tarball_pattern))
    if not archives:
        raise FileNotFoundError(f"No files matching '{tarball_pattern}' in {directory}")

    records: list[pd.Series] = []
    index: list[str] = []
    tarball_labels: list[str] = []

    args = [(str(a), n_bins) for a in archives]
    with ProcessPoolExecutor(max_workers=threads) as ex:
        for batch in ex.map(_process_archive, args):
            for transcript_id, tarball_stem, stats in batch:
                records.append(stats)
                index.append(transcript_id)
                tarball_labels.append(tarball_stem)

    if not records:
        raise ValueError(
            f"No '{_FP_SUFFIX}' files found inside tarballs in {directory}"
        )

    result = pd.DataFrame(records, index=pd.Index(index, name="transcript_id"))
    result.insert(0, "tarball", tarball_labels)
    return result
