# notebooks/g4_analysis_lib.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

# ─── Column schemas ───────────────────────────────────────────────────────────
INTERSECT_COLS = [
    "transcript_id",
    "start",
    "end",
    "pqsfinder_score",
    "motif_length",
    "strand",
    "g4HunterScore",
    "rg4d_chr",
    "rg4d_start",
    "rg4d_end",
    "rg4d_name",
    "rg4d_score",
    "rg4d_strand",
]
PEAKS_COLS = [
    "transcript_id",
    "rg4d_start",
    "rg4d_end",
    "rg4d_name",
    "rg4d_score",
    "rg4d_strand",
]
FILL_VALUES: dict = {
    "start": -1,
    "end": -1,
    "pqsfinder_score": -1.0,
    "motif_length": -1,
    "strand": ".",
    "g4HunterScore": -1.0,
    "rg4d_chr": ".",
    "rg4d_start": -1,
    "rg4d_end": -1,
    "rg4d_name": ".",
    "rg4d_score": -1.0,
    "rg4d_strand": ".",
    "G4D_name": ".",
}


# ─── Loaders ──────────────────────────────────────────────────────────────────
def load_g4discovery(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "transcript_id",
            "start",
            "end",
            "pqsfinder_score",
            "motif_length",
            "strand",
            "g4HunterScore",
        ],
        compression="gzip",
    )
    simple_id = df["transcript_id"].str.split("|").str[0]
    df["G4D_name"] = (
        simple_id + "_" + df["start"].astype(str) + "_" + df["end"].astype(str) + "_G4D"
    )
    return df


def load_intersect(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=INTERSECT_COLS)
    df["overlaps"] = df["rg4d_chr"] != "."
    simple_id = df["transcript_id"].str.split("|").str[0]
    df["G4D_name"] = (
        simple_id + "_" + df["start"].astype(str) + "_" + df["end"].astype(str) + "_G4D"
    )
    df["rg4d_name"] = (
        simple_id
        + "_"
        + df["rg4d_start"].astype(str)
        + "_"
        + df["rg4d_end"].astype(str)
        + "_rg4d"
    )
    df.loc[df["rg4d_chr"] == ".", "rg4d_name"] = "."
    return df


def load_peaks(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=PEAKS_COLS)
    simple_id = df["transcript_id"].str.split("|").str[0]
    df["rg4d_name"] = (
        simple_id
        + "_"
        + df["rg4d_start"].astype(str)
        + "_"
        + df["rg4d_end"].astype(str)
        + "_rg4d"
    )
    return df


def load_coding_annotation(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def load_length_annotation(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


# ─── Union builder ─────────────────────────────────────────────────────────────
def build_union_df(intersect: pd.DataFrame, peaks: pd.DataFrame) -> pd.DataFrame:
    missing = peaks[~peaks["rg4d_name"].isin(intersect["rg4d_name"])]
    union = pd.concat([intersect, missing], ignore_index=True)
    union = union.fillna(FILL_VALUES)
    union["overlaps"] = union["overlaps"].fillna(False).astype(bool)
    union["type"] = np.where(
        union["overlaps"],
        "both",
        np.where(union["rg4d_name"] == ".", "G4D only", "RG4D only"),
    )
    return union


def annotate_union(
    union: pd.DataFrame,
    coding_ann: pd.DataFrame,
    length_ann: pd.DataFrame,
) -> pd.DataFrame:
    union = union.copy()
    union["short_id"] = union["transcript_id"].str.split("|").str[0]
    union = union.join(coding_ann.set_index("transcript_id"), on="short_id", how="left")
    union = union.join(
        length_ann.set_index("transcript_id"), on="transcript_id", how="left"
    )
    return union


# ─── Subset definitions ────────────────────────────────────────────────────────
@dataclass
class SubsetSpec:
    name: str  # short slug for dict keys
    label: str  # human-readable label for plot titles


_PQS_THRESHOLDS = [30, 40]
_G4H_THRESHOLDS = [1.0, 1.2, 1.5]


def _threshold_filter(pqs_th: float, g4h_th: float) -> Callable:
    def _filter(df: pd.DataFrame) -> pd.DataFrame:
        return df[
            (df["type"] != "RG4D only")
            & (df["pqsfinder_score"] >= pqs_th)
            & (df["g4HunterScore"] >= g4h_th)
        ].copy()

    return _filter


def build_subsets(union_df: pd.DataFrame) -> list[tuple[SubsetSpec, pd.DataFrame]]:
    """Return list of (SubsetSpec, filtered_dataframe) for all 8 subsets."""
    result: list[tuple[SubsetSpec, pd.DataFrame]] = [
        (SubsetSpec("union", "Union (all motifs)"), union_df.copy()),
        (
            SubsetSpec("rg4d_only", "rG4D-only motifs"),
            union_df[union_df["type"] == "RG4D only"].copy(),
        ),
    ]
    for pqs_th in _PQS_THRESHOLDS:
        for g4h_th in _G4H_THRESHOLDS:
            name = f"pqs{pqs_th}_g4h{g4h_th}"
            label = f"G4D pqs≥{pqs_th} & G4H≥{g4h_th}"
            result.append(
                (SubsetSpec(name, label), _threshold_filter(pqs_th, g4h_th)(union_df))
            )
    return result
