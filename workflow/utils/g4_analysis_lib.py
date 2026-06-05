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
        g4d_mask = (df["pqsfinder_score"] >= pqs_th) & (df["g4HunterScore"] >= g4h_th)
        rg4d_only_mask = df["type"] == "RG4D only"
        return df[g4d_mask | rg4d_only_mask].copy()

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


# ─── Analysis functions ────────────────────────────────────────────────────────
def compute_aggregate_stats(df: pd.DataFrame) -> dict:
    """Return counts and overlap fractions for the motifs in df."""
    n_both = int((df["type"] == "both").sum())
    n_g4d_only = int((df["type"] == "G4D only").sum())
    n_rg4d_only = int((df["type"] == "RG4D only").sum())
    n_g4d = n_both + n_g4d_only
    n_rg4d = n_both + n_rg4d_only
    return {
        "total": len(df),
        "both": n_both,
        "g4d_only": n_g4d_only,
        "rg4d_only": n_rg4d_only,
        "g4d_total": n_g4d,
        "rg4d_total": n_rg4d,
        "g4d_overlap_frac": n_both / n_g4d if n_g4d > 0 else 0.0,
        "rg4d_overlap_frac": n_both / n_rg4d if n_rg4d > 0 else 0.0,
    }


def compute_coding_breakdown(df: pd.DataFrame) -> pd.DataFrame:
    """Class-level summary: total G4D motifs and rG4D-overlap fraction by coding class."""
    w_class = df.dropna(subset=["real"])
    if w_class.empty:
        return pd.DataFrame(
            columns=["Total motifs", "With rG4D overlap", "Overlap (%)"]
        )
    return (
        w_class.groupby("real")
        .agg(total_motifs=("overlaps", "count"), with_overlap=("overlaps", "sum"))
        .assign(overlap_frac=lambda d: d["with_overlap"] / d["total_motifs"] * 100)
        .rename(index={True: "Protein-coding", False: "Non-coding"})
        .rename(
            columns={
                "total_motifs": "Total motifs",
                "with_overlap": "With rG4D overlap",
                "overlap_frac": "Overlap (%)",
            }
        )
    )


def _vda(u_stat: float, n1: int, n2: int) -> float:
    return u_stat / (n1 * n2)


def _cramers_v(chi2: float, n: int, dof: int) -> float:
    return float(np.sqrt(chi2 / (n * dof)))


def run_statistical_tests(
    df: pd.DataFrame,
    coding_ann: pd.DataFrame,
    length_ann: pd.DataFrame,
) -> dict:
    """Chi-square (motif presence/absence) + Mann-Whitney U (motifs per kb).

    df must already contain only the motifs to be tested (no internal subsetting).
    """
    stats_df = coding_ann.copy()
    stats_df["has_g4d"] = stats_df["transcript_id"].isin(df["short_id"])
    stats_df["g4d_count"] = (
        stats_df["transcript_id"]
        .map(df.groupby("short_id").size())
        .fillna(0)
        .astype(int)
    )

    length_short = length_ann.copy()
    length_short["short_id"] = length_short["transcript_id"].str.split("|").str[0]
    stats_df = stats_df.merge(
        length_short[["short_id", "tx_len_nt"]],
        left_on="transcript_id",
        right_on="short_id",
        how="left",
    )
    stats_df["g4d_per_kb"] = stats_df["g4d_count"] / (stats_df["tx_len_nt"] / 1000)

    # Chi-square: motif presence/absence by class
    presence_table = pd.crosstab(stats_df["real"], stats_df["has_g4d"])
    if presence_table.shape == (2, 2):
        chi2_val, p_chi2, dof, _ = scipy_stats.chi2_contingency(presence_table.values)
        cramers = _cramers_v(chi2_val, len(stats_df), dof)
    else:
        chi2_val = p_chi2 = cramers = np.nan

    # Mann-Whitney U: motif density by class
    coding_kb = stats_df.loc[stats_df["real"] == True, "g4d_per_kb"].dropna()
    lncrna_kb = stats_df.loc[stats_df["real"] == False, "g4d_per_kb"].dropna()
    if len(coding_kb) > 0 and len(lncrna_kb) > 0:
        u_stat, p_mwu = scipy_stats.mannwhitneyu(
            coding_kb, lncrna_kb, alternative="two-sided"
        )
        vda_val = _vda(u_stat, len(coding_kb), len(lncrna_kb))
    else:
        u_stat = p_mwu = vda_val = np.nan

    return {
        "chi2": chi2_val,
        "p_chi2": p_chi2,
        "cramers_v": cramers,
        "u_stat": u_stat,
        "p_mwu": p_mwu,
        "vda": vda_val,
        "median_coding_per_kb": (
            float(coding_kb.median()) if len(coding_kb) > 0 else np.nan
        ),
        "median_lncrna_per_kb": (
            float(lncrna_kb.median()) if len(lncrna_kb) > 0 else np.nan
        ),
        "n_coding": len(coding_kb),
        "n_lncrna": len(lncrna_kb),
    }


def compare_results(all_results: dict[str, dict]) -> pd.DataFrame:
    """Flatten all_results into one row per subset for cross-subset comparison."""
    rows = []
    for name, res in all_results.items():
        spec = res["spec"]
        agg = res["agg_stats"]
        tst = res["test_results"]
        rows.append(
            {
                "subset": spec.label,
                "total_motifs": agg["total"],
                "g4d_total": agg["g4d_total"],
                "rg4d_total": agg["rg4d_total"],
                "both": agg["both"],
                "g4d_overlap_frac": agg["g4d_overlap_frac"],
                "chi2": tst["chi2"],
                "p_chi2": tst["p_chi2"],
                "cramers_v": tst["cramers_v"],
                "u_stat": tst["u_stat"],
                "p_mwu": tst["p_mwu"],
                "vda": tst["vda"],
                "median_coding_per_kb": tst["median_coding_per_kb"],
                "median_lncrna_per_kb": tst["median_lncrna_per_kb"],
            }
        )
    return pd.DataFrame(rows).set_index("subset")
