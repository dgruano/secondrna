"""
rG4detector visualization library.

Provides reusable functions for generating publication-quality static plots
and interactive visualizations from rG4detector summary CSV output.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path
from typing import Optional, List, Dict, Tuple, Union


class RG4SummaryLoader:
    """Load and validate rG4detector summary CSV with flexible column handling."""

    EXPECTED_METRICS = {
        "transcript_id",
        "mean_score",
        "max_score",
        "p95_score",
        "top_peak_score",
        "top_peak_position",
        "top_peak_width",
        "top_peak_prominence",
        "n_peaks_above_low",
        "n_peaks_above_high",
        "rg4_density_per_kb",
        "frac_above_low",
        "frac_above_high",
    }

    @staticmethod
    def load(filepath: Union[str, Path]) -> pd.DataFrame:
        """Load CSV and validate required columns."""
        df = pd.read_csv(filepath)

        # Normalize column names: strip whitespace, lowercase
        df.columns = df.columns.str.strip().str.lower()

        # Check for required columns
        missing = RG4SummaryLoader.EXPECTED_METRICS - set(df.columns)
        if missing:
            print(f"Warning: Missing expected columns: {missing}")

        return df


def plot_score_distribution(
    df: pd.DataFrame,
    metrics: Optional[List[str]] = None,
    format: str = "both",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure]:
    """
    Distribution of rG4 scores across transcripts.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary DataFrame
    metrics : List[str], optional
        Metrics to plot (default: ['mean_score', 'max_score', 'p95_score'])
    format : str
        'static', 'interactive', or 'both'
    output_path : Path, optional
        Save static plot to this path (for PNG)

    Returns
    -------
    Figure or go.Figure
        matplotlib Figure (static) or plotly Figure (interactive)
    """
    if metrics is None:
        metrics = ["mean_score", "max_score", "p95_score"]

    if format in ("static", "both"):
        fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
        if len(metrics) == 1:
            axes = [axes]

        for ax, metric in zip(axes, metrics):
            sns.histplot(
                data=df, x=metric, kde=True, ax=ax, bins=30, color="steelblue"
            )
            ax.set_title(f"Distribution of {metric}")
            ax.set_xlabel(metric)
            ax.set_ylabel("Count")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        fig = go.Figure()
        for metric in metrics:
            fig.add_trace(
                go.Histogram(
                    x=df[metric],
                    name=metric,
                    nbinsx=30,
                    opacity=0.7,
                    showlegend=True,
                )
            )

        fig.update_layout(
            title="Distribution of rG4 Scores",
            xaxis_title="Score",
            yaxis_title="Count",
            hovermode="x unified",
            template="plotly_white",
        )
        if format == "interactive":
            return fig

    return fig


def plot_peak_counts_distribution(
    df: pd.DataFrame,
    format: str = "both",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure]:
    """
    Distribution of peak counts and density across transcripts.

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary DataFrame
    format : str
        'static', 'interactive', or 'both'
    output_path : Path, optional
        Save static plot to this path

    Returns
    -------
    Figure or go.Figure
    """
    metrics = ["n_peaks_above_low", "n_peaks_above_high", "rg4_density_per_kb"]

    if format in ("static", "both"):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        for ax, metric in zip(axes, metrics):
            sns.histplot(data=df, x=metric, kde=True, ax=ax, bins=30, color="coral")
            ax.set_title(f"Distribution of {metric}")
            ax.set_xlabel(metric)
            ax.set_ylabel("Count")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        fig = go.Figure()
        for metric in metrics:
            fig.add_trace(
                go.Histogram(
                    x=df[metric], name=metric, nbinsx=30, opacity=0.7, showlegend=True
                )
            )

        fig.update_layout(
            title="Distribution of Peak Counts and Density",
            xaxis_title="Value",
            yaxis_title="Count",
            hovermode="x unified",
            template="plotly_white",
        )
        if format == "interactive":
            return fig

    return fig


def plot_coverage_distribution(
    df: pd.DataFrame,
    format: str = "both",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure]:
    """
    Distribution of coverage (fraction above thresholds).

    Parameters
    ----------
    df : pd.DataFrame
        Loaded summary DataFrame
    format : str
        'static', 'interactive', or 'both'
    output_path : Path, optional
        Save static plot to this path

    Returns
    -------
    Figure or go.Figure
    """
    metrics = ["frac_above_low", "frac_above_high"]

    if format in ("static", "both"):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))

        for ax, metric in zip(axes, metrics):
            sns.violinplot(data=df, y=metric, ax=ax, color="lightgreen")
            ax.set_title(f"Distribution of {metric}")
            ax.set_ylabel(metric)

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        fig = go.Figure()
        for metric in metrics:
            fig.add_trace(
                go.Box(
                    y=df[metric], name=metric, boxmean="sd", showlegend=True
                )
            )

        fig.update_layout(
            title="Distribution of Coverage Metrics",
            yaxis_title="Fraction",
            hovermode="y unified",
            template="plotly_white",
        )
        if format == "interactive":
            return fig

    return fig


def apply_threshold_filter(
    df: pd.DataFrame, filters: Dict[str, Tuple[str, Optional[float]]]
) -> pd.DataFrame:
    """
    Filter transcripts by threshold criteria.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame
    filters : Dict[str, Tuple[str, Optional[float]]]
        Filter dict: {column: (operator, value)}
        Operators: '>', '<', '>=', '<='

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame
    """
    filtered = df.copy()

    for column, (op_str, threshold) in filters.items():
        if column not in filtered.columns:
            print(f"Warning: Column {column} not found in DataFrame")
            continue

        if op_str == ">":
            filtered = filtered[filtered[column] > threshold]
        elif op_str == "<":
            filtered = filtered[filtered[column] < threshold]
        elif op_str == ">=":
            filtered = filtered[filtered[column] >= threshold]
        elif op_str == "<=":
            filtered = filtered[filtered[column] <= threshold]

    return filtered


def plot_top_transcripts(
    df: pd.DataFrame,
    metric: str = "top_peak_score",
    n: int = 20,
    format: str = "static",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure]:
    """
    Bar chart of top N transcripts by selected metric.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame
    metric : str
        Metric to rank by (default: 'top_peak_score')
    n : int
        Number of top transcripts to show (default: 20)
    format : str
        'static' or 'interactive'
    output_path : Path, optional
        Save plot to this path

    Returns
    -------
    Figure or go.Figure
    """
    top_df = df.nlargest(n, metric)[["transcript_id", metric]].reset_index(drop=True)

    if format == "static":
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.barh(top_df["transcript_id"], top_df[metric], color="steelblue")
        ax.set_xlabel(metric)
        ax.set_title(f"Top {n} Transcripts by {metric}")
        ax.invert_yaxis()
        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        return fig

    else:  # interactive
        fig = go.Figure(
            data=[
                go.Bar(
                    x=top_df[metric],
                    y=top_df["transcript_id"],
                    orientation="h",
                    marker=dict(color="steelblue"),
                )
            ]
        )
        fig.update_layout(
            title=f"Top {n} Transcripts by {metric}",
            xaxis_title=metric,
            yaxis_title="Transcript ID",
            height=max(400, n * 15),
            template="plotly_white",
        )
        return fig


def create_ranking_table(
    df: pd.DataFrame,
    top_n: int = 50,
    key_metrics: Optional[List[str]] = None,
) -> go.Figure:
    """
    Create interactive sortable/filterable table of top candidates.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame
    top_n : int
        Number of top transcripts (default: 50)
    key_metrics : List[str], optional
        Metrics to show in table (default: key ranking metrics)

    Returns
    -------
    go.Figure
        Plotly table figure
    """
    if key_metrics is None:
        key_metrics = [
            "transcript_id",
            "top_peak_score",
            "n_peaks_above_low",
            "rg4_density_per_kb",
            "top_peak_position",
        ]

    # Ensure transcript_id is first
    if "transcript_id" not in key_metrics:
        key_metrics = ["transcript_id"] + key_metrics

    top_df = df.nlargest(top_n, "top_peak_score")[key_metrics].reset_index(drop=True)

    fig = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=list(top_df.columns),
                    fill_color="paleturquoise",
                    align="left",
                    font=dict(color="black"),
                ),
                cells=dict(
                    values=[top_df[col] for col in top_df.columns],
                    fill_color="lavender",
                    align="left",
                    font=dict(color="black"),
                ),
            )
        ]
    )

    fig.update_layout(
        title=f"Top {top_n} rG4 Candidates",
        height=max(400, top_n * 10),
        template="plotly_white",
    )

    return fig


def plot_transcript_profile(
    df: pd.DataFrame,
    transcript_id: str,
    format: str = "both",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure]:
    """
    2x2 panel showing all key metrics for a single transcript.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame
    transcript_id : str
        Transcript ID to plot
    format : str
        'static', 'interactive', or 'both'
    output_path : Path, optional
        Save static plot to this path

    Returns
    -------
    Figure or go.Figure
    """
    row = df[df["transcript_id"] == transcript_id]
    if row.empty:
        print(f"Transcript {transcript_id} not found in data")
        return None

    row = row.iloc[0]

    if format in ("static", "both"):
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle(f"Profile for {transcript_id}", fontsize=14, fontweight="bold")

        # Panel 1: Peak counts
        ax = axes[0, 0]
        peaks_data = [row["n_peaks_above_low"], row["n_peaks_above_high"]]
        ax.bar(["Low threshold", "High threshold"], peaks_data, color=["coral", "red"])
        ax.set_ylabel("Number of peaks")
        ax.set_title("Peak Counts")

        # Panel 2: Top peak properties
        ax = axes[0, 1]
        peak_props = {
            "Score": row["top_peak_score"],
            "Width": row["top_peak_width"] if "top_peak_width" in row else 0,
            "Prominence": row["top_peak_prominence"] if "top_peak_prominence" in row else 0,
        }
        ax.bar(peak_props.keys(), peak_props.values(), color="steelblue")
        ax.set_ylabel("Value")
        ax.set_title("Top Peak Properties")

        # Panel 3: Coverage metrics
        ax = axes[1, 0]
        coverage_data = [row["frac_above_low"], row["frac_above_high"]]
        ax.bar(["Low threshold", "High threshold"], coverage_data, color="lightgreen")
        ax.set_ylabel("Fraction")
        ax.set_ylim([0, 1])
        ax.set_title("Coverage Metrics")

        # Panel 4: Score summary
        ax = axes[1, 1]
        scores = {
            "Mean": row["mean_score"],
            "Max": row["max_score"],
            "P95": row["p95_score"] if "p95_score" in row else 0,
        }
        ax.bar(scores.keys(), scores.values(), color="purple", alpha=0.7)
        ax.set_ylabel("Score")
        ax.set_title("Score Summary")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        fig = go.Figure()

        # Create subplots
        from plotly.subplots import make_subplots

        fig = make_subplots(
            rows=2,
            cols=2,
            subplot_titles=(
                "Peak Counts",
                "Top Peak Properties",
                "Coverage Metrics",
                "Score Summary",
            ),
        )

        # Panel 1: Peak counts
        fig.add_trace(
            go.Bar(
                x=["Low threshold", "High threshold"],
                y=[row["n_peaks_above_low"], row["n_peaks_above_high"]],
                name="Peaks",
                marker_color="coral",
            ),
            row=1,
            col=1,
        )

        # Panel 2: Top peak properties
        peak_props = {
            "Score": row["top_peak_score"],
            "Width": row.get("top_peak_width", 0),
            "Prominence": row.get("top_peak_prominence", 0),
        }
        fig.add_trace(
            go.Bar(
                x=list(peak_props.keys()),
                y=list(peak_props.values()),
                name="Properties",
                marker_color="steelblue",
            ),
            row=1,
            col=2,
        )

        # Panel 3: Coverage
        fig.add_trace(
            go.Bar(
                x=["Low threshold", "High threshold"],
                y=[row["frac_above_low"], row["frac_above_high"]],
                name="Coverage",
                marker_color="lightgreen",
            ),
            row=2,
            col=1,
        )

        # Panel 4: Score summary
        scores = {
            "Mean": row["mean_score"],
            "Max": row["max_score"],
            "P95": row.get("p95_score", 0),
        }
        fig.add_trace(
            go.Bar(
                x=list(scores.keys()),
                y=list(scores.values()),
                name="Scores",
                marker_color="purple",
            ),
            row=2,
            col=2,
        )

        fig.update_layout(
            title_text=f"Profile for {transcript_id}",
            height=700,
            showlegend=False,
            template="plotly_white",
        )
        if format == "interactive":
            return fig

    return fig


def plot_peak_location_summary(
    df: pd.DataFrame,
    transcript_id: str,
    format: str = "static",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure]:
    """
    Show inferred peak positions along transcript.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame
    transcript_id : str
        Transcript ID to plot
    format : str
        'static' or 'interactive'
    output_path : Path, optional
        Save static plot to this path

    Returns
    -------
    Figure or go.Figure
    """
    row = df[df["transcript_id"] == transcript_id]
    if row.empty:
        print(f"Transcript {transcript_id} not found in data")
        return None

    row = row.iloc[0]

    if format == "static":
        fig, ax = plt.subplots(figsize=(10, 4))

        # Show top peak position (normalized)
        position = row.get("top_peak_position", 0)
        ax.bar([position], [1], width=0.1, color="steelblue", alpha=0.7)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1.5])
        ax.set_xlabel("Transcript Position (normalized)")
        ax.set_ylabel("Peak Presence")
        ax.set_title(f"Peak Location for {transcript_id}")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        return fig

    else:  # interactive
        position = row.get("top_peak_position", 0)
        fig = go.Figure(
            data=[
                go.Bar(
                    x=[position],
                    y=[1],
                    marker=dict(color="steelblue"),
                    name="Top Peak",
                    width=0.1,
                )
            ]
        )
        fig.update_layout(
            title=f"Peak Location for {transcript_id}",
            xaxis_title="Transcript Position (normalized)",
            yaxis_title="Peak Presence",
            xaxis=dict(range=[0, 1]),
            yaxis=dict(range=[0, 1.5]),
            template="plotly_white",
        )
        return fig
