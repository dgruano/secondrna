"""
rG4detector visualization library.

Provides reusable functions for generating publication-quality static plots
and interactive visualizations from rG4detector summary CSV output.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns


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

    COLUMN_ALIASES = {
        "index": "transcript_id",
        "top_peak_position_nt": "top_peak_position",
        "top_peak_width_nt": "top_peak_width",
    }

    NUMERIC_COLUMNS = {
        "mean_score",
        "max_score",
        "p95_score",
        "p90_score",
        "p99_score",
        "top_peak_score",
        "top_peak_position",
        "top_peak_width",
        "top_peak_prominence",
        "top_peak_area",
        "n_peaks_above_low",
        "n_peaks_above_high",
        "rg4_density_per_kb",
        "frac_above_low",
        "frac_above_high",
        "transcript_length_nt",
    }

    @staticmethod
    def load(filepath: Union[str, Path]) -> pd.DataFrame:
        """Load CSV and validate required columns."""
        df = pd.read_csv(filepath)

        # Normalize column names: strip whitespace, lowercase
        df.columns = df.columns.str.strip().str.lower()

        # Apply column aliases to handle naming variations in the CSV
        df.rename(columns=RG4SummaryLoader.COLUMN_ALIASES, inplace=True)

        # Coerce expected numeric columns to avoid dtype object errors
        for col in RG4SummaryLoader.NUMERIC_COLUMNS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

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

    valid_metrics = [m for m in metrics if m in df.columns]
    static_fig = None
    interactive_fig = None

    if format in ("static", "both"):
        fig, axes = plt.subplots(
            1, len(valid_metrics), figsize=(5 * len(valid_metrics), 4)
        )
        if len(valid_metrics) == 1:
            axes = [axes]

        for ax, metric in zip(axes, valid_metrics):
            sns.histplot(data=df, x=metric, kde=True, ax=ax, bins=30, color="steelblue")
            ax.set_title(f"Distribution of {metric}")
            ax.set_xlabel(metric)
            ax.set_ylabel("Count")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        static_fig = fig
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        interactive_fig = go.Figure()
        for metric in valid_metrics:
            interactive_fig.add_trace(
                go.Histogram(
                    x=df[metric],
                    name=metric,
                    nbinsx=30,
                    opacity=0.7,
                    showlegend=True,
                )
            )

        interactive_fig.update_layout(
            title="Distribution of rG4 Scores",
            xaxis_title="Score",
            yaxis_title="Count",
            hovermode="x unified",
            template="plotly_white",
        )
        if format == "interactive":
            return interactive_fig

    return interactive_fig if format == "both" else static_fig


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
    valid_metrics = [m for m in metrics if m in df.columns]
    static_fig = None
    interactive_fig = None

    if format in ("static", "both"):
        fig, axes = plt.subplots(
            1, len(valid_metrics), figsize=(5 * len(valid_metrics), 4)
        )
        if len(valid_metrics) == 1:
            axes = [axes]

        for ax, metric in zip(axes, valid_metrics):
            sns.histplot(data=df, x=metric, kde=True, ax=ax, bins=30, color="coral")
            ax.set_title(f"Distribution of {metric}")
            ax.set_xlabel(metric)
            ax.set_ylabel("Count")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        static_fig = fig
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        interactive_fig = go.Figure()
        for metric in valid_metrics:
            interactive_fig.add_trace(
                go.Histogram(
                    x=df[metric], name=metric, nbinsx=30, opacity=0.7, showlegend=True
                )
            )

        interactive_fig.update_layout(
            title="Distribution of Peak Counts and Density",
            xaxis_title="Value",
            yaxis_title="Count",
            hovermode="x unified",
            template="plotly_white",
        )
        if format == "interactive":
            return interactive_fig

    return interactive_fig if format == "both" else static_fig


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
    valid_metrics = [m for m in metrics if m in df.columns]
    static_fig = None
    interactive_fig = None

    if format in ("static", "both"):
        fig, axes = plt.subplots(
            1, len(valid_metrics), figsize=(5 * len(valid_metrics), 4)
        )
        if len(valid_metrics) == 1:
            axes = [axes]

        for ax, metric in zip(axes, valid_metrics):
            sns.violinplot(data=df, y=metric, ax=ax, color="lightgreen")
            ax.set_title(f"Distribution of {metric}")
            ax.set_ylabel(metric)

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        static_fig = fig
        if format == "static":
            return fig

    if format in ("interactive", "both"):
        interactive_fig = go.Figure()
        for metric in valid_metrics:
            interactive_fig.add_trace(
                go.Box(y=df[metric], name=metric, boxmean="sd", showlegend=True)
            )

        interactive_fig.update_layout(
            title="Distribution of Coverage Metrics",
            yaxis_title="Fraction",
            hovermode="y unified",
            template="plotly_white",
        )
        if format == "interactive":
            return interactive_fig

    return interactive_fig if format == "both" else static_fig


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

    # Build combined mask for all filters at once
    mask = pd.Series(True, index=filtered.index)

    for column, (op_str, threshold) in filters.items():
        if column not in filtered.columns:
            print(f"Warning: Column {column} not found in DataFrame")
            continue

        if op_str == ">":
            mask &= filtered[column] > threshold
        elif op_str == "<":
            mask &= filtered[column] < threshold
        elif op_str == ">=":
            mask &= filtered[column] >= threshold
        elif op_str == "<=":
            mask &= filtered[column] <= threshold

    return filtered[mask]


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

    # Filter to valid metrics (handle missing columns)
    key_metrics = [m for m in key_metrics if m in df.columns]

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

    # Pre-compute panel data once for both formats
    peaks_data = [row["n_peaks_above_low"], row["n_peaks_above_high"]]
    peak_props = {
        "Score": row["top_peak_score"],
        "Width": row.get("top_peak_width", 0),
        "Prominence": row.get("top_peak_prominence", 0),
    }
    coverage_data = [row["frac_above_low"], row["frac_above_high"]]
    scores = {
        "Mean": row["mean_score"],
        "Max": row["max_score"],
        "P95": row.get("p95_score", 0),
    }

    static_fig = None
    interactive_fig = None

    if format in ("static", "both"):
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle(f"Profile for {transcript_id}", fontsize=14, fontweight="bold")

        # Panel 1: Peak counts
        ax = axes[0, 0]
        ax.bar(["Low threshold", "High threshold"], peaks_data, color=["coral", "red"])
        ax.set_ylabel("Number of peaks")
        ax.set_title("Peak Counts")

        # Panel 2: Top peak properties
        ax = axes[0, 1]
        ax.bar(peak_props.keys(), peak_props.values(), color="steelblue")
        ax.set_ylabel("Value")
        ax.set_title("Top Peak Properties")

        # Panel 3: Coverage metrics
        ax = axes[1, 0]
        ax.bar(["Low threshold", "High threshold"], coverage_data, color="lightgreen")
        ax.set_ylabel("Fraction")
        ax.set_ylim([0, 1])
        ax.set_title("Coverage Metrics")

        # Panel 4: Score summary
        ax = axes[1, 1]
        ax.bar(scores.keys(), scores.values(), color="purple", alpha=0.7)
        ax.set_ylabel("Score")
        ax.set_title("Score Summary")

        plt.tight_layout()
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
        static_fig = fig
        if format == "static":
            return fig

    if format in ("interactive", "both"):
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
                y=peaks_data,
                name="Peaks",
                marker_color="coral",
            ),
            row=1,
            col=1,
        )

        # Panel 2: Top peak properties
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
                y=coverage_data,
                name="Coverage",
                marker_color="lightgreen",
            ),
            row=2,
            col=1,
        )

        # Panel 4: Score summary
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
        interactive_fig = fig
        if format == "interactive":
            return fig

    return interactive_fig if format == "both" else static_fig


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


def plot_class_comparison(
    df: pd.DataFrame,
    metrics: Optional[List[str]] = None,
    class_col: str = "real",
    format: str = "both",
    output_path: Optional[Path] = None,
) -> Optional[Union[plt.Figure, go.Figure]]:
    """
    Compare rG4 metric distributions between transcript classes.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame containing a class column (e.g., 'real') with
        boolean values (True = Protein-coding, False = lncRNA).
    metrics : List[str], optional
        Metrics to compare (default: 8 core rG4 metrics).
    class_col : str
        Column name containing class labels (default: 'real').
    format : str
        'static', 'interactive', or 'both'.
    output_path : Path, optional
        Directory to save output file(s). Static PNG saved as
        'class_comparison.png'; interactive HTML as
        'class_comparison_interactive.html'.

    Returns
    -------
    plt.Figure, go.Figure, or None
        matplotlib Figure for 'static', plotly Figure for 'interactive' or
        'both'. Returns None if class_col is missing from df.
    """
    if class_col not in df.columns:
        print(f"Warning: Column '{class_col}' not found in DataFrame")
        return None

    if metrics is None:
        metrics = [
            "mean_score",
            "max_score",
            "top_peak_score",
            "n_peaks_above_low",
            "n_peaks_above_high",
            "rg4_density_per_kb",
            "frac_above_low",
            "frac_above_high",
        ]

    # Map boolean class labels to human-readable strings
    label_map = {True: "Protein-coding", False: "lncRNA"}
    plot_df = df.copy()
    plot_df["_class_label"] = plot_df[class_col].map(label_map)

    static_fig = None
    interactive_fig = None

    # --- Static: seaborn violin grid (2 rows x 4 cols) ---
    if format in ("static", "both"):
        n_cols = 4
        n_rows = 2
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
        axes_flat = axes.flatten()

        for i, metric in enumerate(metrics):
            ax = axes_flat[i]
            if metric in plot_df.columns:
                sns.violinplot(
                    data=plot_df,
                    x="_class_label",
                    y=metric,
                    hue="_class_label",
                    ax=ax,
                    palette={"Protein-coding": "steelblue", "lncRNA": "coral"},
                    legend=False,
                )
                ax.set_title(metric)
                ax.set_xlabel("")
                ax.set_ylabel(metric)
            else:
                ax.set_visible(False)

        # Hide any unused axes
        for j in range(len(metrics), len(axes_flat)):
            axes_flat[j].set_visible(False)

        plt.suptitle(
            "rG4 Metric Distributions by Transcript Class",
            fontsize=14,
            fontweight="bold",
        )
        plt.tight_layout()

        if output_path is not None:
            out = Path(output_path)
            save_path = out / "class_comparison.png" if out.is_dir() else out
            plt.savefig(save_path, dpi=300, bbox_inches="tight")

        static_fig = fig

    # --- Interactive: plotly grouped violin (vectorized) ---
    if format in ("interactive", "both"):
        interactive_fig = go.Figure()

        # Melt data once for efficient grouping
        available_metrics = [m for m in metrics if m in plot_df.columns]
        melted_df = plot_df[["_class_label"] + available_metrics].melt(
            id_vars=["_class_label"], var_name="metric", value_name="value"
        )

        classes = sorted(plot_df["_class_label"].dropna().unique())
        colors = {"Protein-coding": "steelblue", "lncRNA": "coral"}

        for cls in classes:
            cls_data = melted_df[melted_df["_class_label"] == cls]
            interactive_fig.add_trace(
                go.Violin(
                    x=cls_data["metric"],
                    y=cls_data["value"],
                    name=cls,
                    legendgroup=cls,
                    scalegroup=cls,
                    showlegend=True,
                    line_color=colors.get(cls, None),
                    box_visible=True,
                    meanline_visible=True,
                )
            )

        interactive_fig.update_layout(
            title="rG4 Metric Distributions by Transcript Class",
            xaxis_title="Metric",
            yaxis_title="Value",
            violinmode="group",
            template="plotly_white",
            legend_title="Class",
        )

        if output_path is not None:
            out = Path(output_path)
            save_path = (
                out / "class_comparison_interactive.html" if out.is_dir() else out
            )
            interactive_fig.write_html(str(save_path))

    if format == "static":
        return static_fig
    return interactive_fig


def plot_correlation_heatmap(
    df: pd.DataFrame,
    metrics: Optional[List[str]] = None,
    class_col: Optional[str] = None,
    format: str = "both",
    output_path: Optional[Path] = None,
) -> Union[plt.Figure, go.Figure, List[Union[plt.Figure, go.Figure]]]:
    """
    Show Pearson correlations between rG4 metrics.

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame.
    metrics : List[str], optional
        Metrics to correlate. If None, auto-selects all numeric columns.
    class_col : str, optional
        Column name for class labels. If provided and present, produces
        N+1 heatmaps (global + one per class value).
    format : str
        'static', 'interactive', or 'both'.
    output_path : Path, optional
        Directory to save output file(s). Files saved with suffixes:
        '_global', '_pc', '_lncrna' for class-specific heatmaps.

    Returns
    -------
    plt.Figure, go.Figure, or List of Figures
        Single figure if class_col is None, list of figures if class_col
        is provided and present.
    """
    # Auto-detect numeric columns if metrics not specified
    if metrics is None:
        metrics = df.select_dtypes(include="number").columns.tolist()

    # Filter to available metrics
    available_metrics = [m for m in metrics if m in df.columns]
    if not available_metrics:
        print("Warning: No valid numeric metrics found")
        return None

    # Determine if we're doing class-specific heatmaps
    do_class_split = class_col is not None and class_col in df.columns

    figures = []

    # Helper function to create one heatmap
    def create_heatmap(data: pd.DataFrame, title_suffix: str, file_suffix: str):
        corr_matrix = data[available_metrics].corr(method="pearson")

        static_fig = None
        interactive_fig = None

        # --- Static: seaborn annotated heatmap ---
        if format in ("static", "both"):
            fig, ax = plt.subplots(figsize=(10, 8))
            sns.heatmap(
                corr_matrix,
                annot=True,
                fmt=".2f",
                cmap="RdBu_r",
                vmin=-1,
                vmax=1,
                center=0,
                square=True,
                linewidths=0.5,
                cbar_kws={"shrink": 0.8},
                ax=ax,
            )
            ax.set_title(f"Metric Correlation Heatmap{title_suffix}")
            plt.tight_layout()

            if output_path is not None:
                out = Path(output_path)
                save_path = (
                    out / f"correlation_heatmap{file_suffix}.png"
                    if out.is_dir()
                    else out.parent / f"correlation_heatmap{file_suffix}.png"
                )
                plt.savefig(save_path, dpi=300, bbox_inches="tight")

            static_fig = fig

        # --- Interactive: plotly heatmap ---
        if format in ("interactive", "both"):
            interactive_fig = go.Figure(
                data=go.Heatmap(
                    z=corr_matrix.values,
                    x=corr_matrix.columns,
                    y=corr_matrix.index,
                    colorscale="RdBu_r",
                    zmid=0,
                    zmin=-1,
                    zmax=1,
                    text=corr_matrix.values,
                    texttemplate="%{text:.2f}",
                    textfont={"size": 10},
                    colorbar=dict(title="Correlation"),
                    hovertemplate="<b>%{x}</b> vs <b>%{y}</b><br>Correlation: %{z:.3f}<extra></extra>",
                )
            )
            interactive_fig.update_layout(
                title=f"Metric Correlation Heatmap{title_suffix}",
                xaxis_title="Metric",
                yaxis_title="Metric",
                template="plotly_white",
                width=800,
                height=700,
            )

            if output_path is not None:
                out = Path(output_path)
                save_path = (
                    out / f"correlation_heatmap{file_suffix}_interactive.html"
                    if out.is_dir()
                    else out.parent
                    / f"correlation_heatmap{file_suffix}_interactive.html"
                )
                interactive_fig.write_html(str(save_path))

        if format == "static":
            return static_fig
        return interactive_fig

    # Create global heatmap
    fig = create_heatmap(df, "", "" if not do_class_split else "_global")
    figures.append(fig)

    # Create per-class heatmaps if requested
    if do_class_split:
        label_map = {True: "Protein-coding", False: "lncRNA"}
        suffix_map = {True: "_pc", False: "_lncrna"}

        for class_value in sorted(df[class_col].dropna().unique()):
            class_data = df[df[class_col] == class_value]
            class_label = label_map.get(class_value, str(class_value))
            file_suffix = suffix_map.get(class_value, f"_{class_value}")

            fig = create_heatmap(class_data, f" — {class_label}", file_suffix)
            figures.append(fig)

    # Return single figure or list
    if len(figures) == 1:
        return figures[0]
    return figures


def plot_transcript_dashboard(
    df: pd.DataFrame,
    transcript_id: str,
    class_col: str = "real",
    output_path: Optional[Path] = None,
) -> Optional[go.Figure]:
    """
    Full interactive profile for a single transcript with percentile ranks.

    Always interactive (returns go.Figure, exports .html).

    Parameters
    ----------
    df : pd.DataFrame
        Summary DataFrame.
    transcript_id : str
        Transcript ID to plot.
    class_col : str
        Column name for class labels (default: 'real').
    output_path : Path, optional
        Directory to save output HTML file as
        'transcript_{transcript_id}_dashboard.html'.

    Returns
    -------
    go.Figure or None
        Plotly figure with 3-row × 2-col subplot grid. Returns None if
        transcript not found.
    """
    from plotly.subplots import make_subplots

    # Find transcript
    row = df[df["transcript_id"] == transcript_id]
    if row.empty:
        print(f"Warning: Transcript {transcript_id} not found in data")
        return None

    row = row.iloc[0]

    # Create 3×2 subplot grid
    fig = make_subplots(
        rows=3,
        cols=2,
        subplot_titles=(
            "Peak Counts",
            "Top Peak Properties",
            "Coverage Metrics",
            "Score Summary",
            "Percentile Ranks (Full Dataset & Within Class)",
        ),
        specs=[
            [{"type": "bar"}, {"type": "bar"}],
            [{"type": "bar"}, {"type": "bar"}],
            [{"type": "bar", "colspan": 2}, None],
        ],
        row_heights=[0.3, 0.3, 0.4],
        vertical_spacing=0.12,
        horizontal_spacing=0.1,
    )

    # Row 1, Col 1: Peak counts
    fig.add_trace(
        go.Bar(
            x=["Low threshold", "High threshold"],
            y=[row["n_peaks_above_low"], row["n_peaks_above_high"]],
            name="Peaks",
            marker_color="coral",
            showlegend=False,
        ),
        row=1,
        col=1,
    )

    # Row 1, Col 2: Top peak properties
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
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    # Row 2, Col 1: Coverage
    fig.add_trace(
        go.Bar(
            x=["Low threshold", "High threshold"],
            y=[row["frac_above_low"], row["frac_above_high"]],
            name="Coverage",
            marker_color="lightgreen",
            showlegend=False,
        ),
        row=2,
        col=1,
    )

    # Row 2, Col 2: Score summary
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
            showlegend=False,
        ),
        row=2,
        col=2,
    )

    # Row 3 (full width): Percentile ranks
    key_metrics = [
        "top_peak_score",
        "n_peaks_above_high",
        "rg4_density_per_kb",
        "max_score",
        "frac_above_high",
    ]

    # Filter to valid metrics
    valid_key_metrics = [m for m in key_metrics if m in df.columns]
    metric_labels = valid_key_metrics

    # Compute percentiles (vectorized)
    has_class = class_col in df.columns and pd.notna(row.get(class_col))

    # Full dataset percentiles (vectorized)
    percentiles_full = []
    for metric in valid_key_metrics:
        value = row[metric]
        percentile_full = (df[metric] <= value).mean() * 100
        percentiles_full.append(percentile_full)

    # Within-class percentiles (vectorized)
    percentiles_class = []
    if has_class:
        class_value = row[class_col]
        class_data = df[df[class_col] == class_value]
        for metric in valid_key_metrics:
            value = row[metric]
            percentile_class = (class_data[metric] <= value).mean() * 100
            percentiles_class.append(percentile_class)

    # Add full dataset percentile bars
    fig.add_trace(
        go.Bar(
            y=metric_labels,
            x=percentiles_full,
            name="Full Dataset",
            orientation="h",
            marker_color="steelblue",
            text=[f"Top {100-p:.0f}%" for p in percentiles_full],
            textposition="outside",
        ),
        row=3,
        col=1,
    )

    # Add within-class percentile bars if available
    if has_class:
        class_label = "Protein-coding" if row[class_col] else "lncRNA"
        fig.add_trace(
            go.Bar(
                y=metric_labels,
                x=percentiles_class,
                name=f"Within {class_label}",
                orientation="h",
                marker_color="coral",
                text=[f"Top {100-p:.0f}%" for p in percentiles_class],
                textposition="outside",
            ),
            row=3,
            col=1,
        )

    # Update layout
    fig.update_xaxes(range=[0, 100], row=3, col=1, title_text="Percentile")
    fig.update_layout(
        title_text=f"Dashboard for {transcript_id}",
        height=1000,
        showlegend=True,
        template="plotly_white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.05,
            xanchor="center",
            x=0.5,
        ),
    )

    # Save if output path provided
    if output_path is not None:
        out = Path(output_path)
        save_path = (
            out / f"transcript_{transcript_id}_dashboard.html" if out.is_dir() else out
        )
        fig.write_html(str(save_path))

    return fig
