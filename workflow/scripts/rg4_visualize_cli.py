#!/usr/bin/env python3
"""
CLI for rG4detector visualization and analysis.

Generates publication-quality static plots and interactive HTML dashboards
from rG4detector summary CSV output.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, List

# Add parent directories to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils.rg4_viz import (
    RG4SummaryLoader,
    plot_score_distribution,
    plot_peak_counts_distribution,
    plot_coverage_distribution,
    plot_top_transcripts,
    create_ranking_table,
    apply_threshold_filter,
    plot_transcript_profile,
    plot_peak_location_summary,
)


def setup_output_dirs(output_dir: Path) -> tuple:
    """Create output directory structure."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    viz_static = output_dir / "viz_static"
    viz_interactive = output_dir / "viz_interactive"
    viz_reports = output_dir / "viz_reports"

    viz_static.mkdir(exist_ok=True)
    viz_interactive.mkdir(exist_ok=True)
    viz_reports.mkdir(exist_ok=True)

    return output_dir, viz_static, viz_interactive, viz_reports


def parse_filters(filter_args: Optional[List[str]]) -> dict:
    """Parse filter arguments into a dictionary."""
    filters = {}
    if not filter_args:
        return filters

    for filter_str in filter_args:
        # Parse "column > value" or "column < value"
        for op in [">=", "<=", ">", "<"]:
            if op in filter_str:
                parts = filter_str.split(op)
                if len(parts) == 2:
                    col = parts[0].strip()
                    val = parts[1].strip()
                    try:
                        filters[col] = (op, float(val))
                    except ValueError:
                        print(f"Warning: Could not parse filter {filter_str}")
                break

    return filters


def main():
    parser = argparse.ArgumentParser(
        description="Visualize and analyze rG4detector summary results"
    )

    parser.add_argument(
        "--input",
        required=True,
        type=str,
        help="Path to rg4_summary.csv",
    )

    parser.add_argument(
        "--output",
        default="results/viz",
        type=str,
        help="Output directory (default: results/viz)",
    )

    parser.add_argument(
        "--plots",
        nargs="+",
        default=["distributions", "ranking", "details"],
        choices=["distributions", "ranking", "details", "all"],
        help="Which visualizations to generate (default: all)",
    )

    parser.add_argument(
        "--format",
        default="both",
        choices=["static", "interactive", "both"],
        help="Output format (default: both)",
    )

    parser.add_argument(
        "--top-n-plot",
        type=int,
        default=20,
        help="Number of top transcripts for ranking plot (default: 20)",
    )

    parser.add_argument(
        "--metric",
        default="top_peak_score",
        choices=[
            "top_peak_score",
            "max_score",
            "n_peaks_above_high",
            "rg4_density_per_kb",
        ],
        help="Metric for ranking (default: top_peak_score)",
    )

    parser.add_argument(
        "--transcript-detail",
        type=str,
        help="Generate detail view for single transcript (optional)",
    )

    parser.add_argument(
        "--filter",
        nargs="*",
        help='Apply filter (e.g., "n_peaks_above_high > 0"), repeatable',
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for PNG files (default: 300)",
    )

    parser.add_argument(
        "--figsize",
        type=str,
        default="10x6",
        help="Figure size for static plots (default: 10x6)",
    )

    args = parser.parse_args()

    # Setup
    print(f"Loading summary from {args.input}")
    df = RG4SummaryLoader.load(args.input)
    print(f"Loaded {len(df)} transcripts")

    output_dir, viz_static, viz_interactive, viz_reports = setup_output_dirs(
        args.output
    )
    print(f"Output directory: {output_dir}")

    # Apply filters if specified
    if args.filter:
        filters = parse_filters(args.filter)
        if filters:
            df = apply_threshold_filter(df, filters)
            print(f"After filtering: {len(df)} transcripts")

    # Determine which plots to generate
    plot_types = args.plots
    if "all" in plot_types:
        plot_types = ["distributions", "ranking", "details"]

    # Generate distribution plots
    if "distributions" in plot_types:
        print("\nGenerating distribution plots...")

        if args.format in ("static", "both"):
            plot_score_distribution(
                df,
                format="static",
                output_path=viz_static / "score_distribution.png",
            )
            print(f"  ✓ {viz_static / 'score_distribution.png'}")

            plot_peak_counts_distribution(
                df, format="static", output_path=viz_static / "peak_counts_distribution.png"
            )
            print(f"  ✓ {viz_static / 'peak_counts_distribution.png'}")

            plot_coverage_distribution(
                df, format="static", output_path=viz_static / "coverage_distribution.png"
            )
            print(f"  ✓ {viz_static / 'coverage_distribution.png'}")

        if args.format in ("interactive", "both"):
            fig = plot_score_distribution(df, format="interactive")
            fig.write_html(viz_interactive / "score_distribution_interactive.html")
            print(f"  ✓ {viz_interactive / 'score_distribution_interactive.html'}")

            fig = plot_peak_counts_distribution(df, format="interactive")
            fig.write_html(viz_interactive / "peak_counts_distribution_interactive.html")
            print(f"  ✓ {viz_interactive / 'peak_counts_distribution_interactive.html'}")

            fig = plot_coverage_distribution(df, format="interactive")
            fig.write_html(viz_interactive / "coverage_distribution_interactive.html")
            print(f"  ✓ {viz_interactive / 'coverage_distribution_interactive.html'}")

    # Generate ranking plots
    if "ranking" in plot_types:
        print("\nGenerating ranking visualizations...")

        if args.format in ("static", "both"):
            fig = plot_top_transcripts(
                df, metric=args.metric, n=args.top_n_plot, format="static"
            )
            output_file = (
                viz_static / f"top_{args.top_n_plot}_by_{args.metric}.png"
            )
            fig.savefig(output_file, dpi=args.dpi, bbox_inches="tight")
            print(f"  ✓ {output_file}")

        if args.format in ("interactive", "both"):
            fig = plot_top_transcripts(
                df, metric=args.metric, n=args.top_n_plot, format="interactive"
            )
            output_file = (
                viz_interactive / f"top_{args.top_n_plot}_by_{args.metric}.html"
            )
            fig.write_html(output_file)
            print(f"  ✓ {output_file}")

        # Ranking table (interactive only)
        fig = create_ranking_table(df, top_n=min(50, len(df)))
        output_file = viz_interactive / "ranking_table.html"
        fig.write_html(output_file)
        print(f"  ✓ {output_file}")

    # Generate per-transcript detail views
    if "details" in plot_types:
        print("\nGenerating transcript detail views...")

        if args.transcript_detail:
            # Specific transcript requested
            transcripts = [args.transcript_detail]
        else:
            # Top 3 transcripts
            transcripts = df.nlargest(3, "top_peak_score")["transcript_id"].tolist()

        for transcript_id in transcripts:
            if args.format in ("static", "both"):
                fig = plot_transcript_profile(
                    df, transcript_id, format="static"
                )
                if fig:
                    output_file = viz_static / f"transcript_{transcript_id}_profile.png"
                    fig.savefig(output_file, dpi=args.dpi, bbox_inches="tight")
                    print(f"  ✓ {output_file}")

            if args.format in ("interactive", "both"):
                fig = plot_transcript_profile(
                    df, transcript_id, format="interactive"
                )
                if fig:
                    output_file = viz_interactive / f"transcript_{transcript_id}_profile.html"
                    fig.write_html(output_file)
                    print(f"  ✓ {output_file}")

    print("\n✓ All visualizations generated successfully!")
    print(f"Static plots: {viz_static}")
    print(f"Interactive plots: {viz_interactive}")
    print(f"Reports: {viz_reports}")


if __name__ == "__main__":
    main()
