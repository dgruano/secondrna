# rG4detector Summary and Visualization Rules
# Integrates rg4detector_summarize.py and rg4_visualize_cli.py

# Configuration defaults for summarization
THRESHOLD_LOW = config.get("threshold_low", 1.0)
THRESHOLD_HIGH = config.get("threshold_high", 2.0)
MIN_PEAK_PROMINENCE = config.get("min_peak_prominence", 0.3)
MIN_PEAK_WIDTH = config.get("min_peak_width", 20)
ANNOTATION_FILE = config.get("annotation_file", None)

# Configuration defaults for group annotation
GROUP_ANNOTATION_FILE = config.get("group_annotation_file", None)
GROUP_ANNOTATION_ID_COLUMN = config.get("group_annotation_id_column", None)
GROUP_ANNOTATION_COLUMNS = config.get("group_annotation_columns", None)
GROUP_ANNOTATION_SEPARATOR = config.get("group_annotation_separator", "\t")
GROUP_ANNOTATION_KEEP_UNMATCHED = config.get("group_annotation_keep_unmatched", False)

# Configuration defaults for visualization
VIZ_OUTPUT_DIR = config.get("viz_output_dir", "results/{sample}/viz")
VIZ_PLOTS = config.get("viz_plots", ["distributions", "ranking", "details"])
VIZ_FORMAT = config.get("viz_format", "both")
VIZ_TOP_N = config.get("viz_top_n", 20)
VIZ_METRIC = config.get("viz_metric", "top_peak_score")
VIZ_DPI = config.get("viz_dpi", 300)
VIZ_FIGSIZE = config.get("viz_figsize", "10x6")
VIZ_CLASS_COL = config.get("viz_class_col", "real")
VIZ_TRANSCRIPT_DETAIL = config.get("viz_transcript_detail", None)
VIZ_FILTER = config.get("viz_filter", None)


rule summarize_rg4_detection:
    """
    Summarize rG4detector detection output per transcript.
    Computes peak-based, distribution, and density metrics.

    Parameters (configurable via --config):
      - threshold_low: Score threshold for candidate rG4 peaks (default: 1.0)
      - threshold_high: Score threshold for high-confidence rG4 peaks (default: 2.0)
      - min_peak_prominence: Minimum peak prominence (default: 0.3)
      - min_peak_width: Minimum peak width in nucleotides (default: 20)
      - annotation_file: Optional BED file with feature annotations
    """
    input:
        detection_csv = "results/{sample}/detection.csv"
    output:
        summary_csv = "results/{sample}/rg4_summary.csv"
    params:
        threshold_low = THRESHOLD_LOW,
        threshold_high = THRESHOLD_HIGH,
        min_peak_prominence = MIN_PEAK_PROMINENCE,
        min_peak_width = MIN_PEAK_WIDTH,
        annotation_file = ANNOTATION_FILE
    log:
        "logs/{sample}/summarize_rg4_detection.log"
    benchmark:
        "benchmarks/{sample}/summarize_rg4_detection.tsv"
    resources:
        runtime=30,
        mem_mb=1024*20,
        cpus_per_task=1
    shell:
        """
        {{
            CMD="python workflow/scripts/rg4detector_summarize.py \
                --input {input.detection_csv} \
                --output {output.summary_csv} \
                --threshold_low {params.threshold_low} \
                --threshold_high {params.threshold_high} \
                --min_peak_prominence {params.min_peak_prominence} \
                --min_peak_width {params.min_peak_width}"

            if [ ! -z "{params.annotation_file}" ] && [ "{params.annotation_file}" != "None" ]; then
                CMD="$CMD --annotation {params.annotation_file}"
            fi

            eval "$CMD"
        }} 2>&1 | tee {log}
        """


rule add_group_annotations:
    """
    Add group annotations (e.g., protein-coding vs lncRNA) to rG4detector summary.
    Enriches summary with external group/class labels.

    Parameters (configurable via --config):
      - group_annotation_file: Annotation file path (TSV, CSV with ID column)
      - group_annotation_id_column: Name of ID column in annotation file (auto-detect if not set)
      - group_annotation_columns: Comma-separated columns to add (default: all)
      - group_annotation_separator: Delimiter for annotation file (default: tab)
      - group_annotation_keep_unmatched: Keep unmatched summary rows (default: False)

    Example annotation file (TSV):
        seq_ID      group       biotype
        TX001       protein     protein_coding
        TX002       lncRNA      lncRNA
    """
    input:
        summary_csv = "results/{sample}/rg4_summary.csv",
        annotation_file = GROUP_ANNOTATION_FILE if GROUP_ANNOTATION_FILE else []
    output:
        annotated_csv = "results/{sample}/rg4_summary_annotated.csv"
    params:
        annotation_file = GROUP_ANNOTATION_FILE,
        id_column = GROUP_ANNOTATION_ID_COLUMN,
        columns = GROUP_ANNOTATION_COLUMNS,
        separator = GROUP_ANNOTATION_SEPARATOR,
        keep_unmatched = GROUP_ANNOTATION_KEEP_UNMATCHED
    log:
        "logs/{sample}/add_group_annotations.log"
    benchmark:
        "benchmarks/{sample}/add_group_annotations.tsv"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=1
    shell:
        """
        {{
            if [ -z "{params.annotation_file}" ] || [ "{params.annotation_file}" = "None" ]; then
                echo "[add_annotations] No annotation file specified, copying summary as-is" >&2
                cp {input.summary_csv} {output.annotated_csv}
            else
                CMD="python workflow/scripts/add_group_annotations.py \
                    --summary {input.summary_csv} \
                    --annotations {params.annotation_file} \
                    --output {output.annotated_csv} \
                    --separator '{params.separator}'"

                if [ ! -z "{params.id_column}" ] && [ "{params.id_column}" != "None" ]; then
                    CMD="$CMD --annotation-id-column {params.id_column}"
                fi

                if [ ! -z "{params.columns}" ] && [ "{params.columns}" != "None" ]; then
                    CMD="$CMD --annotation-columns {params.columns}"
                fi

                if [ "{params.keep_unmatched}" = "True" ]; then
                    CMD="$CMD --keep-unmatched"
                fi

                eval "$CMD"
            fi
        }} 2>&1 | tee {log}
        """


rule visualize_rg4_summary:
    """
    Generate publication-quality visualizations from rG4detector summary.
    Creates static PNG plots and interactive HTML dashboards.
    Uses annotated summary if available, otherwise falls back to plain summary.

    Parameters (configurable via --config):
      - viz_output_dir: Output directory root (default: results/{sample}/viz)
      - viz_plots: Plot types to generate (default: ["distributions", "ranking", "details"])
      - viz_format: Output format - "static", "interactive", or "both" (default: both)
      - viz_top_n: Number of top transcripts for ranking (default: 20)
      - viz_metric: Metric for ranking (default: top_peak_score)
        choices: top_peak_score, max_score, n_peaks_above_high, rg4_density_per_kb
      - viz_dpi: DPI for PNG files (default: 300)
      - viz_figsize: Figure size as "WxH" (default: 10x6)
      - viz_class_col: Column name for class labels (default: real)
      - viz_transcript_detail: Single transcript ID for detail view (optional)
      - viz_filter: Filter string e.g. "n_peaks_above_high > 0" (optional)
    """
    input:
        summary_csv = "results/{sample}/rg4_summary_annotated.csv"
    output:
        viz_marker = "results/{sample}/.viz_done"
    conda:
        "rg4_visualization"
    params:
        output_dir = VIZ_OUTPUT_DIR,
        plots = VIZ_PLOTS,
        format = VIZ_FORMAT,
        top_n = VIZ_TOP_N,
        metric = VIZ_METRIC,
        dpi = VIZ_DPI,
        figsize = VIZ_FIGSIZE,
        class_col = VIZ_CLASS_COL,
        transcript_detail = VIZ_TRANSCRIPT_DETAIL,
        filter_str = VIZ_FILTER
    log:
        "logs/{sample}/visualize_rg4_summary.log"
    benchmark:
        "benchmarks/{sample}/visualize_rg4_summary.tsv"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=2
    shell:
        """
        {{
            CMD="python workflow/scripts/rg4_visualize_cli.py \
                --input {input.summary_csv} \
                --output {params.output_dir} \
                --format {params.format} \
                --top-n-plot {params.top_n} \
                --metric {params.metric} \
                --dpi {params.dpi} \
                --figsize {params.figsize} \
                --class-col {params.class_col} \
                --plots {params.plots}"

            if [ ! -z "{params.transcript_detail}" ] && [ "{params.transcript_detail}" != "None" ]; then
                CMD="$CMD --transcript-detail {params.transcript_detail}"
            fi

            if [ ! -z "{params.filter_str}" ] && [ "{params.filter_str}" != "None" ]; then
                CMD="$CMD --filter {params.filter_str}"
            fi
            echo "dummy"
            eval "$CMD"
            touch {output.viz_marker}
        }} 2>&1 | tee {log}
        """
