# rG4 Peak Count and Statistics Rules

PEAK_THRESHOLD = config.get("peak_threshold", 1.56)
PEAK_LABELS_FILE = config.get("group_annotation_file", None)


rule rg4_count_peaks:
    """
    Count rG4 peaks per transcript and compute summary statistics.
    Streams detection.csv, calls peaks above a threshold, and writes:
      - peak_counts.csv: per-transcript peak counts and lengths
      - peak_stats.tsv:  aggregate stats; class-association tests if labels provided

    Parameters (configurable via --config):
      - peak_threshold:       Peak height threshold (default: 1.56)
      - group_annotation_file: TSV with seq_ID and real columns (optional;
                               enables chi-square and Mann-Whitney tests)
    """
    input:
        detection = "results/{sample}/detection.csv",
        labels    = PEAK_LABELS_FILE if PEAK_LABELS_FILE else []
    output:
        counts = "results/{sample}/peak_counts.csv",
        stats  = "results/{sample}/peak_stats.tsv"
    params:
        threshold = PEAK_THRESHOLD
    log:
        "logs/{sample}/rg4_count_peaks.log"
    benchmark:
        "benchmarks/{sample}/rg4_count_peaks.tsv"
    resources:
        runtime    = 30,
        mem_mb     = 1024 * 8,
        cpus_per_task = 1
    script:
        "../scripts/rg4_count_peaks.py"
