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
        detection="results/{sample}/rg4detector/detection.csv",
        labels=PEAK_LABELS_FILE if PEAK_LABELS_FILE else [],
    output:
        counts="results/{sample}/rg4detector/peak_counts.csv",
        stats="results/{sample}/rg4detector/peak_stats.tsv",
    log:
        "logs/{sample}/rg4detector/rg4_count_peaks.log",
    benchmark:
        "benchmarks/{sample}/rg4detector/rg4_count_peaks.tsv"
    conda:
        "../envs/rg4_visualization.yaml"
    resources:
        runtime=30,
        mem_mb=1024 * 10,
        cpus_per_task=1,
    params:
        threshold=PEAK_THRESHOLD,
    script:
        "../scripts/rg4_count_peaks.py"


rule rg4_peaks_to_bed:
    """
Convert detected rG4 peaks to BED format for visualization.
"""
    input:
        detection="results/{sample}/rg4detector/detection.csv",
    output:
        bed="results/{sample}/rg4detector/{sample}.rG4detector.peaks.bed",
    log:
        "logs/{sample}/rg4detector/rg4_peaks_to_bed.log",
    resources:
        runtime=10,
        mem_mb=1024 * 5,
    params:
        d=config.get("peak_extension", 0),
    script:
        "../scripts/rg4_peaks_to_bed.py"
