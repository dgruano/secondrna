"""
secondrna: Snakemake pipeline for RNA structure prediction using multiple tools.
Integrates rG4detector for G-quadruplex prediction on secondary and tertiary structures.
"""

from pathlib import Path

configfile: "config/config.yaml"

# Include rule modules
include: "workflow/rules/input.smk"
include: "workflow/rules/rg4detector_batch.smk"
#include: "workflow/rules/rg4detector.smk"
include: "workflow/rules/utils.smk"
include: "workflow/rules/rg4detector_summary.smk"
include: "workflow/rules/rg4_summary_viz.smk"

rule all:
    input:
        expand(rules.merge_rg4_predictions.output.csv, sample="gencode.v47"),
        expand(rules.merge_rg4_detections.output.csv, sample="gencode.v47"),
        expand(rules.summarize_rg4_detection.output.summary_csv, sample="gencode.v47"),
        expand(rules.add_group_annotations.output.annotated_csv, sample="gencode.v47"),
        expand(rules.visualize_rg4_summary.output.viz_marker, sample="gencode.v47"),
        expand(rules.rg4_count_peaks.output.counts, sample="gencode.v47"),
        expand(rules.rg4_count_peaks.output.stats, sample="gencode.v47")

# Define default target
rule all_tests:
    input:
        expand("results/{sample}/rG4detector_prediction.csv", sample=[100, 500, 1000, 2000]),
        expand("results/{sample}/detection.csv", sample=[100, 500, 1000, 2000])

# Ensure output directory exists
Path("results").mkdir(exist_ok=True)
