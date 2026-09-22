"""
secondrna: Snakemake pipeline for RNA structure prediction using multiple tools.
Integrates rG4detector for G-quadruplex prediction on secondary and tertiary structures.
"""

from pathlib import Path


configfile: "config/config.yaml"


configured_samples = config["samples"]
samples = (
    configured_samples if isinstance(configured_samples, list) else [configured_samples]
)


# Include rule modules
include: "workflow/rules/rg4detector_batch.smk"
include: "workflow/rules/rg4detector.smk"
include: "workflow/rules/utils.smk"
include: "workflow/rules/rg4detector_summary.smk"
include: "workflow/rules/rg4_summary_viz.smk"
include: "workflow/rules/batch.smk"
include: "workflow/rules/g4Discovery.smk"
include: "workflow/rules/scanfold_core.smk"


rule all:
    input:
        expand(rules.rg4_detect.output, sample=samples),
        expand(rules.merge_rg4_predictions.output.csv, sample=samples),
        expand(rules.merge_rg4_detections.output.csv, sample=samples),
        expand(rules.summarize_rg4_detection.output.summary_csv, sample=samples),
        expand(rules.add_group_annotations.output.annotated_csv, sample=samples),
        expand(rules.visualize_rg4_summary.output.viz_marker, sample=samples),
        expand(rules.rg4_count_peaks.output.counts, sample=samples),
        expand(rules.rg4_count_peaks.output.stats, sample=samples),


# Ensure output directory exists
Path("results").mkdir(exist_ok=True)
