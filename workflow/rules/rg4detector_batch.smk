# Global wildcard constraints
wildcard_constraints:
    sample="[a-zA-Z0-9_\\-\\.]+",  # Disallow forward slashes


rule rg4_predict_batch:
    """
Run rG4detector prediction on individual batch FASTA.
Memory is dynamically allocated based on batch_size.
"""
    input:
        batch="results/{sample}/batches/batch_{batch_id}.fa",
        setup="results/.rg4_setup_done",
    output:
        csv="results/{sample}/predictions/{batch_id}/rG4detector_prediction.csv",
    log:
        "logs/{sample}/predict_batch_{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/predict_batch_{batch_id}.tsv"
    conda:
        "test_rg4"
    resources:
        runtime=120,
        mem_mb=lambda wc: calculate_memory_for_batch(BATCH_SIZE),
        cpus_per_task=1,
    shell:
        """
        {{
            mkdir -p $(dirname {output.csv})
            IN=$(realpath {input.batch})
            OUT=$(dirname $(realpath {output.csv}))
            echo "Running rG4detector prediction on batch {wildcards.batch_id}"
            cd software/rG4detector/code
            python predict_fasta.py -f $IN -o $OUT
            cd ../../../
        }} 2>&1 | tee {log}
        """


rule rg4_detect_batch:
    """
Run rG4detector detection on individual batch FASTA.
Memory is dynamically allocated based on batch_size.
"""
    input:
        batch="results/{sample}/batches/batch_{batch_id}.fa",
        setup="results/.rg4_setup_done",
    output:
        csv="results/{sample}/detections/{batch_id}/detection.csv",
    log:
        "logs/{sample}/detect_batch_{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/detect_batch_{batch_id}.tsv"
    conda:
        "test_rg4"
    resources:
        runtime=120,
        mem_mb=lambda wc: calculate_memory_for_batch(BATCH_SIZE),
        cpus_per_task=1,
    shell:
        """
        {{
            mkdir -p $(dirname {output.csv})
            IN=$(realpath {input.batch})
            OUT=$(dirname $(realpath {output.csv}))
            echo "Running rG4detector detection on batch {wildcards.batch_id}"
            cd software/rG4detector/code
            python predict_fasta.py -d -f $IN -o $OUT
            cd ../../../
        }} 2>&1 | tee {log}
        """


def get_prediction_batches(wc):
    """Collect prediction CSV files from checkpoint."""
    checkpoints.split_fasta_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample)
    return expand(
        "results/{sample}/predictions/{batch_id}/rG4detector_prediction.csv",
        sample=wc.sample,
        batch_id=batch_ids,
    )


rule merge_rg4_predictions:
    """
Merge per-batch predictions into single output, preserving original sequence order.
Uses sequence ID matching for robustness.
"""
    input:
        manifest=lambda wc: checkpoints.split_fasta_batches.get(
            sample=wc.sample
        ).output.manifest,
        predictions=get_prediction_batches,
    output:
        csv="results/{sample}/rG4detector_prediction.csv",
    log:
        "logs/{sample}/merge_predictions.log",
    benchmark:
        "benchmarks/{sample}/merge_predictions.tsv"
    resources:
        runtime=10,
        mem_mb=1024 * 2,
        cpus_per_task=1,
    shell:
        """
        {{
            python workflow/scripts/merge_rg4_results.py {input.manifest} {output.csv} {input.predictions}
        }} 2>&1 | tee {log}
        """


def get_detection_batches(wc):
    """Collect detection CSV files from checkpoint."""
    checkpoints.split_fasta_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample)
    return expand(
        "results/{sample}/detections/{batch_id}/detection.csv",
        sample=wc.sample,
        batch_id=batch_ids,
    )


rule merge_rg4_detections:
    """
Merge per-batch detections into single output, preserving original sequence order.
Uses sequence ID matching for robustness.
"""
    input:
        manifest=lambda wc: checkpoints.split_fasta_batches.get(
            sample=wc.sample
        ).output.manifest,
        detections=get_detection_batches,
    output:
        csv="results/{sample}/detection.csv",
    log:
        "logs/{sample}/merge_detections.log",
    benchmark:
        "benchmarks/{sample}/merge_detections.tsv"
    resources:
        runtime=20,
        mem_mb=1024 * 30,  # Detection files are larger, so we allow more memory for merging
        cpus_per_task=1,
    shell:
        """
        {{
            python workflow/scripts/merge_rg4_results.py {input.manifest} {output.csv} {input.detections}
        }} 2>&1 | tee {log}
        """
