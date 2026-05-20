# Default batch size (can be overridden via --config batch_size=N)
BATCH_SIZE = config.get("batch_size", 1000)  # 50 for testing, 1000 for production benchmarking

# Global wildcard constraints
wildcard_constraints:
    sample="[a-zA-Z0-9_\-\.]+"  # Disallow forward slashes

def calculate_memory_for_batch(batch_size):
    """
    Calculate memory needed for a batch based on actual failures.
    Batch size of 1000 takes:
      - ~20 min
      - ~10 GB of RAM

    Minimum RAM is around 2 GB due to model overhead
    So we can estimate:
      - Base memory: 2 GB (2000 MB)
      - Additional memory per sequence: (10 GB - 2 GB) / 1000 = 8 MB per sequence, rounded up to 10 MB for safety
    """
    memory_per_seq_mb = 10
    total_mb = int(batch_size * memory_per_seq_mb) + 2000  # Base memory + per-sequence memory

    # Cap at reasonable maximum (e.g., 250 GB)
    return min(total_mb, 250000)


def calculate_runtime_for_batch(batch_size):
    """
    Estimate runtime for a batch
    Based on observed performance:
      - Batch size of 1000 takes ~20 minutes
    So we can estimate:
      - Base time: 5 minutes (for setup and overhead)
    """
    base_time = 5  # minutes
    time_per_seq = 20 / 1000  # minutes per sequence
    total_time = int(base_time + (batch_size * time_per_seq))
    return total_time


rule benchmark_batch_sizes:
    """
    Benchmark different batch sizes to determine optimal configuration.
    Tests batch sizes (100, 500, 1000, 2000) and generates performance analysis.
    """
    input:
        fasta = "resources/input_{sample}.fa"
    output:
        analysis = "results/{sample}/batch_performance_analysis.txt"
    log:
        "logs/{sample}/benchmark_batch_sizes.log"
    conda:
        "test_rg4"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=1
    benchmark:
        "benchmarks/{sample}/benchmark_batch_sizes.tsv"
    shell:
        """
        {{
            python workflow/scripts/benchmark_batch_sizes.py {input.fasta} {output.analysis}
        }} 2>&1 | tee {log}
        """


def get_num_batches(input_fasta):
    """Calculate number of batches needed."""
    num_sequences = 0
    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                num_sequences += 1
    #Num seqs = grep -c "^>" input_fasta
    return (num_sequences + BATCH_SIZE - 1) // BATCH_SIZE


def get_batch_ids(sample):
    """
    Checkpoint function: determine batch IDs after split completes.
    Reads the manifest to find all batch IDs.

    NOTE: IDs are determined by us, so we could infer them directly
    instead of reading the manifest.
    """
    manifest_path = f"results/{sample}/batches/batch_manifest.txt"
    batch_ids = []
    try:
        with open(manifest_path) as f:
            for c, line in enumerate(f):
                if c == 0:
                    continue  # Skip header
                # Manifest format: batch_id\tseq_id (one per line)
                parts = line.strip().split("\t")
                if len(parts) >= 1:
                    batch_id = parts[0]
                    if batch_id not in batch_ids:
                        batch_ids.append(batch_id)
        return batch_ids
    except FileNotFoundError:
        return []


checkpoint split_fasta_batches:
    """
    Split input FASTA into fixed-size batches, preserving sequence order.
    Creates batch FASTA files and manifest with sequence ID mapping.
    """
    input:
        fasta = "resources/input_{sample}.fa"
    output:
        marker = "results/{sample}/.batches_created",
        manifest = "results/{sample}/batches/batch_manifest.txt",
        batches = directory("results/{sample}/batches/")
    params:
        batch_size = BATCH_SIZE
    log:
        "logs/{sample}/split_fasta_batches.log"
    conda:
        "test_rg4"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=1
    benchmark:
        "benchmarks/{sample}/split_fasta_batches.tsv"
    shell:
        """
        {{
            python workflow/scripts/split_fasta.py {input.fasta} {output.batches} {params.batch_size}
            touch {output.marker}
        }} 2>&1 | tee {log}
        """


rule rg4_predict_batch:
    """
    Run rG4detector prediction on individual batch FASTA.
    Memory is dynamically allocated based on batch_size.
    """
    input:
        batch = "results/{sample}/batches/batch_{batch_id}.fa",
        setup = "results/.rg4_setup_done"
    output:
        csv = "results/{sample}/predictions/{batch_id}/rG4detector_prediction.csv"
    log:
        "logs/{sample}/predict_batch_{batch_id}.log"
    conda:
        "test_rg4"
    resources:
        runtime=lambda wc: calculate_runtime_for_batch(BATCH_SIZE),  # Scale runtime with batch size (e.g., 20 min per 1000 sequences)
        mem_mb=lambda wc: calculate_memory_for_batch(BATCH_SIZE),
        cpus_per_task=2
    benchmark:
        "benchmarks/{sample}/predict_batch_{batch_id}.tsv"
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
        batch = "results/{sample}/batches/batch_{batch_id}.fa",
        setup = "results/.rg4_setup_done"
    output:
        csv = "results/{sample}/detections/{batch_id}/detection.csv"
    log:
        "logs/{sample}/detect_batch_{batch_id}.log"
    conda:
        "test_rg4"
    resources:
        runtime=lambda wc: calculate_runtime_for_batch(BATCH_SIZE),  # Scale runtime with batch size (e.g., 20 min per 1000 sequences)
        mem_mb=lambda wc: calculate_memory_for_batch(BATCH_SIZE),
        cpus_per_task=2
    benchmark:
        "benchmarks/{sample}/detect_batch_{batch_id}.tsv"
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
        manifest = lambda wc: checkpoints.split_fasta_batches.get(sample=wc.sample).output.manifest,
        predictions = get_prediction_batches
    output:
        csv = "results/{sample}/rG4detector_prediction.csv"
    log:
        "logs/{sample}/merge_predictions.log"
    resources:
        runtime=10,
        mem_mb=4096,
        cpus_per_task=1
    benchmark:
        "benchmarks/{sample}/merge_predictions.tsv"
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
        manifest = lambda wc: checkpoints.split_fasta_batches.get(sample=wc.sample).output.manifest,
        detections = get_detection_batches
    output:
        csv = "results/{sample}/detection.csv"
    log:
        "logs/{sample}/merge_detections.log"
    resources:
        runtime=10,
        mem_mb=4096,
        cpus_per_task=1
    benchmark:
        "benchmarks/{sample}/merge_detections.tsv"
    shell:
        """
        {{
            python workflow/scripts/merge_rg4_results.py {input.manifest} {output.csv} {input.detections}
        }} 2>&1 | tee {log}
        """


rule cleanup_batch_files:
    """
    Optional: Delete intermediate batch files, keep only manifest and final results.
    """
    input:
        checkpoint = lambda wc: checkpoints.split_fasta_batches.get(sample=wc.sample).output.marker,
        predictions = "results/{sample}/rG4detector_prediction.csv",
        detections = "results/{sample}/detection.csv"
    output:
        cleanup_marker = "results/{sample}/.cleanup_done"
    log:
        "logs/{sample}/cleanup_batch_files.log"
    benchmark:
        "benchmarks/{sample}/cleanup_batch_files.tsv"
    shell:
        """
        {{
            echo "Cleaning up intermediate batch files..."
            rm -f results/{wildcards.sample}/batches/batch_*.fa
            rm -f results/{wildcards.sample}/predictions/{batch_id}/rG4detector_prediction.csv
            rm -f results/{wildcards.sample}/detections/{batch_id}/detection.csv
            echo "Kept: results/{wildcards.sample}/batches/batch_manifest.txt (audit trail)"
            touch {output.cleanup_marker}
        }} 2>&1 | tee {log}
        """
