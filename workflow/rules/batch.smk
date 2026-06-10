# Default batch size (can be overridden via --config batch_size=N)
BATCH_SIZE = config.get(
    "batch_size", 1000
)  # 50 for testing, 1000 for production benchmarking


def get_num_batches(input_fasta):
    """Calculate number of batches needed."""
    num_sequences = 0
    with open(input_fasta) as f:
        for line in f:
            if line.startswith(">"):
                num_sequences += 1
    # Num seqs = grep -c "^>" input_fasta
    return (num_sequences + BATCH_SIZE - 1) // BATCH_SIZE


def get_batch_ids(sample, subdir="batches"):
    """
    Checkpoint function: determine batch IDs after split completes.
    Reads the manifest to find all batch IDs.

    NOTE: IDs are determined by us, so we could infer them directly
    instead of reading the manifest.
    """
    manifest_path = f"results/{sample}/{subdir}/batch_manifest.txt"
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
        fasta="resources/input_{sample}.fa",
    output:
        marker="results/{sample}/.batches_created",
        manifest="results/{sample}/batches/batch_manifest.txt",
        batches=directory("results/{sample}/batches/"),
    log:
        "logs/{sample}/split_fasta_batches.log",
    benchmark:
        "benchmarks/{sample}/split_fasta_batches.tsv"
    conda:
        "test_rg4"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=1,
    params:
        batch_size=BATCH_SIZE,
    shell:
        """
        {{
            python workflow/scripts/split_fasta.py {input.fasta} {output.batches} {params.batch_size}
            touch {output.marker}
        }} 2>&1 | tee {log}
        """


checkpoint scanfold_split_fasta_batches:
    input:
        fasta="resources/{sample}.fa",
    output:
        marker="results/{sample}/.scanfold_batches_created",
        manifest="results/{sample}/scanfold_batches/batch_manifest.txt",
        batches=directory("results/{sample}/scanfold_batches/"),
    log:
        "logs/{sample}/scanfold_split_fasta_batches.log",
    benchmark:
        "benchmarks/{sample}/scanfold_split_fasta_batches.tsv"
    conda:
        "test_rg4"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=1,
    params:
        batch_size=BATCH_SIZE,
    shell:
        """
        {{
            python workflow/scripts/split_fasta.py {input.fasta} {output.batches} {params.batch_size}
            touch {output.marker}
        }} 2>&1 | tee {log}
        """


rule compress_batch:
    """
Compress a batch FASTA to .fa.gz.
Downstream rules that need compressed batches declare batch_{batch_id}.fa.gz as input;
Snakemake resolves this rule automatically.
"""
    input:
        fa="results/{sample}/batches/batch_{batch_id}.fa",
    output:
        gz="results/{sample}/batches/batch_{batch_id}.fa.gz",
    resources:
        runtime=5,
        mem_mb=256,
        cpus_per_task=1,
    shell:
        "gzip -k {input.fa}"


rule cleanup_batch_files:
    """
Optional: Delete intermediate batch files, keep only manifest and final results.
"""
    input:
        checkpoint=lambda wc: checkpoints.split_fasta_batches.get(
            sample=wc.sample
        ).output.marker,
        predictions="results/{sample}/rG4detector_prediction.csv",
        detections="results/{sample}/detection.csv",
    output:
        cleanup_marker="results/{sample}/.cleanup_done",
    log:
        "logs/{sample}/cleanup_batch_files.log",
    benchmark:
        "benchmarks/{sample}/cleanup_batch_files.tsv"
    shell:
        """
        {{
            echo "Cleaning up intermediate batch files..."
            rm -f results/{wildcards.sample}/batches/batch_*.fa results/{wildcards.sample}/batches/batch_*.fa.gz
            rm -f results/{wildcards.sample}/predictions/{batch_id}/rG4detector_prediction.csv
            rm -f results/{wildcards.sample}/detections/{batch_id}/detection.csv
            echo "Kept: results/{wildcards.sample}/batches/batch_manifest.txt (audit trail)"
            touch {output.cleanup_marker}
        }} 2>&1 | tee {log}
        """
