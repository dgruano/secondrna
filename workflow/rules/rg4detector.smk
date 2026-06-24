# Configuration defaults
from pathlib import Path


rule rg4_setup:
    """
Clone rG4detector repository, create conda environment, and verify installation.
"""
    output:
        marker="results/.rg4_setup_done",
    log:
        "logs/rg4detector/rg4_setup.log",
    resources:
        runtime=10,
        mem_mb=4096,
        cpus_per_task=1,
    shell:
        """
        exec 2> {log}

        # Clone repository if not already present
        if [ ! -d "rG4detector" ]; then
            git clone https://github.com/OrensteinLab/rG4detector.git
        fi

        # Create conda environment
        conda env create -f workflow/envs/rg4detector.yaml --force --quiet 2>&1 || true

        # Verify installation by testing import
        conda init bash
        conda activate rg4detector
        python -c "import tensorflow; import keras; import numpy; import pandas; print('rG4detector environment verified')" 2>&1

        # Create marker file
        touch {output.marker}
        """


rule rg4_predict:
    """
Run rG4detector in prediction mode on input FASTA.
"""
    input:
        #validated = "results/{sample}/.input_validated",
        marker="results/.rg4_setup_done",
        input_fasta="resources/{sample}.fa",
    output:
        "results/{sample}/rg4detector/rG4detector_prediction.csv",
    log:
        "logs/{sample}/rg4detector/rg4_predict.log",
    benchmark:
        "benchmarks/{sample}/rg4detector/rg4_predict.tsv"
    conda:
        #"../envs/rg4detector.yaml"
        "test_rg4"
    resources:
        runtime=120,
        mem_mb=1024 * 12,
        cpus_per_task=2,
    shell:
        """
        {{
            IN=$(realpath resources/{wildcards.sample}.fa)
            OUT=$(realpath results/{wildcards.sample}/rg4detector/)
            mkdir -p $OUT
            echo "Running rG4detector prediction on $IN, outputting to $OUT"
            cd software/rG4detector/code
            python predict_fasta.py -f $IN -o $OUT
            cd ..
        }} 2>&1 | tee {log}
        """


rule rg4_detect:
    """
Run rG4detector in detection mode on input FASTA.
"""
    input:
        marker="results/.rg4_setup_done",
        input_fasta="resources/{sample}.fa",
    output:
        "results/{sample}/rg4detector/detection.csv",
    log:
        "logs/{sample}/rg4detector/rg4_detect.log",
    benchmark:
        "benchmarks/{sample}/rg4detector/rg4_detect.tsv"
    conda:
        #"../envs/rg4detector.yaml"
        "test_rg4"
    resources:
        runtime=120,
        mem_mb=1024 * 12,
        cpus_per_task=2,
    shell:
        """
        {{
            IN=$(realpath resources/{wildcards.sample}.fa)
            OUT=$(realpath results/{wildcards.sample}/rg4detector/)
            mkdir -p $OUT
            echo "Running rG4detector detection on $IN, outputting to $OUT"
            cd software/rG4detector/code
            python predict_fasta.py -d -f $IN -o $OUT
            cd ..
        }} 2>&1 | tee {log}
        """
