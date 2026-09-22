import os


rule scanfold2_install:
    output:
        "software/ScanFold2.0/installed.txt",
    log:
        "logs/scanfold2/install.log",
    localrule: True
    shell:
        """
        {{
            cd software
            if [ -d "ScanFold2.0/.git" ]; then
                echo "ScanFold2.0 repo already exists, skipping clone"
            else
                git clone https://github.com/moss-lab/ScanFold2.0.git
            fi
            if conda env create -f ScanFold2.0/environment.yml; then
                echo "Conda environment created successfully"
            else
                echo "Conda environment may already exist, attempting update"
                conda env update -f ScanFold2.0/environment.yml --prune
            fi
            touch "ScanFold2.0/installed.txt"
        }} >{log} 2>&1
        """


rule strip_gencode_headers:
    input:
        "resources/{sample}.fa",
    output:
        "resources/{sample}.simple.fa",
    log:
        "logs/{sample}/strip_headers.log",
    threads: 1
    resources:
        mem_mb=2048,
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    header = line[1:].split("|")[0].strip()
                    outfile.write(f">{header}\n")
                else:
                    outfile.write(line)


rule scanfold2_run:
    input:
        flag="software/ScanFold2.0/installed.txt",
        fasta="resources/{sample}.fa",
    output:
        out="results/{sample}/scanfold2/ScanFold_run.log",
        # TODO: Change to sentinel file
    log:
        "logs/{sample}/scanfold2/run.log",
    benchmark:
        "benchmarks/{sample}/scanfold2/run.txt"
    conda:
        "../envs/scanfold2.yaml"
    threads: 1
    resources:
        runtime=1440,
        mem_mb=2048,
    params:
        window=120,
        step=1,
        temperature=37,
        shuffle="mono",
        folder=lambda wc, output: os.path.dirname(output.out),
    shell:
        """
        {{
        input_fasta="$(realpath {input.fasta})"
        output_folder=$(realpath {params.folder})
        cd software/ScanFold2.0
        echo "Running ScanFold2.0 with the following parameters:"
        echo "------------------------------------------------------------"
        echo "Input FASTA: $input_fasta"
        echo "Window size: {params.window}"
        echo "Step size: {params.step}"
        echo "Temperature: {params.temperature}°C"
        echo "Shuffle method: {params.shuffle}"
        echo "Output folder: $output_folder"
        echo "------------------------------------------------------------"

        python ScanFold2.0.py $input_fasta\
        -w {params.window}\
        -s {params.step}\
        -t {params.temperature}\
        --shuffle {params.shuffle}\
        --folder $output_folder
        echo "ScanFold2.0 run completed. Results saved to $output_folder"
        }} >{log} 2>&1
        """


def get_scanfold_batches(wc):
    """Collect per-batch ScanFold log files after batch splitting completes."""
    checkpoints.scanfold_split_fasta_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample, subdir="scanfold_batches")
    return expand(
        "results/{sample}/scanfold2/batch_{batch_id}/ScanFold_run.log",
        sample=wc.sample,
        batch_id=batch_ids,
    )


use rule scanfold2_run as scanfold2_run_batch with:
    input:
        flag="software/ScanFold2.0/installed.txt",
        fasta="results/{sample}/scanfold_batches/batch_{batch_id}.fa",
    output:
        out="results/{sample}/scanfold2/batch_{batch_id}/ScanFold_run.log",
        # TODO: Change to sentinel file
    log:
        "logs/{sample}/scanfold2/batch_{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/scanfold2/batch_{batch_id}.tsv"
    wildcard_constraints:
        batch_id=r"\d+",
