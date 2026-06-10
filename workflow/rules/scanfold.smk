def reduce_fasta_headers(input_file, output_file):
    """Read FASTA and reduce headers to first element before pipe."""


rule scanfold2_all:
    input:
        expand(
            "results/{sample}/scanfold2/{sample}.no_filter.ct",
            sample=["gencode.v47"],
        ),


rule scanfold2_summarize_all:
    input:
        expand("results/{sample}/scanfold2/scanfold_stats.tsv", sample=["gencode.v47"]),


rule scanfold2_summarize:
    """Compute per-transcript ScanFold2 stats and annotate with coding class."""
    input:
        parquet="results/scanfold_all.parquet",
        coding_ann="results/{sample}/{sample}.coding_annotation.tsv",
    output:
        "results/{sample}/scanfold2/scanfold_stats.tsv",
    log:
        "logs/{sample}/scanfold2/summarize.log",
    benchmark:
        "benchmarks/{sample}/scanfold2/summarize.tsv"
    resources:
        runtime=30,
        mem_mb=200000,
        nodelist="apollo01",
    script:
        "../scripts/scanfold_summarize.py"


rule scanfold2_extract_csv_tarball:
    input:
        folder="/mnt/cbib/LNClassifier/RNA_ScanFold2.0",
    output:
        tarball="ScanFold2.0.csv.tar.gz",
    log:
        "logs/scanfold2/extract_csv.log",
    threads: 46
    resources:
        mem_mb=8192,
    shell:
        """
        {{
        set -euo pipefail
        tmpdir=$(mktemp -d -p /scratch/dgarcia)
        trap 'rm -rf "$tmpdir"' EXIT

        extract_archive() {{
            local archive="$1"
            local dest="$2/$(basename "$archive" .tar.gz)"
            mkdir -p "$dest"
            tar -tzf "$archive" | grep -E '\\.csv$' | \
                tar -xzf "$archive" -C "$dest" --transform='s|.*/||' --files-from=-
        }}
        export -f extract_archive
        export DEST="$tmpdir"

        find {input.folder} -name '*.tar.gz' | \
            xargs -P {threads} -I % bash -c 'extract_archive "%" "$DEST"'

        tar -czf {output.tarball} -C "$tmpdir" .
        }} >{log} 2>&1
        """


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
        no_filter_ct="results/{sample}/scanfold2/{sample}.no_filter.ct",
        minus1_ct="results/{sample}/scanfold2/{sample}.minus_1.ct",
        minus2_ct="results/{sample}/scanfold2/{sample}.minus_2.ct",
        no_filter_dbn="results/{sample}/scanfold2/{sample}.no_filter.dbn",
        minus1_dbn="results/{sample}/scanfold2/{sample}.minus_1.dbn",
        minus2_dbn="results/{sample}/scanfold2/{sample}.minus_2.dbn",
    log:
        "logs/{sample}/scanfold2/run.log",
    benchmark:
        "benchmarks/{sample}/scanfold2/run.txt"
    conda:
        "scanfold2"
    threads: 4
    resources:
        runtime=1440,
        mem_mb=2048,
    params:
        window=120,
        step=1,
        temperature=37,
        shuffle="mono",
        folder=lambda wc, output: os.path.dirname(output.no_filter_ct),
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


rule scanfold2_run_gpu:
    input:
        flag="software/ScanFold2.0/installed.txt",
        fasta="resources/{sample}.fa",
    output:
        out="results/{sample}/scanfold2_gpu/ScanFold_run.log",
    log:
        "logs/{sample}/scanfold2_gpu/run.log",
    benchmark:
        "benchmarks/{sample}/scanfold2_gpu/run.txt"
    conda:
        "../envs/scanfold2_gpu.yaml"
    resources:
        runtime=1440,
        mem_mb=6000,
        slurm_partition="gpu",
        gres="gpu:nvidia_h100_nvl_1g.12gb",
        ntasks_per_gpu=0,
        cpus_per_gpu=4,
    params:
        window=120,
        step=1,
        temperature=37,
        shuffle="mono",
        folder=lambda wc, output: os.path.dirname(output.out),
    shell:
        """
        {{
        # GPU-related setup
        echo "Starting lncRNA-BERT training"
        echo "GPU Info:"
        echo "Available physical GPUs:"
        nvidia-smi --query-gpu=index,platform.module_id,name,driver_version,memory.total,compute_cap,mig.mode.current --format=csv
        echo "Available MIG devices:"
        nvidia-smi -L | grep -i mig || echo "No MIG devices found"

        input_fasta="$(realpath {input.fasta})"
        output_folder=$(realpath {params.folder})
        cd software/ScanFold2.0
        echo "Running ScanFold2.0 (GPU) with the following parameters:"
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
        echo "ScanFold2.0 (GPU) run completed. Results saved to $output_folder"
        }} >{log} 2>&1
        """


def get_scanfold_gpu_batches(wc):
    """Collect per-batch ScanFold GPU log files from checkpoint."""
    checkpoints.scanfold_split_fasta_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample, subdir="scanfold_batches")
    return expand(
        "results/{sample}/scanfold2_gpu/batch_{batch_id}/ScanFold_run.log",
        sample=wc.sample,
        batch_id=batch_ids,
    )


use rule scanfold2_run_gpu as scanfold2_run_gpu_batch with:
    input:
        flag="software/ScanFold2.0/installed.txt",
        fasta="results/{sample}/scanfold_batches/batch_{batch_id}.fa",
    output:
        out="results/{sample}/scanfold2_gpu/batch_{batch_id}/ScanFold_run.log",
    log:
        "logs/{sample}/scanfold2_gpu/batch_{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/scanfold2_gpu/batch_{batch_id}.tsv"


rule scanfold2_gpu_all:
    input:
        lambda wc: get_scanfold_gpu_batches(
            type("WC", (), {"sample": "gencode.v47.repeat.simple"})()
        ),
