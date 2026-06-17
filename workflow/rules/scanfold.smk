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
        runtime="2d",
        mem_mb=10000,
        slurm_partition="gpu",
        gres="gpu:nvidia_h100_nvl_1g.12gb",
        ntasks_per_gpu=0,
        cpus_per_gpu=1,
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


SCANFOLD_RETRY_SUBBATCH_SIZE = config.get("scanfold_retry_subbatch_size", 100)
SCANFOLD_OVERSIZED_MAX_LEN = config.get("scanfold_oversized_max_len", 20_000)


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


# ---------------------------------------------------------------------------
# Retry workflow: auto-detect incomplete batches and re-run in sub-batches
# ---------------------------------------------------------------------------


checkpoint scanfold_split_retry_batches:
    """Scan all original batches for incomplete sequences and split into sub-batches."""
    input:
        "results/{sample}/scanfold_batches/batch_manifest.txt",
    output:
        manifest="results/{sample}/scanfold_retry_batches/batch_manifest.txt",
        batches=directory("results/{sample}/scanfold_retry_batches/"),
        oversized_manifest="results/{sample}/scanfold_oversized_batches/batch_manifest.txt",
        oversized_batches=directory("results/{sample}/scanfold_oversized_batches/"),
        skip_tsv="results/{sample}/scanfold_oversized_skip.tsv",
    log:
        "logs/{sample}/scanfold2/split_retry_batches.log",
    benchmark:
        "benchmarks/{sample}/scanfold2/split_retry_batches.tsv"
    resources:
        runtime=10,
        mem_mb=2048,
        cpus_per_task=1,
    params:
        subbatch_size=SCANFOLD_RETRY_SUBBATCH_SIZE,
        max_seq_len=SCANFOLD_OVERSIZED_MAX_LEN,
    shell:
        """
        {{
        python workflow/scripts/check_scanfold_partial.py \
            --sample {wildcards.sample} \
            --subbatch-size {params.subbatch_size} \
            --max-seq-len {params.max_seq_len} \
            --write-fasta \
            --out-dir results/{wildcards.sample}/scanfold_retry_batches \
            --oversized-dir results/{wildcards.sample}/scanfold_oversized_batches \
            --results-dir results
        }} >{log} 2>&1
        """


def get_scanfold_retry_gpu_batches(wc):
    """Collect retry sub-batch log files after the retry split checkpoint resolves."""
    checkpoints.scanfold_split_retry_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample, subdir="scanfold_retry_batches")
    return expand(
        "results/{sample}/scanfold2_gpu/batch_{batch_id}/ScanFold_run.log",
        sample=wc.sample,
        batch_id=batch_ids,
    )


use rule scanfold2_run_gpu_batch as scanfold2_run_gpu_retry_batch with:
    input:
        flag="software/ScanFold2.0/installed.txt",
        fasta="results/{sample}/scanfold_retry_batches/batch_{batch_id}.fa",
    output:
        out="results/{sample}/scanfold2_gpu/batch_{batch_id}/ScanFold_run.log",
    log:
        "logs/{sample}/scanfold2_gpu/retry_{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/scanfold2_gpu/retry_{batch_id}.tsv"


rule scanfold2_retry_all:
    input:
        lambda wc: get_scanfold_retry_gpu_batches(
            type("WC", (), {"sample": "gencode.v47.repeat.simple"})()
        ),


# ---------------------------------------------------------------------------
# Oversized fallback: single-sequence highmem CPU jobs
# ---------------------------------------------------------------------------


rule scanfold2_run_oversized:
    """Process a single ultra-long sequence on a high-memory CPU node."""
    input:
        flag="software/ScanFold2.0/installed.txt",
        fasta="results/{sample}/scanfold_oversized_batches/batch_{batch_id}.fa",
    output:
        out="results/{sample}/scanfold2_oversized/batch_{batch_id}/ScanFold_run.log",
    log:
        "logs/{sample}/scanfold2_oversized/{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/scanfold2_oversized/{batch_id}.tsv"
    conda:
        "scanfold2"
    resources:
        runtime="7d",
        mem_mb=50_000,  # Got OOM with a 370k nt sequence and 10 GB RAM 12 GB VRAM. 50 should do
        cpus_per_task=1,
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
        python ScanFold2.0.py $input_fasta \
            -w {params.window} -s {params.step} \
            -t {params.temperature} --shuffle {params.shuffle} \
            --folder $output_folder
        }} >{log} 2>&1
        """


def get_scanfold_oversized_batches(wc):
    checkpoints.scanfold_split_retry_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample, subdir="scanfold_oversized_batches")
    return expand(
        "results/{sample}/scanfold2_oversized/batch_{batch_id}/ScanFold_run.log",
        sample=wc.sample,
        batch_id=batch_ids,
    )


rule scanfold2_oversized_all:
    input:
        lambda wc: get_scanfold_oversized_batches(
            type("WC", (), {"sample": "gencode.v47.repeat.simple"})()
        ),


# ---------------------------------------------------------------------------
# Completeness validation: confirm every sequence has a .no_filter.ct output
# ---------------------------------------------------------------------------


def get_all_scanfold_batch_outputs(wc):
    """Aggregate log files from all three ScanFold batch stages."""
    return (
        get_scanfold_gpu_batches(wc)
        + get_scanfold_retry_gpu_batches(wc)
        + get_scanfold_oversized_batches(wc)
    )


rule scanfold2_check_complete:
    """Validate all original sequences have .no_filter.ct outputs across all batch stages."""
    input:
        batch_logs=get_all_scanfold_batch_outputs,
        manifest="results/{sample}/scanfold_batches/batch_manifest.txt",
    output:
        report="results/{sample}/scanfold2/completeness_check.tsv",
    log:
        "logs/{sample}/scanfold2/check_complete.log",
    benchmark:
        "benchmarks/{sample}/scanfold2/check_complete.tsv"
    resources:
        runtime=10,
        mem_mb=2048,
    params:
        results_dir="results",
    script:
        "../scripts/check_scanfold_complete.py"


rule scanfold2_validate_all:
    input:
        expand(
            "results/{sample}/scanfold2/completeness_check.tsv",
            sample=["gencode.v47.repeat.simple"],
        ),
