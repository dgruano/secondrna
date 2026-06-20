# GPU-accelerated ScanFold2.0 rules — kept for reference.
# The main pipeline now uses CPU rules in scanfold.smk.


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
        echo "Starting ScanFold2.0 (GPU) run for sample: {wildcards.sample}"
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
