rule scanfold2_all:
    input:
        expand(
            "results/{sample}/scanfold2/{sample}.no_filter.ct",
            sample=["input_100"],
        ),


rule scanfold2_gpu_all:
    input:
        expand(
            "results/{sample}/scanfold2_gpu/{sample}.no_filter.ct",
            sample=["input_100"],
        ),


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
    threads: 52
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
        no_filter_ct="results/{sample}/scanfold2_gpu/{sample}.no_filter.ct",
        minus1_ct="results/{sample}/scanfold2_gpu/{sample}.minus_1.ct",
        minus2_ct="results/{sample}/scanfold2_gpu/{sample}.minus_2.ct",
        no_filter_dbn="results/{sample}/scanfold2_gpu/{sample}.no_filter.dbn",
        minus1_dbn="results/{sample}/scanfold2_gpu/{sample}.minus_1.dbn",
        minus2_dbn="results/{sample}/scanfold2_gpu/{sample}.minus_2.dbn",
    log:
        "logs/{sample}/scanfold2_gpu/run.log",
    benchmark:
        "benchmarks/{sample}/scanfold2_gpu/run.txt"
    conda:
        "../envs/scanfold2_gpu.yaml"
    resources:
        runtime=240,
        mem_mb=8192,
        slurm_partition="gpu",
        gres="gpu:1",
        ntasks_per_gpu=0,
        cpus_per_gpu=4,
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
