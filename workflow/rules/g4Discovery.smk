_ASSEMBLY_PREFIX = "GCA_000001405.15_GRCh38_no_alt_analysis_set"


rule install_g4discovery:
    """
Clone the g4Discovery.PanSN Python/R tool.
No compilation needed; produces the main Python script.
"""
    output:
        script="software/g4Discovery.PanSN/src/g4Discovery.py",
    log:
        "logs/upstream/install_g4discovery.log",
    conda:
        "../envs/g4discovery.yaml"
    shell:
        """
        mkdir -p software
        if [ ! -d software/g4Discovery.PanSN/.git ]; then
            git clone https://github.com/saswat-km/g4Discovery.PanSN.git \
                software/g4Discovery.PanSN 2> {log}
        fi
        echo "g4Discovery cloned at software/g4Discovery.PanSN" >> {log}
        """


use rule compress_batch as compress_resources_fasta with:
    input:
        fa="resources/{sample}.fa",
    output:
        fa_gz="resources/{sample}.fa.gz",
    log:
        "logs/compress_resources_fasta_{sample}.log",


rule run_g4discovery:
    """
Run g4Discovery on the GRCh38 assembly with default settings.

Outputs a gzipped BED file containing G4 positions on both strands,
with both pqsfinder and G4Hunter scores.

Default thresholds: pqsfinder >= 40, |G4Hunter| >= 1.5, tetrads >= 3.
"""
    input:
        fa_gz="resources/{sample}.fa.gz",
        script="software/g4Discovery.PanSN/src/g4Discovery.py",
    output:
        gz="results/{sample}/g4discovery/{sample}.g4Discovery.bed.gz",
    log:
        "logs/{sample}/run_g4discovery.log",
    conda:
        "../envs/g4discovery.yaml"
    params:
        output_dir = lambda wc, output: subpath(output.gz, parent=True)
    threads: 1
    resources:
        mem_mb=32768,
        runtime=2880,  # 48 h — whole-genome G4 prediction is very slow
    shell:
        """
        {{
            mkdir -p {params.output_dir}
            FA=$(realpath {input.fa_gz})
            OUT=$(realpath -m {output.gz})
            SCRIPT=$(realpath {input.script})

            cd software
            python3 $SCRIPT \
                -fa "$FA" \
                -o "$OUT"
        }} 2>&1 | tee {log}
        """


use rule run_g4discovery as run_g4discovery_parallel with:
    input:
        fa_gz="results/{sample}/batches/batch_{batch_id}.fa.gz",
        script="software/g4Discovery.PanSN/src/g4Discovery.py",
    output:
        gz="results/{sample}/g4discovery/batch_results/batch_{batch_id}.bed.gz",
    log:
        "logs/{sample}/run_g4discovery_batch_{batch_id}.log",


def get_g4discovery_batches(wc):
    checkpoints.split_fasta_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample)
    return expand(
        "results/{sample}/g4discovery/batch_results/batch_{batch_id}.bed.gz",
        sample=wc.sample,
        batch_id=batch_ids,
    )


rule merge_g4discovery_batches:
    """
Merge batch-level g4Discovery BED files into a single BED file per sample.
Preserves original sequence order and strand information.
"""
    input:
        gz=get_g4discovery_batches,
    output:
        gz="results/{sample}/g4discovery/{sample}.g4Discovery.bed_merged.gz",
    log:
        "logs/{sample}/merge_g4discovery_batches.log",
    conda:
        "../envs/g4discovery.yaml"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    shell:
        """
        {{
            zcat {input.gz} | gzip -c > {output.gz}
        }} 2>&1 | tee {log}
        """
