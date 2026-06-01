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
            git clone -b transcripts https://github.com/dgruano/g4Discovery.PanSN.git \
                software/g4Discovery.PanSN 2> {log}
        fi
        echo "g4Discovery cloned at software/g4Discovery.PanSN" >> {log}
        """


rule run_g4discovery_transcripts:
    """
Run g4DiscoveryTranscripts on a batch FASTA with multiple sequences in a single call.
This rule requires a new version of g4Discovery wrapper script and a new pqsfinder script. These handle:
    - (1) Non gzipped FASTA input (otherwise it requires us to compress to decompress later)
    - (2) multiple sequences in one FASTA file
Batch level BED-files are produced, and later merged in a separate rule
"""
    input:
        fa_gz="results/{sample}/batches/batch_{batch_id}.fa.gz",
        script="software/g4Discovery.PanSN/src/g4DiscoveryTranscripts.py",  # proxy: repo is cloned from Transcripts version
    output:
        gz="results/{sample}/g4discovery/transcripts_batch_results/batch_{batch_id}.bed.gz",
    log:
        "logs/{sample}/run_g4discovery_transcripts_{batch_id}.log",
    conda:
        "../envs/g4discovery.yaml"
    params:
        output_dir=lambda wc, output: subpath(output.gz, parent=True),
        pqs_score=config.get("pqs_score_threshold", 30),
        g4hunter_score=config.get("g4hunter_score_threshold", 1),

    threads: 1
    resources:
        # OK for batches of 1000 sequences
        mem_mb=2048,
        runtime=10,
    shell:
        """
        {{
            echo "Creating output directory: {params.output_dir}" >&2
            mkdir -p {params.output_dir}
            FA=$(realpath {input.fa_gz})
            OUT=$(realpath -m {output.gz})
            SCRIPT=$(realpath software/g4Discovery.PanSN/src/g4DiscoveryTranscripts.py)

            echo "Input FASTA: $FA" >&2
            echo "Output file: $OUT" >&2
            echo "Script: $SCRIPT" >&2
            echo "Parameters: pqs_score={params.pqs_score}, g4hunter_score={params.g4hunter_score}" >&2
            echo "Starting g4DiscoveryTranscripts..." >&2

            cd software
            python3 "$SCRIPT" \
                -fa "$FA" \
                -o "$OUT" \
                -ps {params.pqs_score} \
                -hs {params.g4hunter_score} \
                -m

            echo "g4DiscoveryTranscripts completed successfully" >&2
        }} 2>&1 | tee {log}
        """


def get_g4discovery_transcripts_batches(wc):
    checkpoints.split_fasta_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample)
    return expand(
        "results/{sample}/g4discovery/transcripts_batch_results/batch_{batch_id}.bed.gz",
        sample=wc.sample,
        batch_id=batch_ids,
    )


rule merge_g4discovery_transcripts:
    """
Merge batch-level g4DiscoveryTranscripts BED files into a single BED file per sample.
"""
    input:
        gz=get_g4discovery_transcripts_batches,
    output:
        gz="results/{sample}/g4discovery/{sample}.g4Discovery.transcripts.bed.gz",
    log:
        "logs/{sample}/merge_g4discovery_transcripts.log",
    conda:
        "../envs/g4discovery.yaml"
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    shell:
        """
        {{
            for f in {input.gz}; do
                zcat "$f" 2>/dev/null || true
            done | gzip -c > {output.gz}
        }} 2>&1 | tee {log}
        """


rule bedtools_intersect:
    input:
        left="results/{sample}/g4discovery/{sample}.g4Discovery.transcripts.bed.gz",  # G4Discovery
        right="results/{sample}/{sample}.rG4detector.peaks.bed"  # rG4detector
    output:
        "results/{sample}/{sample}.g4Discovery.rG4detector.intersect.bed"
    params:
        ## Add optional parameters
        extra="-loj"
    log:
        "logs/intersect/{sample}.log"
    wrapper:
        "v9.4.1/bio/bedtools/intersect"
