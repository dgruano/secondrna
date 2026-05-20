"""Utility rules for common bioinformatics operations."""

rule get_first_n_sequences:
    """Extract the first n sequences from a FASTA file.

    Usage:
        snakemake get_first_n_sequences --config fasta_input=input.fa fasta_output=output.fa n_sequences=100
    """
    output:
        fasta = "results/first_n_sequences.fa"  # Placeholder output path
    params:
        input_fasta = lambda wildcards: config.get("fasta_input", ""),
        n_sequences = lambda wildcards: config.get("n_sequences", 100)
    shell:
        """
        awk 'BEGIN {{count=0}} /^>/ {{if (count >= {params.n_sequences}) exit; count++}} {{print}}' \
            {params.input_fasta} > {output.fasta}
        """
