rule validate_input:
    """
    Validate input FASTA file format and content.
    """
    input:
        "resources/input.fa"
    output:
        marker = "results/.input_validated"
    log:
        "logs/validate_input.log"
    script:
        "../scripts/validate_fasta.py"
