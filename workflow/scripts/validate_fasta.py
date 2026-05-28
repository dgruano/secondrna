#!/usr/bin/env python
"""
Validate input FASTA file for correct format and RNA characters.
Exit codes: 0 = valid, 1 = invalid format, 2 = file not found
"""

import re
import sys

input_file = snakemake.input[0]
marker_file = snakemake.output.marker
log_file = snakemake.log[0]


def validate_fasta(filepath):
    """Check FASTA format, non-empty sequences, valid RNA characters."""
    try:
        with open(filepath, "r") as f:
            content = f.read().strip()
    except FileNotFoundError:
        return False, "File not found"

    if not content:
        return False, "File is empty"

    lines = content.split("\n")
    if not any(line.startswith(">") for line in lines):
        return False, "No FASTA headers found (lines should start with '>')"

    current_header = None
    sequence_count = 0
    valid_chars = set("ACGUTNacgutn")  # Try with T too

    for i, line in enumerate(lines, 1):
        if line.startswith(">"):
            if current_header is None or any(c for c in (current_header or "")):
                current_header = line[1:].strip()
                sequence_count += 1
            if not current_header:
                return False, f"Empty header at line {i}"
        else:
            if not line.strip():
                continue
            if not all(c in valid_chars for c in line):
                invalid = set(line) - valid_chars
                return False, f"Invalid RNA character(s) at line {i}: {invalid}"

    if sequence_count == 0:
        return False, "No sequences found in FASTA file"

    return True, f"Valid FASTA with {sequence_count} sequences"


# Validate
is_valid, message = validate_fasta(input_file)

# Log result
with open(log_file, "w") as f:
    f.write(message + "\n")

if is_valid:
    # Create marker file
    open(marker_file, "w").close()
    print(message)
    sys.exit(0)
else:
    print(f"FASTA validation failed: {message}", file=sys.stderr)
    sys.exit(1)
