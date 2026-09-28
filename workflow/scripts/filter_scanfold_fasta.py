"""Exclude long records before ScanFold and record exactly what was excluded."""

import hashlib
from pathlib import Path

from Bio import SeqIO

source = Path(snakemake.input[0])
limit = int(snakemake.params.max_len)
hasher = hashlib.sha256()
with source.open("rb") as stream:
    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
        hasher.update(chunk)
digest = hasher.hexdigest()

with open(snakemake.output.fasta, "w") as filtered, open(
    snakemake.output.excluded, "w"
) as excluded:
    excluded.write("seq_id\tseq_len\tsource_sha256\tmax_len\n")
    for record in SeqIO.parse(source, "fasta"):
        if len(record.seq) > limit:
            excluded.write(f"{record.id}\t{len(record.seq)}\t{digest}\t{limit}\n")
        else:
            SeqIO.write(record, filtered, "fasta")
