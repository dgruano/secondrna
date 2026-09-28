import hashlib
import runpy
from pathlib import Path
from types import SimpleNamespace


def test_filter_records_and_report(tmp_path):
    source = tmp_path / "input.fa"
    source.write_text(">short\nACG\n>long\nACGT\n")
    fasta = tmp_path / "filtered.fa"
    excluded = tmp_path / "excluded.tsv"
    snakemake = SimpleNamespace(
        input=[str(source)],
        output=SimpleNamespace(fasta=str(fasta), excluded=str(excluded)),
        params=SimpleNamespace(max_len=3),
    )

    runpy.run_path(
        str(Path(__file__).parents[1] / "workflow/scripts/filter_scanfold_fasta.py"),
        init_globals={"snakemake": snakemake},
    )

    assert fasta.read_text() == ">short\nACG\n"
    assert excluded.read_text() == (
        "seq_id\tseq_len\tsource_sha256\tmax_len\n"
        f"long\t4\t{hashlib.sha256(source.read_bytes()).hexdigest()}\t3\n"
    )
