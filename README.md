# Secondarna

This is a Snakemake pipeline that computes and integrates results from several
prediction tools of secondary and tertiary RNA structures.

## Run the pipeline

Place one or more input FASTA files under `resources/`, using the sample name as
the basename (for example, `resources/gencode.v47.fa`). The release defaults
are in `config/config.default.yaml`; copy that file for local overrides and set
the `samples` list, then run:

```bash
conda activate snakemake
snakemake --cores 4 --use-conda
```

The default target runs rG4detector and ScanFold2 under `results/<sample>/`.
Before ScanFold, records longer than `scanfold_oversized_max_len` are removed
from its input and listed in `results/<sample>/scanfold_excluded.tsv`, with the
source FASTA checksum and threshold. Those records have no ScanFold predictions.
For a small local check:

```bash
conda activate snakemake
snakemake --cores 4 --use-conda --config samples=my_sample
```

Set `scanfold_batch_size` independently from the rG4detector `batch_size`.
For example, `--config scanfold_batch_size=1` schedules one ScanFold sequence
per job; larger values reduce scheduler overhead.

For a one-off sample list without editing the config file, use a temporary
config override:

```bash
snakemake --cores 4 --use-conda --config samples=my_sample
```

The ScanFold rules are also available as explicit targets. For example:

```bash
snakemake --cores 4 --use-conda \
	results/my_sample/scanfold2/.complete
```

## Reproducible ScanFold2 dependency

ScanFold runs require upstream commit
`c5cc73291fa5a1c06bf5b3c8367a396a3238141c` plus the tracked
[`vendor/scanfold2/multifasta.patch`](vendor/scanfold2/multifasta.patch).
The patch prevents output overwrites when processing multi-record FASTA files.
Installation pins the upstream revision and applies the patch; execution checks
source and model hashes, including pre-existing installations. FASTA IDs must be
unique and safe as filenames.

See [installation, output naming and tests](vendor/scanfold2/README.md) and
[historical provenance and validation](vendor/scanfold2/PROVENANCE.md). The
historical summary combines multiple result sources whose exact source revisions
are not fully established; the pinned release is the supported version for new
runs.

## Dev notes
After checking different batch sizes, I think that:
- Running time scales O(n), 20 min / 1000 seqs
- Memory requirement scales O(n), with ~ 2GB constant due to model loading
- 1000 sequences takes ~20 minutes for each script (predict/detect) and takes 10Gb of RAM
- I have not controlled for transcript size, but test sequences were 3800 nt on average

The outputs of scripts are:
- predict: <seq_id>,<score>
- detect: <seq_id>,<comma_separated_nucleotides>,<comma_separated_values>

## Notes
- rg4detector is trained on data from HeLa cells, validated agains mouse and Arabidopsis data.
It showed better performance than other methods, but lower than wrt the HeLa validation set.
- We do not know the proportion of lncRNAs (or G4 identified in lncRNAs) with the rG4-seq in the HeLa cells,
so we cannot be sure that what the model learnt is not biased towards mRNAs. This is intrinsically problematic
for my question, if the application is to compare lncRNAs and mRNAs.
