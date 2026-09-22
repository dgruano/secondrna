# Secondarna

This is a Snakemake pipeline that computes and integrates results from several
prediction tools of secondary and tertiary RNA structures.

## Run the pipeline

Add each FASTA basename to `samples` in `config/config.yaml`, with the matching
file under `resources/`, then run:

```bash
conda activate snakemake
snakemake --profile profiles/default
```

For a one-off sample list without editing the config file, use a temporary
config override:

```bash
snakemake --profile profiles/default --config samples=my_sample
```

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
