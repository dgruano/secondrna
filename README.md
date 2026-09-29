# secondrna

Snakemake pipeline for predicting RNA G-quadruplexes with rG4detector and
RNA structure with ScanFold2. The default target runs both tools and produces
per-sample rG4 detections, summaries, and visualizations. It does not produce
a combined rG4/ScanFold feature table.

## Inputs and requirements

Run commands from this directory. Supply one FASTA file per sample at
`resources/<sample>.fa`; the basename must match a name in `samples`.
For the default `gencode.v47` sample, obtain a GENCODE v47 transcript FASTA and
place or link it at `resources/gencode.v47.fa`. FASTA IDs must
be unique and safe as filenames; see the [ScanFold input rules](vendor/scanfold2/README.md#input-and-filenames).

Install Snakemake with Conda support and have `conda` and `git` available.
First use may fetch prediction tools and Conda packages, which requires network
access. Cluster runs also need the Snakemake SLURM executor
plugin and a working SLURM account/partition.

## Run

For a local run using the default sample:

```bash
snakemake --executor local --cores 4 --use-conda
```

For a small DAG check using the included `resources/input_test.fa`:

```bash
snakemake --executor local --cores 4 --use-conda --dry-run --config samples=input_test
```

On SLURM, use [`profiles/default/config.yaml`](profiles/default/config.yaml)
after setting its site-specific partition and, if required, account:

```bash
snakemake --profile profiles/default
```

## Configuration

The Snakefile loads `config/config.default.yaml`. To change its defaults,
copy it to a local file and pass that file explicitly:

```bash
cp config/config.default.yaml config/local.yaml
snakemake --executor local --cores 4 --use-conda --configfile config/local.yaml
```

### Functioning parameters

| Parameter | Default | Meaning |
|---|---|---|
| `samples` | `[gencode.v47]` | FASTA basenames under `resources/`. |
| `batch_size` | `1000` | Sequences per rG4detector batch. |
| `scanfold_batch_size` | `1000` | Sequences per ScanFold batch, independent of `batch_size`. |
| `scanfold_max_len` | `20000` | Maximum sequence length sent to ScanFold, in nucleotides; longer records are reported in `scanfold_excluded.tsv`. |
| `scanfold_retry_subbatch_size` | `100` | Reserved retry batch size; unused by the included Snakefile. |
| `peak_threshold` | `1.56` | Minimum rG4 score used for `peak_counts.csv` and `peak_stats.tsv`. |
| `threshold_low`, `threshold_high` | `1.0`, `2.0` | Score cutoffs for candidate and high-confidence rG4 peaks in `rg4_summary.csv`. |
| `min_peak_prominence`, `min_peak_width` | `0.3`, `20` | Minimum prominence and width (nt) of summarized peaks. |
| `annotation_file` | `null` | Optional BED annotation file for the rG4 summary. |
| `group_annotation_file` | `null` | Optional group-label table for the annotated summary and peak statistics. |
| `final_partners_dirs` | `[]` | Reserved ScanFold FinalPartners input directories; unused by the included Snakefile. |
| `final_partners_output` | `results/gencode.v47/scanfold_final_partners_stats.tsv` | Reserved FinalPartners summary path; unused by the included Snakefile. |
| `final_partners_n_bins` | `10` | Reserved FinalPartners histogram bin count; unused by the included Snakefile. |

For a one-off value, use `--config samples=input_test` or
`--config scanfold_batch_size=1`.

### Visualization customization

| Parameter | Default | Meaning |
|---|---|---|
| `viz_output_dir` | `results/{sample}/viz` | Directory for rG4 visualizations. |
| `viz_plots` | `distributions ranking details comparison` | Plot groups to generate. |
| `viz_format` | `both` | `static`, `interactive`, or both. |
| `viz_top_n`, `viz_metric` | `20`, `top_peak_score` | Number of ranked transcripts and ranking metric. |
| `viz_dpi`, `viz_figsize` | `300`, `10x6` | Resolution and size of static plots. |
| `viz_class_col` | `null` | Optional class column for comparison plots. |
| `viz_transcript_detail` | `null` | Optional transcript ID for a detail plot. |
| `viz_filter` | `null` | Optional expression for filtering visualization rows. |

## Main outputs

| Path under `results/<sample>/` | Meaning |
|---|---|
| `rg4detector/detection.csv` | rG4detector nucleotide-level detection scores. |
| `rg4detector/rg4_summary.csv` | Per-transcript rG4 summary. |
| `rg4detector/rg4_summary_annotated.csv` | Summary with optional group annotations; a copy of the summary when none are configured. |
| `rg4detector/peak_counts.csv`, `rg4detector/peak_stats.tsv` | Per-transcript peak counts and aggregate statistics. |
| `scanfold2/batch_*/` | ScanFold2 results for each FASTA batch. |
| `scanfold_excluded.tsv` | Records omitted from ScanFold, with source FASTA checksum and length threshold. |

The default target also creates `detection.csv` at the sample root and rG4
visualizations under `viz/`. `scanfold2/.complete` is the workflow completion
marker for ScanFold batches; it is not a structure result. For ScanFold alone,
run `snakemake --executor local --cores 4 --use-conda all_scanfold2`.

## Reproducibility and limits

ScanFold2 is pinned to upstream commit
`c5cc73291fa5a1c06bf5b3c8367a396a3238141c` plus the tracked
[`multifasta.patch`](vendor/scanfold2/multifasta.patch). Installation and
execution check source and model hashes. See its [installation and output guide](vendor/scanfold2/README.md)
and [historical provenance](vendor/scanfold2/PROVENANCE.md).
