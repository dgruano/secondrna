# ScanFold Oversized Fallback — Design Spec

**Date:** 2026-06-16
**Branch:** feat/scanfold
**Status:** Approved

---

## Problem

The retry workflow splits failed batches into sub-batches of ≤100 sequences and re-runs them on
CPU. Some sub-batches still fail because they contain at least one pathologically long
transcript that exhausts memory or exceeds the job time limit, blocking all other sequences in
the batch. The two observed failure modes are OOM kill (`batch_00016_retry_00000`, 347k nt) and
O(n²) timeout (`batch_00002_retry_00000`, three ~100k nt sequences).

---

## Solution Overview

Extend the **existing** `scanfold_split_retry_batches` checkpoint to be length-aware. Long
sequences are routed out of the retry sub-batches into a separate `scanfold_oversized_batches/`
directory before the retry jobs run. Each oversized sequence gets its own single-sequence
FASTA and is processed by a dedicated highmem rule rather than a standard retry batch rule.

```
scanfold_batches/           (original, ≤N seqs each)
        │
        ▼
scanfold_split_retry_batches  (checkpoint — modified)
        │
        ├── seq_len ≤ MAX_LEN ──► scanfold_retry_batches/   → scanfold2_run_retry_batch (CPU)
        │
        └── seq_len > MAX_LEN ──► scanfold_oversized_batches/ → scanfold2_run_oversized (CPU)
                                   scanfold_oversized_skip.tsv
```

`MAX_LEN` defaults to 20,000 nt (config key `scanfold_oversized_max_len`).

---

## Files Changed

| File | Change |
|------|--------|
| `workflow/scripts/check_scanfold_partial.py` | Add `--max-seq-len` and `--oversized-dir` flags; length-split pending sequences; write oversized FASTAs, `oversized_manifest.txt`, and `scanfold_oversized_skip.tsv` |
| `workflow/rules/scanfold.smk` | Add `SCANFOLD_OVERSIZED_MAX_LEN` config read; extend `scanfold_split_retry_batches` with oversized outputs; add `rule scanfold2_run_oversized`; add `get_scanfold_oversized_batches` + `rule scanfold2_oversized_all` |

`workflow/scripts/scanfold_summarize.py` is **not changed**.

---

## Detailed Spec

### `check_scanfold_partial.py` additions

Two new CLI flags:

```
--max-seq-len N      Sequences longer than N nt are written to --oversized-dir
                     instead of the normal sub-batch output.
                     Default: no length filtering (all sequences go to normal output).

--oversized-dir DIR  Directory for single-sequence FASTAs of oversized sequences.
                     Also writes oversized_manifest.txt and scanfold_oversized_skip.tsv
                     inside that directory.
                     Required when --max-seq-len is set.
```

Behaviour when `--max-seq-len` is set:

1. After resolving pending sequence IDs for each batch, measure the length of each pending
   sequence via a single-pass scan of the source FASTA (no full load into memory).
2. Sequences ≤ `--max-seq-len` → existing sub-batch logic → `--out-dir`.
3. Sequences > `--max-seq-len`:
   - Write one FASTA per sequence to `--oversized-dir/batch_{seq_id}.fa`
     (or a sanitised batch ID derived from the seq_id and source batch).
   - Append to `--oversized-dir/oversized_manifest.txt` (same TSV format as
     `batch_manifest.txt`: `batch_id`, `seq_count`, `seq_ids`).
   - Append a row to `--oversized-dir/scanfold_oversized_skip.tsv`
     (columns: `seq_id`, `seq_len`, `source_batch`).

The skip TSV is written unconditionally for all oversized sequences, regardless of whether the
oversized rule will actually be run.

### `scanfold.smk` additions

#### Config (top of file, alongside `SCANFOLD_RETRY_SUBBATCH_SIZE`)

```python
SCANFOLD_OVERSIZED_MAX_LEN = config.get("scanfold_oversized_max_len", 20_000)
```

#### Extended `checkpoint scanfold_split_retry_batches`

Add three new outputs and one new param:

```python
output:
    manifest           = "results/{sample}/scanfold_retry_batches/batch_manifest.txt",
    batches            = directory("results/{sample}/scanfold_retry_batches/"),
    oversized_manifest = "results/{sample}/scanfold_oversized_batches/batch_manifest.txt",  # NEW
    oversized_batches  = directory("results/{sample}/scanfold_oversized_batches/"),          # NEW
    skip_tsv           = "results/{sample}/scanfold_oversized_skip.tsv",                     # NEW
params:
    subbatch_size = SCANFOLD_RETRY_SUBBATCH_SIZE,
    max_seq_len   = SCANFOLD_OVERSIZED_MAX_LEN,    # NEW
resources:
    cpus_per_task = 1,                             # NEW
```

Shell call gains `--max-seq-len {params.max_seq_len}` and
`--oversized-dir results/{wildcards.sample}/scanfold_oversized_batches`.

#### New `rule scanfold2_run_oversized`

```python
rule scanfold2_run_oversized:
    """Process a single ultra-long sequence on a high-memory CPU node."""
    input:
        flag  = "software/ScanFold2.0/installed.txt",
        fasta = "results/{sample}/scanfold_oversized_batches/batch_{batch_id}.fa",
    output:
        out = "results/{sample}/scanfold2_oversized/batch_{batch_id}/ScanFold_run.log",
    conda:
        "../envs/scanfold2.yaml"
    resources:
        runtime          = "5d",
        mem_mb           = 262_144,   # 256 GB
        slurm_partition  = "highmem",
        cpus_per_task    = 1,
    params:
        window      = 120,
        step        = 1,
        temperature = 37,
        shuffle     = "mono",
        folder      = lambda wc, output: os.path.dirname(output.out),
    log:
        "logs/{sample}/scanfold2_oversized/{batch_id}.log",
    benchmark:
        "benchmarks/{sample}/scanfold2_oversized/{batch_id}.tsv"
    shell: """
        {{
        input_fasta="$(realpath {input.fasta})"
        output_folder=$(realpath {params.folder})
        cd software/ScanFold2.0
        python ScanFold2.0.py $input_fasta \
            -w {params.window} -s {params.step} \
            -t {params.temperature} --shuffle {params.shuffle} \
            --folder $output_folder
        }} >{log} 2>&1
    """
```

#### New input function and aggregate rule

```python
def get_scanfold_oversized_batches(wc):
    checkpoints.scanfold_split_retry_batches.get(sample=wc.sample)
    batch_ids = get_batch_ids(wc.sample, subdir="scanfold_oversized_batches")
    return expand(
        "results/{sample}/scanfold2_oversized/batch_{batch_id}/ScanFold_run.log",
        sample=wc.sample,
        batch_id=batch_ids,
    )


rule scanfold2_oversized_all:
    input:
        lambda wc: get_scanfold_oversized_batches(
            type("WC", (), {"sample": "gencode.v47.repeat.simple"})()
        ),
```

---

## Architecture Decisions

### Why modify the retry checkpoint rather than add a new one

The retry checkpoint already does the right thing: it scans for pending sequences across all
original batches and writes sub-batch FASTAs. Adding length routing here means oversized
sequences are isolated before any retry job is submitted, so they never block a batch.
A separate tertiary checkpoint would require another round of job submission and waiting.

### Why highmem for oversized sequences

The bottleneck for ultra-long sequences is the ScanFold-Fold partner detection loop (O(n²)),
not TensorFlow inference. A high-memory CPU node covers the RAM requirement (≥256 GB for
100k–200k nt) and avoids tying up a long-running job on a standard partition with its shorter
time limits.

### Why 20,000 nt threshold

The fold loop is empirically O(n²) (R² = 0.9968). At 20k nt the predicted fold time is ~33 min;
a single sequence at 100k nt takes ~14 h. The threshold is exposed as
`config["scanfold_oversized_max_len"]` so it can be tuned without touching rule code.

### Batch ID generation for oversized sequences

Each oversized sequence gets a batch ID derived from its source batch ID and a zero-padded
index (e.g., `00002_retry_00000_oversized_00000`). This keeps the naming consistent with
existing patterns and avoids collisions when multiple batches contribute oversized sequences.
Sequence IDs that contain characters invalid in filenames (e.g., `/`, `|`) must be sanitised
(replace with `_`) before use as a batch ID component.

### Known limitation: summarize does not glob `scanfold2_oversized/`

`scanfold_summarize.py` is not changed in this implementation. Oversized sequences will not
appear in the final `scanfold_stats.tsv` unless the summarize step is updated in a follow-up
to also read CT files from `scanfold2_oversized/`. The `scanfold_oversized_skip.tsv` file
provides a record of which sequences were skipped.

---

## Usage

```bash
# Run retry + oversized in one go (oversized jobs submitted automatically):
snakemake scanfold2_retry_all scanfold2_oversized_all --profile profiles/default

# Override threshold:
snakemake scanfold2_oversized_all --profile profiles/default \
    --config scanfold_oversized_max_len=10000

# Skip oversized jobs (they remain in scanfold_oversized_skip.tsv but are not run):
snakemake scanfold2_retry_all --profile profiles/default
```
