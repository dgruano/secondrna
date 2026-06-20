"""
Evaluate partial ScanFold2 results for interrupted batches.

Compares sequences in each input batch FASTA against completed outputs
(marked by .no_filter.ct) and reports which sequences still need processing.

When --batches is omitted, all batches in the original batch manifest are
scanned automatically and only incomplete ones are processed.

Usage (interactive):
    python check_scanfold_partial.py
    python check_scanfold_partial.py --batches 00002 00016
    python check_scanfold_partial.py --write-fasta --subbatch-size 100

Usage (from Snakemake checkpoint):
    python check_scanfold_partial.py \\
        --sample gencode.v47.repeat.simple \\
        --subbatch-size 100 \\
        --write-fasta \\
        --out-dir results/gencode.v47.repeat.simple/scanfold_retry_batches \\
        --max-seq-len 20000 \\
        --oversized-dir results/gencode.v47.repeat.simple/scanfold_oversized_batches
"""

import argparse
import hashlib
import math
import re
from pathlib import Path

# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------


def get_fasta_ids(fasta_path: Path) -> list[str]:
    ids = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip())
    return ids


def extract_sequences(source_fasta: Path, seq_ids: set[str], out_path: Path) -> None:
    """Write a subset of sequences from source_fasta to out_path."""
    capture = False
    with open(source_fasta) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                capture = line[1:].strip() in seq_ids
            if capture:
                fout.write(line)


def measure_fasta_lengths(source_fasta: Path) -> dict[str, int]:
    """Return {seq_id: length} for all sequences in source_fasta (single pass)."""
    lengths: dict[str, int] = {}
    current_id: str | None = None
    current_len = 0
    with open(source_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].strip()
                current_len = 0
            else:
                current_len += len(line.strip())
    if current_id is not None:
        lengths[current_id] = current_len
    return lengths


def _sanitise_batch_id(s: str) -> str:
    """Replace filename-unsafe characters with underscores."""
    return re.sub(r"[/|\\:*?\"<>]", "_", s)


def write_oversized_batches(
    source_fasta: Path,
    oversized_ids: list[str],
    oversized_dir: Path,
    batch_id: str,
    seq_lengths: dict[str, int],
    skip_tsv_path: Path,
) -> list[tuple[str, list[str], Path]]:
    """Write one single-sequence FASTA per oversized sequence.

    Files are named: batch_{batch_id}_oversized_{i:05d}.fa
    Appends rows to skip_tsv_path and writes batch_manifest.txt.
    Returns list of (oversized_batch_id, [seq_id], path).
    """
    oversized_dir.mkdir(parents=True, exist_ok=True)
    results = []
    safe_batch_id = _sanitise_batch_id(batch_id)

    skip_tsv_exists = skip_tsv_path.exists()
    with open(skip_tsv_path, "a") as skip_f:
        if not skip_tsv_exists:
            skip_f.write("seq_id\tseq_len\tsource_batch\n")
        for i, seq_id in enumerate(oversized_ids):
            oversized_batch_id = f"{safe_batch_id}_oversized_{i:05d}"
            out_path = oversized_dir / f"batch_{oversized_batch_id}.fa"
            extract_sequences(source_fasta, {seq_id}, out_path)
            seq_len = seq_lengths.get(seq_id, 0)
            skip_f.write(f"{seq_id}\t{seq_len}\t{batch_id}\n")
            results.append((oversized_batch_id, [seq_id], out_path))

    return results


# ---------------------------------------------------------------------------
# Batch helpers
# ---------------------------------------------------------------------------


def read_batch_manifest(manifest_path: Path) -> list[str]:
    """Return all batch IDs from a batch_manifest.txt."""
    batch_ids = []
    with open(manifest_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if parts:
                batch_ids.append(parts[0])
    return batch_ids


def write_batch_manifest(out_dir: Path, entries: list[tuple[str, list[str]]]) -> None:
    """Write batch_manifest.txt compatible with get_batch_ids().

    entries: list of (batch_id, seq_ids)
    """
    manifest_path = out_dir / "batch_manifest.txt"
    with open(manifest_path, "w") as f:
        f.write("# batch_id\tseq_count\tseq_ids\n")
        for batch_id, seq_ids in entries:
            f.write(f"{batch_id}\t{len(seq_ids)}\t{','.join(seq_ids)}\n")


def _ct_seq_ids(ct_path: Path) -> set[str]:
    """Parse CT file header lines (length<TAB>seq_id) to extract sequence IDs."""
    ids = set()
    with open(ct_path) as f:
        for line in f:
            parts = line.split("\t", 1)
            if len(parts) == 2 and parts[0].strip().isdigit():
                ids.add(parts[1].strip())
    return ids


def get_completed_ids(result_dir: Path, batch_id: str) -> set[str]:
    prefix = f"batch_{batch_id}."
    suffix = ".no_filter.ct"
    ids = {
        p.name[len(prefix) : -len(suffix)]
        for p in result_dir.glob(f"{prefix}*.no_filter.ct")
    }
    if not ids:
        # ponytail: single-seq batches write only batch_{id}.no_filter.ct, no per-seq files
        batch_ct = result_dir / f"batch_{batch_id}.no_filter.ct"
        if batch_ct.exists():
            ids = _ct_seq_ids(batch_ct)
    return ids


def check_batch(batch_id: str, batches_dir: Path, gpu_dir: Path) -> dict:
    fasta = batches_dir / f"batch_{batch_id}.fa"
    result_dir = gpu_dir / f"batch_{batch_id}"

    if not fasta.exists():
        return {"error": f"FASTA not found: {fasta}"}

    all_ids = get_fasta_ids(fasta)

    # Collect completions from the original dir and any retry/oversized sub-dirs.
    completed = (
        get_completed_ids(result_dir, batch_id) if result_dir.exists() else set()
    )
    for d in gpu_dir.glob(f"batch_{batch_id}_*"):
        if d.is_dir():
            sub_id = d.name[len("batch_") :]
            completed |= get_completed_ids(d, sub_id)

    pending = [seq_id for seq_id in all_ids if seq_id not in completed]

    return {
        "batch_id": batch_id,
        "fasta": fasta,
        "total": len(all_ids),
        "completed": len(completed),
        "pending": len(pending),
        "pending_ids": pending,
    }


# ---------------------------------------------------------------------------
# Sub-batch splitting
# ---------------------------------------------------------------------------


def write_subbatches(
    source_fasta: Path,
    pending_ids: list[str],
    out_dir: Path,
    batch_id: str,
    subbatch_size: int,
) -> list[tuple[str, list[str], Path]]:
    """Split pending sequences into sub-batch FASTA files.

    Files are named: batch_{batch_id}_retry_{i:05d}.fa
    Returns list of (subbatch_id, seq_ids, path).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    n_subbatches = math.ceil(len(pending_ids) / subbatch_size)
    results = []

    for i in range(n_subbatches):
        chunk = pending_ids[i * subbatch_size : (i + 1) * subbatch_size]
        # ponytail: hash of seq IDs so the ID is stable across checkpoint re-runs
        chunk_hash = hashlib.sha1("\n".join(chunk).encode()).hexdigest()[:8]
        subbatch_id = f"{batch_id}_retry_{chunk_hash}"
        out_path = out_dir / f"batch_{subbatch_id}.fa"
        extract_sequences(source_fasta, set(chunk), out_path)
        results.append((subbatch_id, chunk, out_path))

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--batches",
        nargs="+",
        default=None,
        metavar="BATCH_ID",
        help="Batch IDs to check. If omitted, all batches in the original "
        "manifest are scanned and incomplete ones are auto-selected.",
    )
    parser.add_argument(
        "--exclude-batches",
        nargs="+",
        default=None,
        metavar="BATCH_ID",
        help="Batch IDs to skip (e.g. ones still running on SLURM). "
        "Use get_running_scanfold_batches.sh to generate this list.",
    )
    parser.add_argument(
        "--sample",
        default="gencode.v47.repeat.simple",
        help="Sample name (default: gencode.v47.repeat.simple)",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Top-level results directory (default: results/)",
    )
    parser.add_argument(
        "--write-fasta",
        action="store_true",
        help="Write FASTA file(s) of pending sequences",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Directory for retry FASTA output (default: "
        "<results-dir>/scanfold_retry_batches/). "
        "Required when called from Snakemake so the path matches the rule output.",
    )
    parser.add_argument(
        "--subbatch-size",
        type=int,
        default=None,
        metavar="N",
        help="Split pending sequences into sub-batches of N each "
        "(requires --write-fasta; omit to write one file per batch)",
    )
    parser.add_argument(
        "--max-seq-len",
        type=int,
        default=None,
        metavar="N",
        help="Sequences longer than N nt are written to --oversized-dir "
        "instead of the normal sub-batch output. Default: no length filtering.",
    )
    parser.add_argument(
        "--oversized-dir",
        type=Path,
        default=None,
        metavar="DIR",
        help="Directory for single-sequence FASTAs of oversized sequences. "
        "Required when --max-seq-len is set.",
    )
    args = parser.parse_args()

    if args.subbatch_size is not None and not args.write_fasta:
        parser.error("--subbatch-size requires --write-fasta")
    if args.subbatch_size is not None and args.subbatch_size <= 0:
        parser.error("--subbatch-size must be a positive integer")
    if args.max_seq_len is not None and args.oversized_dir is None:
        parser.error("--max-seq-len requires --oversized-dir")
    if args.oversized_dir is not None and args.max_seq_len is None:
        parser.error("--oversized-dir requires --max-seq-len")

    sample_dir = args.results_dir / args.sample
    batches_dir = sample_dir / "scanfold_batches"
    gpu_dir = sample_dir / "scanfold2"
    out_dir = args.out_dir or args.results_dir / "scanfold_retry_batches"

    # --- Resolve which batch IDs to check ---
    if args.batches is not None:
        batch_ids = args.batches
        print(f"Sample : {args.sample}")
        print(f"Batches: {', '.join(batch_ids)} (explicit)")
    else:
        manifest = batches_dir / "batch_manifest.txt"
        if not manifest.exists():
            raise FileNotFoundError(
                f"Batch manifest not found: {manifest}\n"
                "Run the original split checkpoint first, or pass --batches explicitly."
            )
        batch_ids = read_batch_manifest(manifest)
        print(f"Sample : {args.sample}")
        print(f"Batches: auto-detected {len(batch_ids)} from manifest")

    excluded = set(args.exclude_batches) if args.exclude_batches else set()
    if excluded:
        before = len(batch_ids)
        batch_ids = [b for b in batch_ids if b not in excluded]
        print(f"Excluded (still running): {before - len(batch_ids)} batches")
    print()

    # --- Check each batch and collect incomplete ones ---
    manifest_entries: list[tuple[str, list[str]]] = []  # for retry manifest
    oversized_entries: list[tuple[str, list[str]]] = []  # for oversized manifest
    total_pending = 0

    # Path for the oversized skip TSV lives one level above --oversized-dir
    skip_tsv_path: Path | None = None
    if args.oversized_dir is not None:
        skip_tsv_path = args.oversized_dir.parent / "scanfold_oversized_skip.tsv"
        # Remove any stale skip TSV from a previous run so rows don't accumulate
        if skip_tsv_path.exists():
            skip_tsv_path.unlink()

    for batch_id in batch_ids:
        result = check_batch(batch_id, batches_dir, gpu_dir)

        if "error" in result:
            print(f"[batch_{batch_id}] ERROR: {result['error']}")
            continue

        total = result["total"]
        completed = result["completed"]
        pending_count = result["pending"]
        pending_ids = result["pending_ids"]

        # Skip fully-complete batches when auto-detecting (keep output terse)
        if not pending_ids:
            if args.batches is not None:
                print(f"[batch_{batch_id}] complete ({total}/{total})")
            continue

        # --- Length-based routing ---
        normal_ids = pending_ids
        oversized_ids: list[str] = []
        seq_lengths: dict[str, int] = {}

        if args.max_seq_len is not None:
            source_fasta = result["fasta"]
            seq_lengths = measure_fasta_lengths(source_fasta)
            normal_ids = [
                s for s in pending_ids if seq_lengths.get(s, 0) <= args.max_seq_len
            ]
            oversized_ids = [
                s for s in pending_ids if seq_lengths.get(s, 0) > args.max_seq_len
            ]

        total_pending += len(normal_ids) + len(oversized_ids)
        print(f"[batch_{batch_id}]")
        print(f"  Total sequences   : {total}")
        print(f"  Completed         : {completed}  ({completed / total * 100:.1f}%)")
        print(
            f"  Pending (need run): {pending_count}  ({pending_count / total * 100:.1f}%)"
        )
        if oversized_ids:
            print(f"  Oversized (>{args.max_seq_len} nt): {len(oversized_ids)}")
        print()
        print("  Pending sequence IDs:")
        for seq_id in normal_ids:
            print(f"    {seq_id}")
        for seq_id in oversized_ids:
            print(f"    {seq_id}  [OVERSIZED]")
        print()

        if args.write_fasta:
            source_fasta = result["fasta"]

            # Write oversized single-sequence FASTAs
            if oversized_ids and args.oversized_dir is not None:
                oversized = write_oversized_batches(
                    source_fasta,
                    oversized_ids,
                    args.oversized_dir,
                    batch_id,
                    seq_lengths,
                    skip_tsv_path,
                )
                print(f"  Oversized sequences → {args.oversized_dir}:")
                for ob_id, ob_seqs, ob_path in oversized:
                    print(f"    {ob_path}  ({len(ob_seqs)} seq)")
                    oversized_entries.append((ob_id, ob_seqs))
                print()

            # Write normal retry sub-batches (only non-oversized sequences)
            if normal_ids:
                if args.subbatch_size is None:
                    subbatch_id = f"{batch_id}_retry"
                    out_path = out_dir / f"batch_{subbatch_id}.fa"
                    out_dir.mkdir(parents=True, exist_ok=True)
                    extract_sequences(source_fasta, set(normal_ids), out_path)
                    manifest_entries.append((subbatch_id, normal_ids))
                    print(f"  Wrote {len(normal_ids)} sequences to: {out_path}")
                else:
                    subbatches = write_subbatches(
                        source_fasta, normal_ids, out_dir, batch_id, args.subbatch_size
                    )
                    n = len(subbatches)
                    print(
                        f"  Split {len(normal_ids)} pending sequences into {n} sub-batches "
                        f"of ≤{args.subbatch_size} each:"
                    )
                    for subbatch_id, seq_ids, path in subbatches:
                        print(f"    {path}  ({len(seq_ids)} seqs)")
                        manifest_entries.append((subbatch_id, seq_ids))
                print()

    # --- Write manifests for get_batch_ids() compatibility ---
    if args.write_fasta:
        out_dir.mkdir(parents=True, exist_ok=True)
        write_batch_manifest(out_dir, manifest_entries)
        print(f"Manifest written: {out_dir / 'batch_manifest.txt'}")
        print(f"Total retry sub-batches : {len(manifest_entries)}")
        print(f"Total pending sequences : {total_pending}")

        if args.oversized_dir is not None:
            args.oversized_dir.mkdir(parents=True, exist_ok=True)
            write_batch_manifest(args.oversized_dir, oversized_entries)
            print(
                f"Oversized manifest written: {args.oversized_dir / 'batch_manifest.txt'}"
            )
            print(f"Total oversized batches : {len(oversized_entries)}")
            if skip_tsv_path is not None:
                print(f"Skip TSV written: {skip_tsv_path}")

    if not args.write_fasta:
        print(f"Total pending sequences across all batches: {total_pending}")
        print(
            "(Re-run with --write-fasta [--subbatch-size N] to generate retry FASTAs)"
        )


if __name__ == "__main__":
    main()
