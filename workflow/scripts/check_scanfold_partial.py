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
        --out-dir results/gencode.v47.repeat.simple/scanfold_retry_batches
"""

import argparse
import math
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


def get_completed_ids(result_dir: Path, batch_id: str) -> set[str]:
    prefix = f"batch_{batch_id}."
    suffix = ".no_filter.ct"
    return {
        p.name[len(prefix) : -len(suffix)]
        for p in result_dir.glob(f"{prefix}*.no_filter.ct")
    }


def check_batch(batch_id: str, batches_dir: Path, gpu_dir: Path) -> dict:
    fasta = batches_dir / f"batch_{batch_id}.fa"
    result_dir = gpu_dir / f"batch_{batch_id}"

    if not fasta.exists():
        return {"error": f"FASTA not found: {fasta}"}

    all_ids = get_fasta_ids(fasta)

    if not result_dir.exists():
        # Batch never ran at all
        return {
            "batch_id": batch_id,
            "fasta": fasta,
            "total": len(all_ids),
            "completed": 0,
            "pending": len(all_ids),
            "pending_ids": all_ids,
        }

    completed = get_completed_ids(result_dir, batch_id)
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
        subbatch_id = f"{batch_id}_retry_{i:05d}"
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
    args = parser.parse_args()

    if args.subbatch_size is not None and not args.write_fasta:
        parser.error("--subbatch-size requires --write-fasta")
    if args.subbatch_size is not None and args.subbatch_size <= 0:
        parser.error("--subbatch-size must be a positive integer")

    sample_dir = args.results_dir / args.sample
    batches_dir = sample_dir / "scanfold_batches"
    gpu_dir = sample_dir / "scanfold2_gpu"
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
    print()

    # --- Check each batch and collect incomplete ones ---
    manifest_entries: list[tuple[str, list[str]]] = []  # for retry manifest
    total_pending = 0

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

        total_pending += pending_count
        print(f"[batch_{batch_id}]")
        print(f"  Total sequences   : {total}")
        print(f"  Completed         : {completed}  ({completed / total * 100:.1f}%)")
        print(
            f"  Pending (need run): {pending_count}  ({pending_count / total * 100:.1f}%)"
        )
        print()
        print("  Pending sequence IDs:")
        for seq_id in pending_ids:
            print(f"    {seq_id}")
        print()

        if args.write_fasta:
            source_fasta = result["fasta"]

            if args.subbatch_size is None:
                subbatch_id = f"{batch_id}_retry"
                out_path = out_dir / f"batch_{subbatch_id}.fa"
                out_dir.mkdir(parents=True, exist_ok=True)
                extract_sequences(source_fasta, set(pending_ids), out_path)
                manifest_entries.append((subbatch_id, pending_ids))
                print(f"  Wrote {pending_count} sequences to: {out_path}")
            else:
                subbatches = write_subbatches(
                    source_fasta, pending_ids, out_dir, batch_id, args.subbatch_size
                )
                n = len(subbatches)
                print(
                    f"  Split {pending_count} pending sequences into {n} sub-batches "
                    f"of ≤{args.subbatch_size} each:"
                )
                for subbatch_id, seq_ids, path in subbatches:
                    print(f"    {path}  ({len(seq_ids)} seqs)")
                    manifest_entries.append((subbatch_id, seq_ids))
            print()

    # --- Write manifest for get_batch_ids() compatibility ---
    if args.write_fasta:
        out_dir.mkdir(parents=True, exist_ok=True)
        write_batch_manifest(out_dir, manifest_entries)
        print(f"Manifest written: {out_dir / 'batch_manifest.txt'}")
        print(f"Total retry sub-batches : {len(manifest_entries)}")
        print(f"Total pending sequences : {total_pending}")

    if not args.write_fasta:
        print(f"Total pending sequences across all batches: {total_pending}")
        print(
            "(Re-run with --write-fasta [--subbatch-size N] to generate retry FASTAs)"
        )


if __name__ == "__main__":
    main()
