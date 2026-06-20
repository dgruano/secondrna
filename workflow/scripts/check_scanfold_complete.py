#!/usr/bin/env python3
"""
Validate that every sequence in the original ScanFold batch manifest
has a corresponding .no_filter.ct output file.

Scans results/{sample}/scanfold2/, results/{sample}/scanfold2_oversized/, and
results/{sample}/scanfold2_gpu/ (GPU flat output, if present) to find completed
sequences, then cross-references against the original batch manifest.

Writes a TSV report (seq_id, original_batch_id, status) and exits with
code 1 if any sequences are missing.

CLI usage:
    python check_scanfold_complete.py \\
        --sample gencode.v47.repeat.simple \\
        --results-dir results \\
        --output results/gencode.v47.repeat.simple/scanfold2/completeness_check.tsv

Snakemake usage:
    script: "../scripts/check_scanfold_complete.py"
    (expects smk.wildcards.sample, smk.params.results_dir,
             smk.output.report, smk.log[0])
"""

import argparse
import logging
import sys
from pathlib import Path


def setup_logging(log_file=None):
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file, mode="w"))
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
    )


def read_manifest(manifest_path: Path) -> dict[str, str]:
    """Return {seq_id: batch_id} for every sequence in a batch manifest."""
    seq_to_batch: dict[str, str] = {}
    with open(manifest_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            batch_id = parts[0]
            for seq_id in parts[2].split(","):
                seq_id = seq_id.strip()
                if seq_id:
                    seq_to_batch[seq_id] = batch_id
    return seq_to_batch


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
    """Return seq_ids with .no_filter.ct files in a batch result directory."""
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


def collect_all_completed(sample_dir: Path) -> set[str]:
    """Walk scanfold2/, scanfold2_oversized/, and scanfold2_gpu/ to collect all completed seq_ids."""
    completed: set[str] = set()
    for subdir_name in ("scanfold2", "scanfold2_oversized"):
        result_dir = sample_dir / subdir_name
        if not result_dir.exists():
            logging.info("Result directory not found (skipping): %s", result_dir)
            continue
        for batch_dir in sorted(result_dir.glob("batch_*")):
            if not batch_dir.is_dir():
                continue
            batch_id = batch_dir.name[len("batch_") :]
            ids = get_completed_ids(batch_dir, batch_id)
            completed.update(ids)
    # GPU output: flat dir (no batch_* subdirs), parse CT headers directly
    gpu_dir = sample_dir / "scanfold2_gpu"
    if gpu_dir.exists():
        for ct_file in sorted(gpu_dir.rglob("*.no_filter.ct")):
            ids = _ct_seq_ids(ct_file)
            logging.info("GPU CT file %s: %d sequences", ct_file.name, len(ids))
            completed.update(ids)
    return completed


def run(sample: str, results_dir: Path, report_path: Path) -> bool:
    """
    Validate completeness. Returns True if all sequences are present.
    Writes the TSV report to report_path.
    """
    sample_dir = results_dir / sample
    manifest_path = sample_dir / "scanfold_batches" / "batch_manifest.txt"

    logging.info("Sample       : %s", sample)
    logging.info("Manifest     : %s", manifest_path)

    expected = read_manifest(manifest_path)
    logging.info("Expected sequences : %d", len(expected))

    completed = collect_all_completed(sample_dir)
    logging.info("Completed sequences: %d", len(completed))

    missing = {sid: bid for sid, bid in expected.items() if sid not in completed}
    extra = completed - set(expected)

    logging.info("Missing sequences  : %d", len(missing))
    if extra:
        logging.info(
            "Extra sequences    : %d  (completed but not in original manifest)",
            len(extra),
        )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as f:
        f.write("seq_id\toriginal_batch_id\tstatus\n")
        for seq_id, batch_id in sorted(expected.items()):
            status = "COMPLETE" if seq_id in completed else "MISSING"
            f.write(f"{seq_id}\t{batch_id}\t{status}\n")
    logging.info("Report written: %s", report_path)

    if missing:
        logging.error("FAIL: %d sequences not processed:", len(missing))
        for seq_id, batch_id in sorted(missing.items()):
            logging.error("  %s  (original batch: %s)", seq_id, batch_id)
        return False

    logging.info("OK: all %d sequences have .no_filter.ct outputs.", len(expected))
    return True


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--sample", required=True, help="Sample name")
    p.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Top-level results dir",
    )
    p.add_argument("--output", type=Path, required=True, help="Output TSV report path")
    p.add_argument("--log", help="Log file path")
    return p.parse_args()


def main():
    if "snakemake" in globals():
        smk = globals()["snakemake"]
        setup_logging(smk.log[0] if smk.log else None)
        ok = run(
            sample=smk.wildcards.sample,
            results_dir=Path(smk.params.results_dir),
            report_path=Path(smk.output.report),
        )
    else:
        args = parse_args()
        setup_logging(args.log)
        ok = run(
            sample=args.sample,
            results_dir=args.results_dir,
            report_path=args.output,
        )

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
