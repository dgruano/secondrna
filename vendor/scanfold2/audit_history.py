"""Read-only, bounded provenance audit of the local manuscript summary.

Checks all summary IDs/source counts and recomputes a small deterministic sample
from each available loose-file source. Does not claim archive-wide validation.
"""

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.utils.scanfold_lib import (
    compute_final_partners_stats,
    load_final_partners,
)

parser = argparse.ArgumentParser()
parser.add_argument("--summary", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
df = pd.read_csv(args.summary, sep="\t").set_index("versioned_id")
report = {
    "summary": str(args.summary),
    "summary_sha256": hashlib.sha256(args.summary.read_bytes()).hexdigest(),
    "rows": len(df),
    "duplicate_versioned_ids": int(df.index.duplicated().sum()),
    "source_counts": {k: int(v) for k, v in df.source_dir.value_counts().items()},
    "scope": "One available output directory per source; up to 3 summary records per directory. Not an exhaustive completeness audit.",
    "samples": [],
}
for source in sorted(df.source_dir.unique()):
    base = Path(source)
    if not base.is_dir():
        report["samples"].append({"source": source, "missing": True})
        continue
    candidates = [base] + sorted(p for p in base.iterdir() if p.is_dir())
    for directory in candidates:
        files = sorted(directory.glob("*.ScanFold.FinalPartners.txt"))
        if not files or not any(
            f.name.split(".win_")[0] in df.index
            and df.loc[f.name.split(".win_")[0], "source_dir"] == source
            for f in files
        ):
            continue
        log = directory / "ScanFold_run.log"
        text = log.read_text() if log.exists() else ""
        sample = {
            "source": source,
            "directory": str(directory),
            "final_partners_files": len(files),
            "no_filter_ct_files": len(list(directory.glob("*.no_filter.ct"))),
            "run_log_sha256": (
                hashlib.sha256(log.read_bytes()).hexdigest() if log.exists() else None
            ),
            "log_header": text.split("Starting ScanFold analysis")[0],
            "source_hash_in_log": "SHA256 ScanFold2.0.py:" in text,
            "recomputed_records": [],
        }
        for file in files:
            name = file.name.split(".win_")[0]
            if name not in df.index or df.loc[name, "source_dir"] != source:
                continue
            stats = compute_final_partners_stats(load_final_partners(file))
            columns = [c for c in stats.index if c in df.columns]
            expected = df.loc[name, columns].astype(float).to_numpy()
            observed = stats[columns].astype(float).to_numpy()
            sample["recomputed_records"].append(
                {
                    "id": name,
                    "file": str(file),
                    "sha256": hashlib.sha256(file.read_bytes()).hexdigest(),
                    "columns_compared": len(columns),
                    "matches": bool(
                        np.allclose(
                            expected, observed, rtol=1e-10, atol=1e-12, equal_nan=True
                        )
                    ),
                }
            )
            if len(sample["recomputed_records"]) == 3:
                break
        report["samples"].append(sample)
        break
args.output.write_text(json.dumps(report, indent=2) + "\n")
print(
    json.dumps(
        {
            "rows": report["rows"],
            "source_counts": report["source_counts"],
            "sample_matches": [
                [r["matches"] for r in s.get("recomputed_records", [])]
                for s in report["samples"]
            ],
        },
        indent=2,
    )
)
