#!/usr/bin/env python3
"""Install or verify the pinned upstream source plus the secondrna patch.

Never resets an existing checkout. --check is offline and checks source/model
content even when an old installed.txt exists. Runtime caches are excluded.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mismatches(repo, expected):
    return [
        name
        for name, checksum in expected.items()
        if not (repo / name).is_file() or digest(repo / name) != checksum
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--destination", type=Path, default=HERE.parents[1] / "software/ScanFold2.0"
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    repo = args.destination.resolve()
    lock = json.loads((HERE / "source-lock.json").read_text())
    patch = HERE / "multifasta.patch"
    if digest(patch) != lock["patch_sha256"]:
        parser.error("Patch checksum does not match source-lock.json")
    if not repo.exists():
        if args.check:
            parser.error(
                f"Missing ScanFold2 source: {repo}; run this script without --check"
            )
        repo.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", "--no-checkout", lock["upstream_url"], str(repo)],
            check=True,
        )
        subprocess.run(
            ["git", "-C", str(repo), "checkout", "--detach", lock["upstream_commit"]],
            check=True,
        )
    changed = mismatches(repo, lock["patched_files"])
    if changed and not args.check:
        baseline_changes = mismatches(repo, lock["upstream_files"])
        if baseline_changes:
            parser.error(
                "Existing source differs from both supported versions; preserve your edits and use a new destination. "
                "Mismatches: " + ", ".join(baseline_changes[:12])
            )
        subprocess.run(
            ["git", "-C", str(repo), "apply", "--check", str(patch)], check=True
        )
        subprocess.run(["git", "-C", str(repo), "apply", str(patch)], check=True)
        changed = mismatches(repo, lock["patched_files"])
    if changed:
        parser.error(
            "ScanFold2 source/model verification failed: " + ", ".join(changed[:12])
        )
    print(
        "Verified ScanFold2 upstream "
        + lock["upstream_commit"]
        + " + patch "
        + lock["patch_sha256"]
    )


if __name__ == "__main__":
    main()
