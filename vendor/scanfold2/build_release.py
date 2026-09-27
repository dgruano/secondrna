"""Build a deterministic standalone archive from verified ScanFold2 sources."""

import argparse
import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
subprocess.run([sys.executable, str(HERE / "install.py"), "--check"], check=True)
lock = json.loads((HERE / "source-lock.json").read_text())
files = {
    f"software/ScanFold2.0/{name}": ROOT / "software/ScanFold2.0" / name
    for name in lock["patched_files"]
}
for path in HERE.iterdir():
    if path.is_file():
        files[str(path.relative_to(ROOT))] = path
for name in [
    "README.md",
    "workflow/rules/scanfold_core.smk",
    "workflow/rules/scanfold_gpu.smk",
    "workflow/envs/scanfold2.yaml",
    "workflow/envs/scanfold2_gpu.yaml",
    "workflow/scripts/summarize_scanfold_final_partners.py",
    "workflow/utils/scanfold_lib.py",
]:
    files[name] = ROOT / name
manifest = {}
args.output.parent.mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(args.output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for name, path in sorted(files.items()):
        data = path.read_bytes()
        manifest[name] = hashlib.sha256(data).hexdigest()
        info = zipfile.ZipInfo(name)
        info.compress_type = zipfile.ZIP_DEFLATED
        archive.writestr(info, data)
    archive.writestr(
        zipfile.ZipInfo("SHA256SUMS.json"), json.dumps(manifest, indent=2) + "\n"
    )
print(
    f"{args.output}: {len(files)} files; SHA256 {hashlib.sha256(args.output.read_bytes()).hexdigest()}"
)
