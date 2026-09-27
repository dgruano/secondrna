"""Real folding regression: run with the ScanFold2 conda environment's Python."""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2] / "software/ScanFold2.0"


class MultiFastaTest(unittest.TestCase):
    def test_outputs_and_validation(self):
        with tempfile.TemporaryDirectory(prefix="scanfold-regression-") as tmp:
            root = Path(tmp)
            sequences = {
                "alpha": "GCGCAAAATTTTGCGC" * 9,
                "beta": "GGGAAACCCUUUGGGCCC" * 8,
            }

            def run(stem, records, success=True):
                fasta = root / (stem + ".fa")
                fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
                out = root / stem
                result = subprocess.run(
                    [
                        sys.executable,
                        str(REPO / "ScanFold2.0.py"),
                        str(fasta),
                        "-w",
                        "120",
                        "-s",
                        "1",
                        "--folder",
                        str(out),
                    ],
                    cwd=REPO,
                    env={
                        **os.environ,
                        "CUDA_VISIBLE_DEVICES": "-1",
                        "TF_NUM_INTRAOP_THREADS": "1",
                        "TF_NUM_INTEROP_THREADS": "1",
                        "OMP_NUM_THREADS": "1",
                    },
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                self.assertEqual(result.returncode == 0, success, result.stdout[-6000:])
                return out, result.stdout

            multi, _ = run("multi", list(sequences.items()))
            for name, seq in sequences.items():
                single, _ = run("single_" + name, [(name, seq)])
                prefix = name + ".win_120.stp_1.csv"
                for suffix in [
                    "",
                    ".ScanFold.FinalPartners.txt",
                    ".bp",
                    ".Zavg.wig",
                    ".scan-MFE.wig",
                    ".scan-ED.wig",
                    ".scan-zscores.wig",
                ]:
                    a, b = multi / (prefix + suffix), single / (prefix + suffix)
                    self.assertTrue(a.is_file(), str(a))
                    self.assertEqual(a.read_bytes(), b.read_bytes(), suffix)
                for suffix in [
                    "no_filter.ct",
                    "minus_1.ct",
                    "minus_2.ct",
                    "no_filter.dbn",
                    "minus_1.dbn",
                    "minus_2.dbn",
                    "ALL.bp",
                    "README.md",
                ]:
                    self.assertGreater(
                        (multi / f"multi.{name}.{suffix}").stat().st_size, 0
                    )
                    self.assertGreater(
                        (single / f"single_{name}.{suffix}").stat().st_size, 0
                    )
                session = json.loads(
                    (multi / f"multi.{name}.igv_session.json").read_text()
                )
                self.assertTrue(session["locus"].startswith(name + ":"))
                self.assertTrue(
                    all(Path(track["path"]).exists() for track in session["tracks"])
                )
                self.assertTrue((single / "igv_session.json").is_file())
            self.assertFalse((multi / "igv_session.json").exists())
            self.assertIn(
                "SHA256 ScanFold2.0.py:", (multi / "ScanFold_run.log").read_text()
            )
            for stem, records, error in [
                (
                    "duplicate",
                    [("same", sequences["alpha"])] * 2,
                    "Duplicate FASTA record ID",
                ),
                (
                    "unsafe",
                    [("../escape", sequences["alpha"])],
                    "Unsafe FASTA record ID",
                ),
                ("empty", [], "contains no records"),
            ]:
                out, log = run(stem, records, success=False)
                self.assertIn(error, log)
                self.assertFalse(out.exists())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=REPO)
    args, rest = parser.parse_known_args()
    REPO = args.repo.resolve()
    unittest.main(argv=[sys.argv[0], *rest])
