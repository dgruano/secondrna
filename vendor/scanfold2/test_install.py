"""Offline installation test using a local Git remote at the upstream baseline."""

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent


class InstallTest(unittest.TestCase):
    def test_clean_install_and_modification_rejection(self):
        local_source = HERE.parents[1] / "software/ScanFold2.0"
        with tempfile.TemporaryDirectory(prefix="scanfold-install-") as tmp:
            root = Path(tmp)
            vendor = root / "vendor"
            vendor.mkdir()
            for name in ["install.py", "multifasta.patch", "source-lock.json"]:
                shutil.copy2(HERE / name, vendor / name)
            lock = json.loads((vendor / "source-lock.json").read_text())
            lock["upstream_url"] = str(local_source)
            (vendor / "source-lock.json").write_text(json.dumps(lock))
            dest = root / "installed"

            def run(*args):
                return subprocess.run(
                    [
                        sys.executable,
                        str(vendor / "install.py"),
                        "--destination",
                        str(dest),
                        *args,
                    ],
                    capture_output=True,
                    text=True,
                )

            self.assertNotEqual(run("--check").returncode, 0)
            result = run()
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(run("--check").returncode, 0)
            self.assertEqual(run().returncode, 0)  # idempotent
            source = dest / "ScanFold2.0.py"
            source.write_text(source.read_text() + "\n# local modification\n")
            self.assertNotEqual(run("--check").returncode, 0)
            self.assertIn("preserve your edits", run().stderr)
            self.assertTrue(source.read_text().endswith("# local modification\n"))
            with (vendor / "multifasta.patch").open("a") as f:
                f.write("\n")
            self.assertIn("Patch checksum", run().stderr)


if __name__ == "__main__":
    unittest.main()
