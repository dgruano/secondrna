import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_rg4detector_summary_handles_block_detection_csv(tmp_path):
    detection = tmp_path / "detection.csv"
    output = tmp_path / "summary.csv"

    detection.write_text(
        "ENST00000001|gene1,0,1,2,3,4\n"
        ",0.8,0.1,1.7,0.3,0.9\n"
        "ENST00000002|gene2,0,1,2,3\n"
        ",0.2,0.1,0.5,0.6\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow" / "scripts" / "rg4detector_summarize.py"),
            "--input",
            str(detection),
            "--output",
            str(output),
            "--threshold_low",
            "1.0",
            "--threshold_high",
            "2.0",
            "--min_peak_prominence",
            "0.3",
            "--min_peak_width",
            "1",
        ],
        cwd=str(ROOT),
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert output.exists()

    rows = output.read_text(encoding="utf-8").strip().splitlines()
    assert len(rows) >= 2
    assert rows[0].startswith("transcript_id")
