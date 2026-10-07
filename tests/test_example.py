import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "example"


def test_example_matches_expected_output(tmp_path):
    outdir = tmp_path / "arcsv_out"
    cmd = [
        sys.executable,
        "-m",
        "arcsv",
        "call",
        "-i",
        str(EXAMPLE / "input.bam"),
        "-r",
        "20:0-250000",
        "-R",
        str(EXAMPLE / "reference.fa"),
        "-G",
        str(EXAMPLE / "gaps.bed"),
        "-o",
        str(outdir),
    ]
    env = dict(os.environ, PYTHONPATH=str(REPO))
    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]

    expected = (EXAMPLE / "expected_output.tab").read_text()
    assert (outdir / "arcsv_out.tab").read_text() == expected

    # the VCF header records the run date and the arcsv version, so skip those
    def vcf_without_date_and_version(text):
        skip = ("##fileDate=", "##source=")
        return [line for line in text.splitlines() if not line.startswith(skip)]

    expected_vcf = (EXAMPLE / "expected_output.vcf").read_text()
    assert vcf_without_date_and_version(
        (outdir / "arcsv_out.vcf").read_text()
    ) == vcf_without_date_and_version(expected_vcf)
