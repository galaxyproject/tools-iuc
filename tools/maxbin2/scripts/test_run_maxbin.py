"""Regression checks for precise failure classification without expensive binning."""

import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("run_maxbin", Path(__file__).with_name("run_maxbin.py"))
WRAPPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WRAPPER)


class MaxBinFailureTests(unittest.TestCase):
    def test_unrelated_failures_remain_fatal(self):
        for code, message in ((1, WRAPPER.NO_BINS), (255, "hmmsearch failed"), (-9, "Killed")):
            with self.subTest(code=code), patch.object(
                WRAPPER.subprocess, "run", return_value=subprocess.CompletedProcess([], code, message.encode())
            ):
                self.assertEqual(WRAPPER.main([]), code)

    def test_message_followed_by_real_error_remains_fatal(self):
        message = WRAPPER.NO_BINS + "\nError: disk full"
        with patch.object(
            WRAPPER.subprocess, "run", return_value=subprocess.CompletedProcess([], 255, message.encode())
        ):
            self.assertEqual(WRAPPER.main([]), 255)

    def test_marker_outcome_returns_unclassified_contigs(self):
        with tempfile.TemporaryDirectory() as directory:
            contigs = Path(directory) / "input.fasta"
            contigs.write_text(">long\nACGTACGT\n>short\nAC\n")
            prefix = str(Path(directory) / "out")
            with patch.object(
                WRAPPER.subprocess, "run",
                return_value=subprocess.CompletedProcess([], 255, (WRAPPER.NO_BINS + "\n").encode()),
            ), patch.object(WRAPPER.Path, "glob", return_value=[]):
                self.assertEqual(WRAPPER.main([
                    "-contig", str(contigs), "-out", prefix, "-min_contig_length", "4"
                ]), 0)
            self.assertEqual(Path(prefix + ".noclass").read_text(), ">long\nACGTACGT\n")
            self.assertEqual(Path(prefix + ".tooshort").read_text(), ">short\nAC\n")
            self.assertEqual(Path(prefix + ".summary").stat().st_size, 0)

    def test_partial_bins_remain_fatal(self):
        with patch.object(
            WRAPPER.subprocess, "run",
            return_value=subprocess.CompletedProcess([], 255, (WRAPPER.NO_BINS + "\n").encode()),
        ), patch.object(WRAPPER.Path, "glob", return_value=[Path("out.001.fasta")]):
            self.assertEqual(WRAPPER.main(["-out", "out"]), 255)

    def test_stderr_error_remains_fatal(self):
        with patch.object(
            WRAPPER.subprocess, "run",
            return_value=subprocess.CompletedProcess([], 255, WRAPPER.NO_BINS.encode(), b"hmmsearch failed"),
        ):
            self.assertEqual(WRAPPER.main([]), 255)


if __name__ == "__main__":
    unittest.main()
