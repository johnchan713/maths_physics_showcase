#!/usr/bin/env python3
"""Actual CLI checks: frozen replay, historical endpoint, cadence and failed gradient."""
import argparse
import csv
from pathlib import Path
import subprocess
import tempfile
import unittest

parser = argparse.ArgumentParser()
parser.add_argument("--optimizer", required=True)
parser.add_argument("--legacy-optimizer")
options, unittest_args = parser.parse_known_args()
binary = str(Path(options.optimizer).resolve(strict=True))
source = Path(__file__).resolve().parents[1] / "candidates/wave_k3_amplification_t004_screening.csv"


class LateGrowthCli(unittest.TestCase):
    def test_replay_clock_cadence_and_gradient_gate(self):
        with tempfile.TemporaryDirectory(prefix="late-growth-test-", dir=Path.cwd()) as temporary:
            folder = Path(temporary)
            original = source.read_bytes()

            def run(name, extra, executable=binary):
                prefix = folder / name
                command = [executable, "--state-input", str(source), "--grid", "16", "--fine-grid", "32",
                           "--seed-bandwidth", "3", "--energy", "10", "--final-time", "0.004",
                           "--dt", "0.0005", "--fine-max-dt", "0.0005", "--iterations", "0",
                           "--search-track", "amplification", "--output", str(prefix) + "-trace.csv",
                           "--state-output", str(prefix) + "-state.csv", "--evidence-output", str(prefix) + "-evidence.csv"]
                completed = subprocess.run(command + extra, capture_output=True, text=True, timeout=180)
                self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
                with Path(str(prefix) + "-trace.csv").open(newline="") as stream:
                    rows = list(csv.DictReader(stream))
                self.assertTrue(rows and all(None not in r and None not in r.values() for r in rows))
                evidence_path = Path(str(prefix) + "-evidence.csv")
                with evidence_path.open(newline="") as stream:
                    evidence_rows = list(csv.DictReader(stream))
                self.assertEqual(len(evidence_rows), 18,
                                 name + " reported success without complete evidence: " + completed.stdout)
                return rows

            run("endpoint", [])
            for cadence in (1, 97):
                name = "late-" + str(cadence)
                run(name, ["--growth-objective", "late-rate", "--diagnostic-every", str(cadence),
                           "--late-rate-output", str(folder / (name + "-rates.csv"))])
                self.assertEqual((folder / (name + "-state.csv")).read_bytes(), original)
                self.assertEqual((folder / (name + "-evidence.csv")).read_bytes(),
                                 (folder / "endpoint-evidence.csv").read_bytes())
            self.assertEqual((folder / "late-1-rates.csv").read_bytes(),
                             (folder / "late-97-rates.csv").read_bytes())
            with (folder / "late-1-rates.csv").open(newline="") as stream:
                rates = list(csv.DictReader(stream))
            self.assertEqual(len(rates), 10)
            for grid in ("coarse", "fine"):
                times = [float(row["time"]) for row in rates if row["resolution"] == grid]
                self.assertEqual(times, [0.002, 0.0025, 0.003, 0.0035, 0.004])
            failed = run("failed-gradient", ["--growth-objective", "late-rate", "--iterations", "1",
                                              "--gradient-tolerance", "1e-15"])
            self.assertFalse(any(row["stage"] == "accepted" for row in failed))
            self.assertTrue(any(row["stage"] == "gradient-check" for row in failed))
            self.assertEqual((folder / "failed-gradient-state.csv").read_bytes(), original)
            if options.legacy_optimizer:
                run("legacy", [], str(Path(options.legacy_optimizer).resolve(strict=True)))
                for suffix in ("state", "evidence"):
                    self.assertEqual((folder / ("legacy-" + suffix + ".csv")).read_bytes(),
                                     (folder / ("endpoint-" + suffix + ".csv")).read_bytes())
            self.assertEqual(source.read_bytes(), original)


if __name__ == "__main__":
    unittest.main(argv=[__file__, *unittest_args])
