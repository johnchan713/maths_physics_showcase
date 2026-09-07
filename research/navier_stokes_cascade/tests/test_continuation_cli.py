#!/usr/bin/env python3
"""Exercise actual saved-field continuation, restart and evidence gates."""
import argparse
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

research = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("continuation_driver", research / "scripts/continue_candidates.py")
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
parser = argparse.ArgumentParser()
parser.add_argument("--solver", required=True)
options, unittest_args = parser.parse_known_args()
solver = str(Path(options.solver).resolve())


class ContinuationCliTests(unittest.TestCase):
    def test_saved_field_restart_and_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source = research / "candidates/wave_k3_amplification_t004_screening.csv"
            seed_sha = driver.sha(source)
            common = [solver, "--state-input", str(source), "--grid", "16", "--sampling-grid", "32",
                      "--dense-sampling-grid", "64", "--observation-interval", "0.002", "--dt", "0.0005"]

            def run(command, success=True):
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                self.assertEqual(result.returncode == 0, success, result.stdout)
                return result.stdout

            run(common + ["--final-time", "0.002", "--output", str(folder / "a.csv"),
                          "--checkpoint-output", str(folder / "split.chk")])
            checkpoint_before = driver.sha(folder / "split.chk")
            resumed = [solver, "--restart", str(folder / "split.chk"), "--final-time", "0.004",
                       "--output", str(folder / "b.csv"), "--checkpoint-output", str(folder / "split.chk")]
            run(resumed + ["--dt", "0.0001"], False)
            self.assertEqual(driver.sha(folder / "split.chk"), checkpoint_before)
            self.assertFalse((folder / "b.csv").exists())
            run(resumed)
            run(common + ["--final-time", "0.004", "--output", str(folder / "full.csv"),
                          "--checkpoint-output", str(folder / "full.chk")])
            self.assertEqual((folder / "full.chk").read_bytes(), (folder / "split.chk").read_bytes())
            a, b, full = (driver.read_evidence(folder / (name + ".csv")) for name in ("a", "b", "full"))
            self.assertEqual(a[-1], b[0])
            self.assertEqual(a + b[1:], full)
            self.assertEqual(seed_sha, driver.sha(source))
            run(common + ["--final-time", "0.003", "--output", str(folder / "invalid.csv"),
                          "--checkpoint-output", str(folder / "invalid.chk")], False)
            self.assertFalse((folder / "invalid.chk").exists())
            # Output/input aliases must fail before truncating either input.
            run(common + ["--final-time", "0.004", "--output", str(source),
                          "--checkpoint-output", str(folder / "invalid.chk")], False)
            self.assertEqual(seed_sha, driver.sha(source))
            with (folder / "full.chk").open("ab") as corrupt:
                corrupt.write(b"x")
            run([solver, "--restart", str(folder / "full.chk"), "--final-time", "0.006",
                 "--output", str(folder / "bad.csv"), "--checkpoint-output", str(folder / "full.chk")], False)
            # A real weakly growing run can still fail the late-growth gate.
            assessment = driver.assess_pair(full, full, 0.004, interval=0.002)
            self.assertEqual(assessment["numerical_failures"], [])
            bad = [dict(row) for row in full]
            bad[-1]["vorticity_dense_ratio"] *= 1.2
            failed = driver.assess_pair(full, bad, 0.004, interval=0.002)
            self.assertEqual(failed["status"], "resolution-or-run-failure")
            self.assertIn("cross-resolution vorticity_dense_ratio", failed["failures"])
            incomplete = driver.assess_pair(full[:-1], full, 0.004, interval=0.002)
            self.assertEqual(incomplete["status"], "resolution-or-run-failure")


if __name__ == "__main__":
    unittest.main(argv=[__file__, *unittest_args])
