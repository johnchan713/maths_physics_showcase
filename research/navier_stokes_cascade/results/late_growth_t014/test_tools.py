#!/usr/bin/env python3
"""Small regression checks for evidence stitching; no simulation is performed."""
import unittest
import hashlib
from pathlib import Path
import tempfile

from extend import stitch, inherit_history
from analyze import scientific_status
from archive import verify_parts


class ObservationClockTests(unittest.TestCase):
    def setUp(self):
        self.rows = [{"time": i * .01, "step": i * 10, "value": 1 + i} for i in range(4)]

    def test_valid_restart_preserves_all_rows(self):
        self.assertEqual(stitch(self.rows[:2], [self.rows[1:3], self.rows[2:]], .03), self.rows)

    def test_changed_boundary_is_rejected(self):
        segment = [dict(self.rows[1], value=99), self.rows[2]]
        with self.assertRaisesRegex(ValueError, "boundary"):
            stitch(self.rows[:2], [segment], .02)

    def test_missing_observation_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            stitch(self.rows[:2], [[self.rows[1], self.rows[3]]], .03)

    def test_shifted_clock_is_rejected(self):
        segment = [self.rows[1], dict(self.rows[2], time=.021)]
        with self.assertRaisesRegex(ValueError, "clock changed"):
            stitch(self.rows[:2], [segment], .02)

    def test_nonadvancing_step_is_rejected(self):
        segment = [self.rows[1], dict(self.rows[2], step=10)]
        with self.assertRaisesRegex(ValueError, "step clock"):
            stitch(self.rows[:2], [segment], .02)

    def test_duplicate_endpoint_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "boundary"):
            stitch(self.rows[:2], [[self.rows[1]]], .01)


class InheritedHistoryTests(unittest.TestCase):
    def setUp(self):
        self.rows = [{"time": t} for t in (.00, .01, .02, .03)]
        self.budgets = [{"time": t, "rate": 1 + t} for t in (.01, .02, .03)]
        self.parent = {"initial_state_sha256": "same-initial", "settings": {"grid": 64, "start_time": .02},
                       "budgets": self.budgets[-1:], "checkpoints": [{"time": .03, "identity": "new-state"}]}
        self.joined = {"initial_state_sha256": "same-initial", "budgets": self.budgets,
                       "checkpoints": [{"time": .02, "identity": "old-state"}]}

    def test_keeps_full_budgets_and_immediate_parent_restart(self):
        inherited, rows = inherit_history(self.joined, self.rows, self.parent, 64, .03)
        self.assertEqual(rows, self.rows)
        self.assertEqual(inherited["budgets"], self.budgets)
        self.assertEqual(inherited["checkpoints"], self.parent["checkpoints"])
        self.assertEqual(inherited["settings"], self.parent["settings"])

    def test_wrong_grid_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Wrong source"):
            inherit_history(self.joined, self.rows, self.parent, 128, .03)

    def test_truncated_history_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Wrong source"):
            inherit_history(self.joined, self.rows[:-1], self.parent, 64, .03)

    def test_mixed_initial_fields_are_rejected(self):
        self.joined["initial_state_sha256"] = "another-initial"
        with self.assertRaisesRegex(ValueError, "initial identity"):
            inherit_history(self.joined, self.rows, self.parent, 64, .03)

    def test_duplicate_budget_is_rejected(self):
        self.joined["budgets"] = self.budgets + self.budgets[-1:]
        with self.assertRaisesRegex(ValueError, "budget clock"):
            inherit_history(self.joined, self.rows, self.parent, 64, .03)

    def test_changed_parent_budget_is_rejected(self):
        self.parent["budgets"] = [{"time": .03, "rate": 999}]
        with self.assertRaisesRegex(ValueError, "Joined budgets differ"):
            inherit_history(self.joined, self.rows, self.parent, 64, .03)


class ScientificStatusTests(unittest.TestCase):
    def test_missing_comparison_is_not_a_pass(self):
        with self.assertRaisesRegex(ValueError, "Missing paired"):
            scientific_status({}, [])

    def test_clean_checks_pass(self):
        self.assertEqual(scientific_status({".12": {"failures": [], "numerical_failures": []}}, []),
                         "passed-preliminary-finite-amplification-continuation-gates")

    def test_stalled_growth_does_not_pass(self):
        self.assertEqual(scientific_status({".12": {"failures": ["growth"], "numerical_failures": []}}, []),
                         "finite-growth-stalled")

    def test_resolution_failure_takes_precedence(self):
        self.assertEqual(scientific_status({".12": {"failures": ["growth", "resolution"],
                                                   "numerical_failures": ["resolution"]}}, []),
                         "resolution-or-sampling-gate-failed")

    def test_sampling_failure_is_not_a_pass(self):
        self.assertEqual(scientific_status({".12": {"failures": [], "numerical_failures": []}}, ["sampling"]),
                         "resolution-or-sampling-gate-failed")


class ArchiveIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="ns-archive-test-", dir=Path(__file__).parent)
        self.addCleanup(self.temporary.cleanup)
        self.folder = Path(self.temporary.name)
        self.payload = b"small independent test payload"
        (self.folder / "part").write_bytes(self.payload)
        digest = hashlib.sha256(self.payload).hexdigest()
        self.index = {"parts": [{"path": "part", "bytes": len(self.payload), "sha256": digest}],
                      "compressed_sha256": digest, "compressed_bytes": len(self.payload)}

    def test_valid_parts_pass(self):
        verify_parts(self.folder, self.index)

    def test_corrupted_part_is_rejected(self):
        (self.folder / "part").write_bytes(b"x" + self.payload[1:])
        with self.assertRaisesRegex(ValueError, "part identity"):
            verify_parts(self.folder, self.index)

    def test_wrong_combined_identity_is_rejected(self):
        self.index["compressed_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "stream changed"):
            verify_parts(self.folder, self.index)

    def test_part_cannot_escape_archive_directory(self):
        self.index["parts"][0]["path"] = "../outside"
        with self.assertRaises(ValueError):
            verify_parts(self.folder, self.index)


if __name__ == "__main__":
    unittest.main()
