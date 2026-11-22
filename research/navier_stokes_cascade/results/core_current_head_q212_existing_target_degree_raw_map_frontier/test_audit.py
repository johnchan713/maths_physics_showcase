import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q212


class Q212Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q212.statuses().values()))

    def test_finite_degree_max(self):
        self.assertEqual(q212.finite_degree_max([0, Q(3, 2), 7]), 7)
        self.assertEqual(q212.finite_degree_max([Q(1, 3)]), Q(1, 3))

    def test_invalid_degree_requests(self):
        for values in ([], [-1], [True], ["not-rational"]):
            with self.assertRaises(ValueError):
                q212.finite_degree_max(values)

    def test_target_degree_join(self):
        row = q212.target_degree_join()
        self.assertTrue(row["conditional_degree_array_exists"])
        self.assertFalse(row["Q160_displayed_degree_zero_claimed"])
        self.assertFalse(row["degree_array_numerically_evaluated"])
        self.assertIn("S^(d_(j,I))", row["expanded_Q160_target"])

    def test_complete_pressure_memory(self):
        rows = q212.pressure_memory()
        self.assertEqual(len(rows), 8)
        self.assertTrue(any("six covariance" in row for row in rows))
        self.assertTrue(any("moving pressure" in row for row in rows))

    def test_corrected_frontier(self):
        row = q212.corrected_frontier()
        self.assertTrue(row["Q211_next_step_repeats_existing_conditional_target_theorem"])
        self.assertTrue(row["conditional_target_primitive_degree_obligation_closed"])
        self.assertFalse(row["raw_selected_map_majorant_independently_proved"])
        self.assertFalse(row["common_gap_and_inverse_archive_independently_proved"])

    def test_review_scope(self):
        row = q212.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["conditional_target_primitive_degree_existence_joined"])
        self.assertFalse(row["scope"]["actual_degree_array_evaluated"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_identity_guards(self):
        for call in (q212.statuses, q212.target_degree_join, q212.pressure_memory,
                     q212.corrected_frontier, q212.review, q212.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q212.finite_degree_max([0], identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 460)
        self.assertGreaterEqual(record["closure"]["edges"], 1690)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6190)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
