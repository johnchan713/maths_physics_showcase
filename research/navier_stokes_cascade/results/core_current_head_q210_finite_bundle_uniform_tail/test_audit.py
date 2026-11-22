import json
import tempfile
import unittest
from pathlib import Path
import audit
import q210


class Q210Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q210.statuses().values()))

    def test_finite_bundle(self):
        row = q210.finite_bundle()
        self.assertEqual(row["order"], 36479350301)
        self.assertEqual(row["stage"], 602)
        self.assertTrue(row["jet_index_set_finite"])
        self.assertTrue(row["stage_prefix_finite"])
        self.assertFalse(row["finite_maximum_needs_numerical_materialization"])
        self.assertTrue(row["conditional_on_named_continuum_contracts"])
        self.assertFalse(row["contracts_independently_reproved_here"])

    def test_countermodel(self):
        for d in (0, 1, 100, 10000):
            row = q210.all_band_obstruction(d)
            self.assertEqual(row["witness_n"], 8 * (d + 1))
            self.assertEqual(row["witness_ell"], row["witness_n"] ** 4)
            self.assertTrue(row["lower_bound_positive"])
            self.assertTrue(row["absolute_source_tends_to_zero_for_positive_B"])
            self.assertFalse(row["normalized_polynomial_supremum_finite"])
            self.assertFalse(row["counterexample_claimed_for_actual_selected_source"])

    def test_corrected_frontier(self):
        row = q210.corrected_frontier()
        self.assertFalse(row["Q209_numeric_A_R_is_first_obligation"])
        self.assertTrue(row["conditional_finite_bundle_follows_from_existing_contracts"])
        self.assertIn("sup_(ell>=ell*)", row["first_independent_inequality"])
        self.assertFalse(row["per_band_smoothness_suffices"])
        self.assertFalse(row["unconditional_residual_or_blowup"])

    def test_review_scope(self):
        row = q210.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["stage_602_bundle_finite_relative_to_named_contracts"])
        self.assertFalse(row["scope"]["numerical_stage_602_bundle_evaluated"])
        self.assertFalse(row["scope"]["uniform_actual_source_tail_proved_independently"])
        self.assertFalse(row["scope"]["unconditional_pulse_residual_force_exclusion_or_blowup"])

    def test_invalid_inputs(self):
        for call in (q210.statuses, q210.finite_bundle, q210.all_band_obstruction,
                     q210.corrected_frontier, q210.review, q210.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        for bad in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                q210.finite_bundle(order=bad)
            with self.assertRaises(ValueError):
                q210.all_band_obstruction(bad)
        with self.assertRaises(ValueError):
            q210.finite_bundle(stage=-1)

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 459)
        self.assertGreaterEqual(record["closure"]["edges"], 1677)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6168)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
