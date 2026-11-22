import json
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path
import audit
import q217


class Q217Tests(unittest.TestCase):
    def test_provider_statuses(self):
        self.assertTrue(all(q217.statuses().values()))

    def test_localization_requests(self):
        row = q217.localization_requests(5)
        self.assertEqual(row["direct"], {"potential": 7, "direct_swirl": 6, "pressure": 5})
        self.assertEqual(row["cauchy"], {"potential": 8, "direct_swirl": 7, "pressure": 6})
        self.assertEqual(row["supplied"], {"potential": 9, "direct_swirl": 8, "pressure": 7})
        self.assertTrue(all(value == 1 for value in row["spare_over_cauchy"].values()))
        with self.assertRaises(ValueError):
            q217.localization_requests(-1)

    def test_energy_scope(self):
        row = q217.energy()
        self.assertEqual(row["zero_safe_regularization"], "sqrt(E+epsilon^2)")
        self.assertFalse(row["energy_and_dissipation_extrema_added_under_one_rhs"])
        self.assertFalse(row["terminal_field_value_used"])

    def test_comparison_and_limit_order(self):
        row = q217.comparison()
        self.assertTrue(row["individual_pressures_keep_minus_div_f"])
        self.assertTrue(row["common_force_cancels_exactly_once"])
        self.assertEqual(row["expanding_ball_error_rate"], Q(2))
        self.assertEqual(row["limit_order"], ["fix T<1", "send R to infinity", "send tau to zero"])

    def test_protected_core(self):
        row = q217.protected_core()
        self.assertEqual(row["support_hole_margin"], Q(3, 101))
        self.assertEqual(row["velocity_exponent"], Q(101, 200))
        self.assertEqual(row["background_remainder_exponent"], Q(1, 100))
        self.assertTrue(row["annular_corrections_vanish_exactly"])
        self.assertFalse(row["background_vanishes"])

    def test_finite_energy_obstruction(self):
        row = q217.first_new_conditional_premise()
        example = row["counterexample"]
        self.assertFalse(row["implied_by_local_smoothness_zero_data_and_compact_force"])
        self.assertTrue(example["divergence_free"])
        self.assertTrue(example["zero_initial_data"])
        self.assertTrue(example["navier_stokes_residual_zero"])
        self.assertFalse(example["finite_L2_at_positive_time"])

    def test_review(self):
        row = q217.review()
        self.assertTrue(all(row["checks"].values()))
        self.assertTrue(row["scope"]["finite_energy_competitor_exclusion_implication"])
        self.assertFalse(row["scope"]["all_locally_smooth_competitors_excluded"])
        self.assertFalse(row["scope"]["unconditional_force_residual_or_blowup"])

    def test_identity_guards(self):
        for call in (q217.statuses, q217.energy, q217.comparison, q217.protected_core,
                     q217.first_new_conditional_premise, q217.review, q217.record):
            with self.assertRaises(ValueError):
                call(identity="wrong")
        with self.assertRaises(ValueError):
            q217.localization_requests(1, identity="wrong")

    def test_frozen_evidence(self):
        record = audit.generate()
        self.assertEqual(record, json.loads((audit.HERE / "evidence.json").read_text()))
        self.assertGreaterEqual(record["closure"]["nodes"], 467)
        self.assertGreaterEqual(record["closure"]["edges"], 1750)
        self.assertGreaterEqual(record["closure"]["source_bindings"], 6280)

    def test_fresh_output(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as folder:
            path = Path(folder) / "fresh.json"
            path.write_text(json.dumps(audit.generate(), indent=2, sort_keys=True) + "\n")
            self.assertEqual(json.loads(path.read_text()), audit.generate())


if __name__ == "__main__":
    unittest.main()
