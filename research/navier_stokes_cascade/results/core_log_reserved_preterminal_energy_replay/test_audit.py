from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from preterminal_energy import (
    combined_budget_counterexample, fixed_time_budget, parent,
    regularized_step, review,
)


class LogReservedPreterminalEnergyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_parent_localization_and_force_are_same_field(self):
        scope = parent.review()["scope"]
        self.assertTrue(scope["same_field_axisymmetric_localization_replayed_conditionally"])
        self.assertTrue(scope["one_smooth_compact_force_extension_constructed_conditionally"])

    def test_fixed_time_budgets_keep_factor_two(self):
        row = fixed_time_budget(Q(3, 2))
        self.assertEqual(row["joint_budget"], Q(9, 4))
        self.assertEqual(row["energy_budget"], Q(9, 4))
        self.assertEqual(row["dissipation_budget"], Q(9, 8))

    def test_zero_force_budget(self):
        self.assertEqual(fixed_time_budget(0)["joint_budget"], 0)

    def test_regularized_square_root(self):
        row = regularized_step(Q(1, 1000))
        self.assertIn("delta", row["quantity"])
        self.assertEqual(row["limit"], "sqrt(E(t))<=F(t)")

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):
            fixed_time_budget(-1)
        with self.assertRaises(ValueError):
            regularized_step(0)
        for value in (0, 1, Q(3, 2)):
            with self.assertRaises(ValueError):
                combined_budget_counterexample(value)

    def test_exact_counterexample_values(self):
        row = combined_budget_counterexample()
        self.assertEqual(row["force_L1"], Q(201, 200))
        self.assertEqual(row["force_budget_squared"], Q(40401, 40000))
        self.assertEqual(row["lower_bound_for_sup_y2_plus_twice_integral_y2"], Q(37349, 22350))
        self.assertGreater(row["strict_gap"], 0)

    def test_review_keeps_separate_endpoint_bounds(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        self.assertIn("sup_", result["endpoint"]["energy"])
        self.assertIn("/2", result["endpoint"]["dissipation"])
        self.assertIn("sup_", result["invalid_strengthening"]["rejected"])

    def test_scope(self):
        scope = review()["scope"]
        self.assertTrue(scope["same_field_preterminal_energy_identity_replayed_conditionally"])
        self.assertFalse(scope["fixed_gauge_pressure_recovery_replayed_here"])
        self.assertFalse(scope["expanding_ball_comparison_or_competitor_exclusion_proved"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["fixed_gauge_pressure_recovery_replayed_here"] = True
        with self.assertRaises(ValueError):
            audit.validate(bad)

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
