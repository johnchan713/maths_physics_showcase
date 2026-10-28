from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from preterminal_energy import (
    combined_budget_counterexample, fixed_time_budget, q8,
    regularized_step, review,
)


class CurrentHeadQ9PreterminalEnergyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_exact_Q8_parent_field_and_force(self):
        scope = q8.review()["scope"]
        self.assertTrue(scope["exact_Q1_Q7_field_localized_conditionally"])
        self.assertTrue(
            scope["one_smooth_compact_force_extension_constructed_conditionally"]
        )
        self.assertFalse(scope["Q7_endpoint_uniform_all_order_bound_proved_here"])

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
        self.assertEqual(
            row["lower_bound_for_sup_y2_plus_twice_integral_y2"],
            Q(37349, 22350),
        )
        self.assertGreater(row["strict_gap"], 0)

    def test_exact_identity_and_separate_limits(self):
        self.assertEqual(len(self.result["checks"]), 12)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertIn("+D(t)", self.result["identity"]["differential"])
        self.assertIn("2 integral", self.result["fixed_time"])
        self.assertIn("sup_", self.result["endpoint"]["energy"])
        self.assertIn("/2", self.result["endpoint"]["dissipation"])

    def test_invalid_combined_terminal_strengthening_is_rejected(self):
        self.assertIn(
            "sup_", self.result["invalid_strengthening"]["rejected"]
        )
        self.assertTrue(
            self.result["invalid_strengthening"]["historical_value_matches"]
        )
        self.assertTrue(
            self.result["invalid_strengthening"][
                "diagnostic_is_not_a_Navier_Stokes_solution"
            ]
        )

    def test_endpoint_obstruction_is_unchanged(self):
        self.assertIn(
            "sup_", self.result["blocking_inequality"]
        )
        self.assertTrue(
            self.result["checks"]["Q8_endpoint_premise_remains_conditional"]
        )

    def test_scope_does_not_promote_comparison(self):
        scope = self.result["scope"]
        self.assertTrue(
            scope[
                "exact_Q1_Q8_preterminal_energy_identity_replayed_conditionally"
            ]
        )
        self.assertFalse(scope["fixed_gauge_pressure_recovery_replayed_here"])
        self.assertFalse(
            scope["expanding_ball_comparison_or_competitor_exclusion_proved"]
        )

    def test_dependency_closure_and_controls(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertGreaterEqual(self.record["graph"]["nodes"], 247)

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["review"]["scope"][
            "terminal_velocity_pressure_or_energy_trace_proved"
        ] = True
        with self.assertRaises(ValueError):
            audit.validate(changed)

    def test_frozen_record_and_no_overwrite(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
