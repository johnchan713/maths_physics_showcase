from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from pressure_comparison import (
    current, historical, limit_review, pressure_review, q9,
    radius_review, review,
)


class CurrentHeadQ10PressureComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = review()
        cls.record = audit.generate()

    def test_exact_Q9_energy_and_force_parent(self):
        scope = q9.review()["scope"]
        self.assertTrue(
            scope[
                "exact_Q1_Q8_preterminal_energy_identity_replayed_conditionally"
            ]
        )
        self.assertTrue(
            scope[
                "separate_finite_energy_and_total_dissipation_bounds_obtained_conditionally"
            ]
        )
        self.assertFalse(scope["Q8_endpoint_uniform_all_order_bound_proved_here"])

    def test_forced_pressure_equations_and_single_cancellation(self):
        row = pressure_review()
        self.assertIn("-div f", row["individual_pressures"])
        self.assertTrue(row["checks"]["Q9_common_force_cancels_once"])

    def test_finite_energy_removes_harmonic_gradient(self):
        row = pressure_review()
        self.assertTrue(row["checks"]["harmonic_Hminus_three_gradient_vanishes"])
        self.assertIn("v=t e1", row["harmonic_counterexample"])

    def test_pressure_threshold_is_strict(self):
        self.assertTrue(
            pressure_review()["checks"][
                "L1_stress_requires_Hminus_s_strictly_above_three_halves"
            ]
        )

    def test_commutator_and_radius_exponents(self):
        row = radius_review()
        self.assertEqual(row["commutator_exponent"], Q(4, 3))
        self.assertEqual(row["commutator_radial_coefficient"], 4)
        self.assertEqual(row["commutator_norm_decay"], Q(3, 4))
        self.assertEqual(
            [item["Young_decay"] for item in row["flux_rows"]],
            [4, Q(14, 5)],
        )
        self.assertEqual(row["worst_decay"], 2)

    def test_support_separation_validation(self):
        self.assertEqual(historical.separated_radius(2, 1), 1)
        with self.assertRaises(ValueError):
            historical.separated_radius(1, 1)

    def test_limit_order_and_fixed_interval_constants(self):
        row = limit_review()
        self.assertTrue(all(row["checks"].values()))
        self.assertIn("R->infinity", row["limit_order"])
        self.assertIn("fix T<1", row["Q9_limit_order"])

    def test_current_core_is_same_conditional_chain(self):
        row = current.current_core_review()
        self.assertTrue(all(row["checks"].values()))
        self.assertIn("tau^-A", row["bound"])

    def test_review_conclusion_and_scope(self):
        self.assertEqual(len(self.result["checks"]), 13)
        self.assertTrue(all(self.result["checks"].values()))
        self.assertIn(
            "v=u", self.result["conditional_conclusion"]["comparison"]
        )
        scope = self.result["scope"]
        self.assertTrue(scope["Q9_expanding_ball_comparison_replayed_conditionally"])
        self.assertFalse(scope["finite_energy_hypothesis_removed"])
        self.assertFalse(
            scope["unconditional_admissible_force_exclusion_or_blowup_proved"]
        )

    def test_endpoint_obstruction_remains(self):
        self.assertIn("sup_", self.result["blocking_inequality"])
        self.assertFalse(
            self.result["scope"]["Q8_endpoint_uniform_all_order_bound_proved_here"]
        )

    def test_dependency_closure_and_controls(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 6)
        self.assertTrue(all(self.record["gates"].values()))
        self.assertGreaterEqual(self.record["graph"]["nodes"], 248)

    def test_scope_tampering_is_rejected(self):
        changed = deepcopy(self.record)
        changed["review"]["scope"][
            "unconditional_admissible_force_exclusion_or_blowup_proved"
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
