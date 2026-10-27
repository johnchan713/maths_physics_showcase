from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from pressure_comparison import (
    core, energy, historical, limit_review, pressure_review, radius_review, review,
)


class LogReservedPressureComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_same_energy_and_force_parent(self):
        scope = energy.review()["scope"]
        self.assertTrue(scope["same_field_preterminal_energy_identity_replayed_conditionally"])
        self.assertTrue(energy.parent.review()["scope"]["one_smooth_compact_force_extension_constructed_conditionally"])

    def test_forced_pressure_equations_and_single_cancellation(self):
        row = pressure_review()
        self.assertIn("-div f", row["individual_pressures"])
        self.assertTrue(row["checks"]["common_force_cancels_exactly_once"])

    def test_finite_energy_removes_harmonic_gradient(self):
        row = pressure_review()
        self.assertTrue(row["checks"]["harmonic_Hminus_three_gradient_vanishes"])
        self.assertIn("v=t e1", row["harmonic_counterexample"])

    def test_pressure_threshold_is_strict(self):
        self.assertTrue(pressure_review()["checks"]["L1_stress_requires_Hminus_s_strictly_above_three_halves"])

    def test_commutator_and_radius_exponents(self):
        row = radius_review()
        self.assertEqual(row["commutator_exponent"], Q(4, 3))
        self.assertEqual(row["commutator_radial_coefficient"], 4)
        self.assertEqual(row["commutator_norm_decay"], Q(3, 4))
        self.assertEqual([x["Young_decay"] for x in row["flux_rows"]], [4, Q(14, 5)])
        self.assertEqual(row["worst_decay"], 2)

    def test_support_separation_validation(self):
        self.assertEqual(historical.separated_radius(2, 1), 1)
        with self.assertRaises(ValueError):
            historical.separated_radius(1, 1)

    def test_fixed_time_radius_core_limit_order(self):
        row = limit_review()
        self.assertIn("fix T<1", row["limit_order"])
        self.assertTrue(all(row["checks"].values()))

    def test_same_field_core_lower_bound_is_only_conditional(self):
        self.assertTrue(core.review()["analytic_conclusions"]["same_complete_diagonal_has_positive_core_lower_bound"])
        self.assertFalse(core.review()["scope"]["unconditional_provider_chain_or_blowup_claimed"])

    def test_review_and_scope(self):
        result = review()
        self.assertTrue(all(result["checks"].values()))
        scope = result["scope"]
        self.assertTrue(scope["same_field_expanding_ball_comparison_replayed_conditionally"])
        self.assertTrue(scope["global_smooth_finite_energy_competitor_exclusion_follows_conditionally"])
        self.assertFalse(scope["unconditional_admissible_force_or_blowup_proved"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_admissible_force_or_blowup_proved"] = True
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
