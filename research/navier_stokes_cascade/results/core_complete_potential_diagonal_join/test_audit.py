from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from diagonal import (
    KAPPA, current, parameters, representative_orders, tuple_loss, cut_loss,
    total_loss, residual_loss, gain, comparison, choose_log_scale,
    scale_review, closure_statement, universal_controls,
)


class CorrectedCompletePotentialDiagonalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_corrected_seed_and_gain(self):
        self.assertEqual(current.rates(0)["regular"], Q(3, 2) - 2 * KAPPA)
        self.assertEqual(current.parameters()["gain"], Q(1, 2) - 3 * KAPPA)

    def test_representative_slope_and_margin(self):
        delta = current.parameters()["gain"]
        first, second = representative_orders(1), representative_orders(2)
        self.assertEqual(second, {key: value + delta for key, value in first.items()})
        self.assertEqual(min(first.values()) - delta, 1)

    def test_complete_fixed_gauge_pressure_and_direct_swirl_present(self):
        rows = representative_orders(1)
        self.assertIn("complete_fixed_gauge_pressure_difference", rows)
        self.assertIn("direct_axisymmetric_swirl", rows)

    def test_losses_reserve_curl_and_total_state(self):
        self.assertEqual(cut_loss(3), tuple_loss(4))
        self.assertEqual(total_loss(3), cut_loss(3) + 1)
        self.assertEqual(residual_loss(0), Q(251, 100))

    def test_default_fixed_comparison_stage(self):
        row = comparison()
        self.assertEqual(row["comparison_stage"], 12816)
        self.assertTrue(all(row["checks"].values()))

    def test_multiple_target_comparisons(self):
        for h in (Q(1, 200), Q(9, 1000), Q(1, 10000)):
            for m, target in ((0, 0), (1, 3), (5, 11), (20, 2)):
                self.assertTrue(all(comparison(m, target, h)["checks"].values()))

    def test_scale_majorant_absorbs_extreme_finite_constants(self):
        gamma = gain(1) / 2
        chosen = choose_log_scale(1, gamma, [(10**60, 10**30)])
        margin = Q(10**60) + Q(10**30)**2 / (2 * gamma) - gamma * chosen / 2 + 1
        self.assertLess(margin, 0)
        self.assertTrue(all(scale_review()["checks"].values()))

    def test_invalid_inputs_rejected(self):
        for h in (0, -1, Q(1, 100)):
            with self.assertRaises(ValueError):
                parameters(h)
        for constraints in ([], [(-1, 0)], [(0, -1)]):
            with self.assertRaises(ValueError):
                choose_log_scale(1, 1, constraints)

    def test_complete_cut_curl_and_residual_fixture(self):
        review = universal_controls().polynomial_review()
        self.assertTrue(review["checks"]["complete_cut_velocity_exactly_divergence_free"])
        self.assertEqual(review["nonzero_individually_necessary_terms"],
                         ["new_old", "new_square", "old_new", "pressure", "time", "viscosity"])

    def test_scope_is_conditional_and_global_obligations_open(self):
        scope = closure_statement()
        self.assertTrue(scope["one_locally_finite_divergence_free_field_follows_conditionally"])
        self.assertTrue(scope["all_fixed_jets_of_complete_local_residual_are_flat_conditionally"])
        self.assertTrue(scope["protected_inner_velocity_and_heat_exterior_preserved_conditionally"])
        self.assertFalse(scope["primitive_analytic_contracts_independently_proved"])
        self.assertFalse(scope["closed_time_localization_forcing_energy_global_comparison_or_blowup_proved"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["closure"]["primitive_analytic_contracts_independently_proved"] = True
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
