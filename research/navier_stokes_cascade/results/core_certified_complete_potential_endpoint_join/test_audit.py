from copy import deepcopy
from fractions import Fraction as Q
from math import factorial
import json
from pathlib import Path
import tempfile
import unittest
import audit
from endpoint import (
    diagonal, band_window, active_stage_bound, closed_centers, requirements,
    heat_A, heat_S, heat_terminal, pressure_terminal, atlas_review,
    representative_review, punctured_cover_review, exterior_review,
    closure_statement,
)


class CertifiedCompletePotentialEndpointJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_same_certified_complete_diagonal_seed_and_gain(self):
        self.assertEqual(diagonal.current.rates(0)["regular"], Q(3, 2) - 2 * Q(1, 100000))
        self.assertEqual(diagonal.current.parameters()["gain"], Q(1, 2) - 3 * Q(1, 100000))
        self.assertTrue(representative_review()["checks"]["certified_fixed_stage_domains_are_the_diagonal_inputs"])

    def test_full_enlarged_closed_band_window(self):
        self.assertEqual(band_window(Q(2) ** -24, Q(2) ** -20),
                         dict(first=18, last=26, count=9))
        self.assertTrue(all(atlas_review()["checks"].values()))

    def test_closed_terminal_mesh_and_finite_stage_bound(self):
        self.assertEqual(closed_centers(0, 0, Q(1, 10**36)), (-1, 1))
        self.assertEqual(active_stage_bound(Q(2) ** -24, Q(2) ** -16), 7)

    def test_endpoint_derivative_requirements(self):
        for m in (0, 1, 5, 10000):
            row = requirements(m)
            self.assertEqual(row["potential_for_velocity_Cauchy"], m + 2)
            self.assertEqual(row["direct_swirl_for_velocity_Cauchy"], m + 1)
            self.assertEqual(row["pressure_Cauchy"], m + 1)

    def test_invalid_endpoint_orders_and_domains_rejected(self):
        for m in (-1, True, Q(1, 2)):
            with self.assertRaises(ValueError):
                requirements(m)
        for args in ((0, Q(1, 2)), (Q(1, 2), Q(1, 4))):
            with self.assertRaises(ValueError):
                band_window(*args)

    def test_representative_blocker_and_punctured_cover(self):
        self.assertTrue(all(representative_review()["checks"].values()))
        self.assertIn("sup_", representative_review()["blocking_inequality"])
        self.assertTrue(all(punctured_cover_review()["checks"].values()))

    def test_complete_mixed_exterior_jets(self):
        review = exterior_review()
        self.assertEqual(len(review["mixed_terminal_rows"]), 36)
        self.assertTrue(all(review["checks"].values()))

    def test_pressure_keeps_full_second_leibniz_cross(self):
        h = Q(1, 200)
        self.assertEqual(heat_S(h, 2) - 2 * heat_A(h, 2), 2 * heat_A(h, 1) ** 2)
        self.assertLess(pressure_terminal(h, 0, 2), 0)

    def test_zero_radius_factorial_obstruction(self):
        h = Q(1, 200)
        for j in (1, 2, 40):
            self.assertGreaterEqual(heat_A(h, j), h * factorial(j - 1) * factorial(j))
            self.assertGreaterEqual(heat_S(h, j), 2 * heat_A(h, j))
        for j in (1, 1000):
            ratio = 4 * heat_A(h, j + 1) / heat_A(h, j) / (j + 1)
            self.assertEqual(ratio - 4 * j, 4 * h * (2 * j + 1 + h) / (j + 1))

    def test_conditional_scope_and_origin_limit(self):
        scope = closure_statement()
        self.assertTrue(scope["same_certified_complete_diagonal_field_used"])
        self.assertTrue(scope["compatible_one_sided_C_infinity_potential_direct_swirl_velocity_and_pressure_jets_on_punctured_terminal_compacts_follow_conditionally"])
        self.assertFalse(scope["terminal_origin_regularized"])
        self.assertFalse(scope["primitive_all_order_endpoint_contracts_independently_proved"])
        self.assertFalse(scope["localization_initialization_admissible_force_energy_global_comparison_or_blowup_proved"])

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["closure"]["terminal_origin_regularized"] = True
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
