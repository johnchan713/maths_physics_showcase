from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from localization import (
    TERMS, endpoint, order, force_requirements, derivative_factor,
    first_controlled, domain_review, residual_review, extension_review,
    energy_review, closure_statement,
)


class CertifiedCompletePotentialForceJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_same_certified_endpoint_is_used(self):
        self.assertTrue(endpoint.closure_statement()["same_certified_complete_diagonal_field_used"])
        self.assertTrue(domain_review()["checks"]["certified_endpoint_parent_used"])

    def test_enlarged_chart_margin(self):
        self.assertIn("16q<=3qstar/8", domain_review()["enlarged_chart"])
        self.assertTrue(domain_review()["checks"]["enlarged_chart_margin_strict"])

    def test_force_derivative_requirements(self):
        self.assertEqual(force_requirements(5), dict(
            force_order=5, potential_order=8, direct_swirl_order=7,
            pressure_order=6, Cauchy_potential_order=9,
            Cauchy_direct_swirl_order=8, Cauchy_pressure_order=7))

    def test_invalid_orders_rejected(self):
        for value in (-1, True, Q(1, 2)):
            with self.assertRaises(ValueError): order(value)
            with self.assertRaises(ValueError): force_requirements(value)

    def test_all_twelve_localization_terms_retained(self):
        self.assertEqual(len(TERMS), 12)
        self.assertEqual(len(set(TERMS)), 12)
        self.assertTrue(residual_review()["checks"]["time_viscosity_pressure_and_all_nonlinear_terms_present"])

    def test_complete_pressure_and_force_divergence_retained(self):
        review = residual_review()
        self.assertTrue(review["checks"]["complete_pressure_cutoff_term_retained"])
        self.assertIn("-div f_minus", review["pressure_Poisson"])

    def test_precise_corner_and_punctured_blockers_recorded(self):
        review = residual_review()
        self.assertIn("sup_corner", review["corner_blocker"])
        self.assertIn("punctured terminal compact", review["punctured_blocker"])

    def test_Borel_Leibniz_factor_and_invalid_index(self):
        self.assertEqual(derivative_factor(5, 2, 1), Q(1, 12))
        self.assertEqual(derivative_factor(1, 3, 1), 0)

    def test_every_fixed_mixed_jet_enters_one_diagonal(self):
        for a, m in ((0, 0), (1, 3), (20, 4), (6, 70), (1000, 1000)):
            j = first_controlled(a, m)
            self.assertLessEqual(a + m, j // 2)
            self.assertLess(m, j)

    def test_future_force_does_not_assert_future_solution(self):
        scope = closure_statement()
        self.assertTrue(scope["compact_smooth_force_and_zero_initial_data_follow_conditionally"])
        self.assertFalse(scope["future_velocity_or_pressure_constructed"])

    def test_energy_extrema_correction(self):
        review = energy_review()
        values = review["counterexample"]
        self.assertGreater(values["lower_combined"], values["force_norm"] ** 2)
        self.assertTrue(review["checks"]["separate_extrema_cannot_be_added_with_same_budget"])

    def test_conditional_scope_and_remaining_comparison(self):
        scope = closure_statement()
        self.assertFalse(scope["primitive_corner_and_punctured_bounds_independently_proved"])
        self.assertFalse(scope["global_competitor_comparison_or_blowup_proved"])

    def test_frozen_dependencies_and_sources(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 4)
        self.assertGreater(len(self.record["inherited_source_checks"]), 39)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["global_competitor_comparison_or_blowup_proved"] = True
        with self.assertRaises(ValueError): audit.validate(bad)

    def test_existing_output_not_overwritten(self):
        with tempfile.TemporaryDirectory(dir=audit.HERE) as directory:
            path = Path(directory) / "kept.json"
            path.write_text("keep")
            with self.assertRaises(SystemExit):
                audit.main(["--output", str(path)])
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
