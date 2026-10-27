from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from comparison import (
    force, young_decay, separated_radius, interface_review, pressure_review,
    commutator_review, radius_review, limit_review, closure_statement,
)


class CertifiedCompletePotentialComparisonJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = audit.generate()

    def test_same_certified_force_is_used(self):
        self.assertTrue(force.closure_statement()["same_certified_complete_diagonal_field_localized"])
        self.assertTrue(interface_review()["checks"]["certified_force_parent_used"])

    def test_finite_energy_is_essential(self):
        review = pressure_review()
        self.assertTrue(interface_review()["checks"]["finite_energy_is_essential_to_remove_harmonic_pressure_gradients"])
        self.assertIn("v=t e1", review["harmonic_counterexample"])

    def test_pressure_recovery_keeps_forced_individual_equations(self):
        review = pressure_review()
        self.assertIn("-div f", review["forced_individual"])
        self.assertTrue(review["checks"]["force_divergence_cancels_once_and_is_not_subtracted_again"])

    def test_pressure_threshold_is_strict(self):
        self.assertTrue(pressure_review()["checks"]["L1_stress_embeds_in_H_minus_s_only_for_s_strictly_above_three_halves"])

    def test_commutator_near_and_far_exponents(self):
        review = commutator_review()
        self.assertEqual(review["near_radial_power"], Q(-2, 3))
        self.assertEqual(review["far_radial_power"], Q(-2))
        self.assertEqual(review["radial_coefficient"], 4)
        self.assertEqual(review["norm_radius_decay"], Q(3, 4))

    def test_cutoff_power_failure_control(self):
        checks = commutator_review()["checks"]
        self.assertTrue(checks["eighth_cutoff_power_is_sufficient_for_L_four_thirds_pairing"])
        self.assertTrue(checks["seventh_power_is_insufficient_for_this_pairing"])

    def test_young_exponents(self):
        self.assertEqual(young_decay(Q(3, 2), 1), 4)
        self.assertEqual(young_decay(Q(3, 4), Q(7, 4)), Q(14, 5))

    def test_invalid_young_exponents_rejected(self):
        for args in ((0, 1), (2, 1), (1, 0)):
            with self.assertRaises(ValueError): young_decay(*args)

    def test_support_separation(self):
        self.assertEqual(separated_radius(2, 1), 1)
        for args in ((Q(1, 2), 0), (1, 1)):
            with self.assertRaises(ValueError): separated_radius(*args)

    def test_compact_mixed_moment_not_discarded(self):
        review = radius_review()
        self.assertTrue(review["checks"]["mixed_moment_is_not_assumed_isotropic"])
        self.assertIn("R^-3", review["compact_pressure"])

    def test_complete_remainder_is_R_minus_two(self):
        review = radius_review()
        self.assertEqual(review["worst_decay"], 2)
        self.assertEqual([row["Young_decay"] for row in review["rows"]], [4, Q(14, 5)])

    def test_limit_order_and_scope(self):
        checks = limit_review()["checks"]
        self.assertTrue(checks["constants_need_not_be_uniform_as_T_approaches_one"])
        self.assertTrue(checks["equality_on_every_fixed_interval_gives_v_equals_u_for_all_t_below_one"])

    def test_conditional_not_unconditional_exclusion(self):
        scope = closure_statement()
        self.assertTrue(scope["global_smooth_bounded_energy_competitor_exclusion_follows_conditionally_from_inherited_core_growth"])
        self.assertFalse(scope["primitive_field_force_and_positive_core_contracts_independently_proved"])
        self.assertFalse(scope["unconditional_admissible_force_or_blowup_proved"])

    def test_frozen_dependencies_and_sources(self):
        self.assertEqual(len(self.record["dependency_evidence_sha256"]), 5)
        self.assertGreater(len(self.record["inherited_source_checks"]), 35)
        self.assertTrue(all(self.record["inherited_source_checks"].values()))

    def test_frozen_record_and_tampering(self):
        frozen = json.loads((audit.HERE / "evidence.json").read_text())
        audit.validate(frozen)
        bad = deepcopy(frozen)
        bad["review"]["scope"]["unconditional_admissible_force_or_blowup_proved"] = True
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
