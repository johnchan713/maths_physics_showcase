#!/usr/bin/env python3
"""Finite-jet diagnostics, failure controls and scientific-scope bindings."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import Jet,log_one,majorant_powers,review,tame_bound


class FixedProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_log_coordinate_majorants_at_poles_and_tiny_h(self):
        self.controls("composition","coordinate_log_majorant")

    def test_direct_power_and_Euler_composition_match(self):
        self.controls("composition","direct_power_equals_Euler_composition")

    def test_mixed_radial_and_density_majorants(self):
        self.controls("composition","coefficient_majorant")

    def test_omitting_density_or_radial_generator_shift_is_detected(self):
        self.controls("composition","omitted_density_radial_shift_detected")

    def test_exact_majorant_convolution(self):
        self.controls("composition","constant_term_is_one")
        self.controls("composition","linear_majorant_coefficients")
        self.controls("composition","quadratic_mixed_convolution")
        self.assertEqual(majorant_powers(1)[1].c[(1,0)],Q(32,3))

    def test_weighted_radial_order_zero_factor(self):
        self.assertEqual(tame_bound(1,2,0,[Q(5)]),40)

    def test_all_six_pressure_source_terms_and_mass(self):
        for key in ("exact_six_term_commutator","each_moving_product_term_nonzero","each_term_has_zero_radial_mass"):
            self.controls("commutator",key)

    def test_compact_commutator_pressure_and_derivative(self):
        self.controls("commutator","commutator_pressure_compact")
        self.controls("commutator","commutator_pressure_radial_identity")

    def test_commutator_norm_budget_at_three_orders(self):
        self.controls("commutator","all_six_terms_budget")

    def test_each_pressure_commutator_omission_is_detected(self):
        self.controls("commutator","all_six_term_omissions_detected")

    def test_positive_margins_have_a_finite_common_threshold(self):
        self.controls("margins","finite_common_threshold_is_positive")
        self.controls("margins","all_positive_integer_margins_absorb_fixed_constants")

    def test_zero_margin_and_growing_target_constants_are_not_absorbed(self):
        self.controls("margins","zero_margin_not_absorbed")
        self.controls("margins","unbounded_target_constant_can_cancel_gain")
        self.controls("margins","infinite_orders_need_not_share_positive_threshold")

    def test_invalid_orders_and_inexact_constants_rejected(self):
        for args in ((0,0,-1,[]),(2,0,0,[1]),(0,0,1,[1]),(0,0,0,[-1])):
            with self.assertRaises(ValueError):tame_bound(*args)
        with self.assertRaises(TypeError):tame_bound(0,0,0,[1.0])
        with self.assertRaises(ValueError):log_one(Jet(2,{(0,0):1}))

    def test_all_six_covariances_and_axial_velocity_retained(self):
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_both_complete_background_pressure_tails_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))

    def test_historical_records_and_sources_match(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),4)

    def test_correction_preserves_actual_target_and_full_residual_limits(self):
        self.assertTrue(self.record["flags"]["numerical_profile_required_for_asymptotic_existence_corrected"])
        for key in ("actual_profile_or_accepted_tail_changed","uniformity_over_all_profiles_or_infinite_orders_proved","actual_current_target_higher_jets_verified","complete_residual_improvement_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_tampering_is_rejected(self):
        record=deepcopy(self.record);record["flags"]["actual_current_target_higher_jets_verified"]=True
        with self.assertRaises(ValueError):audit.validate(record)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_existing_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"evidence.json";path.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(path)])
            self.assertEqual(path.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
