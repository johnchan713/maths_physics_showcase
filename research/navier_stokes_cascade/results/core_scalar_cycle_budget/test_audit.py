#!/usr/bin/env python3
"""Scalar-cycle norm, finite-stage quantifier and pressure scale-loss tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import ScalePolynomial,cycle_polynomials,exponent_margins,review,tail_upper


class ScalarCycleBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_preslow_norm_polynomial_includes_both_wave_steps_and_fast(self):
        self.controls("cycle_polynomials","preslow_expanded_budget")

    def test_postslow_norm_polynomial_keeps_old_Haar_crosses_and_square(self):
        self.controls("cycle_polynomials","postslow_expanded_budget")

    def test_accumulated_Haar_bound_is_a_finite_sum(self):
        self.controls("cycle_polynomials","Haar_update_is_finite_sum")

    def test_raw_exponent_ledgers_are_bounded_without_new_tail(self):
        self.controls("cycle_polynomials","raw_exponent_budget_below_Pi_without_tail_change")
        self.controls("cycle_polynomials","raw_slow_budget_below_Theta")

    def test_whole_half_line_margins_and_prefix_rates(self):
        self.assertTrue(all(self.result["parts"]["exponent_margins"]["checks"].values()))
        self.assertEqual(exponent_margins(2-Q(2,100000))["next_input"],Q(1,2)+Q(1,100000))

    def test_fresh_linear_row_and_current_scalar_recomputation(self):
        for key in ("nonzero_current_angular_linear_row_divisor","fresh_current_row_cancels_entire_incoming_scalar","independently_recomputed_postslow_scalar"):
            self.controls("current_slow",key)

    def test_old_Haar_crosses_and_nonzero_postslow_scalar(self):
        for key in ("old_Haar_means_are_not_zero","old_Haar_axial_cross_retained","old_Haar_angular_cross_retained","postslow_scalar_not_zero"):
            self.controls("current_slow",key)

    def test_slow_divergence_and_velocity_moments(self):
        self.controls("current_slow","slow_increment_divergence_exact")
        self.controls("current_slow","slow_increment_both_velocity_moments_zero")

    def test_stale_row_mismatch_even_with_exact_divergence(self):
        self.controls("current_slow","stale_angular_row_leaves_exact_scalar_mismatch")
        self.controls("current_slow","stale_row_control_still_has_exact_divergence")

    def test_scalar_and_integrated_residual_norm_budgets(self):
        self.controls("current_slow","postslow_nonlinear_norm_budget")
        self.controls("current_slow","integrated_residual_jet_budget")

    def test_existing_tail_scale_loss_bound_and_tiny_h(self):
        self.controls("scale_tail","positive_bound")
        self.controls("scale_tail","monotone_whole_tail_bound")
        self.controls("scale_tail","given_tail_passes_diagnostic_tolerance")

    def test_no_zero_margin_polynomial_absorption(self):
        self.controls("scale_tail","zero_margin_positive_S_power_grows")
        with self.assertRaises(ValueError):tail_upper(1,2,Q(1,200),0,10)

    def test_invalid_order_and_scale_data_rejected(self):
        for kwargs in ({-1:1},{0:-1}):
            with self.assertRaises(ValueError):ScalePolynomial(kwargs)
        with self.assertRaises(TypeError):ScalePolynomial({0:1.0})
        with self.assertRaises(ValueError):ScalePolynomial(1).evaluate(0)
        with self.assertRaises(ValueError):exponent_margins(1)
        with self.assertRaises(ValueError):tail_upper(1,-1,Q(1,200),1,10)
        with self.assertRaises(ValueError):tail_upper(1,0,Q(1,200),1,0)

    def test_finite_polynomial_multiplication_preserves_scale_powers(self):
        p=ScalePolynomial({1:2});q=ScalePolynomial({2:3})
        self.assertEqual((p*q).record(),{"3":"6"})

    def test_complete_covariance_and_both_pressure_memories_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_historical_bindings_match(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)

    def test_new_budget_not_promoted_to_input_or_infinite_cycle_proof(self):
        self.assertTrue(self.record["flags"]["existing_postslow_order_not_presented_as_new"])
        for key in ("actual_seed_wave_mean_inverse_and_other_defect_inputs_independently_verified","uniform_infinite_stage_or_all_order_constant_proved","actual_first_band_or_full_residual_certified","actual_profile_or_accepted_tail_changed","unconditional_pulse_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_tampering_is_rejected(self):
        r=deepcopy(self.record);r["flags"]["actual_first_band_or_full_residual_certified"]=True
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_record_reproduces(self):
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_existing_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
