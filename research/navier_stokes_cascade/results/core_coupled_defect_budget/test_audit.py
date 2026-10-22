#!/usr/bin/env python3
"""Complete-current moments, coupled norm budgets and honest-scope tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import geometry,order,review,wave_budget,fast_budget,slow_budget


class CoupledDefectBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_complete_wave_changes_recompute_both_current_moments(self):
        for fragment in ("current_wave", "all_six_covariance", "cross_plus_square"):
            self.controls("complete_wave",fragment)

    def test_axial_covariance_and_extra_radial_axial_derivative_retained(self):
        for fragment in ("axial_velocity_covariance", "radial_axial_derivative"):
            self.controls("complete_wave",fragment)

    def test_wave_norm_budgets_at_higher_orders(self):
        for fragment in ("wave_pressure_budget", "wave_axial_scalar_budget"):
            self.controls("complete_wave",fragment)

    def test_fast_zero_Haar_only_removes_linear_terms(self):
        for fragment in ("zero_Haar", "complete_nonlinear_change", "correlation_nonzero", "products_of_Haar"):
            self.controls("fast_correlations",fragment)

    def test_fast_plain_weighted_and_axial_moments_independently_recomputed(self):
        for fragment in ("plain_moment_exact", "weighted_source_moment_exact", "current_axial_scalar_exact"):
            self.controls("fast_correlations",fragment)

    def test_fast_compact_pressure_keeps_whole_moving_moment(self):
        for fragment in ("pressure_compact", "moving_pressure_moment", "compensator_derivative_nonzero"):
            self.controls("fast_correlations",fragment)

    def test_fast_coupled_norm_budgets(self):
        for fragment in ("pressure_norm_bound", "axial_norm_bound"):
            self.controls("fast_correlations",fragment)

    def test_slow_current_remainders_follow_full_source_difference(self):
        for fragment in ("full_source_before", "remainder_moment", "difference_after_linear_row"):
            self.controls("slow_remainder",fragment)

    def test_slow_divergence_and_both_velocity_moments(self):
        for fragment in ("divergence_exact", "both_velocity_moments_zero"):
            self.controls("slow_remainder",fragment)

    def test_radial_viscosity_boundary_cancellation_preserves_axial_viscosity(self):
        for fragment in ("radial_viscosity", "axial_viscosity_retained_nonzero"):
            self.controls("slow_remainder",fragment)

    def test_slow_pressure_compact_and_entire_moving_product_kept(self):
        for fragment in ("pressure_compact", "whole_slow_pressure_moment", "moving_product_derivative_required"):
            self.controls("slow_remainder",fragment)

    def test_reduced_scalar_and_full_axial_moment_norm_budgets(self):
        for fragment in ("reduced_postslow", "whole_reduced_axial_residual_bound"):
            self.controls("slow_remainder",fragment)

    def test_reduced_polynomials_dominate_unreduced_epsilon_budgets(self):
        self.controls("polynomial_budget","unreduced_below_polynomial")
        self.controls("polynomial_budget","whole_half_line_exponent")

    def test_both_current_wave_steps_and_fast_couple_preslow_inputs(self):
        self.controls("polynomial_budget","both_wave_and_fast")
        self.controls("polynomial_budget","all_coupled_outputs_supply_next_cycle_order")

    def test_four_extra_incoming_jets_are_distinct_from_field_input_order(self):
        self.controls("polynomial_budget","scalar_input_target_order")
        self.controls("polynomial_budget","axial_moment_input_target_order")
        self.assertTrue(self.record["flags"]["axial_moment_requires_up_to_four_extra_incoming_target_jets"])

    def test_invalid_geometry_and_derivative_orders_rejected(self):
        for a,d in ((0,2),(2,1),(1,1)):
            with self.assertRaises(ValueError):geometry(a,d)
        for m in (-1,True,Q(1,2)):
            for fn in (order,lambda x:wave_budget(x,{},0,0),lambda x:fast_budget(x,{},0,0,0,0),lambda x:slow_budget(x,{},0,0,0,0,0,0)):
                with self.assertRaises(ValueError):fn(m)

    def test_both_pressure_tails_and_complete_covariance_remain(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_historical_bindings_and_frozen_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_conditional_budget_not_promoted_to_full_residual_proof(self):
        self.assertTrue(self.record["flags"]["existing_postslow_order_not_presented_as_new"])
        for key in ("actual_seed_and_inverse_inputs_independently_verified","uniform_infinite_stage_or_all_order_constant_proved","complete_pointwise_residual_certified","actual_profile_or_accepted_tail_changed","unconditional_pulse_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])
        altered=deepcopy(self.record);altered["flags"]["complete_pointwise_residual_certified"]=True
        with self.assertRaises(ValueError):audit.validate(altered)

    def test_existing_evidence_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
