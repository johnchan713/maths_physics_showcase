#!/usr/bin/env python3
"""Physical row scaling, compact response and all-fixed-order budget tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import background_positive_bound,composition_bound,inverse_budget,reciprocal_bound,responses,review


class SlowInverseBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_three_response_profiles_solve_all_five_signed_rows(self):
        self.controls("matrix_scaling","signed_five_rows")

    def test_fixed_axial_gap_cost_retained(self):
        self.controls("matrix_scaling","axial_gap_cost")
        self.controls("matrix_scaling","angular_limit_nonsingular")
        self.controls("matrix_scaling","axial_limit_singular")
        with self.assertRaises(ValueError):responses(2,1,(1,1,1,1,1))

    def test_all_three_nonzero_physical_rows_have_correct_scale(self):
        self.controls("matrix_scaling","three_physical_rows_have_zero_scale_power")

    def test_omitted_swirl_and_common_target_scale_are_detected(self):
        self.controls("matrix_scaling","omitting_swirl_scale_detected")
        self.controls("matrix_scaling","common_P_Jz_scale_detected")

    def test_potential_carries_its_half_power(self):
        self.controls("matrix_scaling","potential_has_extra_half_power")

    def test_actual_compact_bumps_realize_angular_rows(self):
        for key in ("angular_velocity_moment_zero","pressure_response","axial_response"):
            self.controls("compact_profiles",key)

    def test_actual_compact_bumps_realize_axial_rows(self):
        self.controls("compact_profiles","axial_velocity_moment_zero")
        self.controls("compact_profiles","angular_response")

    def test_radial_potential_satisfies_curl_identity_piecewise(self):
        self.controls("compact_profiles","potential_radial_identity")
        self.controls("compact_profiles","potential_continuity")

    def test_potential_compactness_requires_whole_hull(self):
        for key in ("potential_zero_at_left_hull","potential_zero_at_right_hull","potential_fills_gap_between_bumps","entire_axial_potential_hull_used"):
            self.controls("compact_profiles",key)

    def test_fractional_euler_generator_matches_direct_chart_powers(self):
        self.controls("finite_jets","fractional_generator_equals_direct_power")
        self.controls("finite_jets","wrong_unshifted_generator_detected")

    def test_generalized_profile_norms_include_fractional_and_radial_costs(self):
        self.controls("finite_jets","generalized_profile_majorant")
        self.assertEqual(composition_bound(Q(1,2),2,0,[3]),6)
        self.assertEqual(composition_bound(Q(-3,4),1,0,[3]),24)

    def test_background_composition_supplies_actual_reciprocal_input_bound(self):
        self.controls("finite_jets","actual_background_composition_bound")
        self.assertEqual(background_positive_bound(0,[]),0)

    def test_reciprocal_jets_independently_multiply_to_one(self):
        self.controls("finite_jets","reciprocal_identity")
        self.controls("finite_jets","finite_reciprocal_majorant")

    def test_no_small_high_norm_requirement_for_finite_reciprocal(self):
        self.controls("finite_jets","no_high_norm_smallness_needed")
        self.assertEqual(reciprocal_bound(0,Q(1,3),10**9),3)
        self.assertEqual(reciprocal_bound(2,2,4),Q(7,2))

    def test_explicit_inverse_budget_preserves_scale_polynomials(self):
        self.controls("coupled_substitution","explicit_inverse_polynomial_substitution")
        self.controls("coupled_substitution","reciprocal_polynomial_evaluation")

    def test_explicit_constants_feed_both_coupled_scalar_channels(self):
        self.controls("coupled_substitution","closed_coupled_scalar_polynomials")

    def test_radial_curl_extra_jet_and_moment_target_cost_remain(self):
        self.controls("coupled_substitution","field_vs_target_jet_ledger")
        r=inverse_budget(2,1,2,3,4,5,6,7,8,9,10)
        self.assertEqual(r,dict(U=53,Z=45,X=360))

    def test_invalid_majorant_and_order_inputs_are_rejected(self):
        for args in ((0,0,0),(1,-1,2),(1,1,-2)):
            with self.assertRaises(ValueError):reciprocal_bound(*args)
        with self.assertRaises(TypeError):composition_bound(.5,0,0,[1])
        with self.assertRaises(ValueError):composition_bound(1,0,2,[1])
        with self.assertRaises(ValueError):composition_bound(1,0,0,[-1])
        with self.assertRaises(ValueError):background_positive_bound(2,[1])
        with self.assertRaises(ValueError):background_positive_bound(1,[-1])
        with self.assertRaises(ValueError):inverse_budget(-1,*([1]*10))

    def test_full_pressure_memory_and_covariance_are_not_simplified(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_historical_evidence_and_source_bindings_match(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_scope_preserves_retained_seed_background_and_source_obligations(self):
        self.assertTrue(self.record["flags"]["qualitative_inverse_existence_not_presented_as_new"])
        for key in ("actual_core_background_floor_or_canonical_seed_independently_verified","actual_sourced_and_fast_inverse_inputs_independently_verified","uniform_infinite_stage_or_all_order_constant_proved","actual_profile_tilt_width_or_accepted_tail_changed","complete_residual_unconditional_pulse_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])
        changed=deepcopy(self.record);changed["flags"]["complete_residual_unconditional_pulse_or_blowup_proved"]=True
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_existing_frozen_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
