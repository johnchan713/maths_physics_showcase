#!/usr/bin/env python3
"""Current fast-time error, actual axial velocity and derivative ledger tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest

import audit
from proof_checks import fast_majorants,flat_budget_data,map_constants,particular,review,slow_count,temporal_constant


class FastInverseBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=review();cls.record=audit.generate()

    def controls(self,part,fragment):
        rows={k:v for k,v in self.result["parts"][part]["checks"].items() if fragment in k}
        self.assertTrue(rows);self.assertTrue(all(rows.values()))

    def test_actual_divisors_have_nonzero_integer_norm(self):
        self.controls("actual_divisors","nonzero_integer_norm")
        self.controls("actual_divisors","conjugate_upper_bound")
        self.assertEqual(self.result["parts"]["actual_divisors"]["divisors"],240)

    def test_explicit_temporal_constant_and_whole_lattice_sum(self):
        self.controls("actual_divisors","explicit_fourier_lattice_constant")
        self.controls("actual_divisors","whole_lattice_sum_bound")
        self.assertEqual(temporal_constant(0),384)
        self.assertEqual(temporal_constant(3),3072)

    def test_raw_to_normalized_slow_norm_count_is_retained(self):
        self.assertEqual([slow_count(m) for m in (0,1,3)],[1,3,10])

    def test_compact_oscillatory_inverse_is_independently_recomputed(self):
        self.controls("formal_current_fast","particular_polynomial_solves")
        self.controls("formal_current_fast","whole_line_integral_is_graph_homogeneous")
        self.controls("formal_current_fast","full_compact_radial_identity")

    def test_actual_A1_velocity_is_not_discarded(self):
        self.controls("formal_current_fast","actual_A1_velocity_nonzero")
        self.controls("formal_current_fast","omitting_actual_A1_velocity_detected")

    def test_fast_time_error_is_A1_of_current_residual(self):
        self.controls("formal_current_fast","actual_fast_time_error_is_A1_of_current_residual")
        self.controls("formal_current_fast","current_residual_A1_error_nonzero")
        self.controls("formal_current_fast","wrong_fast_time_error_sign_detected")

    def test_commutation_removes_temporal_inverse_factor(self):
        self.controls("formal_current_fast","commutation_removes_temporal_inverse_factor")
        self.controls("explicit_budgets","flat_error_uses_current_residual_directly")
        self.controls("explicit_budgets","flat_constant_has_no_c_inverse")

    def test_full_moving_cutoff_divergence_and_omission_control(self):
        self.controls("formal_current_fast","full_moving_cutoff_divergence_exact")
        self.controls("formal_current_fast","moving_cutoff_Z_term_is_required")

    def test_both_compact_endpoints_cancel(self):
        self.controls("formal_current_fast","left_compact_endpoint_cancels")
        self.controls("formal_current_fast","right_compact_endpoint_cancels")
        self.controls("formal_current_fast","cutoff_zero_one_and_nonconstant_slow_part")

    def test_nonzero_zero_mode_moment_is_not_flat(self):
        self.controls("formal_current_fast","zero_mode_nonzero_moment_not_flat")

    def test_profile_constants_keep_cutoff_and_positive_hull_costs(self):
        self.assertEqual(map_constants(0,1,Q(1,2),2,2,3,4,5),dict(T=7,A=5))
        self.assertEqual(map_constants(1,1,Q(1,2),2,2,3,4,5),dict(T=56,A=40))

    def test_actual_fast_norms_keep_S_and_A1_factors(self):
        self.controls("explicit_budgets","explicit_actual_fast_norm_substitution")
        self.controls("explicit_budgets","both_scalar_channels_keep_actual_A1")
        self.assertEqual(fast_majorants(0,2,3,4,5,6,7,1),dict(Ft=26880,Fr=32256))

    def test_fast_scalar_channels_retain_five_and_six_extra_residual_jets(self):
        self.controls("explicit_budgets","fast_scalar_input_orders")

    def test_joined_slow_scalar_and_whole_axial_moment_costs(self):
        self.controls("explicit_budgets","joined_slow_scalar_input_orders")
        self.controls("explicit_budgets","joined_axial_moment_input_order")

    def test_flat_error_input_derivatives_are_explicit(self):
        r=flat_budget_data(2,3,1,2,3)
        self.assertEqual(r["input_order"],11)
        self.assertGreater(r["constant"],0)

    def test_strict_exponent_margin_required_for_flatness(self):
        self.controls("explicit_budgets","strict_exponent_margin")
        self.controls("explicit_budgets","smoothness_alone_not_a_uniform_band_bound")

    def test_invalid_orders_floors_and_profile_majorants_rejected(self):
        for m in (-1,True,Q(1,2)):
            with self.assertRaises(ValueError):temporal_constant(m)
        with self.assertRaises(ValueError):flat_budget_data(0,0,1,1,1)
        with self.assertRaises(ValueError):flat_budget_data(0,1,0,1,1)
        with self.assertRaises(ValueError):map_constants(0,0,1,1,1,1,1,1)
        with self.assertRaises(ValueError):map_constants(0,1,-1,1,1,1,1,1)
        with self.assertRaises(ValueError):fast_majorants(0,0,1,1,1,1,1,1)
        with self.assertRaises(ValueError):particular(1,0)

    def test_full_pressure_memory_and_all_covariances_retained(self):
        self.assertTrue(all(self.result["parts"]["pressure_tails"]["checks"].values()))
        self.controls("complete_covariance","all_six_covariances_nonzero")
        self.controls("complete_covariance","omitted_axial_velocity_detected")

    def test_historical_bindings_and_frozen_record_reproduce(self):
        self.assertTrue(all(self.record["gates"].values()))
        self.assertEqual(len(self.record["dependency_evidence_sha256"]),6)
        self.assertEqual(json.loads(Path(__file__).with_name("evidence.json").read_text()),self.record)

    def test_operator_bounds_do_not_certify_actual_inputs_or_complete_residual(self):
        self.assertTrue(self.record["flags"]["qualitative_inverse_existence_and_flatness_not_presented_as_new"])
        for key in ("actual_canonical_seed_and_current_residual_input_norms_verified","actual_sourced_inverse_inputs_independently_verified","uniform_infinite_stage_or_all_order_constant_proved","actual_profile_tilt_width_or_accepted_tail_changed","complete_residual_unconditional_pulse_or_blowup_proved"):
            self.assertFalse(self.record["flags"][key])

    def test_scope_promotion_is_rejected(self):
        changed=deepcopy(self.record);changed["flags"]["complete_residual_unconditional_pulse_or_blowup_proved"]=True
        with self.assertRaises(ValueError):audit.validate(changed)

    def test_frozen_evidence_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"evidence.json";p.write_text("sentinel")
            with self.assertRaises(SystemExit):audit.main(["--output",str(p)])
            self.assertEqual(p.read_text(),"sentinel")


if __name__=="__main__":unittest.main()
