#!/usr/bin/env python3
"""Exact compact compensation, fresh-map rates and deferred input checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import rates,oscillatory_recipe,remaining_row_envelope
from controls import review


class FinalMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        values=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(values,fragment);self.assertTrue(all(values),fragment)

    def test_moving_unit_mass_compensation(self):self.control('full_moving_compensator_has_exact_unit_mass')
    def test_integrated_zero_mode_is_not_pointwise_source_zero(self):self.control('compensated_Haar_source_has_zero_whole_radial_moment')
    def test_full_Haar_pressure_retained(self):self.control('Haar_pressure_field_not_zero_or_order_promoted')
    def test_exact_radial_Haar_residual_and_A0_zero_mode(self):
        self.control('Haar_radial_residual_is_exactly_minus_rho_P');self.control('compensated_A0_zero_mode_exactly_zero')
    def test_compact_nonzero_mode_identity(self):self.control('nonzero_mode_compact_identity_exact')
    def test_pointwise_A0_is_not_discarded(self):
        self.control('pointwise_A0_nonzero_despite_zero_Haar');self.control('full_pointwise_radial_nonzero_residual_is_minus_A0')
    def test_full_compact_pressure_has_both_endpoints(self):self.control('compact_pressure_')
    def test_moving_cutoff_derivatives_remain(self):self.control('moving_cutoff_derivatives_commute_with_graph')
    def test_moving_pressure_compensation_derivatives_remain(self):
        self.control('full_pressure_compensation_stays_zero_after_Z_derivative');self.control('pointwise_pressure_source_needs_rho_Z_P')
    def test_nonunit_or_missing_compensator_is_a_real_obstruction(self):
        self.control('nonunit_compensator_leaves_nonflat_A0_zero_mode');self.control('omitting_compensation_leaves_nonflat_A0_zero_mode')
    def test_actual_radial_moment_includes_U_to_R_Jacobian(self):self.control('e0_U_to_R_jacobian_restores_exact_radial_moment')
    def test_full_A0_input_source_B_includes_covariance_graph_terms(self):
        self.control('full_radial_source_B_derived');self.control('independent_full_radial_source_keeps_covariance_and_graph')
    def test_fresh_pressure_and_viscosity_difference_orders(self):self.control('both_fresh_difference_orders_stronger_than_C')
    def test_independent_physical_fast_size_substitution(self):self.control('independent_fresh_fast_')
    def test_independent_physical_slow_size_substitution(self):self.control('independent_fresh_slow_')
    def test_pressure_is_not_given_extra_epsilon_without_its_axial_derivative(self):self.control('pressure_not_falsely_promoted')
    def test_two_fresh_A1_integrations_suffice(self):self.control('only_two_fresh_A1_integrations_needed')
    def test_A0_integration_order_is_minimal_for_strict_C_margin(self):self.control('strict_A0_margin_and_minimal_integer')
    def test_final_all_three_regular_mean_rows_C(self):self.control('all_three_final_regular_mean_rows_C')
    def test_combined_regular_order_is_D_and_physical_tails_stay_open(self):self.control('combined_regular_still_D_not_full_physical')
    def test_stage_exact_cancellation_does_not_reuse_postwave_row(self):self.control('final_stage_does_not_reuse_postwave_row_as_final')
    def test_sharp_fast_and_slow_recipes_come_from_original_maps(self):
        self.control('sharp_fast_D_is_derived_not_old_B_relabelled');self.control('fresh_slow_C_uses_same_scalar_derived_recipe')
    def test_nested_high_jet_A1_not_truncated(self):
        self.control('nested_original_A1_high_jet_not_discarded');self.control('fast_difference_radial_input_retains_old_source_200009')
    def test_A0_current_full_radial_source_high_jets_not_fabricated(self):
        self.control('A0_source_is_current_final_full_graph_gr');self.control('A0_source_high_jets_deferred_not_sampled_or_allocated')
    def test_A0_uses_its_own_original_raw_kernel_profiles(self):self.control('A0_uses_e0_raw_profiles_not_A1_operator_norm')
    def test_poles_original_width_gap_and_weight_floor_retained(self):
        self.control('source_pole_12p_and_original_floor_width_gap_retained');self.control('current_gr_lower_pole_embedding_precedes_products')
    def test_actual_A1_velocity_and_postslow_P_retained(self):
        self.control('actual_A1_velocities_still_in_fast_increment');self.control('radial_compensator_is_actual_postslow_P_recipe')
    def test_original_moderate_Fourier_constant_matches(self):self.control('moderate_original_Fourier_constant_matches')
    def test_actual_input_and_tail_certificates_not_inferred(self):
        self.control('deferred_recipe_cannot_certify_actual_inputs');self.control('exact_additive_tails_excluded_not_forwarded')
    def test_invalid_range_or_state_rejected(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):remaining_row_envelope(0,{}, {},{}, {})
        with self.assertRaises(ValueError):oscillatory_recipe(0,1,e=0,source='test',source_rate=1,target_rate=2,kappa=Q(1,10),geometry={},pole=1)
    def test_seven_dependencies_no_new_assumptions_and_frozen_evidence(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),7)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_obligation_and_provenance_mutations_rejected(self):
        for k,v in self.record['flags'].items():
            bad=deepcopy(self.record);bad['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['dependency_evidence_sha256'][next(iter(bad['dependency_evidence_sha256']))]='wrong'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_frozen_record_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
