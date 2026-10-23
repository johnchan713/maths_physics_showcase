#!/usr/bin/env python3
"""Full mean rows, current scalars, compact high jets and scope tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import rates,strict_integrations,dependency_trace,operator_ledger,mean_feedback,P
from controls import review


class PostSignedMeanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        values=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(values,fragment);self.assertTrue(all(values),fragment)

    def test_both_complete_current_radial_source_differences(self):self.control('complete_current_radial_source_difference')
    def test_full_tangential_difference_including_pressure_viscosity_and_time(self):self.control('whole_tangential_difference')
    def test_actual_three_component_divergence_and_velocity_moments(self):
        self.control('actual_three_component_divergence');self.control('both_weighted_Haar_velocity_moments')
    def test_zero_Haar_fast_and_zero_slow_temporal_principal(self):
        self.control('all_new_fast_components_zero_Haar');self.control('new_slow_has_zero_temporal_principal_term')
    def test_full_old_fast_correlations_are_not_discarded(self):self.control('old_new_fast_correlations_survive')
    def test_previous_slow_and_nonzero_reserved_G_contract_not_forgotten(self):
        self.control('new_slow_Haar_old_field_is_old_Haar_plus_previous_slow');self.control('nonzero_reserved_axial_background_can_leave_linear_scalar_leaks')
    def test_full_recomputed_pressure_and_all_six_covariance_entries(self):
        self.control('full_recomputed_pressure_increment_derivatives_not_erased');self.control('all_six_complete_covariances_remain_fixed_through_means')
    def test_complete_two_step_identity_is_full_not_Qcirc_only(self):self.control('complete_two_step_identity')
    def test_original_A1_actual_velocity_positive_residual_sign_and_moving_cutoff(self):
        self.control('original_exact_A1_time_error_has_positive_residual_sign');self.control('original_actual_A1_velocity_and_moving_cutoff_are_retained')
    def test_independent_physical_full_mean_row_budget(self):self.control('full_regular_mean_row_bound')
    def test_full_pressure_not_promoted_with_row_gain(self):self.control('full_pressure_not_promoted_with_row_gain')
    def test_three_complete_postfast_scalar_functionals(self):self.control('complete_postfast_scalar_polynomial_substitution')
    def test_preserved_Haar_recipe_keeps_the_previous_slow_increment(self):self.control('preserved_Haar_recipe_keeps_previous_slow')
    def test_new_fast_mean_feedback_uses_derived_D_not_coarse_size(self):self.control('derived_new_fast_mean_feedback_order_D')
    def test_both_full_row_differences_strictly_above_C0(self):self.control('fast_and_slow_full_difference_strictly_above_C0')
    def test_all_three_regular_mean_rows_and_fresh_scalar_gain(self):
        self.control('all_regular_mean_rows_close_at_C0');self.control('fresh_current_postfast_scalars_stay_C');self.control('fresh_Haar_scalar_reset_C_plus_one')
    def test_compact_integration_counts_strict_and_minimal_for_rule(self):
        self.control('both_compact_counts_have_strict_margin');self.control('preceding_compact_integers_do_not_meet_this_strict_rule')
    def test_fast_envelope_keeps_raw_factorials_A1_and_input_orders(self):
        self.control('fast_envelope_reads_tangential_n4_and_radial_n5');self.control('fast_envelope_keeps_A1_and_raw_factorial_conversion')
    def test_new_4q_domain_and_old_current_mean_poles_kept(self):
        self.control('new_postfast_4q_domain_exceeds_signed_pressure_2q');self.control('wave_only_pole_does_not_forget_old_mean_domain')
    def test_original_kernel_recipes_cover_all_new_proof_domains(self):
        self.control('original_shift_proves_all_2q_4q_domains');self.control('moderate_new_kernel_constants_are_dominated')
    def test_geometry_and_width_cannot_be_changed_by_estimate_shift(self):
        self.control('same_original_width_hull_and_gap_in_operator_recipe');self.control('different_actual_width_or_hull_is_rejected')
    def test_dependency_graph_traces_new_covariance_and_target_requests(self):
        self.control('new_next_target_input_is_A0_source_plus_14');self.control('new_prefast_row_is_A0_source_plus_10');self.control('order_zero_A1_A0_and_next_target_counts')
    def test_whole_request_pole_not_low_output_pole(self):
        self.control('one_common_pole_uses_whole_huge_request_graph');self.control('huge_operator_constants_deferred_instead_of_allocated')
    def test_old_DAGs_are_not_collapsed_or_actual_jets_evaluated(self):self.control('older_DAGs_and_actual_profiles_not_evaluated')
    def test_new_source_recipes_have_current_stages_and_full_radial_source(self):
        self.control('new_rows_derived_from_old_final_recipe_and_all_six_covariance_entries');self.control('new_A1_source_is_post_signed_row_recipe');self.control('new_A0_source_is_FULL_final_gr_not_stale_gr')
    def test_actual_A1_whole_axial_pressure_and_raw_kernel_factors_kept(self):
        self.control('A1_actual_velocity_kept_in_fresh_fast_recipe');self.control('whole_new_axial_pressure_product_stays_current');self.control('all_high_jet_source_profiles_floor_and_raw_kernel_recipes_retained')
    def test_full_difference_gets_increment_jets_through_n_plus_three(self):self.control('fast_difference_reads_increment_jets_through_n_plus_three')
    def test_no_actual_input_certificate_from_deferred_norm_program(self):self.control('no_actual_input_or_full_cycle_certificate_from_finite_recipes')
    def test_invalid_exact_parameters_and_dependency_counts_rejected(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):strict_integrations(1,2,0)
        with self.assertRaises(ValueError):dependency_trace(0,0,1)
        with self.assertRaises(ValueError):operator_ledger(-1,1)
        with self.assertRaises(ValueError):mean_feedback(0,P(1),{'u':P(1)},graph=P(1),reciprocal=P(1),normal=P(1))
    def test_eight_frozen_dependencies_and_record_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
        self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_obligation_or_dependency_mutation_rejected(self):
        for k,v in self.record['flags'].items():
            bad=deepcopy(self.record);bad['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['dependency_evidence_sha256'][next(iter(bad['dependency_evidence_sha256']))]='wrong'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_frozen_record_never_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
