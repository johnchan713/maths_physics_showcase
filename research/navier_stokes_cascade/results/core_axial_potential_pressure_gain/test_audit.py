"""Axial potential, full-pressure viscosity, moving cutoff and scope checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import (stage_rates,A1_request,direct_A1_recipe,complement_envelope,
    potential_envelope,tail_envelope,COMPLEMENT_KEYS)
from controls import review

class AxialPotentialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_complete_axial_row_and_compact_potential(self):self.control('complete_axial_row_has_exact_compact_potential_factorization')
    def test_actual_potential_preserves_divergence_and_commuting_derivatives(self):
        self.control('actual_potential_preserves_exact_divergence');self.control('fixed_label_time_and_slow_derivatives_commute_with_D1')
    def test_old_coarse_source_temporal_inverse_gives_exact_Q_potential(self):self.control('coarse_source_gives_exact_nonzero_Haar_potential_after_Nt_inverse')
    def test_dropping_time_or_radial_viscosity_is_detected(self):
        self.control('potential_time_term_is_not_dropped');self.control('potential_radial_viscosity_is_not_dropped')
    def test_full_and_reduced_pressure_forms_have_correct_distinct_viscous_factors(self):self.control('full_pressure_and_reduced_pressure_forms_require_different_viscous_factor')
    def test_missing_second_viscous_pressure_contribution_is_detected(self):
        self.control('removing_local_viscous_pressure_without_doubling_is_detected');self.control('local_viscous_pressure_still_has_nonzero_axial_force')
    def test_variable_G_projection_and_compact_product_identity(self):
        self.control('full_A1_source_reduction_retains_variable_G_radial_derivative');self.control('compact_product_identity_for_G_gamma_is_exact')
    def test_omitted_G_R_is_detected(self):self.control('wrongly_omitting_G_R_Psi_changes_the_axial_source')
    def test_auxiliary_independent_stress_removed_only_by_Q(self):self.control('Q_allows_auxiliary_independent_stress_only_after_projection')
    def test_all_pressure_covariance_and_axial_velocity_square_retained(self):
        self.control('all_current_pressure_and_gamma_square_stay_in_Phi');self.control('dropping_full_pressure_from_Phi_changes_pointwise_axial_force');self.control('dropping_axial_velocity_square_changes_projected_flux')
    def test_compensator_removal_is_scoped_to_Q_channel(self):self.control('Q_removes_only_the_Haar_compensator_not_its_full_axial_force')
    def test_finer_Haar_potential_not_inferred_zero(self):self.control('Haar_potential_and_Haar_axial_velocity_are_nonzero')
    def test_whole_compact_potential_traces(self):
        self.control('whole_potential_Lambda_trace_vanishes_at_both_ends');self.control('whole_G_Psi_potential_trace_vanishes_at_both_ends')
    def test_local_viscous_pressure_drops_only_inside_exact_projection(self):self.control('actual_Viscous_pressure_may_only_drop_after_A1_Dz')
    def test_independent_T1_compact_inverse_and_weighted_moments(self):
        self.control('T1_D1_compact_potential_identity_is_exact');self.control('whole_weighted_axial_characteristic_moment_is_zero')
    def test_independent_characteristic_G_gamma_projection(self):
        self.control('G_gamma_characteristic_moment_keeps_G_R_Psi');self.control('G_R_term_cannot_be_omitted_from_A1')
    def test_Haar_potential_uses_primitive_not_temporal_inverse(self):self.control('Haar_potential_uses_weighted_primitive_not_zero_inverse')
    def test_moving_T1_commutator_is_kept(self):
        self.control('nonzero_characteristic_source_moment_exercises_moving_cutoff');self.control('T1_Dz_commutator_contains_full_cutoff_derivative');self.control('commuting_T1_through_Dz_for_general_X_is_invalid')
    def test_moving_A1_commutator_is_kept(self):
        self.control('A1_Dz_commutator_contains_full_cutoff_derivative');self.control('commuting_A1_through_Dz_for_general_Phi_is_invalid')
    def test_Haar_kernel_witness_obeys_moment_and_has_fine_row_response(self):
        self.control('Haar_kernel_comparison_preserves_weighted_axial_velocity_moment');self.control('Haar_kernel_comparison_has_nonzero_potential_and_fine_viscous_row_response');self.control('Haar_kernel_response_is_finer_than_both_existing_rows')
    def test_coarse_Haar_rescaling_still_has_exponential_gap(self):self.control('coarse_Haar_potential_rescaling_can_grow_exponentially')
    def test_stage_source_orders_potential_beta_and_regular_frontier(self):
        self.control('coarse_RHS_source_is_one_with_no_fine_pressure_assumption');self.control('potential_one_beta_two_gamma_and_regular_unchanged');self.control('angular_pressure_extension_uses_only_row_at_least_E')
    def test_direct_A1_margin_and_distinct_physical_flux_requests(self):
        self.control('direct_A1_order_two_has_double_kappa_target_margin');self.control('same_physical_target_and_two_lower_distinct_flux_jets');self.control('physical_Phi_raw60000006_and_previous_H1_raw60000008_are_distinct_sources')
    def test_original_larger_A0_graph_and_both_stage_domains_remain(self):
        self.control('new_local_coarse_pressure_field60000014_is_below_old_final_A0_field60200007');self.control('postwave6q_and_final12q_not_interchanged')
    def test_positive_polynomial_composition_matches_rational_data(self):self.control('composed_positive_polynomials_match_exact_rational_inputs')
    def test_independent_full_epsilon_norm_budget(self):self.control('independent_full_epsilon_coarse_to_fine_potential_pressure_and_Phi_budget')
    def test_Q_route_has_no_additional_finer_Haar_or_moment_callback(self):
        self.control('Q_tail_ignores_finer_Haar_norms_and_Haar_pressure_moment');self.control('Q_only_program_requires_no_Haar_or_moment_callback')
    def test_tail_and_full_beta_have_different_input_obligations(self):self.control('tail_uses_one_extra_Phi_jet_and_beta_separately_needs_Haar_history')
    def test_normal_form_uses_old_coarse_inputs_without_circularity(self):self.control('potential_source_uses_OLD_coarse_pressure_and_Psi_to_avoid_circularity')
    def test_postwave_current_pressure_and_whole_hull_identity_kept(self):self.control('same_POSTWAVE_pressure_and_whole_Psi_hull_are_explicit')
    def test_complete_pressure_scope_and_frontier_are_retained(self):
        self.control('complete_pressure_mean_variance_and_Pi_memory_remain');self.control('no_actual_input_full_regular_or_blowup_promotion')
    def test_final_state_cannot_replace_postwave_tail(self):
        with self.assertRaises(ValueError):tail_envelope(0,{},stage='final')
        with self.assertRaises(ValueError):stage_rates('stale')
    def test_invalid_orders_kappa_and_missing_sources_rejected(self):
        for k in (0,Q(1,4),Q(-1)):
            with self.assertRaises(ValueError):stage_rates(kappa=k)
        for fn in (A1_request,direct_A1_recipe):
            with self.assertRaises(ValueError):fn(-1)
        with self.assertRaises(ValueError):A1_request(0,physical_power=-1)
        with self.assertRaises(ValueError):potential_envelope(0,{})
        with self.assertRaises(ValueError):complement_envelope(0,{})
    def test_large_requests_and_raw_direct_conversion_remain_deferred(self):
        r=A1_request(10,physical_power=100,h=Q(1,1000));self.assertGreater(r['raw_Phi_order'],1000000)
        self.assertFalse(r['giant_factorials_or_actual_source_norms_evaluated'])
        a=direct_A1_recipe(3);self.assertEqual((a['scalar'],a['six_power'],a['source_factorial']),(2,3,4));self.assertFalse(a['materialized'])
    def test_review_cache_does_not_alias_frozen_evidence(self):
        r=review();r['checks'].clear();self.assertTrue(review()['checks'])
    def test_six_historical_records_and_new_evidence_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertTrue(all(self.record['gates'].values()))
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_stage_pressure_scope_and_provenance_mutations_rejected(self):
        for field in ('flags','protocol','dependency_evidence_sha256','remaining_obligations'):
            r=deepcopy(self.record)
            if field=='flags':r[field]['actual_all_scale_high_jet_Haar_profile_or_support_inputs_validated']=True
            elif field=='protocol':r[field]['pressure']='use final pressure in postwave source and delete viscosity'
            elif field=='dependency_evidence_sha256':r[field][next(iter(r[field]))]='wrong'
            else:r[field].pop()
            with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')
if __name__=='__main__':unittest.main()
