"""Full current angular/radial pressure identities and scoped norm checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import rates,request_counts,pressure_envelope,temporal_recipe
from controls import review

class AngularPressureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_full_angular_equation_and_exact_swirl_normal_form(self):self.control('full_current_angular_row_yields_exact_swirl_normal_form')
    def test_current_angular_equation_has_no_pressure_or_A0(self):self.control('no_current_pressure_or_A0_enters_angular_equation')
    def test_full_vector_viscosity_connection_identity(self):self.control('full_vector_angular_viscosity_has_correct_cylindrical_connection')
    def test_wrong_viscosity_sign_is_detected(self):self.control('wrong_viscosity_sign_in_F_destroys_exact_swirl_identity')
    def test_missing_slow_time_is_detected(self):self.control('missing_slow_time_in_F_destroys_exact_swirl_identity')
    def test_selected_stress_removed_only_after_Haar_projection(self):self.control('auxiliary_independent_selected_stress_drops_only_after_Q')
    def test_variable_V_multiplier_and_complete_graph_product_rule(self):
        self.control('varying_V_requires_the_full_radial_multiplier_derivative');self.control('weighted_radial_product_identity_includes_full_graph_drift')
    def test_omitted_V_derivative_is_detected(self):self.control('omitting_V_derivative_is_detected')
    def test_old_radial_source_with_both_viscosities(self):self.control('old_radial_source_reconstructed_with_both_viscosities')
    def test_new_radial_source_and_local_pressure_reconstructed_independently(self):
        self.control('full_new_radial_source_and_local_pressure_are_exact');self.control('pressure_source_rewrite_is_exact_before_projection')
    def test_H_K_zero_Haar_does_not_mean_zero_fields(self):self.control('H_K_have_zero_Haar_but_are_not_zero')
    def test_Haar_pressure_and_full_fast_variance_remain(self):
        self.control('all_Haar_pressure_density_and_moving_P_are_unchanged');self.control('full_fast_variance_is_not_Haar_v_squared');self.control('mean_pressure_keeps_full_variance')
    def test_omitted_new_local_pressure_changes_axial_force(self):self.control('omitting_new_local_pressure_changes_axial_pressure_force')
    def test_axial_velocity_flux_and_whole_endpoint_traces_remain(self):
        self.control('axial_velocity_and_pressure_flux_stay_in_radial_remainder');self.control('entire_K0_compact_trace_vanishes');self.control('new_local_potential_has_zero_both_endpoint_traces')
    def test_independent_exact_characteristic_moments(self):
        self.control('same_full_characteristic_moment_kills_exact_radial_derivative');self.control('weighted_angular_flux_reduces_without_graph_derivative')
    def test_pressure_primitive_is_the_nonzero_local_potential(self):
        self.control('pressure_primitive_of_derivative_is_the_potential_not_zero');self.control('dropping_local_pressure_would_leave_nonzero_radial_force')
    def test_actual_new_kernel_factors_keep_V_U_and_reciprocal_profiles(self):
        self.control('actual_new_J_kernel_has_the_extra_background_radial_jet');self.control('actual_F_kernel_is_the_retained_reciprocal_U_times_two_V')
    def test_conditional_flux_and_pressure_exponents_for_original_kappa_range(self):
        self.control('J_is_one_F_is_E_and_both_strictly_controlled');self.control('mean_pressure_improves_conditionally_but_regular_E_A1_do_not')
    def test_all_four_distinct_compact_counts_have_strict_margins(self):
        self.control('all_four_same_target_strict_margins');self.control('normalized_F_requires_two_integrations_quadratic_one');self.control('physical_request_keeps_every_distinct_source_and_profile')
    def test_normalized_and_physical_exact_raw_indices(self):
        self.control('normalized_raw_orders_are_J200000_F6_square5_rest200001');self.control('physical_raw_orders_are_J60200004_F60000010_square60000009_rest60200005')
    def test_previous_local_field_max_and_whole_old_graph_retained(self):self.control('retained_older_local_field_max_is_not_declared_reduced')
    def test_positive_polynomial_pressure_program_with_exact_rational_substitution(self):self.control('positive_polynomial_pressure_program_agrees_with_exact_rational')
    def test_independent_epsilon_budget_including_two_graph_viscous_losses(self):self.control('independent_epsilon_J1_FE_and_full_pressure1_budget')
    def test_four_temporal_plus_two_viscous_jets_kept(self):self.control('four_temporal_jets_plus_two_viscous_jets_retained')
    def test_existing_angular_input_not_a_circular_radial_assumption(self):self.control('explicit_non_circular_input_is_the_existing_angular_component_only')
    def test_full_local_pressure_and_actual_frontier_not_erased(self):
        self.control('complete_Haar_variance_local_K1_and_old_pressure_memory_kept');self.control('actual_input_and_complete_residual_frontier_not_promoted')
    def test_invalid_orders_inputs_and_kappa_rejected(self):
        for k in (Q(0),Q(1,4),Q(-1)):
            with self.assertRaises(ValueError):rates(k)
        with self.assertRaises(ValueError):request_counts(-1)
        with self.assertRaises(ValueError):request_counts(0,physical_power=-1)
        with self.assertRaises(ValueError):pressure_envelope(0,{})
    def test_large_requests_remain_integer_recipes(self):
        r=request_counts(50,physical_power=100,h=Q(1,1000));self.assertGreater(r['J']['raw_order'],1000000)
        self.assertFalse(r['actual_input_or_giant_constant_evaluated'])
    def test_temporal_constant_keeps_projection_four_jets_and_raw_factorial(self):
        r=temporal_recipe(3)
        self.assertEqual((r['scalar'],r['two_power'],r['six_power'],r['source_factorial']),(48,7,3,7))
        self.assertEqual(temporal_recipe(60200004)['source_factorial'],60200008)
        self.assertFalse(r['materialized'])
    def test_review_cache_does_not_mutate_saved_record(self):
        r=review();r['checks'].clear();self.assertTrue(review()['checks'])
    def test_five_frozen_records_and_new_evidence_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),5)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertTrue(all(self.record['gates'].values()))
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_pressure_input_scope_and_provenance_mutations_rejected(self):
        for field in ('flags','protocol','dependency_evidence_sha256','remaining_obligations'):
            r=deepcopy(self.record)
            if field=='flags':r[field]['actual_all_band_high_jet_Haar_profile_or_support_norms_validated']=True
            elif field=='protocol':r[field]['pressure']='delete K1 and fast variance'
            elif field=='dependency_evidence_sha256':r[field][next(iter(r[field]))]='wrong'
            else:r[field].pop()
            with self.assertRaises(ValueError):audit.validate(r)
    def test_existing_frozen_output_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')
if __name__=='__main__':unittest.main()
