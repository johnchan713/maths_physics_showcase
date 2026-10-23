"""Joint divisor and current nonlinear centrifugal/pressure input tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import divisor_data,kernel_recipe,sharp_boundary,radial_counts,split_envelope,radial_profile_recipe
from controls import review

class MixedCompactTimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_exact_original_joint_divisor_and_quadratic_conjugacy(self):
        self.control('exact_same_mode_joint_integer_form');self.control('temporal_divisor_is_b_times_radial_quadratic_conjugate')
    def test_all_224_integer_modes_have_nonzero_joint_norm(self):
        self.control('all_224_integer_grid_joint_norms_nonzero');self.control('all_224_joint_products_have_nonzero_integer_floor')
    def test_joint_divisor_constant_and_all_Pell_modes(self):
        self.control('joint_constant_b_inverse_is_less_than_three');self.control('actual_joint_divisor_product_is_exact_plus_minus_b')
    def test_temporal_inverse_gains_one_frequency_on_radial_resonance(self):self.control('temporal_inverse_gains_one_over_P_on_radial_resonance')
    def test_different_input_output_divisor_mutation_rejected(self):self.control('different_input_output_divisor_mutation_detected')
    def test_raw_derivative_counts_for_every_output_split(self):self.control('all_derivative_splits_read_exact_m_plus_two_p_plus_two')
    def test_joint_multiplier_has_summable_fourier_tail(self):
        self.control('joint_divisor_leaves_summable_r_minus_three');self.control('four_r_modes_and_inverse_square_tail_bound')
    def test_two_pi_direction_and_five_derivative_sequential_comparison(self):
        self.control('two_pi_power_is_negative_and_not_dropped_in_wrong_direction');self.control('five_derivative_saving_is_only_vs_sequential_upper_rules')
    def test_original_raw_profiles_and_source_factorial_not_erased(self):self.control('raw_original_profiles_and_source_factorial_retained')
    def test_mixed_finite_jet_lower_boundary_and_unresolved_equality(self):
        self.control('one_below_mixed_lower_boundary_is_obstructed');self.control('equality_left_open_and_upper_gap_three')
    def test_independent_rational_fast_compliant_ratio_grows(self):self.control('independent_rational_mixed_ratio_grows_like_P')
    def test_exact_polynomial_radial_integration_by_parts(self):
        self.control('complete_radial_integration_with_no_endpoint_terms');self.control('fixture_has_every_raw_source_jet_at_endpoints');self.control('nonzero_real_compact_integral_is_retained')
    def test_real_temporal_inverse_normalization_and_sign(self):
        self.control('real_Nt_inverse_has_correct_two_pi_and_sign');self.control('dropping_temporal_divisor_changes_exact_real_inverse')
    def test_actual_radial_profile_reduces_to_reciprocal_U(self):self.control('actual_f0_over_R_is_exact_reciprocal_U_profile')
    def test_all_nine_actual_profile_derivatives_and_through_jet_bounds(self):
        self.control('exact_raw_derivative_of_actual_new_kernel_factor');self.control('entire_through_jet_upper_recipe')
    def test_giant_actual_profile_recipe_keeps_factorial_deferred(self):self.control('huge_actual_profile_recipe_keeps_factorial_symbolic_and_positive_hull')
    def test_complete_nonlinear_Q_split_c_normalization_and_sign(self):self.control('complete_nonlinear_centrifugal_Q_split_with_correct_c_and_sign')
    def test_Fsum_derived_from_full_current_swirl(self):self.control('Fsum_is_full_current_time_derivative_not_new_fine_row')
    def test_slow_fast_fields_and_quadratic_factor_two_are_essential(self):
        self.control('full_slow_and_all_fast_parts_remain_inside_Gtheta');self.control('missing_factor_two_in_Nt_v_square_is_detected')
    def test_Haar_pressure_keeps_fast_variance_and_old_new_correlations(self):
        self.control('Haar_centrifugal_pressure_keeps_all_fast_variance_and_old_new_cross');self.control('Haar_v_square_is_not_square_of_Haar_v')
    def test_Q_split_does_not_replace_full_source_or_erase_pressure(self):
        self.control('using_Q_split_for_full_C_would_erase_pressure_moment');self.control('radial_Q_split_preserves_complete_K0_compensated_pressure')
    def test_auxiliary_multiplier_inverse_and_product_commutators(self):
        self.control('auxiliary_dependent_multiplier_inverse_commutator_nonzero');self.control('auxiliary_multiplier_Nt_product_rule_keeps_derivative_term')
    def test_fast_compliant_comparison_moments_and_fundamental_source(self):
        self.control('comparison_obeys_actual_fast_formula_and_both_velocity_moments');self.control('fast_compliant_centrifugal_square_does_not_cancel_fundamental');self.control('fast_compliant_pressure_variance_remains_nonzero')
    def test_full_positive_polynomial_and_independent_epsilon_budgets(self):
        self.control('positive_polynomial_and_exact_rational_split_envelopes_agree');self.control('independent_full_epsilon_budget_Gtheta_B_rest_one')
    def test_current_Fsum_and_rest_local_derivative_costs(self):self.control('generic_Fsum_route_keeps_one_v_and_two_beta_derivatives')
    def test_normalized_and_physical_strict_split_request_counts(self):
        self.control('normalized_strict_split_counts_and_local_input_max');self.control('same_physical_loss_and_strict_split_recipe')
    def test_different_source_functions_and_older_graph_retained(self):
        self.control('normalized_G200002_rest200001_are_different_sources_not_gr');self.control('physical_G60200006_rest60200005_and_original_graph_retained')
    def test_mixed_epsilon_route_still_has_huge_finite_row_threshold(self):self.control('mixed_epsilon_boundary_still_requires_huge_row_order')
    def test_complete_pressure_and_full_regular_and_A1_orders_unchanged(self):
        self.control('complete_pressure_variance_and_all_original_shifts_remain');self.control('whole_regular_E_A1_gain_and_actual_input_gate_unchanged')
    def test_large_orders_only_integer_recipes_not_allocated_factorials(self):
        r=kernel_recipe(2,30000000);self.assertEqual(r['source_factorial'],60000004)
        self.assertFalse(r['actual_input_or_factorial_evaluated'])
        self.assertFalse(radial_profile_recipe(60000004,1,1)['numerical_factorial_or_actual_source_materialized'])
    def test_invalid_zero_modes_counts_hulls_and_sources_rejected(self):
        for k in ((0,0),(1.0,2),(1,)):
            with self.assertRaises(ValueError):divisor_data(k)
        for fn,args in ((kernel_recipe,(0,0)),(sharp_boundary,(0,0,-1)),(split_envelope,(0,{})),(radial_profile_recipe,(0,0,1)),(radial_profile_recipe,(0,1,0))):
            with self.assertRaises(ValueError):fn(*args)
        with self.assertRaises(ValueError):radial_counts(0,physical_power=-1)
    def test_review_cache_does_not_alias_saved_evidence(self):
        r=review();r['checks'].clear();self.assertTrue(review()['checks'])
    def test_five_historical_records_and_frozen_evidence_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),5)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertTrue(all(self.record['gates'].values()))
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_provenance_pressure_and_remaining_mutations_rejected(self):
        for field in ('flags','protocol','dependency_evidence_sha256','remaining_obligations'):
            r=deepcopy(self.record)
            if field=='flags':r[field]['actual_all_scale_current_source_high_jet_or_support_inputs_validated']=True
            elif field=='protocol':r[field]['Haar_pressure']='drop fast variance'
            elif field=='dependency_evidence_sha256':r[field][next(iter(r[field]))]='wrong'
            else:r[field].pop()
            with self.assertRaises(ValueError):audit.validate(r)
    def test_no_existing_frozen_output_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')
if __name__=='__main__':unittest.main()
