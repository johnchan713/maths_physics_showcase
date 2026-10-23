"""Current fresh mean/radial and complete physical request regression checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import rates,request_trace,compact_physical_orders,regular_physical,SourcedMeanTailJoin,current_old_Haar,P
from controls import review

class NextMeanTailTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_full_tangential_identity_and_Haar_differences(self):
        self.control('full_tangential_identity_keeps_entire_fast_and_slow_differences');self.control('extra_Qcirc_on_differences_would_erase_Haar_terms')
    def test_A1_positive_residual_sign(self):self.control('exact_A1_residual_sign_is_positive')
    def test_compensated_A0_zero_Haar_and_nonzero_pointwise_remainder(self):
        self.control('compensated_A0_zero_Haar_piece_cancels');self.control('nonzero_radial_A0_does_not_follow_from_Haar_cancellation')
    def test_current_pressure_and_new_A0_not_stale_prewave(self):
        self.control('current_radial_compensator_uses_new_pressure');self.control('stale_prewave_A0_residual_is_detected');self.control('dropping_compensation_adds_spurious_Haar_A0_pressure')
    def test_whole_axial_pressure_and_full_old_fast_correlations(self):
        self.control('whole_axial_pressure_derivative_keeps_moving_compensator');self.control('FULL_old_fast_correlations_survive_Haar')
    def test_actual_A1_velocity_not_erased_with_tail(self):self.control('mean_feedback_retains_actual_A1_velocity')
    def test_independent_physical_full_tangential_C0_budget(self):self.control('full_current_tangential_C0')
    def test_increment_pressure_not_promoted_with_row_gain(self):self.control('increment_pressure_not_promoted_to_full_row_order')
    def test_all_unreduced_postslow_scalar_budgets_C0_plus_one(self):self.control('all_complete_postslow_scalars_C0_plus_one')
    def test_H4_feedback_has_derived_actual_fast_E_gain(self):self.control('H4_mean_feedback_has_actual_E_gain')
    def test_current_final_full_gr_includes_old_new_means_and_covariance(self):self.control('full_final_gr_B_with_all_means_and_covariance')
    def test_complete_postfast_three_channel_polynomial(self):self.control('postfast_three_channel_positive_polynomial_substitution')
    def test_current_Haar_sum_keeps_all_historical_slow_terms(self):self.control('preserved_Haar_is_all_three_historical_slow_summands')
    def test_postwave_Qcirc_rows_are_derived_with_full_pressure(self):
        self.control('postwave_Qcirc_rows_derived_from_old_C0_and_complete_covariance');self.control('postwave_axial_row_keeps_full_reconstructed_pressure')
    def test_actual_fast_input_order_and_A1_factor(self):self.control('actual_fast_recipe_reads_derived_n4_n5_and_keeps_A1')
    def test_fresh_postfast_scalar_slow_and_postslow_orders(self):
        self.control('derived_postfast_scalar_and_slow_orders_C0');self.control('postslow_full_scalars_and_axial_pressure_moment')
    def test_full_mean_difference_margins_and_all_three_final_rows(self):
        self.control('full_fast_and_slow_remainders_strictly_above_C0');self.control('all_three_regular_mean_rows_C0')
    def test_A1_two_integration_strict_rule_and_A0_source_bound(self):
        self.control('A1_two_integrations_but_not_one_for_strict_rule');self.control('A0_strict_count_and_B_source_derived');self.control('normalized_raw_A1_7_and_A0_200003_requests')
    def test_whole_regular_source_remains_E_not_mean_C0(self):self.control('whole_regular_E_not_mean_C0')
    def test_normalized_entire_new_source_normal_and_older_graph(self):
        self.control('normalized_trace_reads_current_old_source_A0_J_plus14');self.control('normalized_master_q_traces_older_graph_not_low_output')
    def test_physical_request_counts_use_current_E_A1_source(self):self.control('physical_current_A0_and_A1_requests_use_current_E_source')
    def test_physical_old_source_normal_and_earlier_signed_target_requests(self):
        self.control('physical_trace_reaches_current_source_and_normal_jets');self.control('physical_q_covers_older_target60400038')
    def test_original_6q_12q_and_full_nonzero9q_domains(self):
        self.control('original_6q_12q_domains_cover_actual_fast_and_final_gr');self.control('full_nonzero_H4_feedback_source_pole9q')
    def test_archived_shifts_and_oriented_Ie_remain(self):self.control('all_archived_shifts_and_separate_Ie6q_retained')
    def test_huge_inputs_not_materialized_or_geometry_changed(self):self.control('huge_source_graphs_not_expanded_evaluated_or_geometry_changed')
    def test_raw_current_profile_factorial_and_source_poles(self):self.control('current_full_source_raw_profiles_factorials_and_poles_retained')
    def test_current_scalar_history_full_pressure_and_radial_compensation(self):
        self.control('current_scalar_recipe_uses_entire_Haar_slow_history');self.control('whole_pressure_derivative_and_new_radial_compensation_kept')
    def test_exact_current_tail_and_assigned_base_allocation(self):
        self.control('exact_current_compact_tail_and_assigned_base_ledger');self.control('mean_updates_do_not_create_or_forward_Gaussian_tails')
    def test_physical_margins_and_negative_regular_exponent(self):
        self.control('strict_physical_margin_and_count');self.control('regular_E_physical_exponent_negative_not_residual_lower_bound')
    def test_conditional_gain_not_actual_all_scale_theorem(self):
        self.control('conditional_full_regular_gain_E_but_not_actual_complete_theorem');self.control('actual_all_scale_input_validity_not_certified_by_recipes')
    def test_invalid_parameters_and_current_prefix_rejected(self):
        with self.assertRaises(TypeError):SourcedMeanTailJoin(object())
        with self.assertRaises(ValueError):request_trace(0,0,1)
        with self.assertRaises(ValueError):compact_physical_orders(2,0,0,Q(1,200))
        with self.assertRaises(ValueError):regular_physical(0,Q(1,200),reserve=0)
        with self.assertRaises(ValueError):current_old_Haar({'x':P(1)}, {}, {})
        for k in (0,Q(1,4),0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
    def test_eight_frozen_dependencies_and_record_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
        self.assertTrue(all(self.record['gates'].values()));self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_remaining_and_history_mutations_rejected(self):
        for k,v in self.record['flags'].items():
            bad=deepcopy(self.record);bad['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['dependency_evidence_sha256'][next(iter(bad['dependency_evidence_sha256']))]='wrong'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_frozen_record_never_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')

if __name__=='__main__':unittest.main()
