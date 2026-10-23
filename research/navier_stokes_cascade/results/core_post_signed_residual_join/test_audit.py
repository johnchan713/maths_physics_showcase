"""Current full signed source and physical request/domain regression checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import rates,physical_request,whole_domain,regular_source_envelope,P
from controls import review,fixture

class CurrentResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_full_homogeneous_pressure_normal_projection(self):
        self.control('full_moving_normal_pressure_cancels_normal_principal_residual');self.control('homogeneous_pressure_is_not_zero')
    def test_moving_normal_pressure_cannot_be_omitted(self):self.control('omitting_moving_normal_pressure_is_detected')
    def test_signed_scalar_must_be_auxiliary_time_independent(self):
        self.control('auxiliary_independent_current_target_scalar_preserves_homogeneous_scaling');self.control('auxiliary_dependent_scalar_leaves_a_principal_velocity_leak')
    def test_transverse_pulse_solve_is_not_certified_by_pressure_formula(self):self.control('transverse_equation_still_requires_original_projected_solve')
    def test_both_K_and_moving_normal_pressure_norm_terms(self):self.control('both_K_and_moving_normal_pressure_terms_retained')
    def test_exact_current_signed_plus_two_mean_nonzero_difference(self):self.control('exact_current_signed_then_two_mean_nonzero_difference')
    def test_full_zero_mode_stays_in_mean_system(self):self.control('zero_mode_not_discarded_from_actual_full_system')
    def test_both_old_signed_wave_crosses_retained(self):self.control('both_signed_wave_cross_orientations_required')
    def test_signed_square_not_cancelled_by_homogeneous_linear_equation(self):self.control('signed_square_is_not_homogeneous_linear_cancellation')
    def test_both_new_mean_wave_orientations_retained(self):self.control('both_mean_feedback_orientations_required')
    def test_third_new_and_fourth_inherited_residual_harmonics_retained(self):
        self.control('signed_addition_creates_third_harmonics');self.control('inherited_fourth_harmonic_still_required');self.control('signed_update_does_not_introduce_fifth_or_higher_modes')
    def test_complete_velocities_and_conjugate_reality(self):
        self.control('updated_complete_velocities_have_only_first_and_second_harmonics');self.control('conjugate_output_reality_preserved')
    def test_actual_A1_repair_velocity_not_erased_with_residual_tail(self):self.control('deleting_actual_repair_velocity_changes_nonzero_feedback')
    def test_positive_polynomial_and_independent_physical_source_bounds(self):
        self.control('full_positive_source_polynomial_substitution');self.control('independent_physical_complete_source_bound_D')
    def test_old_regular_source_added_not_cancelled(self):self.control('old_regular_source_is_added_not_cancelled')
    def test_curl_pressure_and_both_mean_omissions_detected(self):
        for k in ('R','normal_time','fastM','slowM'):self.control('omission_of_'+k+'_detected')
    def test_actual_cut_velocity_and_derived_full_signed_pressure(self):
        self.control('actual_cut_velocity_uses_cutoff_once');self.control('full_signed_pressure_derived_not_independent_callback')
    def test_source_remains_D_while_mean_rows_are_C0(self):
        self.control('all_current_nonzero_terms_have_D_or_better');self.control('new_mean_C0_does_not_promote_nonzero_source');self.control('quantitative_seed_relation')
    def test_original_root_and_curl_recurrence_really_executed(self):self.control('signed_wave_uses_original_root_and_complete_curl_recurrence')
    def test_physical_requests_larger_than_normalized_requests(self):
        self.control('physical_counts_exceed_normalized_mean_counts');self.control('current_physical_target_request_is_A0_source_plus_14');self.control('whole_physical_request_needs_larger_common_pole')
    def test_actual_current_2q_3q_4q_domains(self):
        self.control('physical_A0_and_A1_use_current_4q_and_2q_sources');self.control('entire_nonzero_source_uses_three_q_domain')
    def test_separate_oriented_primitive_constant_not_lost(self):
        self.control('separate_oriented_primitive_constant_kept');self.control('moderate_sum_shift_dominates_both_distinct_operator_constants')
    def test_archived_shifts_retained_without_pole_monotonicity(self):
        self.control('archived_lower_shifts_explicitly_retained_without_pole_monotonicity');self.control('omitting_archived_shift_need_not_dominate_old_input_polynomial')
    def test_raw_kernel_and_current_tail_allocation(self):
        self.control('raw_profiles_source_factorials_and_graph_floor_retained');self.control('current_tail_signs_and_pressure_compensation_kept');self.control('new_signed_tail_has_no_inhomogeneous_source_term');self.control('base_stress_not_counted_twice')
    def test_giant_operators_and_actual_inputs_not_materialized(self):
        self.control('operators_and_huge_actual_fields_not_materialized');self.control('actual_input_certification_is_refused')
    def test_physical_margin_trace_and_negative_regular_exponent(self):
        self.control('physical_margins_strict');self.control('source_jet_trace_not_low_output_pole');self.control('regular_source_still_does_not_certify_absolute_decay');self.control('regular_grammar_only_not_actual_complete_theorem')
    def test_invalid_parameters_or_missing_full_data_rejected(self):
        for k in (0,Q(1,4),0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):physical_request(0,0,Q(1,100))
        with self.assertRaises(ValueError):whole_domain({'maxima':{}},-1,1)
        with self.assertRaises(ValueError):regular_source_envelope(0,P(1),{})
    def test_frozen_dependencies_and_evidence_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
        self.assertTrue(all(self.record['gates'].values()));self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_remaining_or_provenance_mutations_rejected(self):
        for k,v in self.record['flags'].items():
            bad=deepcopy(self.record);bad['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['dependency_evidence_sha256'][next(iter(bad['dependency_evidence_sha256']))]='wrong'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_frozen_record_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')

if __name__=='__main__':unittest.main()
