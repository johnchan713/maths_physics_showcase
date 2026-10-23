#!/usr/bin/env python3
"""Original-root signed covariance, current pressure, domains and scope tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import rates,pressure_domains,edge,RootSignedEnvelope,P
from controls import review


class GenericSignedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        values=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(values,fragment);self.assertTrue(all(values),fragment)

    def test_localized_original_root_linear_cross_matches_current_target(self):self.control('original_root_linear_selected_assembly_exact')
    def test_original_cutoffs_angular_half_and_Haar_integral_kept(self):
        self.control('original_complete_column_contains_both_cutoffs');self.control('wrong_angular_half_or_missing_factor_two_detected')
    def test_partition_identity_including_nonzero_derivatives(self):
        self.control('square_partition_counted_once');self.control('localized_cross_identity_including_partition_derivatives');self.control('individual_partition_derivatives_nonzero')
    def test_updated_root_denominator_and_signed_square_shortcuts_detected(self):
        self.control('updated_root_divisor_breaks_original_cross_identity');self.control('signed_square_not_part_of_cancelled_linear_term')
    def test_other_four_principal_entries_retained(self):self.control('all_four_unselected_principal_entries_nonzero')
    def test_full_generic_tensor_identity_and_three_error_groups(self):
        self.control('complete_generic_tensor_identity_all_six');self.control('all_three_signed_error_groups_are_nonzero')
    def test_complex_quadratures_and_signed_square_are_not_dropped(self):
        self.control('complex_old_wave_has_both_quadratures');self.control('deleting_signed_square_is_detected')
    def test_multiharmonic_zero_mode_is_actual_convolution(self):
        self.control('multiharmonic_complete_generic_identity_all_six');self.control('second_harmonic_old_wave_is_nonzero');self.control('no_m2_m1_cross_reclassified_as_angular_mean')
    def test_both_current_targets_cancel_full_principal_cross(self):self.control('current_signed_target_cancels_complete_principal_cross')
    def test_new_Haar_rows_reconstruct_with_current_moment_channel(self):self.control('whole_new_Haar_row_reconstructs_from_target_and_current_moment')
    def test_moving_projection_commutator_exact_and_nonzero(self):
        self.control('moving_projection_commutator_exact');self.control('moving_projection_commutator_nonzero')
    def test_new_pressure_same_unit_mass_compensator_and_full_moment(self):
        self.control('new_pressure_reconstructed_with_same_unit_mass_rho');self.control('whole_axial_pressure_moment_retains_moving_product');self.control('freezing_c_rho_derivative_detected')
    def test_target_cancellation_does_not_cancel_pressure_scalar_defects(self):
        self.control('projected_target_cancellation_does_not_cancel_scalar_pressure_defect');self.control('other_four_covariance_entries_cannot_be_discarded')
    def test_total_wave_norm_not_used_as_original_error_bound(self):self.control('original_principal_error_not_inferred_from_total_norm')
    def test_positive_tensor_polynomials_match_independent_substitution(self):self.control('covariance_positive_polynomials_equal_direct_substitution')
    def test_raw_physical_covariance_and_pressure_budgets(self):
        self.control('complete_covariance_C_and_error_C0_physical_bounds');self.control('Haar_pressure_unreduced_bound_keeps_epsilon');self.control('post_signed_target_C0_not_pressure_order_promotion')
    def test_all_three_generic_errors_tie_at_C0(self):self.control('all_three_generic_errors_tie_at_C0')
    def test_pressure_scalar_radial_mean_orders_not_promoted(self):
        self.control('axial_pressure_response_stronger_than_target');self.control('pressure_scalar_and_radial_mean_not_promoted')
    def test_original_root_recurrence_and_full_curl_extra_jet(self):
        self.control('same_original_root_signed_principal_recurrence');self.control('complete_signed_curl_uses_same_operator_and_one_extra_target_jet')
    def test_both_target_components_and_whole_axial_moment_reach_original_root(self):
        self.control('both_distinct_current_targets_and_whole_axial_moment_feed_original_root');self.control('both_target_edge_allowances_are_covered')
        with self.assertRaises(ValueError):RootSignedEnvelope(edge.root.fixture(),edge.fixture(),lambda n:P(1))
    def test_old_operator_domains_do_not_cover_new_pressure_call(self):self.control('new_pressure_domain_not_covered_by_old_6p_12p_proof')
    def test_pressure_and_target_source_orders_retained(self):self.control('source_orders_keep_pressure_and_curl_extra_jets')
    def test_new_pressure_and_moment_proof_domains_from_original_shift(self):
        self.control('replacement_exact_shift_covers_new_pressure_and_moment_domains');self.control('new_fixed_shift_dominates_actual_raw_kernel_constants')
    def test_oriented_target_primitive_has_separate_correct_domain(self):self.control('oriented_new_target_primitive_has_its_own_correct_source_pole')
    def test_shift_changes_no_actual_field_width_hull_or_band(self):self.control('estimate_shift_preserves_actual_width_hull_and_band')
    def test_example_pressure_pole60_is_sufficient_not_actual_input(self):self.control('zero_order_example_q30_new_pressure60')
    def test_current_final_rows_whole_moments_and_prefix_identity(self):
        self.control('current_target_uses_final_row_recipe_not_new_fine_callback');self.control('current_moment_uses_whole_parent_pressure_scalar_recipe');self.control('different_current_prefix_is_rejected')
    def test_actual_input_certificate_not_inferred_from_recipes(self):self.control('deferred_recipes_do_not_certify_actual_inputs')
    def test_invalid_kappa_order_and_pole_rejected(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):pressure_domains(-1,1,edge.fixture())
        with self.assertRaises(ValueError):pressure_domains(0,-1,edge.fixture())
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
    def test_frozen_evidence_never_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
