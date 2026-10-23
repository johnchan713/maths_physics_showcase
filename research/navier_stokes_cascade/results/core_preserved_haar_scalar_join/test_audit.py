#!/usr/bin/env python3
"""Complete pressure/moment, strict rate, finite recipe and scope checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import rates,haar_slow_reset,P
from controls import review,fixture


class PreservedHaarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        values=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(values,fragment);self.assertTrue(all(values),fragment)

    def test_auxiliary_independent_slow_products_use_exact_Haar_factor(self):self.control('slow_products_use_old_Haar_exactly')
    def test_old_Haar_is_kept_and_old_products_not_factorized(self):
        self.control('old_Haar_is_nonzero');self.control('full_old_correlations_are_not_products_of_Haars')
    def test_fast_zero_Haar_does_not_remove_old_new_correlations(self):
        self.control('fast_all_components_zero_Haar');self.control('fast_old_new_correlation_survives');self.control('using_Haar_only_for_fast_would_erase_actual_cross')
    def test_complete_plain_pressure_and_weighted_pressure_moments(self):
        self.control('exact_complete_pressure_remainder');self.control('exact_weighted_pressure_remainder')
    def test_both_complete_tangential_scalar_remainders_after_linear_rows(self):
        self.control('exact_angular_scalar_after_linear_row');self.control('exact_axial_scalar_after_linear_row')
    def test_pressure_field_is_not_discarded_and_compact(self):self.control('pressure_field_not_discarded')
    def test_whole_moving_pressure_derivative_not_frozen(self):
        self.control('whole_moving_axial_pressure_identity');self.control('freezing_pressure_compensator_derivative_detected')
    def test_axial_velocity_square_and_old_Haar_crosses_not_omitted(self):
        self.control('omitting_axial_velocity_square_detected');self.control('old_Haar_crosses_cannot_be_reset')
    def test_complete_covariance_and_velocity_moments_retained(self):
        self.control('complete_covariance_unchanged_all_six_nonzero');self.control('reserved_velocity_moments_exact')
    def test_radial_viscosity_only_vanishes_in_moments_and_axial_stays(self):
        self.control('radial_vector_viscosity_has_zero_plain_and_weighted_moments');self.control('axial_viscosity_is_retained')
    def test_positive_polynomials_match_independent_exact_substitution(self):self.control('positive_polynomials_match_direct_substitution')
    def test_unreduced_physical_terms_are_bounded_at_C_plus_one(self):self.control('all_four_unreduced_physical_budgets_below_C_plus_one')
    def test_complete_complex_covariance_estimate_keeps_mixed_terms(self):
        self.control('all_six_complex_covariance_changes_below_cross_plus_square');self.control('complex_mixed_terms_not_deleted')
    def test_updated_wave_norm_cannot_replace_original_primary_bound(self):
        self.control('updated_wave_norm_not_used_as_original_primary_bound');self.control('first_covariance_uses_separate_original_primary_and_signed_recipes')
    def test_first_postfast_rate_comes_from_three_complete_changes(self):self.control('first_slow_order_derived_from_three_complete_changes')
    def test_finite_Haar_sum_preserves_sharper_rates(self):
        self.control('both_accumulated_slow_orders_above_one');self.control('Haar_preservation_is_exact_finite_slow_sum')
    def test_all_three_scalars_gain_kappa_and_whole_moment_has_extra_epsilon(self):
        self.control('all_complete_postslow_scalars_C_plus_one');self.control('strict_kappa_improvement_over_full_mean_budget');self.control('whole_axial_pressure_residual_has_extra_epsilon')
    def test_coarse_full_initial_callback_cannot_be_relabelled(self):
        self.control('coarse_B_cannot_be_relabelled_Haar_one');self.control('coarse_bound_allows_unbounded_rescaled_Haar_norm')
    def test_initial_defects_from_original_seed_and_complete_first_update(self):
        self.control('initial_defects_from_original_SeedBudget_not_coarse_initial_callback');self.control('first_postfast_all_three_complete_changes')
    def test_first_fast_keeps_full_old_means_and_latest_slow_is_fresh(self):
        self.control('first_fast_uses_FULL_initial_means');self.control('fresh_slow_uses_current_scalar_derived_C_envelope')
    def test_current_radial_rho_and_whole_pressure_moment_keep_same_state(self):
        self.control('radial_compensator_keeps_current_rho');self.control('whole_axial_moment_keeps_pressure')
    def test_nested_original_seed_scalar_jet_counts_are_traced(self):
        self.control('postslow_initial_primary_scalar_channel_uses_n_plus_13');self.control('whole_axial_initial_primary_scalar_channel_uses_n_plus_14')
    def test_compact_tail_high_jets_and_regular_D_not_promoted(self):
        self.control('same_previous_final_rows_and_compact_high_jet_orders');self.control('regular_nonzero_source_not_promoted')
    def test_actual_weighted_input_gap_and_certificate_refusal(self):
        self.control('original_weighted_input_gap_not_hidden');self.control('finite_recipe_cannot_certify_actual_inputs')
        with self.assertRaises(NotImplementedError):fixture().certify_actual_inputs()
    def test_order_boundary_and_type_validation(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):fixture().old_Haar(-1)
        with self.assertRaises(ValueError):haar_slow_reset(0,{'u':P(1)},{'x':P(1),'u':P(1),'z':P(1)},P(1),P(1))
    def test_axial_output_uses_next_jet_and_whole_multiplier(self):
        j=fixture();self.assertEqual(j.axial_residual(1,P(7)).c,(2*j.whole_axial_moment(2,P(7))).c)
    def test_seven_frozen_dependencies_and_new_record_reproduce(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),7)
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
    def test_frozen_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
