#!/usr/bin/env python3
"""Actual fast-row gain controls; no finite diagnostic proves the original inputs."""
from copy import deepcopy
from fractions import Fraction as Q
from math import factorial
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from budgets import P,rates,flat_recipe,materialize_flat,difference


class FastRowGainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        a=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(a);self.assertTrue(all(a))

    def test_full_radial_source_subtraction(self):self.control('complete_radial_source_difference')
    def test_full_compensated_pressure_subtraction(self):self.control('full_linear_pressure_reconstruction_difference')
    def test_only_fast_time_removed_from_theta(self):self.control('tangential_difference_excludes_only_fast_time')
    def test_only_fast_time_removed_from_axial(self):self.control('axial_difference_excludes_only_fast_time')
    def test_pressure_viscosity_and_cylindrical_omissions_detected(self):
        for s in ('pressure_gradient_cannot_be_omitted','viscosity_cannot_be_omitted','cylindrical_connection_cannot_be_omitted'):self.control(s)

    def test_old_crosses_and_whole_pressure_compensator_retained(self):
        self.control('old_new_mean_crosses_are_nonzero');self.control('full_moving_pressure_subtraction_keeps_rho')

    def test_both_complete_mean_curls_divergence_free(self):self.control('both_mean_curls_exactly_divergence_free')
    def test_full_covariance_row_change(self):
        for s in ('complete_covariance_radial_difference','complete_covariance_theta_difference','complete_covariance_axial_and_pressure_difference','all_six_covariance_entries_present'):self.control(s)

    def test_independent_unreduced_physical_norm_bounds(self):self.control('independent_unreduced_physical_norm_')
    def test_all_old_remainder_terms_exceed_D(self):self.control('every_old_difference_term_exceeds_D')
    def test_wave_change_sets_postwave_D(self):self.control('wave_sets_exact_postwave_D')
    def test_exact_finite_A1_margin_and_no_extra_time_loss(self):
        self.control('strict_finite_A1_gain');self.control('A1_has_no_temporal_inverse_plus_four')

    def test_raw_A1_constant_and_factorial_conversions(self):self.control('independent_raw_A1_constant_and_all_conversions')
    def test_two_viscous_graph_losses_are_necessary(self):self.control('viscosity_two_graph_losses_cannot_be_erased')
    def test_current_feedback_and_regular_source_gain_are_conditional(self):
        self.control('fresh_actual_means_close_feedback_at_D');self.control('next_regular_source_at_D_only_conditionally')

    def test_large_default_dependency_is_exposed(self):
        self.control('actual_small_kappa_high_jet_dependency_is_not_hidden')
        r=flat_recipe(0);self.assertEqual((r['integration_order'],r['raw_source_order']),(99999,200001))
        self.assertEqual(flat_recipe(8)['raw_source_order']+3,200012)

    def test_exact_actual_A1_time_error_and_sign(self):
        for s in ('actual_fast_time_error_is_A1_of_current_residual','current_residual_A1_error_nonzero','wrong_fast_time_error_sign_detected'):self.control(s)

    def test_actual_A1_velocity_and_moving_cutoff_divergence(self):
        self.control('actual_A1_velocity_nonzero');self.control('full_moving_cutoff_divergence_exact');self.control('moving_cutoff_Z_term_is_required')

    def test_same_old_fast_slow_and_complete_wave_recipe_sum(self):self.control('same_old_fast_slow_complete_wave_sum_')
    def test_every_row_contribution_and_pressure_is_nonzero(self):
        for s in ('old_fast_difference_nonzero_','old_slow_difference_nonzero_','complete_wave_difference_nonzero_','full_pressure_in_both_mean_changes_and_covariance'):self.control(s)

    def test_same_scalar_constructor_and_pole_floor_contract(self):
        self.control('scalar_C_join_uses_same_existing_mean_object');self.control('original_higher_pole_coverage_and_floor_retained')

    def test_actual_inputs_and_uninstantiated_inverse_not_promoted(self):
        self.control('actual_input_enclosures_not_promoted');self.control('uninstantiated_high_jet_inverse_callback_is_blocked')

    def test_raw_to_factorial_through_conversion_bound(self):
        for n in range(13):self.assertLessEqual(sum(Q(5**j,factorial(j)) for j in range(n+1)),6**n)

    def test_materializer_guard_is_explicit_execution_limit(self):
        with self.assertRaises(ValueError):materialize_flat(flat_recipe(0),L=1,Bn=P(1),AJ=P(1),source=P(1),inverse_floor=P(1),rho_g=Q(1),weight_factors=P(1))

    def test_bad_ranges_and_incomplete_fields_rejected(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):difference(0,{}, {},{})
        with self.assertRaises((ValueError,TypeError)):flat_recipe(-1)

    def test_eight_dependencies_no_new_gain_hypothesis_and_frozen_record(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))

    def test_scope_and_dependency_mutations_rejected(self):
        for key,v in self.record['flags'].items():
            r=deepcopy(self.record);r['flags'][key]=not v
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['dependency_evidence_sha256'][next(iter(r['dependency_evidence_sha256']))]='bad'
        with self.assertRaises(ValueError):audit.validate(r)

    def test_existing_record_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
