#!/usr/bin/env python3
"""Haar stress and radial-mean controls; actual fields remain analytical inputs."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from budgets import rates,haar_rows


class HaarStressJoinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        a=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(a);self.assertTrue(all(a))

    def test_both_covariances_symmetric_and_transverse(self):
        self.control('all_tensors_symmetric');self.control('both_covariances_transverse_to_same_normal')

    def test_positive_rank_one_weights_prove_PSD(self):self.control('positive_rank_one_weights_give_PSD')
    def test_energy_trace_unchanged(self):self.control('full_energy_trace_unchanged')
    def test_full_radial_source_and_pressure_unchanged(self):self.control('complete_radial_source_and_pressure_unchanged')
    def test_all_three_scalar_defects_unchanged(self):self.control('all_three_scalar_defects_unchanged')
    def test_selected_fluxes_change_despite_same_moments(self):
        self.control('changed_selected_theta_flux');self.control('changed_selected_axial_flux')

    def test_both_spatial_Haar_rows_are_nonzero(self):
        self.control('nonzero_full_Haar_theta_row');self.control('nonzero_full_Haar_axial_row')

    def test_global_residual_moments_miss_spatial_error(self):self.control('both_residual_global_moments_zero')
    def test_Q1_also_misses_auxiliary_independent_error(self):self.control('both_rows_are_auxiliary_independent')
    def test_compact_endpoints_and_full_pressure_retained(self):
        self.control('every_covariance_entry_smoothly_compact_to_fixture_order');self.control('full_pressure_compact_at_both_ends')

    def test_moving_density_unit_mass_and_axial_pressure_memory(self):
        self.control('moving_density_has_exact_unit_mass');self.control('reference_axial_full_moving_pressure_memory_nonzero')

    def test_fine_flux_and_other_Haar_terms_have_sufficient_orders(self):
        self.control('fine_flux_minimum_is_C');self.control('all_nonflux_rows_stronger_than_C')

    def test_mixed_postslow_scalar_bound_against_independent_dynamic_budget(self):self.control('independent_mixed_postslow_')
    def test_radial_Haar_gain_is_B_plus_C(self):self.control('radial_postslow_gain_is_B_plus_C')
    def test_normalized_principal_stress_requirement_has_positive_gap(self):self.control('normalized_principal_gain_is_positive')
    def test_close_positive_weights_preserve_energy_but_not_flux(self):self.control('positive_close_weights_preserve_energy_not_flux')
    def test_actual_principal_identity_is_not_asserted(self):self.control('actual_matching_not_promoted')
    def test_every_fresh_postslow_scalar_and_radial_compensator_is_kept(self):
        self.control('fresh_postslow_complete_');self.control('radial_Haar_obstruction_uses_fresh_postslow_P')

    def test_whole_pressure_product_and_unpromoted_pointwise_pressure(self):
        self.control('whole_pressure_product_not_frozen');self.control('current_Haar_pressure_full_and_not_order_promoted')

    def test_leading_gap_cannot_be_silently_deleted(self):self.control('unmatched_leading_stress_retained_in_both_rows')
    def test_no_graph_or_fast_time_loss_after_Haar(self):self.control('no_graph_or_fast_time_coefficient_in_Haar_grammar')
    def test_same_current_mean_wave_and_covariance_ledger(self):
        self.control('fine_row_uses_all_current_means_and_wave');self.control('all_complete_covariance_differences_accumulate')

    def test_common_operator_poles_unchanged(self):self.control('common_higher_pole_recipes_unchanged')
    def test_no_Boolean_principal_matching_certificate(self):self.control('actual_principal_matching_cannot_be_asserted_by_a_flag')
    def test_finite_original_mean_input_dependency(self):self.assertEqual(self.r['parts']['composition']['highest_original_initial_jet'],31)
    def test_invalid_range_and_missing_mean_components(self):
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        with self.assertRaises(ValueError):haar_rows(0,{}, {},0,0,0)

    def test_eight_dependencies_and_frozen_record(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))

    def test_all_scope_and_remaining_obligation_mutations_rejected(self):
        for key,v in self.record['flags'].items():
            r=deepcopy(self.record);r['flags'][key]=not v
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_dependency_digest_mutation_rejected(self):
        r=deepcopy(self.record);r['dependency_evidence_sha256'][next(iter(r['dependency_evidence_sha256']))]='bad'
        with self.assertRaises(ValueError):audit.validate(r)

    def test_existing_frozen_record_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
