#!/usr/bin/env python3
"""Multiharmonic source and sharp mean-feedback gap regression checks."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review
from budgets import P,harmonic_constants,increment_recipe,rates,source_pole,NextResidualBudget


class MultiharmonicResidualTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        rows=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(rows);self.assertTrue(all(rows))

    def test_complete_curl_divergence_for_all_four_signed_inputs(self):self.control('complete_curl_divergence')
    def test_all_48_phase_components_use_input_ratio(self):
        self.assertEqual(self.r['parts']['harmonics']['exact_phase_components'],48);self.control('exact_input_ratio')

    def test_old_ratio_one_fails_for_a_nonzero_second_mode_pair(self):self.control('single_frequency_ratio_one_underestimates')
    def test_independent_full_fourier_square_norm_majorant(self):self.control('independent_full_fourier_square_has_l1_majorant')
    def test_zero_fourth_and_negative_outputs_are_kept(self):self.control('zero_and_highest_fourth_harmonic_are_retained')
    def test_conjugate_reality_of_complete_transport(self):self.control('conjugate_complete_transport_reality')
    def test_no_output_frequency_divisor(self):self.control('all_three_output_sets_use_inputs_not_output_division')
    def test_exact_mean_phase_advection_is_nonzero(self):self.control('mean_feedback_is_exact_nonzero_phase_advection')
    def test_reverse_term_vanishes_in_sharp_witness(self):self.control('reverse_mean_transport_is_zero_in_witness')
    def test_witness_wave_is_full_divergence_free(self):self.control('wave_normal_constraint_and_full_divergence_exact')
    def test_witness_mean_has_both_exact_velocity_moments(self):self.control('mean_is_divergence_free_and_has_both_zero_velocity_moments')
    def test_witness_uses_compact_overlapping_potentials(self):
        self.control('compact_potential_endpoints_preserve_axis_support');self.control('f_and_wave_overlap_nontrivially')

    def test_polynomial_cost_cannot_absorb_positive_order_gap(self):self.control('positive_gap_exponential_cannot_be_absorbed_by_fixed_polynomial')
    def test_positive_recipe_substitutes_exactly_at_rational_scales(self):self.control('positive_polynomial_rational_substitution')
    def test_single_frequency_linear_limit_matches_historical_operator(self):self.control('single_frequency_linear_recipe_matches_historical_operator')
    def test_full_pressure_uncut_curl_and_fresh_mean_omissions_detected(self):
        for k in ('Pi','R','dM'):self.control('omission_of_'+k+'_detected')

    def test_second_mode_mixed_viscosity_and_damping_factors(self):self.control('second_mode_increases_mixed_and_damping_costs')
    def test_coarse_size_only_reopens_old_order_source(self):self.control('coarse_fresh_means_do_not_supply_improved_residual_order')
    def test_candidate_sharper_mean_targets_are_only_sufficient(self):self.control('desired_fast_and_slow_input_orders_are_sufficient_only')
    def test_updated_mean_wave_feedback_requires_nine_pole(self):self.control('new_mean_feedback_sets_source_pole_nine_p')
    def test_generated_full_pressure_and_current_mean_rows_survive(self):self.control('full_current_pressure_and_all_scalars_retained')
    def test_exact_raw_inverse_coefficient_transfer(self):self.control('raw_regular_source_converts_exactly_to_inverse_interface')
    def test_feedback_is_added_to_the_wave_stage(self):self.control('mean_feedback_not_omitted_or_promoted_to_wave_stage')
    def test_full_source_pressure_and_uncut_coefficient_curl(self):self.control('full_unprojected_source_pressure_and_uncut_curl_positive')
    def test_uncut_coefficients_are_not_already_temporally_cut(self):self.control('uncut_principal_bound_does_not_contain_temporal_cutoff')
    def test_finite_composed_jet_loss_and_same_band(self):
        self.control('finite_original_mean_jet_loss_is_31');self.control('same_parent_shift_width_and_band_preserved')

    def test_only_coarse_order_is_advertised_at_inverse_interface(self):self.control('inverse_interface_explicitly_carries_only_coarse_order_and_nine_pole')
    def test_nontrivial_factorial_conversion(self):
        b=object.__new__(NextResidualBudget);b.regular_source=lambda n:{'F':P({2:7})}
        self.assertEqual(b.raw_source(3).c,{2:Q(42)})

    def test_invalid_harmonics_orders_and_inputs_rejected(self):
        for h in ((),(0,),(-1,1,1),(1,),(-1.0,1),(-True,True)):
            with self.assertRaises(ValueError):harmonic_constants(h)
        for k in (0,Q(1,4),-1,0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
        for p in (-1,0.5):
            with self.assertRaises((ValueError,TypeError)):source_pole(p)
        with self.assertRaises(ValueError):increment_recipe(0,{})
        with self.assertRaises((ValueError,TypeError)):increment_recipe(-1,{})

    def test_six_records_and_no_new_proved_gain_hypothesis(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),6)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        self.assertTrue(all(self.record['gates'].values()));self.control('actual_primitive_enclosures_not_promoted')

    def test_frozen_record_reproduces(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_all_scope_flags_and_remaining_obligations_are_bound(self):
        for k,v in self.record['flags'].items():
            r=deepcopy(self.record);r['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)

    def test_frozen_output_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
