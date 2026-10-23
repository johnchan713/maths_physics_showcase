"""Next current-source eight-mode constrained-wave regression tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from budgets import inverse,P,FiniteHarmonicWave,next_domains,NextSourcedRecipes,rates
from controls import review

class NextSourcedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        values=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(values,key);self.assertTrue(all(values),key)
    def test_all_eight_full_forced_equations_and_constraints(self):self.control('full_forced_equation_and_constraint')
    def test_full_moving_normal_and_unprojected_source_pressure(self):self.control('full_source_normal_moving_pressure')
    def test_zero_entrance_and_parameter_entrance_jets(self):self.control('zero_entrance_and_parameter_entrance_jets')
    def test_whole_box_mixed_velocity_and_pressure_majorants(self):
        self.control('full_box_mixed_velocity_majorant');self.control('full_box_mixed_pressure_majorant')
    def test_fourth_damping_and_same_Gaussian_ambient(self):
        self.control('fourth_mode_generator_retains_sixteen_times_damping');self.control('third_and_fourth_modes_use_same_Gaussian_ambient_not_short_residence')
    def test_pure_normal_pressure_not_deleted_with_projected_source(self):self.control('pure_normal_pressure_not_bounded_from_projected_force')
    def test_direct_source_pressure_not_inverse_square_in_harmonic(self):self.control('pressure_direct_source_is_not_inverse_square_in_m')
    def test_ambient_normal_cannot_use_high_mode_residence(self):self.control('normal_ambient_cannot_use_high_mode_residence')
    def test_manufactured_forcing_not_flat_entrance_enclosure(self):self.control('manufactured_source_is_not_flat_entrance_certificate')
    def test_all_eight_cut_modes_and_full_uncut_pressure_retained(self):
        self.control('all_eight_signed_input_modes_are_summed');self.control('full_cut_pressure_and_coefficient_curl_are_positive')
    def test_fourth_cut_pressure_keeps_every_cutoff_jet(self):self.control('fourth_mode_cut_pressure_retains_all_cutoff_jets')
    def test_generated_complete_covariance_crosses_and_squares(self):self.control('full_covariance_cross_and_signed_squares_are_retained')
    def test_original_normal_floor_frame_recipe_recomputed(self):self.control('original_normal_floor_and_frame_certificate_recomputed')
    def test_first_generated_constructor_insufficient_for_modes_three_four(self):self.control('old_single_constructor_omits_actual_third_and_fourth_modes')
    def test_actual_source_and_moving_normal_jet_requests(self):
        self.control('order_two_complete_wave_reads_source_through_three');self.control('order_two_complete_wave_reads_moving_normal_through_three')
    def test_complete_curl_divergence_and_all_sixty_four_input_pairs(self):
        self.control('complete_curl_divergence');self.control('input_ratio_and_zero_output_safe')
    def test_eighth_residual_output_and_H4_ratio4(self):
        self.control('new_complete_wave_has_nonzero_eighth_harmonic_square');self.control('fourth_input_requires_H4_ratio4_and_outputs_through8')
    def test_all_six_covariances_include_new_high_mode_squares(self):
        self.control('full_covariance_has_all_six_nonzero_entries');self.control('absent_old_high_modes_still_change_all_six_covariances')
    def test_positive_wave_stage_polynomial_and_H4_constants(self):
        self.control('wave_stage_positive_polynomial_substitution');self.control('old_H2_constants_underestimate_H4_linear_and_transport')
    def test_stage_before_fresh_mean_feedback(self):self.control('no_fresh_mean_feedback_yet')
    def test_old_scalar_and_full_pressure_memory_kept(self):
        self.control('old_pressure_and_scalar_memory_not_reset');self.control('full_graph_pressure_source_and_compensation_retained')
    def test_sourced_velocity_not_confused_with_pressure_order(self):self.control('sourced_D_velocity_not_promoted_to_C0_pressure_order')
    def test_half_minus_kappa_wave_stage_gain_and_current_scalar_orders(self):
        self.control('wave_stage_regular_gain_is_half_minus_kappa');self.control('full_tangential_E_but_Haar_target_and_scalars_C0')
    def test_same_derived_current_source_not_new_fine_callback(self):self.control('same_current_source_program_not_new_fine_source_callback')
    def test_source_normal_and_padded_old_jet_graph(self):
        self.control('complete_wave_source_jet_and_full_row_jet_ladders');self.control('padded_old_graph_covers_all_new_requests_without_collapsing_old_DAGs')
    def test_old_4q_domain_does_not_cover_new_6q_source(self):self.control('new_6q_pressure_domain_not_old_4q')
    def test_original_new_kernel_and_oriented_primitive_recipes(self):
        self.control('new_6q_12q_original_kernel_recipes_and_Ie6q_retained');self.control('moderate_original_shift_covers_new_6q_12q_pressure')
    def test_archived_shifts_and_unchanged_fields_geometry(self):
        self.control('all_archived_shifts_and_separate_oriented_constants_kept');self.control('no_new_actual_band_fields_or_huge_constants_materialized')
    def test_current_regular_source_only_in_tail_ledger(self):self.control('new_source_tail_uses_CURRENT_regular_source_only')
    def test_new_A0_not_inferred_from_pressure_size_or_stale_tail(self):
        self.control('recomputed_pressure_cannot_reuse_stale_A0_residual');self.control('new_pointwise_radial_A0_not_claimed_from_tangential_pressure_size')
    def test_support_and_actual_source_materialization_limits(self):
        self.control('whole_path_and_flat_entrance_support_still_required');self.control('actual_high_jet_source_materialization_refused');self.control('no_fresh_mean_complete_cycle_or_actual_validity_certificate')
    def test_invalid_current_source_modes_and_poles_rejected(self):
        with self.assertRaises(TypeError):NextSourcedRecipes(object())
        with self.assertRaises(ValueError):FiniteHarmonicWave({1:inverse.fixture(harmonic=1)},lambda n:P(1),lambda n:P(1),lambda n:P(1))
        modes={m:inverse.fixture(harmonic=1) for m in (1,2,3,4)}
        with self.assertRaises(ValueError):FiniteHarmonicWave(modes,lambda n:P(1),lambda n:P(1),lambda n:P(1))
        with self.assertRaises(ValueError):next_domains(-1)
        for k in (0,Q(1,4),0.1):
            with self.assertRaises((ValueError,TypeError)):rates(k)
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
