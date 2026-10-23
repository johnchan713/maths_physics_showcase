#!/usr/bin/env python3
"""Physical exponent, exact tail-allocation and uniform-domain diagnostics."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from bounds import scales,compact_orders,absorption_constant,gaussian_recipe,regular_transfer,hypothetical_steps
from controls import review


class PhysicalTailTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()

    def control(self,fragment):
        values=[v for k,v in self.r['checks'].items() if fragment in k]
        self.assertTrue(values,fragment);self.assertTrue(all(values),fragment)

    def test_exact_current_residual_allocation(self):
        self.control('exact_fresh_cancellation');self.control('full_residual_is_regular_plus_compact_Gaussian_and_assigned_base')
    def test_full_radial_pressure_compensator_stays_regular(self):self.control('pressure_compensator_remains_in_radial_regular_row')
    def test_compact_tail_signs_and_omissions(self):
        self.control('A0_and_A1_have_opposite_derived_signs');self.control('omitting_A0_or_A1_is_detected')
    def test_old_Q1_rows_not_double_counted(self):self.control('old_Q1_rows_not_added_again')
    def test_residual_tail_is_not_a_velocity(self):self.control('residual_tail_not_a_velocity_or_nonlinear_source')
    def test_actual_A1_velocities_remain_in_regular_interactions(self):
        self.control('actual_A1_velocity_stays_in_both_regular_interactions');self.control('compact_residual_and_actual_A1_velocity_are_distinct_symbols')
    def test_angular_projection_keeps_compact_tail(self):
        self.control('angular_projection_does_not_discard_compact_mean_tail');self.control('temporal_wave_tail_has_zero_angular_projection')
    def test_spatial_pressure_cutoff_not_called_temporal_Gaussian(self):self.control('spatial_pressure_cutoff_not_erased_as_temporal_Gaussian')
    def test_compact_physical_counts_have_strict_margin(self):self.control('strict_physical_margin')
    def test_counts_are_minimal_for_displayed_sufficient_rule(self):self.control('previous_integration_integer_insufficient_for_this_rule')
    def test_epsilon_and_physical_q_powers_not_confused(self):self.control('physical_power_not_epsilon_power')
    def test_finite_current_source_jet_accounting(self):self.control('finite_current_source_jet')
    def test_current_regular_upper_bound_does_not_certify_physical_decay(self):self.control('first_regular_native_physical_bound_has_negative_power')
    def test_hypothetical_future_stage_is_not_certified_induction(self):self.control('hypothetical_stage_threshold_is_only_conditional')
    def test_illustrative_physical_counts_exceed_normalized_count(self):
        self.control('illustrative_A0_physical_zero_power');self.control('illustrative_A1_physical_zero_power');self.control('fixed_normalized_A0_C_count_does_not_alone_prove_physical_flatness')
    def test_polynomial_absorption_without_new_band(self):self.control('exponential_series_proves_polynomial_absorption_without_new_band')
    def test_log_chart_enlargement_constants(self):
        self.control('natural_log_chart_upper_from_enlargement');self.control('natural_log_chart_lower_from_enlargement')
    def test_band_Gaussian_transfers_with_safe_constant(self):
        self.control('band_Gaussian_dominates_chosen_physical_Gaussian');self.control('Gaussian_recipe_retains_explicit_band_to_natural_log_constant')
    def test_Gaussian_square_completion(self):
        self.control('Gaussian_recipe_keeps_square_completion_constant');self.control('Gaussian_square_completion_exact')
    def test_current_embedded_full_gr_recipe_remains_deferred(self):self.control('A0_uses_current_full_embedded_gr_without_evaluating_huge_jets')
    def test_nested_original_A1_recipe_not_truncated(self):self.control('A1_keeps_nested_original_high_jet_recipe_deferred')
    def test_actual_raw_profile_and_edge_geometry_kept(self):
        self.control('physical_tail_uses_raw_e0_and_e1_profiles_separately');self.control('same_original_width_gap_and_pole_domains')
    def test_tail_recipes_are_not_actual_input_certificates(self):self.control('tails_are_not_materialized_actual_input_certificates')
    def test_later_cutoffs_and_assigned_base_stress_not_double_counted(self):
        self.control('spatial_cutoffs_and_later_diagonals_not_called_Gaussian');self.control('base_ledger_does_not_readd_the_cancelled_bare_stress')
    def test_illustrative_stage_count_and_no_actual_lower_bound(self):
        self.assertEqual(hypothetical_steps(0,0,Q(1,200))['additional_steps'],1002)
        self.assertFalse(regular_transfer(0,Q(1,200))['actual_residual_lower_bound'])
    def test_invalid_parameters_and_nonpositive_absorption_margin_rejected(self):
        for h in (0,Q(1,100),-1,0.005):
            with self.assertRaises((ValueError,TypeError)):scales(h)
        with self.assertRaises(ValueError):compact_orders(2,0,0,Q(1,200))
        with self.assertRaises(ValueError):absorption_constant({0:Q(1)},1,0)
        with self.assertRaises(ValueError):gaussian_recipe(0,1,1,1)
        with self.assertRaises(ValueError):regular_transfer(0,Q(1,200),reserve=0)
    def test_seven_dependencies_frozen_record_and_no_new_assumptions(self):
        self.assertEqual(len(self.record['dependency_evidence_sha256']),7)
        self.assertEqual(self.record['protocol']['new_analytical_assumptions'],[])
        audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_obligation_and_dependency_mutations_rejected(self):
        for k,v in self.record['flags'].items():
            bad=deepcopy(self.record);bad['flags'][k]=not v
            with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(bad)
        bad=deepcopy(self.record);bad['dependency_evidence_sha256'][next(iter(bad['dependency_evidence_sha256']))]='wrong'
        with self.assertRaises(ValueError):audit.validate(bad)
    def test_frozen_record_not_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'keep.json';p.write_text('sentinel')
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
            self.assertEqual(p.read_text(),'sentinel')


if __name__=='__main__':unittest.main()
