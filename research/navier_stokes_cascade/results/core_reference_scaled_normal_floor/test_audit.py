"""Original reference-relative normal floor and unchanged inverse tests."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
import audit
from controls import review,fixture
from bounds import normal,relative_error_floor,certificate,coefficientwise_leq,P

class ReferenceNormalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=review();cls.record=audit.generate()
    def control(self,key):
        vals=[v for k,v in self.r['checks'].items() if key in k]
        self.assertTrue(vals,key);self.assertTrue(all(vals),key)
    def test_relative_triangle_information_retained(self):self.control('relative_triangle_slack_is_theta_times_reference_gap');self.control('reference_relative_triangle_preserved')
    def test_full_vector_norm_bounds_include_error_direction(self):self.control('exact_squared_normal_floor')
    def test_31_over_32_is_sharp_at_full_antiparallel_error(self):self.control('relative_31_over_32_bound_is_attained')
    def test_coarse_floor_alone_shortcut_is_invalid(self):self.control('coarse_uniform_floor_and_calibration_alone_do_not_prove_phase_floor');self.control('coarse_countermodel_violates_retained_full_error_bound')
    def test_existing_Riccati_slack_factor(self):self.control('original_Riccati_factor_leaves_one_over_32');self.control('original_Riccati_slack_closes_relative_reference_floor')
    def test_original_reference_calibration_with_three_halves(self):self.control('original_three_halves_dominates_phase_beta')
    def test_original_ceiling_carrier_factor_four(self):self.control('same_ceiling_k_retains_both_calibration_bounds');self.control('calibration_factor_four_cannot_be_dropped');self.control('original_rounding_upper_gap_factorization')
    def test_recover_phase_inverse_without_optional_native_margin(self):self.control('requested_phase_inverse_needs_no_extra_native_margin');self.control('coarse_bmin_square_can_be_below_calibrated_beta')
    def test_existing_bands_lengths_and_original_inputs_retained(self):self.control('existing_band_length_and_primitives_retained');self.control('uniform_profile_constant_not_S_power')
    def test_unknown_band_and_unlinked_length_rejected(self):self.control('unknown_early_band_rejected_without_reselection');self.control('length_semantics_lost_rejected_without_reselection')
    def test_new_conservative_norm_not_inserted_into_old_length(self):self.control('new_Mz_cannot_replace_old_norm_rejected_without_reselection')
    def test_minimum_preserves_better_original_reciprocal(self):self.control('minimum_retains_old_bound_when_it_is_better')
    def test_optional_native_failure_does_not_block_original_route(self):self.control('optional_native_comparison_can_fail_while_original_reference_route_passes');self.control('optional_failure_does_not_replace_original_Mz_Db_or_length')
    def test_original_tangent_length_and_all_full_callbacks_retained(self):
        self.control('original_tangent_and_length_not_recomputed');self.control('full_source_callback_retained');self.control('all_full_normal_matrix_and_profile_callbacks_retained');self.control('original_tangent_certificate_retains_all_gates')
    def test_complete_inverse_positive_recurrences_monotone(self):self.control('coefficientwise_monotone_at_same_other_inputs');self.control('ambient_majorant_also_monotone')
    def test_same_cutoffs_preserve_inverse_improvement(self):self.control('same_cutoff_preserves_monotonicity')
    def test_pressure_majorant_uses_the_new_valid_denominator(self):self.control('new_reciprocal_strictly_improves_zero_order_pressure_majorant')
    def test_unprojected_source_pressure_does_not_vanish_with_velocity(self):self.control('pure_normal_source_velocity_projection_zero');self.control('pure_normal_source_pressure_not_zero');self.control('projecting_pressure_source_erases_nonzero_term')
    def test_moving_normal_term_retained(self):self.control('moving_normal_pressure_cannot_be_dropped')
    def test_no_actual_primitive_or_historical_band_promotion(self):self.control('actual_primitive_and_historical_band_claim_not_promoted');self.control('calibrated_bound_not_inferred_from_coarse_uniform_floor_alone')
    def test_invalid_reference_and_error_bounds_rejected(self):
        for args in ((0,1,0),(2,1,0),(1,1,-1),(1,1,1)):
            with self.assertRaises(ValueError):relative_error_floor(*args)
        with self.assertRaises(TypeError):relative_error_floor(1,1,0.5)
    def test_incomplete_original_primitives_and_wrong_types_rejected(self):
        p,sp,ell,_=fixture();q=dict(p);q.pop('Mz')
        with self.assertRaises(ValueError):certificate(q,sp,ell)
        q=dict(p,u=0.5)
        with self.assertRaises(TypeError):certificate(q,sp,ell)
        with self.assertRaises(TypeError):coefficientwise_leq(1,P(1))
    def test_no_call_to_selector_in_production_certificate(self):
        p,sp,ell,_=fixture();old=normal.selector.select_tail
        def forbidden(*a,**k):raise AssertionError('Band reselection forbidden')
        normal.selector.select_tail=forbidden
        try:self.assertTrue(all(certificate(p,sp,ell)['gates'].values()))
        finally:normal.selector.select_tail=old
    def test_all_controls_eight_historical_records(self):
        self.assertTrue(all(self.r['checks'].values()));self.assertTrue(all(self.record['gates'].values()))
        self.assertEqual(len(self.record['dependency_evidence_sha256']),8)
    def test_frozen_roundtrip(self):
        audit.validate(json.loads(json.dumps(self.record)))
        p=Path(__file__).with_name('evidence.json')
        if p.exists():audit.validate(json.loads(p.read_text()))
    def test_frontier_and_scope_mutations_rejected(self):
        for key in ('new_analytical_assumptions_or_selector_thresholds_added','actual_primitive_enclosures_unknown_historical_band_or_generated_source_bounds_certified','unconditional_pulse_complete_residual_infinite_field_or_blowup_proved'):
            r=deepcopy(self.record);r['flags'][key]=True
            with self.assertRaises(ValueError):audit.validate(r)
        r=deepcopy(self.record);r['remaining_obligations'].pop()
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_no_frozen_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--output',str(p)])
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
