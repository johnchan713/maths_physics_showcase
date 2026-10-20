"""Actual fifth native propagation, memories, earlier root and strict limits."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import native_bounds as n

class ActualFifthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_cover_sums_all_four_E_regions(self):self.assertEqual(n.regional_cover(2**16,16)['E_C5_cover'],1+21*16+22*16+1000*2**16)
    def test_cover_sums_all_four_U_regions(self):self.assertEqual(n.regional_cover(2**16,16)['U_C5_cover'],36)
    def test_E5_cover_below_A2_for_worst_P(self):
        for A in (2**16,2**17,2**20):self.assertLess(n.regional_cover(A,A)['E_C5_cover'],A*A)
    def test_U5_cover_below_A(self):self.assertLess(n.regional_cover(2**16,16)['U_C5_cover'],2**16)
    def test_relative_annular_three_bumps_and_two_U_bumps(self):
        r=n.regional_cover(2**16,16);self.assertEqual(r['annular_relative_C5'],3*64*Q(572,10**16));self.assertEqual(r['annular_U_C5'],8+2*64*Q(572,10**16))
    def test_cover_invalid_parameters_rejected(self):
        for A,P,d in ((2,16,0),(2**16,15,0),(2**16,2**17,0),(2**16,16,Q(1,100))):
            with self.assertRaises(ValueError):n.regional_cover(A,P,d)
    def test_cover_requires_exact_inputs(self):
        with self.assertRaises(TypeError):n.regional_cover(65536.,16)
    def test_no_FINAL_gate_promotion_from_native_cover(self):self.assertIs(n.regional_cover(2**16,16)['actual_FINAL_field_or_gate_verified'],False)
    def test_two_actual_interval_precisions(self):self.assertEqual([v['interval_digits'] for v in self.r['parameters']],[80,110])
    def test_actual_continuation_J_below_B_over100(self):self.assertTrue(all(v['checks']['prescribed_log_CE_C5_below_J'] and v['checks']['J_below_B_over100'] for v in self.r['parameters']))
    def test_both_g_errors_compared_without_width_exponentiation(self):self.assertTrue(all(v['checks']['native_g_error_below_j_over2'] and v['checks']['width_g_error_below_j_over2'] for v in self.r['parameters']))
    def test_core_regular_scaled_field_not_compact_floor(self):self.assertTrue(all(v['checks']['native_regular_CE_over_sqrtX_C5_below_expB'] for v in self.r['parameters']))
    def test_shape_tame_exponential_fits_Eb(self):self.assertTrue(all(v['checks']['shape_tame_C5_fits_same_exp204B'] and v['checks']['actual_Eb_below_one'] for v in self.r['parameters']))
    def test_integrable_core_pressure_and_log_length(self):self.assertTrue(all(v['checks']['integrable_core_Cp_and_actual_log_length_fit_exp412B'] for v in self.r['parameters']))
    def test_all_five_normalized_target_rows_small(self):self.assertTrue(all(v['checks']['actual_normalized_join_base_C5_below_epsilon_over4'] for v in self.r['parameters']))
    def test_annular_normalizer_C5_is_25(self):self.assertEqual(self.r['algebra']['f_inverse_square_C4_C5'],['25','25'])
    def test_energy_row_new_cost_8125(self):self.assertEqual(self.r['algebra']['annular_energy_row_constant'],8125)
    def test_energy_offset_square_cancellation(self):self.assertTrue(self.r['algebra']['checks']['annular_energy_g_square_cancels_linear_term'])
    def test_same_annular_root_contraction_C5(self):self.assertTrue(all(v['checks']['same_C5_root_contraction'] and v['checks']['same_C5_root_coefficient_below_5p72e14'] for v in self.r['parameters']))
    def test_outer_prefix_finite_exponential_budget(self):self.assertTrue(self.r['algebra']['checks']['outer_finite_tame_exponential_below1000'])
    def test_new_pressure_envelope_not_old_Ahi(self):
        r=self.r['algebra'];self.assertEqual(r['pressure_C4_A_power'],1546);self.assertGreater(r['pressure_C4_A_power'],r['preexisting_native_C4_envelope_A_power'])
    def test_actual_normalized_fifth_and_entry_powers(self):
        r=self.r['algebra'];self.assertEqual((r['native_field_F_E5_A_power'],r['native_field_F_U5_A_power'],r['native_moment_and_entry_C5_A_power']),(2,1,8))
    def test_full_pressure_R0_and_T2_constants_charged(self):self.assertTrue(self.r['algebra']['checks']['fourth_pressure_A1546_controls_complete_R0_and_T2'])
    def test_original_Volterra_majorants_at_minimum_A(self):self.assertTrue(all(v['T1_within_A517'] and v['T2_within_A517'] for v in self.r['volterra']['rows']))
    def test_original_fourth_pressure_majorants_fit_new_envelope(self):self.assertTrue(all(v['p1_within_A1546'] and v['p2_within_A1546'] for v in self.r['volterra']['rows']))
    def test_entering_T2_remains_at_both_poles(self):self.assertTrue(all(v['incoming_T2_positive'] for v in self.r['volterra']['rows']))
    def test_zero_length_keeps_incoming_fifth_memory(self):self.assertTrue(self.r['volterra']['checks']['zero_length_still_keeps_all_entering_memory'])
    def test_fifth_quadratic_exact_at_three_eta_centers(self):self.assertTrue(all(v['exact'] for v in self.r['quadratic']['rows']))
    def test_freezing_eta_dependent_Q_is_detected(self):self.assertTrue(all(v['freezing_Q_eta_detected'] for v in self.r['quadratic']['rows']))
    def test_fourth_Q_derivative_survives_but_fifth_vanishes(self):self.assertTrue(all(v['fourth_Q_coefficient_nonzero'] and v['fifth_Q_coefficient_zero'] for v in self.r['quadratic']['rows']))
    def test_nine_actual_stage_and_eta_cases(self):self.assertEqual(len(self.r['flows']['rows']),9);self.assertEqual({v['stage'] for v in self.r['flows']['rows']},{'activation','axial_taper','shear_blend'})
    def test_both_angular_endpoints_exercised(self):self.assertEqual({v['eta'] for v in self.r['flows']['rows']},{'-1','1/3','1'})
    def test_actual_fifth_shear_identities(self):self.assertTrue(all(v['both_shears_total5_exact'] for v in self.r['flows']['rows']))
    def test_all_actual_five_forward_density_identities(self):self.assertTrue(all(v['five_actual_moment_flows_total5_exact'] for v in self.r['flows']['rows']))
    def test_actual_full_pressure_flow(self):self.assertTrue(all(v['full_pressure_flow_total5_exact'] for v in self.r['flows']['rows']))
    def test_all_six_entering_fifth_memories_nonzero(self):self.assertTrue(all(v['five_entering_fifth_memories_nonzero'] and v['axis_fifth_memory_nonzero'] for v in self.r['flows']['rows']))
    def test_reference_zero_shear_not_copied(self):self.assertTrue(all(v['reference_zero_shear_not_copied'] for v in self.r['flows']['rows']))
    def test_actual_moments_not_reused_reference_moments(self):self.assertTrue(all(v['actual_and_reference_moments_not_reused'] for v in self.r['flows']['rows']))
    def test_shape_fifth_switch_commutes_with_eta(self):self.assertTrue(self.r['flows']['checks']['shape_fifth_angular_is_same_convex_combination'])
    def test_manufactured_not_actual_profile_or_root(self):self.assertIs(self.r['flows']['manufactured_not_selected_profile'],True);self.assertIs(self.r['quadratic']['manufactured_not_actual_root'],True)
    def test_parameters_not_materialized(self):self.assertTrue(all(v['enormous_C_G_H_N_or_actual_profile_materialized'] is False for v in self.r['parameters']))
    def test_invalid_precision_rejected(self):
        with self.assertRaises(ValueError):n.parameters(40)
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'0963f1379a76c213f16db60a5acb52c3cd0cb658');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_false_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_is_not_boolean_false(self):
        r=deepcopy(self.record);r['flags']['independent_foundation_review_or_formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_domain_or_cutoff_scope_not_erased(self):
        for k in ('domain','continuation','moments','annular','cover','pressure','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_fourth_maps_and_variance_obligations_retained(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_evidence_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(path)]);audit.validate(json.loads(path.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(path)])

if __name__=='__main__':unittest.main()
