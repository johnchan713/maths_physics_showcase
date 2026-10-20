"""Exact nonlinear fifth jets, full incoming pressure and conditional gate scope."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import moving5 as b

class FifthLoopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_sixth_q_coefficient(self):self.assertTrue(self.r['q6']['checks']['exact_sixth_q_coefficient'])
    def test_sixth_q_not_fifth_required_for_R5(self):self.assertTrue(self.r['q6']['checks']['sixth_q_required_by_fifth_removable_R'])
    def test_normalized_q_z_derivative_cost6(self):self.assertEqual(self.r['q6']['differentiated_q_operator_cost'],6)
    def test_q_reciprocal_cost64(self):self.assertTrue(self.r['exponents']['checks']['sixth_q_reciprocal_norm_has64'])
    def test_root_complete_nonlinear_fifth_residual(self):self.assertTrue(self.r['scalar_root']['checks']['complete_nonlinear_fifth_implicit_residual_zero'])
    def test_root_keeps_pressure_fifth(self):self.assertTrue(self.r['scalar_root']['checks']['top_p5_coefficient_exact']);self.assertNotEqual(self.r['scalar_root']['mu5_p5_coefficient'],'0')
    def test_root_keeps_target_fifth(self):self.assertTrue(self.r['scalar_root']['checks']['top_s5_coefficient_exact'])
    def test_root_affine_not_high_input_squared(self):self.assertTrue(self.r['scalar_root']['checks']['high_fifth_inputs_affine_not_squared'])
    def test_root_all_lower_cross_terms(self):self.assertTrue(self.r['scalar_root']['checks']['lower_cross_terms_present'])
    def test_signed_root_origin_not_variance_derivative(self):self.assertTrue(self.r['scalar_root']['checks']['signed_origin_derivative_nonzero']);self.assertTrue(self.r['scalar_root']['checks']['signed_origin_p_derivative_zero'])
    def test_moving_inverse_complete_fifth_residual(self):self.assertTrue(self.r['moving']['checks']['complete_fifth_moving_inverse_residual_zero'])
    def test_moving_inverse_nonunit_divisor_retained(self):self.assertTrue(self.r['moving']['checks']['inverse_Jacobian_not_normalized_to_one']);self.assertTrue(self.r['moving']['checks']['fifth_inverse_keeps_actual_base_divisor'])
    def test_full_source_fifth_sensitivity(self):self.assertTrue(self.r['moving']['checks']['full_pulled_fifth_sensitivity_exact'])
    def test_frozen_angle_wrong_top_coefficient(self):self.assertTrue(self.r['moving']['checks']['frozen_angle_misses_nonzero_top_term'])
    def test_moved_source_affine(self):self.assertTrue(self.r['moving']['checks']['full_fifth_moving_source_affine_in_top_input'])
    def test_moving_lower_remainder_nonzero(self):self.assertTrue(self.r['moving']['checks']['nontrivial_lower_moving_remainder_retained'])
    def test_original_pressure_all_three_sixth_moment_rows(self):self.assertTrue(self.r['native_pressure']['checks']['full_pressure_principal_M6_S6_Pi6_exact'])
    def test_original_pressure_axial_U_M6_not_deleted(self):self.assertTrue(self.r['native_pressure']['checks']['nonzero_axial_U_M6_term_required'])
    def test_D_eta_M6_not_multiplied_by_wrong_d(self):self.assertTrue(self.r['native_pressure']['checks']['D_eta_M6_has_no_spurious_one_minus_eta_square'])
    def test_sixth_full_pressure_memory_nonzero(self):self.assertTrue(self.r['native_pressure']['checks']['full_pressure_sixth_memory_coefficient_nonzero']);self.assertIs(self.r['native_pressure']['Pi6_is_full_axis_plus_Cp'],True)
    def test_a5_total_order_six_boundary(self):self.assertTrue(self.r['native_pressure']['checks']['a5_needs_total_order_six_field_row'])
    def test_bs5_total_order_six_boundary(self):self.assertTrue(self.r['native_pressure']['checks']['bs5_needs_total_order_six_field_row'])
    def test_joint_five_not_joint_six(self):self.assertTrue(self.r['native_pressure']['checks']['joint_five_does_not_contain_1_5']);self.assertEqual(self.r['native_pressure']['total_field_order_required'],6)
    def test_zero_mean_primitive_and_fifth_commutation(self):self.assertTrue(self.r['primitive']['checks']['zero_mean_primitive_exact']);self.assertTrue(self.r['primitive']['checks']['fixed_external_phase_fifth_commutes'])
    def test_charge_both_primitives_four_not_two(self):self.assertEqual(self.r['primitive']['both_pure_angular_primitive_operator_cost'],4)
    def test_weighted_loop_axial_and_variance_means(self):
        for k in ('weighted_loop_first_density_exact','weighted_axial_density_exact','mean_one_uses_full_variance','axial_B_source_keeps_E_and_bs'):self.assertTrue(self.r['primitive']['checks'][k])
    def test_all_tame_exponents_strictly_fit(self):self.assertTrue(self.r['exponents']['checks']['all_fifth_tame_packets_have_strict_headroom'])
    def test_fifth_partition_count64(self):self.assertEqual(self.r['exponents']['exact_scalar_count_by_number_inputs'][1],64)
    def test_two_external_colors_not_new_fifth_y(self):self.assertTrue(self.r['exponents']['checks']['two_external_color_inflation_2048']);self.assertIs(self.record['flags']['fifth_y_or_full_joint_three_variable_fifth_loop_bound_derived'],False)
    def test_eight_input_count_bound(self):self.assertLess(max(self.r['exponents']['exact_scalar_count_by_number_inputs']),10**9)
    def test_Lambda_tame_gate(self):self.assertEqual(self.r['exponents']['Lambda5_gate'],'exp(1000000*G)*T5')
    def test_original_native_cutoff_gate(self):self.assertEqual(self.r['exponents']['native_cutoff_gate'],'exp(1100000*G+10*A^20)*Z5*R5')
    def test_gate_accepts_supplied_auxiliary_comparison(self):self.assertIs(b.input_gate(0,0,0,2**32)['auxiliary_B5_comparison'],True)
    def test_gate_rejects_supplied_auxiliary_comparison(self):self.assertIs(b.input_gate(2**100,0,0,2**32)['auxiliary_B5_comparison'],False)
    def test_exact_rational_auxiliary_gate(self):self.assertIs(b.input_gate(Q(1,3),Q(1,7),Q(1,9),2**32)['auxiliary_B5_comparison'],True)
    def test_gate_rejects_inexact_bool_negative_or_small_G(self):
        for a in ((True,0,0,2**32),(0.,0,0,2**32),(-1,0,0,2**32),(0,0,-1,2**32),(0,0,0,2**16)):
            with self.assertRaises(ValueError):b.input_gate(*a)
    def test_gate_does_not_promote_actual_J_or_input(self):
        r=b.input_gate(0,0,0,2**32)
        for k in ('actual_R5_Z5_enclosed','actual_joint_target_J_verified','actual_compact_global_or_physical_pulse_accepted'):self.assertIs(r[k],False)
    def test_two_outward_precisions(self):self.assertEqual([v['interval_digits'] for v in self.r['universal']],[80,110]);self.assertTrue(all(all(v['checks'].values()) for v in self.r['universal']))
    def test_outward_checks_materialize_no_actual_inputs(self):self.assertTrue(all(v['actual_R5_Z5_T5_Lambda5_or_huge_parameters_materialized'] is False for v in self.r['universal']))
    def test_invalid_interval_precision_rejected(self):
        for d in (30,80.,True):
            with self.assertRaises(ValueError):b.universal_review(d)
    def test_reciprocal_zero_guard(self):
        with self.assertRaises(ValueError):b.reciprocal(b.jet(0))
    def test_implicit_zero_divisor_guard(self):
        with self.assertRaises(ValueError):b.fifth_implicit(lambda x:x,0)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'0657b89fddc9b290f204f3d0fd719238949ad61f');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_fixed135_record_manifest(self):h=self.record['inherited_provenance'];self.assertEqual(len(h['frozen_records']),135);self.assertEqual(h['fixed_manifest_sha256'],audit.MANIFEST_SHA)
    def test_all_historical_bindings(self):h=self.record['inherited_provenance'];self.assertIs(h['all_bindings_verified'],True);self.assertGreater(h['binding_count'],2401)
    def test_exact_frozen_evidence(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_all_unverified_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_not_false(self):
        r=deepcopy(self.record);r['flags']['actual_R5_Z5_T5_Lambda5_or_B5_enclosed']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_all_original_input_and_scope_protocol_fields_bound(self):
        for k in ('premises','inputs','signed_root','moving','primitive','cutoff','gate','native_boundary','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_B5_cannot_replace_J_obligation(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_frozen_angle_cannot_replace_full_sensitivity(self):
        r=deepcopy(self.record);r['reviews']['moving']['full_source5_x5_coefficient']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
