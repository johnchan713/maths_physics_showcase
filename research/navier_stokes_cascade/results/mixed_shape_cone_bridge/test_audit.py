"""Mixed colored jets, complete state/pressure and conditional scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import bridge as b

class MixedBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_native_J3_ledger_below_unchanged_Bn(self):self.assertTrue(self.r['exponents']['checks']['native_mixed_inputs_fit_Bn'])
    def test_pressure_third_not_promoted_to_fourth(self):self.assertEqual(self.r['exponents']['native_exponents']['p2_J3'],26644)
    def test_signed_root_three_orders(self):self.assertEqual([r['selected_power'] for r in self.r['exponents']['mu']],[80,240,400]);self.assertTrue(all(r['accepted'] for r in self.r['exponents']['mu']))
    def test_inverse_phase_all_three_external_directions(self):self.assertEqual([r['selected_power'] for r in self.r['exponents']['theta']],[11000,33000,55000])
    def test_color_counts_not_one_direction_only(self):self.assertEqual(self.r['exponents']['color_counts'],[6,45,324])
    def test_full_primitive_phase_operator_is_three(self):self.assertEqual(self.r['exponents']['primitive_joint_operator_bound'],3)
    def test_full_primitive_budget(self):self.assertEqual(self.r['exponents']['primitive_joint3_exponent'],400000)
    def test_three_variable_inverse_identity(self):self.assertTrue(self.r['colored']['checks']['all_three_variable_inverse_coefficients_satisfy_identity'])
    def test_external_phase_derivative_not_frozen(self):self.assertTrue(self.r['colored']['checks']['external_phase_derivative_one'])
    def test_exact_mixed_y_eta_phase_coefficient(self):self.assertEqual(self.r['colored']['mixed_y_eta_phase_coefficient'],'19/504 - sqrt(3)/44')
    def test_frozen_phase_rejected(self):self.assertTrue(self.r['colored']['checks']['freezing_phase_would_give_wrong_zero'])
    def test_colored_l1_product_majorant(self):self.assertTrue(self.r['colored']['checks']['all_colored_product_degrees_bounded'])
    def test_cancellation_control_prevents_majorant_equality_assumption(self):self.assertTrue(self.r['colored']['checks']['cancellation_makes_majorant_an_inequality'])
    def test_shape_positive_majorant(self):r=b.shape_majorant(2**16,7,Q(3,2));self.assertLess(r['raw_product'],r['majorant'])
    def test_shape_width_and_native_powers(self):r=b.shape_majorant(2**16,1,2);self.assertTrue(r['checks']['width_six_kept']);self.assertTrue(r['checks']['native_power_exact'])
    def test_shape_majorant_invalid_inputs(self):
        for a,z,d in ((1,1,1),(2**16,0,1),(2**16,1,0),(2**16,1.,1),(2**16,1,True)):
            with self.assertRaises(ValueError):b.shape_majorant(a,z,d)
    def test_stronger_shape4_implies_mixed3_gap(self):self.assertTrue(self.r['exponents']['checks']['stronger_shape4_implies_mixed3_constant_gap'])
    def test_mixed_gate_does_not_imply_shape4(self):
        A=2**16;z=A**32-7*A**20;self.assertTrue(b.shape_gate(A,z,3)['auxiliary_sufficient_comparison']);self.assertFalse(b.shape_gate(A,z,4)['auxiliary_sufficient_comparison'])
    def test_shape_gate_not_actual_acceptance(self):
        for order in (3,4):
            r=b.shape_gate(2**16,0,order);self.assertTrue(r['auxiliary_sufficient_comparison']);self.assertIs(r['actual_Z_enclosed'],False);self.assertIs(r['actual_S_or_compact_cone_accepted'],False)
    def test_shape_gate_invalid_inputs(self):
        for a,z,o in ((1,0,3),(2**16,-1,3),(2**16,0,2),(2**16,0,3.),(2**16,True,3)):
            with self.assertRaises(ValueError):b.shape_gate(a,z,o)
    def test_fast_radial_a_term_retained(self):self.assertTrue(self.r['state_memory']['checks']['fast_radial_a_shear_exact'])
    def test_b_exponential_denominator_required(self):self.assertTrue(self.r['state_memory']['checks']['exponential_denominator_b_shear_exact']);self.assertTrue(self.r['state_memory']['checks']['dropping_exponential_changes_shear'])
    def test_full_five_row_exact_memory_restoration(self):self.assertTrue(self.r['state_memory']['checks']['all_five_exact_moments_restore_both_pressure_coordinates'])
    def test_Cp_only_defect_still_visible(self):self.assertTrue(self.r['state_memory']['checks']['pressure_only_memory_is_not_erased']);self.assertEqual(self.r['state_memory']['Cp_only_p2_at_equator'],'-X*eps/E')
    def test_nonzero_U_pressure_transport_not_dropped(self):self.assertTrue(self.r['state_memory']['checks']['nonzero_U_DeltaW_memory_required'])
    def test_original_repair_radius_and_swirl(self):self.assertTrue(self.r['auxiliary']['checks']['original_repair_radius_accepts']);self.assertTrue(self.r['auxiliary']['checks']['original_repair_preserves_swirl'])
    def test_original_raw_repair_state_below_Ccorr(self):self.assertTrue(self.r['auxiliary']['checks']['original_raw_repair_state_below_Ccorr'])
    def test_complete_state_error_not_just_modulation(self):self.assertTrue(self.r['auxiliary']['checks']['actual_raw_positive_errors_below_state_bound']);self.assertTrue(self.r['auxiliary']['checks']['full_state_bound_below_epsilon_half'])
    def test_inverse_lambda_five_row_target_kept(self):self.assertTrue(self.r['auxiliary']['checks']['complete_five_row_target_with_inverse_lambda_fits'])
    def test_same_N_contraction_and_cone_gap(self):self.assertTrue(self.r['auxiliary']['checks']['exact_C1_quadratic_contraction']);self.assertTrue(self.r['auxiliary']['checks']['cone_gap_loss_below_half'])
    def test_frequency_fixture_not_actual_H_N(self):self.assertIs(self.r['auxiliary']['auxiliary_exact_rational_not_actual_H_N_or_profile'],True)
    def test_auxiliary_domain_guard(self):
        for A,k in ((1,600),(2**16,599),(True,600)):
            with self.assertRaises(ValueError):b.auxiliary_acceptance(A,k)
    def test_outward_comparisons_both_precisions(self):self.assertEqual([r['interval_digits'] for r in self.r['universal']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['universal']))
    def test_huge_actual_scales_not_materialized(self):self.assertTrue(all(r['actual_G_H_N_or_cutoff_materialized'] is False for r in self.r['universal']))
    def test_low_interval_precision_rejected(self):
        with self.assertRaises(ValueError):b.universal_review(40)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'66ae6450fb32aad2720462f45ca29c47eb7cd80f');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_history_manifest_fixed_before_future_checkpoints(self):self.assertEqual(len(audit.HISTORICAL_RECORDS),130);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),130);self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
    def test_historical_source_bindings(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreater(self.record['inherited_provenance']['binding_count'],2242)
    def test_exact_frozen_evidence(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_unconditional_or_FINAL_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_cutoff_budget_or_M_or_Z_gate_verified']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_shape_moving_state_and_memory_protocol_cannot_be_erased(self):
        for k in ('native','mixed_shape','stronger_gate','moving','state','cone','memory','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_actual_shape_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
