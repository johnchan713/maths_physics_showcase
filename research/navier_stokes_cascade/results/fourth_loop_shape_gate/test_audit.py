"""Fourth angular scalar/moving estimates with strict conditional shape guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import loop_fourth as f

class FourthLoopTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_mu_all_four_orders_have_headroom(self):self.assertTrue(all(v['accepted'] for v in self.r['exponents']['mu']))
    def test_mu_normalized_fourth_max_input_is480(self):self.assertEqual(self.r['exponents']['mu'][3]['largest_input_power'],480)
    def test_mu_fourth_selected560_not_old_third500(self):self.assertEqual(self.r['exponents']['mu'][3]['selected_power'],560)
    def test_mu_exact_Taylor_partition_counts(self):self.assertEqual([v['lower_composition_terms'] for v in self.r['exponents']['mu']],[2,5,12,28])
    def test_insufficient_mu_budget_is_rejected(self):self.assertTrue(any(not v['accepted'] for v in f.implicit_budget(64,10,[1,1,1,1])))
    def test_implicit_budget_invalid_inputs(self):
        for p,i,j in ((-1,1,[1]*4),(1.,1,[1]*4),(1,1,[1]*3),(1,1,[True]*4)):
            with self.assertRaises(ValueError):f.implicit_budget(p,i,j)
    def test_theta_all_four_orders_have_headroom(self):self.assertTrue(all(v['accepted'] for v in self.r['exponents']['theta']))
    def test_theta_normalized_fourth_max_input66000(self):self.assertEqual(self.r['exponents']['theta'][3]['largest_input_power'],66000)
    def test_theta_fourth_selected77000(self):self.assertEqual(self.r['exponents']['theta'][3]['selected_power'],77000)
    def test_q5_not_q4_budget(self):self.assertEqual(self.r['exponents']['constants']['q_joint5'],4);self.assertTrue(self.r['q5']['checks']['exact_fifth_q_coefficient'])
    def test_q5_expansion_has_all_three_terms(self):self.assertEqual(self.r['q5']['q5'],'h**5/120 - h**3/24 + 3*h/64')
    def test_R2_regular_fourth_coefficient(self):self.assertTrue(self.r['q5']['checks']['exact_removable_R2_zero_and_fourth'])
    def test_dropping_q5_changes_R2_fourth(self):self.assertTrue(self.r['q5']['checks']['omitting_q5_changes_fourth_R2'])
    def test_signed_Fmu_at_zero_nonzero(self):self.assertTrue(self.r['q5']['checks']['signed_Fmu_at_zero'])
    def test_Fp_zero_at_zero_variance(self):self.assertTrue(self.r['q5']['checks']['signed_Fp_at_zero'])
    def test_unsigned_Vmu_zero_not_inverted(self):self.assertTrue(self.r['q5']['checks']['unsigned_variance_mu_derivative_vanishes_at_zero'])
    def test_complete_variance_implies_mean_one_lift(self):self.assertTrue(self.r['q5']['checks']['exact_variance_implies_full_mean_one_lift'])
    def test_three_moving_centers_include_both_endpoints(self):self.assertEqual([r['eta'] for r in self.r['moving']['rows']],[-1,0,1])
    def test_all_four_inverse_coefficients_exact(self):self.assertTrue(all(v['implicit_residual_zero'] for v in self.r['moving']['rows']))
    def test_moving_fourth_135_over4096_not_explicit_map_control(self):self.assertEqual(self.r['moving']['rows'][1]['fourth'],'135/4096');self.assertNotEqual(self.r['moving']['rows'][1]['fourth'],'3/4096')
    def test_pole_fourth_controls_nonzero(self):self.assertEqual([v['fourth'] for v in self.r['moving']['rows']],['7/81','135/4096','114/15625'])
    def test_frozen_phase_fails_at_all_centers(self):self.assertTrue(all(v['frozen_phase_gives_wrong_zero'] for v in self.r['moving']['rows']))
    def test_positive_phase_inverse_at_both_poles(self):self.assertTrue(all(v['positive_w'] for v in self.r['moving']['rows']))
    def test_independent_nine_solved_root_comparison(self):self.assertTrue(all(v['independent_solved_root_interpolation'] and v['solved_values']==9 for v in self.r['root']['rows']))
    def test_root_tests_zero_pressure_and_zero_variance(self):self.assertEqual([v['case'] for v in self.r['root']['rows']],['zero_pressure','zero_variance'])
    def test_numerical_root_tolerance_not_frozen_roundoff(self):self.assertTrue(all(v['relative_error_limit']=='1e-18' and 'relative_error' not in v for v in self.r['root']['rows']))
    def test_numerical_root_not_interval_or_actual_profile(self):self.assertIs(self.r['root']['manufactured_floating_not_interval_or_actual_profile'],True)
    def test_universal_two_precisions(self):self.assertEqual([v['interval_digits'] for v in self.r['universal']],[80,110])
    def test_polynomial_absorption_all_larger_G(self):self.assertTrue(all(v['checks']['polynomial_absorption_at_minimum'] and v['checks']['absorption_log_ratio_decreases'] for v in self.r['universal']))
    def test_native_A10000_below_G_all_larger_A(self):self.assertTrue(all(v['checks']['A10000_below_exp_A32'] and v['checks']['A_power_comparison_ratio_decreases'] for v in self.r['universal']))
    def test_primitive_exp400000G_below_H(self):self.assertTrue(all(v['checks']['primitive_400000G_below_logH_G4'] for v in self.r['universal']))
    def test_all_fixed_angle_and_moving_composition_costs(self):self.assertTrue(self.r['exponents']['checks']['fixed_angle_removable_t_composition_fits4000']);self.assertTrue(self.r['exponents']['checks']['complete_fourth_moving_composition_fits350000'])
    def test_no_actual_hierarchy_materialization(self):self.assertTrue(all(v['enormous_actual_G_H_N_materialized'] is False for v in self.r['universal']))
    def test_auxiliary_universal_inputs_checked(self):
        for A,G,d in ((2,2**16,80),(2**16,2,80),(2**16,2**16,40)):
            with self.assertRaises(ValueError):f.universal_gate(A,G,d)
    def test_original_five_normalizers_bounded(self):self.assertTrue(self.r['target']['checks']['all_five_original_C4_normalizers_fit_A5'])
    def test_original_complete_target_affines_bounded(self):self.assertTrue(self.r['target']['checks']['all_three_actual_target_affine_coefficients_fit_A522'])
    def test_inverse_lambda_row_kept(self):self.assertTrue(self.r['target']['checks']['inverse_lambda_transformed_row_is_retained'])
    def test_original_conditional_compact_gate_and_majorant_agree(self):self.assertTrue(self.r['target']['checks']['original_substituted_compact_inequality_accepts_auxiliary_bound']);self.assertTrue(self.r['target']['checks']['original_positive_compact_majorant_also_accepts'])
    def test_target_gate_auxiliary_not_actual_HN(self):self.assertIs(self.r['target']['auxiliary_algebraic_parameters_not_actual_H_N'],True)
    def test_small_Z_shape_gate_is_only_auxiliary(self):
        r=f.shape_gate(2**16,10**100);self.assertTrue(r['sufficient']);self.assertIs(r['actual_shape_budget_supplied'],False);self.assertIs(r['actual_compact_gate_accepted'],False)
    def test_shape_log_gate_can_fail(self):self.assertFalse(f.shape_gate(2**16,1,logZ=2**513)['sufficient'])
    def test_shape_gate_invalid_exact_inputs(self):
        for A,Z in ((2,1),(2**16,0),(2**16,1.),(True,1)):
            with self.assertRaises(ValueError):f.shape_gate(A,Z)
    def test_shape_gate_invalid_log_bound(self):
        with self.assertRaises(ValueError):f.shape_gate(2**16,1,logZ=-1)
    def test_parent_SHA_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'66ff8c49955a1b17e04a617d42681633cc19b51b');self.assertEqual(audit.protocol()['manuscript_sha256'],'0e779481c4da40bd28d1e642e1d8ca57447d129610df28dfa5a11e9af8ae228f')
    def test_historical_source_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_record_exact_reproduction(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_actual_shape_and_FINAL_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_not_boolean_false(self):
        r=deepcopy(self.record);r['flags']['actual_variance_cutoff_shape_enclosed_or_s_C4_condition_verified']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_unproved_condition_cannot_be_erased(self):
        for k in ('domain','unproved_condition','scalar','implicit','moving','compact','memory','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_actual_shape_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_evidence_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
