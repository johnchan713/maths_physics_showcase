"""Full joint fourth loop/primitive, conditional Gamma field gate and remaining flux input."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import joint as b

class JointFourthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_native_fourth_pressure_not_old_Bn(self):self.assertEqual(self.r['exponents']['native_shear_pressure_J4_A_powers']['p2'],1547)
    def test_native_fourth_shear_inputs_below_G(self):self.assertTrue(self.r['exponents']['checks']['native_joint4_inputs_below_A1579'])
    def test_all_four_signed_root_orders(self):self.assertEqual([r['selected_power'] for r in self.r['exponents']['mu']],[80,240,400,560]);self.assertTrue(all(r['accepted'] for r in self.r['exponents']['mu']))
    def test_mu_two_color_counts_include_fourth448(self):self.assertEqual(self.r['exponents']['mu_two_color_counts'],[4,20,96,448])
    def test_inverse_three_color_counts_include_fourth2268(self):self.assertEqual(self.r['exponents']['theta_three_color_counts'],[6,45,324,2268])
    def test_fourth_inverse_sufficient_power76012(self):self.assertEqual(self.r['exponents']['theta'][-1]['sufficient_power'],76012)
    def test_all_four_inverse_phase_orders(self):self.assertEqual([r['selected_power'] for r in self.r['exponents']['theta']],[11000,33000,55000,77000]);self.assertTrue(all(r['accepted'] for r in self.r['exponents']['theta']))
    def test_t_fixed_angle_joint_fourth_4000(self):self.assertEqual(self.r['exponents']['exponents']['t_fixed_angle_J4'],4000)
    def test_full_moving_source_exponent350000(self):self.assertEqual(self.r['exponents']['exponents']['moving_sources_J4'],350000)
    def test_full_primitive_exponent400000(self):self.assertEqual(self.r['exponents']['exponents']['primitives_full_jointJ4'],400000)
    def test_three_variable_inverse_identity_all_orders(self):self.assertTrue(self.r['colored']['checks']['all_three_variable_inverse_coefficients_through_four_exact'])
    def test_positive_nonunit_inverse_Jacobian(self):self.assertEqual(self.r['colored']['positive_base_Jacobian'],'394/385')
    def test_external_phase_derivative_retains_inverse_Jacobian(self):self.assertEqual(self.r['colored']['inverse_phase_first'],'385/394')
    def test_independent_diagonal_all_four_coefficients(self):self.assertTrue(self.r['colored']['checks']['all_four_diagonal_coefficients_match_independent_series'])
    def test_mixed_y_eta_phase2_fourth_coefficient(self):self.assertEqual(self.r['colored']['mixed_y_eta_phase2_source_coefficient'],'2266162369625/14242045476336')
    def test_frozen_phase_loses_fourth_mixed_term(self):self.assertTrue(self.r['colored']['checks']['frozen_angle_wrong_fourth_zero'])
    def test_colored_l1_product_bound_not_assumed_equality(self):self.assertTrue(self.r['colored']['checks']['colored_l1_products_through_four_bounded']);self.assertTrue(self.r['colored']['checks']['cancellation_control_prevents_majorant_equality'])
    def test_phase_fixture_is_not_actual_profile(self):self.assertIs(self.r['colored']['manufactured_not_actual_phase'],True)
    def test_all_twenty_normalized_positive_phase_rows(self):self.assertEqual(self.r['primitive']['positive_phase_rows'],20);self.assertTrue(self.r['primitive']['checks']['all_normalized_phase_derivative_rows_exact'])
    def test_phase_free_and_positive_phase_operator_costs(self):self.assertEqual((self.r['primitive']['phase_free_operator_bound'],self.r['primitive']['positive_phase_operator_bound']),(2,1))
    def test_FULL_joint_primitive_bound_three_not_two(self):self.assertEqual(self.r['primitive']['full_joint_operator_bound'],3)
    def test_pure_slow_fourth_primitive_terms_are_not_dropped(self):self.assertTrue(self.r['primitive']['checks']['pure_slow_fourth_primitive_terms_present'])
    def test_exact_primitive_zero_mean(self):self.assertTrue(self.r['primitive']['checks']['primitive_zero_mean_exact'])
    def test_exact_variance_gives_leading_mean_one(self):self.assertTrue(self.r['leading_means']['checks']['exact_variance_gives_mean_one_mass'])
    def test_weighted_aL_density_constant(self):self.assertTrue(self.r['leading_means']['checks']['weighted_aL_density_is_constant'])
    def test_weighted_bL_density_keeps_axial_t(self):self.assertTrue(self.r['leading_means']['checks']['weighted_bL_density_keeps_axial_t'])
    def test_both_leading_sources_have_zero_phase_mean(self):self.assertTrue(self.r['leading_means']['checks']['A_source_zero_mean_after_exact_mass']);self.assertTrue(self.r['leading_means']['checks']['B_source_zero_mean_uses_bs_minus_a_ts'])
    def test_leading_means_not_complete_wave_corrections(self):self.assertIs(self.r['leading_means']['complete_physical_wave_covariance_or_mean_claim'],False)
    def test_conditional_Gamma_bound_is_2H_not_zero(self):self.assertEqual(self.r['exponents']['Gamma_A_and_B_bound'],'2*H')
    def test_same_frequency_field_gate_accepts_Gamma2H(self):
        f=self.r['auxiliary']['flags'];self.assertTrue(all(f[k] for k in ('known_below_half','Gamma_sufficient_comparison','full_majorant_below_budget','repair_J4_below_one')))
    def test_auxiliary_acceptance_not_actual_enclosure(self):self.assertIs(self.r['auxiliary']['flags']['actual_Gamma_enclosed'],False);self.assertIs(self.r['auxiliary']['flags']['actual_FINAL_global_C4_accepted'],False)
    def test_auxiliary_field_domain_guard(self):
        for a,k in ((1,600),(2**16,599),(True,600),(2**16,600.)):
            with self.assertRaises(ValueError):b.auxiliary_field(a,k)
    def test_two_outward_precisions_all_hierarchy_comparisons(self):self.assertEqual([r['interval_digits'] for r in self.r['universal']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['universal']))
    def test_Phi_full_integration_operator_not_phase_free_only(self):self.assertTrue(all(r['checks']['full_Phi_integral_operator_below_eight'] for r in self.r['universal']))
    def test_Gamma2H_gate_has_universal_headroom(self):self.assertTrue(all(r['checks']['Gamma_2H_sufficient_gate_has_headroom'] for r in self.r['universal']))
    def test_huge_actual_hierarchy_not_materialized(self):self.assertTrue(all(r['actual_G_H_N_Z_or_selected_profile_materialized'] is False for r in self.r['universal']))
    def test_low_or_inexact_interval_precision_rejected(self):
        for d in (40,80.,True):
            with self.assertRaises(ValueError):b.universal_review(d)
    def test_fast_radial_shear_and_b_exponential_retained(self):self.assertTrue(self.r['full_memory']['checks']['fast_radial_a_shear_exact']);self.assertTrue(self.r['full_memory']['checks']['exponential_denominator_b_shear_exact'])
    def test_all_five_matching_functions_restore_pressure(self):self.assertTrue(self.r['full_memory']['checks']['all_five_exact_moments_restore_both_pressure_coordinates'])
    def test_Cp_only_pressure_memory_detected(self):self.assertEqual(self.r['full_memory']['Cp_only_p2_at_equator'],'-X*eps/E');self.assertTrue(self.r['full_memory']['checks']['pressure_only_memory_is_not_erased'])
    def test_nonzero_U_DeltaW_memory_required(self):self.assertTrue(self.r['full_memory']['checks']['nonzero_U_DeltaW_memory_required'])
    def test_original_flux_fourth_uses_fifth_angular_memory(self):self.assertTrue(self.r['flux_boundary']['checks']['original_radial_flux_fourth_principal_coefficient_exact'])
    def test_fifth_memory_coefficient_not_zero_interior(self):self.assertTrue(self.r['flux_boundary']['checks']['fifth_angular_memory_coefficient_nonzero_at_equator']);self.assertTrue(self.r['flux_boundary']['checks']['no_endpoint_only_cancellation_in_interior'])
    def test_FINAL_fifth_moment_and_flux_bound_still_open(self):self.assertIs(self.r['flux_boundary']['actual_FINAL_fifth_angular_moment_or_full_flux_C4_enclosed'],False)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'4bca4deafaf306be4fbd87e4d33a99ef7bfa5707');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_fixed_manifest133_before_future_additions(self):self.assertEqual(len(audit.HISTORICAL_RECORDS),133);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),133);self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
    def test_all_inherited_source_bindings_verified(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreaterEqual(self.record['inherited_provenance']['binding_count'],2328)
    def test_frozen_evidence_exactly_reproduced(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_unconditional_or_global_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_Z_or_joint_target_J4_condition_verified']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_condition_phase_memory_field_and_boundary_protocol_preserved(self):
        for k in ('native','condition','scalar','moving','primitive','field','boundary','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_actual_target_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_primitive_operator_or_Gamma_bound_tampering_rejected(self):
        for path,key,value in (('primitive','full_joint_operator_bound',2),('exponents','Gamma_A_and_B_bound','0')):
            r=deepcopy(self.record);r['reviews'][path][key]=value
            with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
