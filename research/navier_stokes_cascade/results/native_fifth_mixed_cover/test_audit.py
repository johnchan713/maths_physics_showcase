"""Native fifth mixed cover, sixth radial step, complete pressure and scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
from math import factorial
import json,tempfile,unittest
from pathlib import Path
import sympy as sp
import audit
import mixed_bounds as b

class MixedFifthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_logistic_recurrence_all_seven_orders(self):self.assertEqual(len(self.r['step']['logistic_rows']),7);self.assertTrue(all(r['accepted'] for r in self.r['step']['logistic_rows']))
    def test_sixth_logistic_coefficient_norm9366(self):self.assertEqual(self.r['step']['logistic_rows'][-1]['coefficient_l1'],9366)
    def test_complete_sixth_Bell_coefficient_sum203(self):self.assertEqual(self.r['step']['Bell6_coefficient_sum'],203)
    def test_Bell_weight_six_not_six_factors_only(self):self.assertTrue(self.r['step']['checks']['Bell6_weighted_order_six']);self.assertTrue(self.r['step']['checks']['Bell6_factor_count_at_most_six'])
    def test_Bell_majorant_coefficients_positive(self):self.assertTrue(self.r['step']['checks']['Bell6_coefficients_positive'])
    def test_edge_derivative_constant10080(self):self.assertEqual(self.r['step']['maximum_edge_derivative_constant'],10080)
    def test_central_and_edge_two_outward_precisions(self):self.assertEqual([r['interval_digits'] for r in self.r['step_intervals']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['step_intervals']))
    def test_edge_geometric_n6_not_n5(self):self.assertTrue(all(r['checks']['edge_n6_geometric_ratio_below_half'] for r in self.r['step_intervals']))
    def test_sixth_step_not_actual_huge_profile(self):self.assertTrue(all(r['actual_huge_C_or_profile_materialized'] is False for r in self.r['step_intervals']))
    def test_precision_inputs_rejected(self):
        for v in (40,80.,True):
            for f in (b.step_comparisons,b.parameters):
                with self.assertRaises(ValueError):f(v)
    def test_relative_bump_sigma6_principal(self):self.assertTrue(self.r['annular']['checks']['relative_E_sigma6_principal_exact'])
    def test_additive_U_sigma6_principal(self):self.assertTrue(self.r['annular']['checks']['additive_U_sigma6_principal_exact'])
    def test_missing_sixth_step_changes_actual_bump_jet(self):self.assertTrue(self.r['annular']['checks']['dropping_sigma6_changes_fifth_radial_jet'])
    def test_no_invented_unit_mass_prefactor(self):self.assertTrue(self.r['annular']['checks']['incorrect_unit_mass_prefactor_rejected'])
    def test_earlier_width_point08_not_later_point1(self):self.assertTrue(self.r['annular']['checks']['width_point08_not_later_point1'])
    def test_radial_bump_cost_sum_fits_one_million(self):self.assertLess(sum(Q(25,2)**j/factorial(j) for j in range(6)),10**6)
    def test_no_sixth_angular_root_input_silently_requested(self):self.assertIs(self.r['annular']['sixth_angular_root_derivative_required'],False)
    def test_all_17_positive_exponent_rows_fit(self):self.assertEqual(len(self.r['exponents']['rows']),17);self.assertTrue(all(r['accepted'] for r in self.r['exponents']['rows']))
    def test_reference_pressure_input_NOT_actual_p2_ODE(self):self.assertIn('nsy=Sn/L-ns',audit.protocol()['reference'])
    def test_original_ns_conversion_full_source_ODE(self):self.assertTrue(self.r['reference_ns']['checks']['original_p2_to_ns_full_source_ODE_exact'])
    def test_ns_conversion_keeps_both_X_and_E(self):self.assertTrue(self.r['reference_ns']['checks']['dropping_X_in_reference_ns_conversion_changes_flow']);self.assertTrue(self.r['reference_ns']['checks']['dropping_E_in_reference_ns_conversion_changes_flow'])
    def test_reference_full_Pi_J4_bound2050(self):self.assertEqual(next(r['selected_C_power'] for r in self.r['exponents']['rows'] if r['name']=='reference_full_Pi_J4'),2050)
    def test_actual_U_prescription_keeps_two_cutoff_costs(self):self.assertEqual(next(r['raw_C_power'] for r in self.r['exponents']['rows'] if r['name']=='actual_prescribed_b_J4_two_cutoffs'),4298)
    def test_value_sensitive_E_not_exp_of_high_norm(self):self.assertEqual(next(r['raw_C_power'] for r in self.r['exponents']['rows'] if r['name']=='actual_continuation_E_J5_value_sensitive'),10241)
    def test_five_region_cover_summed(self):self.assertEqual(self.r['exponents']['region_C_powers'],[16384,16384,128,128,256]);self.assertTrue(self.r['exponents']['checks']['five_regional_norms_SUMMED_not_maximized'])
    def test_native_J5_A5_below_unchanged_Bn(self):self.assertEqual(self.r['exponents']['full_native_field_J5_A_power'],5);self.assertTrue(self.r['exponents']['checks']['new_native_fifth_still_below_unchanged_Bn'])
    def test_shear_a_and_bs_J4_A262(self):self.assertEqual((self.r['exponents']['shear_A_powers']['a'],self.r['exponents']['shear_A_powers']['bs']),(262,262))
    def test_inverse_a_J4_charges_four_products(self):self.assertEqual(self.r['exponents']['shear_A_powers']['inverse_a'],1054)
    def test_ts_and_vs_keep_axial_square(self):self.assertEqual((self.r['exponents']['shear_A_powers']['ts'],self.r['exponents']['shear_A_powers']['vs']),(1316,1579))
    def test_full_native_pressure_A1547_retained(self):self.assertEqual((self.r['exponents']['shear_A_powers']['p1'],self.r['exponents']['shear_A_powers']['p2']),(1547,1547))
    def test_positive_sum_power_not_maximum(self):self.assertEqual(b.add_power(17,17,3),18);self.assertGreater(b.add_power(17,17,3),17)
    def test_power_domain_guard(self):
        for v in ((),(-1,),(1.,),(True,)):
            with self.assertRaises(ValueError):b.add_power(*v)
    def test_prescription_missing_indices_independent_expansion(self):self.assertTrue(self.r['fixtures']['checks']['eight_missing_coefficients_match_independent_polynomial_expansion'])
    def test_activation_axial_taper_and_shear_blend_all_present(self):self.assertEqual([r['stage'] for r in self.r['fixtures']['rows']],['activation','axial_taper','shear_blend'])
    def test_all_prescribed_shear_total4_jets_exact(self):self.assertTrue(self.r['fixtures']['checks']['all_prescribed_shear_total4_jets_exact'])
    def test_nonzero_fifth_radial_E_and_U_before_stop(self):self.assertTrue(self.r['fixtures']['checks']['fifth_radial_E_nonzero_in_every_stage']);self.assertTrue(self.r['fixtures']['checks']['fifth_radial_U_nonzero_before_axial_stop'])
    def test_stopped_axial_slope_has_zero_radial_jet(self):self.assertTrue(self.r['fixtures']['checks']['stopped_U_has_zero_positive_radial_jets'])
    def test_fixture_not_actual_profile(self):self.assertIs(self.r['fixtures']['manufactured_not_actual_parameter_profile'],True)
    def test_original_full_pressure_identities(self):self.assertTrue(self.r['complete_pressure_identity']['checks']['p1_original_forward_density_ODE_exact']);self.assertTrue(self.r['complete_pressure_identity']['checks']['p2_original_forward_density_ODE_exact'])
    def test_axis_and_Cp_cannot_be_erased(self):self.assertTrue(self.r['complete_pressure_identity']['checks']['complete_Pi_has_axis_and_Cp']);self.assertTrue(self.r['complete_pressure_identity']['checks']['dropping_axis_source_changes_original_ODE']);self.assertIs(self.r['complete_pressure_identity']['actual_incoming_axis_or_moments_set_zero'],False)
    def test_actual_parameter_two_precisions_all_fit(self):self.assertEqual([r['interval_digits'] for r in self.r['parameters']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['parameters']))
    def test_actual_parameters_not_materialized(self):self.assertTrue(all(r['actual_C_G_H_N_Z_or_Gamma_materialized'] is False for r in self.r['parameters']))
    def test_joint_target_coefficient168_billion_and_power13685(self):self.assertEqual((self.r['target']['coefficient'],self.r['target']['power']),(168000000000,13685))
    def test_exact_auxiliary_target_majorant(self):self.assertTrue(b.target_majorant(2**16,3,Q(7,2))['raw_below_majorant'])
    def test_target_auxiliary_acceptance_not_actual_Z_or_Gamma(self):
        r=b.target_majorant(2**16,1,1);self.assertIs(r['actual_Z_enclosed'],False);self.assertIs(r['actual_Gamma_or_FINAL_field_accepted'],False)
    def test_target_inexact_or_outside_domain_rejected(self):
        for a,z,d in ((1,1,1),(2**16,0,1),(2**16,1,0),(2**16,1.,1),(True,1,1)):
            with self.assertRaises(ValueError):b.target_majorant(a,z,d)
    def test_same_explicit_Z_gate_has_positive6988logA_gap(self):self.assertEqual(20673-13685,6988);self.assertTrue(all(r['checks']['old_explicit_Z_gate_implies_new_joint_target_gate'] for r in self.r['parameters']))
    def test_bare_angular_bound_does_not_imply_joint_bound(self):
        y,e,n=sp.symbols('y eta n',positive=True)
        self.assertEqual(sp.diff(sp.cos(n*y),e,4),0)
        self.assertEqual(sp.diff(sp.cos(n*y),y,4).subs(y,0),n**4)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'1548fb8b06b5855585d3568980393a973ccc8e95');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_fixed_manifest_132_historical_records(self):self.assertEqual(len(audit.HISTORICAL_RECORDS),132);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),132);self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
    def test_history_all_original_source_bindings(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreaterEqual(self.record['inherited_provenance']['binding_count'],2295)
    def test_frozen_evidence_exact_replay(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_cutoff_Gamma_FINAL_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_raw_variance_cutoff_Z_budget_verified']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_step_reference_memory_and_cover_protocol_not_erasable(self):
        for k in ('native','step','reference','actual','shear','target','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_Gamma_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_tampered_field_power_or_omitted_sixth_step_rejected(self):
        for obj,key,value in (('exponents','full_native_field_J5_A_power',0),('annular','sigma6_coefficient','0')):
            r=deepcopy(self.record);r['reviews'][obj][key]=value
            with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
