"""Tame fifth moment/root budget, complete original maps and high-input scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import audit
import fifth as b

class FinalFifthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_all_seven_fifth_exponential_partitions(self):self.assertTrue(self.r['exponential']['checks']['all_seven_fifth_exponential_partitions'])
    def test_top_fifth_primitive_linear_not_exponentiated(self):self.assertTrue(self.r['exponential']['checks']['top_primitive_jet_enters_linearly'])
    def test_lower_a2_a3_partition_retained(self):self.assertTrue(self.r['exponential']['checks']['lower_a2_a3_partition_cannot_be_dropped'])
    def test_finite_exponential_sum_below_three(self):self.assertTrue(self.r['exponential']['checks']['known_taylor_sum_below_three'])
    def test_exact_full_fifth_product_split(self):self.assertTrue(self.r['product']['checks']['exact_full_fifth_product_split'])
    def test_no_product_of_two_unknown_fifth_jets(self):self.assertTrue(self.r['product']['checks']['no_product_of_two_unknown_fifth_jets'])
    def test_top_fifth_jets_multiply_only_values(self):self.assertTrue(self.r['product']['checks']['top_fields_multiply_only_value_jets'])
    def test_both_high_value_cross_terms_retained(self):self.assertTrue(self.r['product']['checks']['both_high_value_cross_terms_retained'])
    def test_all_five_density_principal_rows(self):self.assertEqual([r['row'] for r in self.r['density']['rows']],['M','I','J','S','Cp']);self.assertTrue(all(r['exact'] for r in self.r['density']['rows']))
    def test_energy_keeps_both_squares(self):self.assertTrue(self.r['density']['checks']['energy_row_keeps_both_U_square_and_E_square'])
    def test_Cp_keeps_E_square_and_not_U_square(self):self.assertTrue(self.r['density']['checks']['Cp_row_has_E_square_not_axial_square'])
    def test_fifth_normalizers_still_five_and25(self):self.assertEqual(self.r['normalizers']['normalizer_polynomial_C5'],['5','25'])
    def test_fifth_normalizer25_not20(self):self.assertTrue(self.r['normalizers']['checks']['inverse_e_square_C5_cost25_not20'])
    def test_inverse_lambda_M_J_transformation_kept(self):self.assertTrue(self.r['normalizers']['checks']['transformed_M_J_row_keeps_two_rows_and_inverse_lambda'])
    def test_exact_fifth_constant_Q_root_equation(self):self.assertTrue(self.r['root']['checks']['exact_fifth_constant_Q_repair_identity'])
    def test_root_fifth_top_multiplier_is_existing_Jacobian(self):self.assertTrue(self.r['root']['checks']['fifth_top_multiplier_is_pointwise_Jacobian'])
    def test_root_fifth_uses_four_lower_pairs(self):self.assertEqual(self.r['root']['lower_pairs'],4)
    def test_no_high_fifth_square_in_root_lower_sum(self):self.assertTrue(self.r['root']['checks']['high_fifth_square_not_in_lower_sum'])
    def test_independent_manufactured_fifth_root_recovered(self):self.assertTrue(self.r['root']['checks']['independent_manufactured_fifth_root_recovered'])
    def test_pointwise_inverse_bound2000_not_new_contraction(self):self.assertEqual(self.r['root']['Jacobian_inverse_bound'],2000)
    def test_known_root_quadratic_tail3p2e15(self):self.assertEqual(self.r['root']['known_quadratic_tail_coefficient'],3200000000000000)
    def test_auxiliary_raw_moments_and_normalized_target_fit(self):self.assertTrue(self.r['auxiliary']['checks']['modulation_moment_cost_below_A515']);self.assertTrue(self.r['auxiliary']['checks']['all_normalized_target_costs_below_A522'])
    def test_auxiliary_same_root_C5_cap(self):self.assertTrue(self.r['auxiliary']['checks']['same_pointwise_Jacobian_C5_cap'])
    def test_auxiliary_complete_field_and_moment_C5_caps(self):self.assertTrue(self.r['auxiliary']['checks']['complete_modulation_plus_repair_field_C5_cap']);self.assertTrue(self.r['auxiliary']['checks']['complete_five_moment_C5_difference_cap'])
    def test_lower_C4_stays_small_without_high_C5_smallness(self):self.assertTrue(self.r['auxiliary']['checks']['lower_C4_change_remains_below_one'])
    def test_auxiliary_actual_Lambda_not_enclosed(self):self.assertIs(self.r['auxiliary']['actual_Lambda5_enclosed'],False);self.assertIs(self.r['auxiliary']['manufactured_not_actual_profile_or_Lambda5'],True)
    def test_auxiliary_domain_guard(self):
        for a,h,l in ((1,600,2),(2**16,599,2),(2**16,600,-1),(True,600,2),(2**16,600,2.)):
            with self.assertRaises(ValueError):b.auxiliary_bounds(a,h,l)
    def test_moment_jointJ5_budget_affine(self):self.assertEqual(self.r['ledger']['moment_J5_form'],'H^260+A^1044*Lambda5/N')
    def test_FINAL_reciprocal_E_charges_four_products(self):self.assertEqual(self.r['ledger']['reciprocal_E_J4_H_power'],514)
    def test_complete_Qs_and_Ns_maps_both_present(self):self.assertEqual((self.r['ledger']['map_H_powers']['Qs'],self.r['ledger']['map_H_powers']['Ns']),(520,132))
    def test_complete_p1_and_p2_ledger(self):self.assertEqual((self.r['ledger']['map_H_powers']['p1'],self.r['ledger']['map_H_powers']['p2']),(526,650))
    def test_velocity_flux_not_only_fields(self):self.assertEqual(self.r['ledger']['map_H_powers']['V'],4)
    def test_all_logarithmic_maps_fit_common660(self):self.assertEqual(self.r['ledger']['common_logarithmic_map_H_power'],660)
    def test_ordinary_conversion_keeps_all_fourth_Stirling_terms(self):self.assertEqual(self.r['ordinary']['rows'][4]['normalized_Dy_to_DX_coefficients'],['0','-1/4','11/12','-3/2','1'])
    def test_ordinary_conversion_complete_sum19_over2(self):self.assertEqual(self.r['ordinary']['absolute_coefficient_sum'],'19/2')
    def test_sufficient_separate_budget1024_not_replacing_K(self):self.assertEqual(self.r['ledger']['sufficient_compact_ordinary_budget_H_power'],1024)
    def test_gate_accept_and_reject_auxiliary_supplied_fifth_inputs(self):
        a=2**16;h=2;n=h**32+1
        self.assertTrue(b.fifth_input_gate(a,h,n,0)['sufficient_fifth_input_comparison'])
        self.assertFalse(b.fifth_input_gate(a,h,n,1)['sufficient_fifth_input_comparison'])
    def test_gate_does_not_promote_actual_input_target_or_global(self):
        r=b.fifth_input_gate(2**16,2,2**32+1,0)
        for k in ('actual_Lambda5_enclosed','actual_joint_target_verified','actual_global_leading_C4_accepted'):self.assertIs(r[k],False)
    def test_gate_rejects_inexact_negative_or_reselected_frequency(self):
        for args in ((2**16,2,2**32,0),(2**16,2,2**33,0),(2**16,2,2**32+1,-1),(2**16,2,2**32+1,1.),(True,2,2**32+1,0)):
            with self.assertRaises(ValueError):b.fifth_input_gate(*args)
    def test_two_precision_outward_hierarchy_checks(self):self.assertEqual([r['interval_digits'] for r in self.r['universal']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['universal']))
    def test_actual_Lambda_Z_and_huge_scales_not_materialized(self):self.assertTrue(all(r['actual_Lambda5_Z_or_huge_parameters_materialized'] is False for r in self.r['universal']))
    def test_invalid_precision_rejected(self):
        for d in (40,80.,True):
            with self.assertRaises(ValueError):b.universal_review(d)
    def test_fifth_memory_countercontrol_linear_growth(self):self.assertTrue(self.r['counter']['checks']['fixed_N_fifth_memory_grows_linearly'])
    def test_full_flux_countercontrol_limit(self):self.assertEqual(self.r['counter']['normalized_flux_fourth_over_n_limit'],'-c/(24*N)');self.assertTrue(self.r['counter']['checks']['complete_flux_limit_keeps_all_lower_terms'])
    def test_countercontrol_not_actual_profile_failure(self):self.assertIs(self.r['counter']['manufactured_functional_family_not_actual_C9_or_cone_profile'],True);self.assertIs(self.r['counter']['actual_profile_failure_inferred'],False)
    def test_all_five_matching_functions_keep_pressure_memory(self):self.assertTrue(self.r['full_memory']['checks']['all_five_exact_moments_restore_both_pressure_coordinates'])
    def test_axis_Cp_and_nonzero_U_DeltaW_retained(self):self.assertTrue(self.r['full_memory']['full_actual_axis_datum_unchanged']);self.assertTrue(self.r['full_memory']['checks']['pressure_only_memory_is_not_erased']);self.assertTrue(self.r['full_memory']['checks']['nonzero_U_DeltaW_memory_required'])
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'827e2c1f9915d28a0246830b9992cdddb2d1a382');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_fixed134_record_history_manifest(self):self.assertEqual(len(audit.HISTORICAL_RECORDS),134);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),134);self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
    def test_all_inherited_source_bindings(self):self.assertIs(self.record['inherited_provenance']['all_bindings_verified'],True);self.assertGreaterEqual(self.record['inherited_provenance']['binding_count'],2365)
    def test_exact_frozen_evidence(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_actual_Lambda_cutoff_or_global_scope_promotion_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_fifth_primitive_Lambda5_or_sufficient_F5_gate_enclosed']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_all_premise_high_input_memory_and_map_protocol_retained(self):
        for k in ('premises','high','modulation','normalization','repair','moments','maps','control','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_actual_fifth_input_obligation_not_removable(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_omitting_high_moment_or_lower_conversion_terms_rejected(self):
        r=deepcopy(self.record);r['reviews']['ledger']['moment_J5_form']='H^260'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
