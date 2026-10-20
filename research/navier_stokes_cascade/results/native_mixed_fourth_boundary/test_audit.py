"""Original pressure identities, fourth jets, precise missing-input and scope checks."""
from copy import deepcopy
from fractions import Fraction as Q
from math import factorial
import json,tempfile,unittest
from pathlib import Path
import sympy as sp
import audit
import boundary as b

class NativeFourthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_original_p1_forward_density_identity(self):self.assertTrue(self.r['pressure']['checks']['p1_original_forward_density_ODE_exact'])
    def test_original_p2_forward_density_identity(self):self.assertTrue(self.r['pressure']['checks']['p2_original_forward_density_ODE_exact'])
    def test_full_axis_pressure_cannot_be_dropped(self):self.assertTrue(self.r['pressure']['checks']['dropping_axis_source_changes_original_ODE'])
    def test_wrong_p2_y_plus_p2_ODE_rejected(self):self.assertTrue(self.r['pressure']['checks']['wrong_p2_y_plus_p2_equation_rejected'])
    def test_energy_density_keeps_both_squares(self):self.assertTrue(self.r['pressure']['checks']['both_U_and_E_square_density_retained'])
    def test_full_pressure_contains_axis_and_Cp(self):self.assertTrue(self.r['pressure']['checks']['complete_Pi_has_axis_and_Cp'])
    def test_no_incoming_memory_reset(self):self.assertIs(self.r['pressure']['actual_incoming_axis_or_moments_set_zero'],False)
    def test_native_not_FINAL_pressure(self):self.assertIs(self.r['pressure']['native_not_FINAL_pressure'],True)
    def test_source_coefficients_34_and_60(self):self.assertEqual((self.r['ledger']['Sq_coefficient'],self.r['ledger']['Sn_coefficient']),(34,60))
    def test_pressure_y_coefficients_107_and_184(self):self.assertEqual((self.r['ledger']['p1_y_coefficient'],self.r['ledger']['p2_y_coefficient']),(107,184))
    def test_radial_mixed_pressure_power771(self):self.assertEqual(self.r['ledger']['radial_mixed_A_power'],771)
    def test_full_native_pressure_power1547_not_Bn(self):self.assertEqual(self.r['ledger']['full_native_J4_A_power'],1547);self.assertGreater(1547,256)
    def test_full_native_ledger_fits_Bn(self):self.assertTrue(self.r['ledger']['checks']['native_fields_pressure_inverse_and_coordinates_fit_Bn'])
    def test_inverse_L_with_positive_order_geometric_sum(self):self.assertTrue(self.r['ledger']['checks']['inverse_L_J3_below_two'])
    def test_radial_pressure_normalization_including_factor_r(self):
        y,e=sp.symbols('y eta');p=3+y**4*e+7*y**2*e**2+11*y*e**3
        for r in range(1,5):
            s=4-r
            lhs=sp.diff(p,y,r,e,s).subs({y:0,e:0})/(factorial(r)*factorial(s))
            rhs=sp.diff(sp.diff(p,y),y,r-1,e,s).subs({y:0,e:0})/(r*factorial(r-1)*factorial(s))
            self.assertEqual(lhs,rhs)
    def test_derivative_J4_to_J3_factor_four(self):
        self.assertEqual(max(a+1 for a in range(4) for c in range(4-a)),4)
    def test_five_native_shear_principal_rows(self):self.assertTrue(all(all(r['checks'].values()) for r in self.r['shear']['rows']))
    def test_shear_r0_uses_existing_one_four_field(self):self.assertEqual(self.r['shear']['rows'][0]['field_index'],[1,4])
    def test_eight_remaining_fifth_mixed_field_inputs(self):self.assertEqual(self.r['shear']['remaining_field_indices'],[[2,3],[3,2],[4,1],[5,0]])
    def test_missing_fifth_mixed_inputs_not_accepted(self):self.assertIs(self.r['shear']['actual_fifth_mixed_field_inputs_enclosed'],False)
    def test_all_five_exact_normalized_frequency_rows(self):self.assertTrue(all(r['exact'] for r in self.r['frequency']['rows']))
    def test_no_extra_binomial_in_normalized_chain(self):
        y,e,phi,N=sp.symbols('y eta phi N')
        p=y**2*e*phi
        actual=sp.diff(p.subs(phi,N*y)/N,y,3,e,1).subs({y:0,e:0})/6
        self.assertEqual(actual,1)
        self.assertNotEqual(actual,sp.binomial(3,1))
    def test_all_fourth_exponential_partitions(self):self.assertTrue(self.r['frequency']['checks']['full_fourth_exponential_partitions'])
    def test_lower_square_partition_cannot_be_omitted(self):self.assertTrue(self.r['frequency']['checks']['missing_lower_square_rejected'])
    def test_Gamma_enters_exponential_fourth_linearly(self):self.assertTrue(self.r['frequency']['checks']['Gamma_enters_linearly'])
    def test_exponential_constant233_over12_below20(self):self.assertEqual(self.r['frequency']['known_exponential_constant'],'233/12')
    def test_primitive_degree_geometric_sum_below_three(self):self.assertTrue(self.r['frequency']['checks']['scaled_primitive_J4_known_constant_fits3'])
    def test_original_width_bump_cost_charged(self):self.assertTrue(self.r['frequency']['checks']['width_point_one_bump_J4_below_Bn'])
    def test_auxiliary_field_gate_accepts_with_actual_flags_false(self):
        r=self.r['auxiliary']['accepted'];self.assertTrue(r['full_majorant_below_budget']);self.assertIs(r['actual_Gamma_enclosed'],False);self.assertIs(r['actual_FINAL_global_C4_accepted'],False)
    def test_large_Gamma_auxiliary_rejected(self):
        r=self.r['auxiliary']['rejected'];self.assertFalse(r['Gamma_sufficient_comparison']);self.assertFalse(r['full_majorant_below_budget'])
    def test_field_gate_keeps_same_frequency_sandwich(self):
        for n in (2**32,2**33):
            with self.assertRaises(ValueError):b.field_gate(2**16,2,n,1)
    def test_field_gate_rejects_inexact_or_negative_inputs(self):
        for args in ((True,2,2**32+1,1),(2**16,2.,2**32+1,1),(2**16,2,2**32+1,-1),(2**16,2,2**32+1,1.),(2**16,2,2**32+1,True)):
            with self.assertRaises(ValueError):b.field_gate(*args)
    def test_auxiliary_hierarchy_not_actual_materialization(self):self.assertIs(self.r['auxiliary']['not_actual_hierarchy'],True)
    def test_manufactured_composed_fourth_exact(self):self.assertTrue(self.r['counter']['checks']['exact_normalized_fourth_composition'])
    def test_fixed_N_linear_limit_is_nonzero(self):self.assertEqual(self.r['counter']['normalized_fourth_limit'],'1/(24*N)');self.assertTrue(self.r['counter']['checks']['fixed_N_limit_linear_in_n'])
    def test_deleting_slow_fourth_loses_growth(self):self.assertTrue(self.r['counter']['checks']['omitting_pure_slow_fourth_loses_linear_term'])
    def test_four_independent_numerical_derivatives(self):self.assertEqual(len(self.r['counter']['numerical_rows']),4);self.assertTrue(all(r['checked'] for r in self.r['counter']['numerical_rows']))
    def test_control_not_actual_C9_or_failure_claim(self):self.assertIs(self.r['counter']['manufactured_functional_counter_not_actual_C9_or_global_pressure_profile'],True);self.assertIs(self.r['counter']['actual_profile_failure_inferred'],False)
    def test_two_outward_precisions_all_comparisons(self):self.assertEqual([r['interval_digits'] for r in self.r['universal']],[80,110]);self.assertTrue(all(all(r['checks'].values()) for r in self.r['universal']))
    def test_actual_huge_scales_not_materialized(self):self.assertTrue(all(r['actual_A_G_H_N_Gamma_materialized'] is False for r in self.r['universal']))
    def test_low_or_inexact_precision_rejected(self):
        for d in (40,80.,True):
            with self.assertRaises(ValueError):b.universal_review(d)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'dca2a2cdbadd9f3e0a56f3a49fc0d4abf6012441');self.assertEqual(json.loads((audit.HERE/'protocol.json').read_text()),audit.protocol())
    def test_history_fixed_parent_manifest(self):self.assertEqual(len(audit.HISTORICAL_RECORDS),131);self.assertEqual(len(set(audit.HISTORICAL_RECORDS)),131);self.assertNotIn(audit.HERE.name,audit.HISTORICAL_RECORDS)
    def test_historical_source_bindings(self):self.assertTrue(self.record['inherited_provenance']['all_bindings_verified']);self.assertGreaterEqual(self.record['inherited_provenance']['binding_count'],2264)
    def test_exact_frozen_evidence(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_unconditional_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_zero_cannot_replace_false(self):
        r=deepcopy(self.record);r['flags']['actual_cutoff_Z_or_Gamma_inequality_verified']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_protocol_memory_and_slow_input_erasure_rejected(self):
        for k in ('native','pressure','shear','frequency','field','control','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_next_obligation_removal_rejected(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_altered_pressure_bound_rejected(self):
        r=deepcopy(self.record);r['reviews']['ledger']['full_native_J4_A_power']=256
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(out)])

if __name__=='__main__':unittest.main()
