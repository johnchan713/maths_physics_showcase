"""Nested domains, raw jets, source parity, cutoff dependence and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import json,tempfile,unittest
from pathlib import Path
import audit
import axis_bounds as b

class AxisSecondTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews'];cls.c=cls.r['nested_core_fixture']
    def test_nested_radial_source_disk(self):self.assertTrue(self.r['checks']['first_profile_radial_disks_remain_inside_five'])
    def test_second_core_Cauchy_disk(self):self.assertTrue(self.r['checks']['second_core_Cauchy_disks_remain_inside_nine_halves'])
    def test_first_eta_raw_disks(self):self.assertTrue(self.r['checks']['first_profile_eta_disks_remain_inside_first_solution'])
    def test_second_source_and_solution_strips(self):self.assertTrue(self.r['checks']['second_source_and_Picard_strips_nested'])
    def test_reference_four_delta_strip(self):self.assertTrue(self.r['checks']['leading_four_delta_inside_joint_reference_strip'])
    def test_reference_three_radial_jets(self):self.assertTrue(self.r['checks']['leading_reference_raw_radial_three_enclosed'])
    def test_raw_second_eta_factorial_retained(self):self.assertTrue(self.r['checks']['first_raw_second_eta_factor_includes_two_factorial'])
    def test_second_radial_Cauchy_factors(self):self.assertEqual(b.nested_core(2,1,3,F(1,200))['radial_factors'],(1,12,288,10368))
    def test_nested_delta_exact(self):self.assertEqual(self.c['delta'],'1/1600');self.assertEqual(self.c['delta2'],'1/6400')
    def test_axis_to_outer_disk_ratio(self):self.assertEqual(self.c['Xa'],'4/3');self.assertEqual(self.c['second_complex_radius'],'3/2')
    def test_first_H1_factor_is_raw_not_factorial_norm(self):
        c=b.nested_core(2,1,3,F(1,200));self.assertEqual(c['H1_M1core_factor'],6*13**3*(1+3200+8*1600**2))
    def test_true_nested_exponential_remains_symbolic(self):self.assertIn('Tcore^2',self.c['M2core']);self.assertIn('M1core',self.c['Tcore'])
    def test_source_and_beta_use_q_and_q_squared(self):
        c=b.nested_core(2,1,3,F(1,200));self.assertEqual(c['z2_factor'],2*c['xi_path_upper']*c['source_q']);self.assertEqual(c['beta2_factor'],6*c['z2_factor']**2/c['delta2'])
    def test_complex_xi_ray_enclosure(self):self.assertTrue(self.r['checks']['complex_xi_ray_length_bound'])
    def test_all_monomial_Cauchy_controls(self):self.assertTrue(self.r['checks']['analytic_monomials_obey_second_core_Cauchy'])
    def test_all_current_matrix_entries_reflect(self):self.assertTrue(self.r['symbolic']['checks']['current_full_matrix_reflects'])
    def test_all_derivative_matrix_entries_reflect(self):self.assertTrue(self.r['symbolic']['checks']['current_full_derivative_matrix_reflects'])
    def test_literal_current_matrix_pattern(self):self.assertEqual(len(self.r['symbolic']['nonzero_A0']),13);self.assertEqual(len(self.r['symbolic']['nonzero_A1']),6)
    def test_full_source_parity(self):self.assertTrue(self.r['symbolic']['checks']['full_second_source_reflects'])
    def test_current_pressure_exponent(self):self.assertEqual(self.r['symbolic']['current_pressure_exponent'],'2*h - 1')
    def test_pressure_columns_not_removed(self):self.assertTrue(self.r['symbolic']['checks']['current_pressure_and_derivative_couplings_retained'])
    def test_diagonal_preserves_parity(self):self.assertTrue(self.r['symbolic']['checks']['diagonal_preserves_parity'])
    def test_variable_diagonal_kernel(self):self.assertTrue(self.r['symbolic']['checks']['variable_diagonal_kernel_blocks_annihilate'])
    def test_derivative_of_right_sparse_block(self):self.assertTrue(self.r['symbolic']['checks']['differentiated_right_block_annihilates'])
    def test_odd_source_negative_control(self):self.assertTrue(self.r['symbolic']['checks']['odd_swirl_source_control_breaks_even_X_descent'])
    def test_nine_complex_full_source_samples(self):
        d=self.r['complex_source_diagnostics'];self.assertEqual(len(d['rows']),9);self.assertTrue(d['checks']['all_nine_complex_full_source_controls']);self.assertTrue(d['manufactured_not_selected_PDE']);self.assertTrue(d['not_interval_certification'])
    def test_exact_axis_pressure_after_removable_quotient_cancellation(self):
        d=self.r['complex_source_diagnostics'];self.assertEqual(d['complete_axis_pressure_slope'],'44/119');self.assertTrue(d['checks']['removable_axis_quotient_matches_original_hierarchy'])
    def test_full_regular_sources_have_no_hidden_radial_poles(self):self.assertTrue(self.r['complex_source_diagnostics']['checks']['all_full_sources_have_no_radial_denominator'])
    def test_core_positive_collar_overlap_all_scales(self):
        for L in (F(1,1000),F(1),F(4),F(1000)):
            j=b.joined_jets(2,1,L,F(1,200),5/L,1+5/L);self.assertTrue(j['overlap']);self.assertLessEqual(j['xi_positive_lower']**2,j['core']['Xa'])
    def test_core_and_collar_M_inputs_remain_distinct(self):self.assertTrue(self.r['checks']['separate_core_collar_source_exponentials_not_identified'])
    def test_common_smaller_eta_strip(self):self.assertEqual(self.r['inner_join_fixture']['common_smaller_eta_strip'],'1/12800')
    def test_join_core_Ceta1_factors(self):
        j=b.joined_jets(2,1,3,F(1,200),2,2);self.assertEqual(j['core_Ceta1_factors'],tuple(F(12801)*factorial(r)*12**r for r in range(4)))
    def test_join_collar_T_powers(self):self.assertEqual(self.r['inner_join_fixture']['collar_T_powers'],[0,0,2,4])
    def test_full_binomial_cut_product(self):
        rows=b.cut_rows((1,2,3,4),((5,6),(7,8),(9,10),(11,12)),2,(1,2,3,5))
        self.assertEqual(rows[3],(45,480,534))
    def test_cut_product_identity(self):
        rows=b.cut_rows((1,2,3,4),((5,6),(7,8),(9,10),(11,12)),2,(1,0,0,0));self.assertEqual(rows,((1,5,6),(2,7,8),(3,36,40),(4,176,192)))
    def test_third_cutoff_jet_not_closed_by_axis_bound(self):self.assertTrue(self.r['checks']['third_cutoff_cost_still_present'])
    def test_invalid_core_inputs(self):
        for S,rho,L,h in ((1,1,1,F(1,200)),(2,0,1,F(1,200)),(2,2,1,F(1,200)),(2,1,0,F(1,200)),(2,1,1,0),(2,1,1,F(1,50))):
            with self.assertRaises(ValueError):b.nested_core(S,rho,L,h)
    def test_float_bool_rejected(self):
        for bad in (2.,True):
            with self.assertRaises(TypeError):b.nested_core(bad,1,1,F(1,200))
    def test_short_collar_and_invalid_enclosing_b_rejected(self):
        for x,bb in ((1,2),(F(4,3),2),(2,1)):
            with self.assertRaises(ValueError):b.joined_jets(2,1,3,F(1,200),x,bb)
    def test_bad_cut_rows_rejected(self):
        for core,aff,T,kk in (((1,)*3,((1,2),)*4,1,(1,0,0,0)),((1,)*4,((1,2),)*4,0,(1,0,0,0)),((1,)*4,((1,2),)*4,1,(1,0,0,-1))):
            with self.assertRaises(ValueError):b.cut_rows(core,aff,T,kk)
    def test_protocol_and_frozen_evidence(self):self.assertEqual(audit.protocol()['parent_commit'],'f2acf72d1abd5d68653ad25cddbe58507aa90040');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_false_scope_promotion_rejected(self):
        for key,value in audit.FLAGS.items():
            if not value:
                r=deepcopy(self.record);r['flags'][key]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_boolean_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_nested_exponent_pressure_and_cutoff_not_erased(self):
        for key in ('exponential','source','core_identity','cutoff'):
            r=deepcopy(self.record);r['protocol'][key]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_remaining_global_H23_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(2)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_nonfinite_axis_diagnostic_rejected(self):
        r=deepcopy(self.record);r['reviews']['complex_source_diagnostics']['complete_axis_pressure_slope']='nan'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
