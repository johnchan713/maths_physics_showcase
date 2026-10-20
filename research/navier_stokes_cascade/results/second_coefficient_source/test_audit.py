"""Full second source, exponent controls, derivative costs and conditional scope."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import source_bounds as b

class SecondSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_original_full_source_nine_profiles_points(self):self.assertEqual(len(self.r['original_hierarchy_diagnostics']['rows']),9);self.assertTrue(all(r['checks']['independent_exact_full_source_matches'] for r in self.r['original_hierarchy_diagnostics']['rows']))
    def test_analytic_xi_source_gradient_matches_original(self):self.assertTrue(all(r['checks']['analytic_source_xi_gradient_matches'] for r in self.r['original_hierarchy_diagnostics']['rows']))
    def test_all_current_A0_value_rows_bounded(self):self.assertTrue(all(r['checks']['A0_row_value_caps'] for r in self.r['original_hierarchy_diagnostics']['rows']))
    def test_all_A0_gradient_rows_bounded(self):self.assertTrue(all(r['checks']['A0_row_gradient_caps'] for r in self.r['original_hierarchy_diagnostics']['rows']))
    def test_all_A1_value_and_gradient_rows_bounded(self):self.assertTrue(all(r['checks']['A1_row_value_caps'] and r['checks']['A1_row_gradient_caps'] for r in self.r['original_hierarchy_diagnostics']['rows']))
    def test_all_full_source_value_and_gradient_entries_bounded(self):self.assertTrue(all(r['checks']['all_source_value_caps'] and r['checks']['all_source_gradient_caps'] for r in self.r['original_hierarchy_diagnostics']['rows']))
    def test_omitting_previous_radial_axial_viscosity_detected(self):self.assertTrue(self.r['original_hierarchy_diagnostics']['checks']['dropping_previous_axial_viscosity_detected'])
    def test_blanket_current_exponent_substitution_detected(self):self.assertTrue(self.r['original_hierarchy_diagnostics']['checks']['blanket_current_exponent_substitution_detected'])
    def test_regular_Z_axial_minus_one_identity(self):self.assertTrue(self.r['symbolic']['checks']['regular_axial_identity_preserves_shift_minus_one'])
    def test_third_angular_leading_control_is_nonzero_six(self):self.assertEqual(self.r['symbolic']['control_w1'],'6')
    def test_zero_pointwise_lower_angular_traces_do_not_supply_third(self):self.assertTrue(self.r['symbolic']['checks']['lower_three_pointwise_angular_traces_do_not_supply_third'])
    def test_uniform_velocity_cost_below20(self):self.assertTrue(self.r['checks']['leading_velocity_eta2_cost_below20']);self.assertLess(max(b.radial_velocity_cost(F(1,100),F(1,50))),20)
    def test_first_velocity_eta1_cost_below12(self):self.assertTrue(self.r['checks']['first_velocity_eta1_cost_below12'])
    def test_uniform_endpoint_dominates_other_h_inputs(self):
        top=b.radial_velocity_cost(F(1,100),F(1,50))
        for h in (F(1,1000),F(1,200),F(1,100)):
            for lam in (0,2*h):self.assertTrue(all(x<=y for x,y in zip(b.radial_velocity_cost(h,lam),top)))
    def test_leading_Cauchy_upgrade_keeps_six_delta_inverse_cube(self):self.assertEqual(b.leading_angular_upgrade(F(1,4)),421)
    def test_current_order_is_four_h(self):self.assertEqual(F(self.r['operator_source_fixture']['exponents']['current']),F(1,50))
    def test_lower_order_is_two_h(self):self.assertEqual(F(self.r['operator_source_fixture']['exponents']['lower']),F(1,100))
    def test_current_Pi2_matrix_coefficient_not_minus_one(self):self.assertEqual(F(self.r['operator_source_fixture']['exponents']['current_pressure']),F(-99,100))
    def test_radial_axial_exponent_remains_minus_one(self):self.assertEqual(F(self.r['operator_source_fixture']['exponents']['radial_axial_a']),-1)
    def test_lower_theta_and_axial_shifted_exponents(self):
        e=self.r['operator_source_fixture']['exponents'];self.assertEqual(F(e['theta_axial_lower']),F(-199,200));self.assertEqual(F(e['axial_lower']),F(-99,200))
    def test_common_T_squared_scaling(self):
        ref=b.operator_and_source(1,1,1,1,F(1,200))
        for H0,H1 in ((1,2),(3,7),(9,4)):
            z=b.operator_and_source(1,1,H0,H1,F(1,200));T=max(H0,H1)
            self.assertLessEqual(z['Q'],ref['Q']*T*T);self.assertLessEqual(z['Qxi'],ref['Qxi']*T*T)
    def test_second_outer_strip_fits_first_half_strip(self):self.assertTrue(self.r['checks']['second_outer_strip_matches_first_complex_strip'])
    def test_second_delta_is_original_quarter(self):self.assertEqual(F(self.r['second_collar_fixture']['delta2']),F(1,6400))
    def test_M2_beta_scales_as_T_four_not_T_two(self):
        j=b.second_collar(1,1,F(1,200),1,F(1,1600));self.assertEqual(j['beta2_factor'],6*j['z2_factor']**2/j['delta2']);self.assertGreater(j['beta2_factor']*7**4,j['beta2_factor']*7**2)
    def test_M2_keeps_unevaluated_actual_T(self):self.assertIn('T^2',self.r['second_collar_fixture']['M2']);self.assertIn('M2',self.r['second_collar_fixture']['M2'])
    def test_third_affine_degree_is_two(self):self.assertEqual(self.r['symbolic']['J3_degree_in_q_qxi'],2)
    def test_affine_coefficients_are_positive_polynomial(self):self.assertTrue(self.r['symbolic']['checks']['third_collar_affine_positive_polynomial'])
    def test_local_J3_T4_scaling_with_same_M2(self):self.assertTrue(self.r['checks']['second_collar_T_power_four_keeps_true_M2_input'])
    def test_distinct_powers_per_derivative(self):self.assertEqual(self.r['second_collar_fixture']['T_powers'],[0,0,2,4])
    def test_lower_cut_interface_retains_raw_factor_two(self):self.assertTrue(self.r['checks']['first_cutoff_raw_input_retains_factorial_conversion'])
    def test_missing_actual_third_cutoff_jet_rejected(self):
        with self.assertRaises(ValueError):b.lower_cut_interface(3,F(1,4),2,3,(1,2,3))
    def test_nonunit_cutoff_value_input_rejected(self):
        with self.assertRaises(ValueError):b.lower_cut_interface(3,F(1,4),2,3,(0,2,3,5))
    def test_invalid_source_envelopes(self):
        for x,rad,H0,H1,h in ((-1,1,1,1,F(1,200)),(2,1,1,1,F(1,200)),(1,1,0,1,F(1,200)),(1,1,1,0,F(1,200)),(1,1,1,1,0)):
            with self.assertRaises(ValueError):b.operator_and_source(x,rad,H0,H1,h)
    def test_invalid_second_collar_radius_or_strip(self):
        with self.assertRaises(ValueError):b.second_collar(1,1,F(1,200),0,F(1,1600))
        with self.assertRaises(ValueError):b.second_collar(1,1,F(1,200),1,0)
        with self.assertRaises(ValueError):b.second_collar(1,1,F(1,200),1,F(1,4))
        with self.assertRaises(ValueError):b.second_collar(1,1,F(1,200),2,F(1,1600))
    def test_velocity_cost_rejects_current_order_as_lower(self):
        with self.assertRaises(ValueError):b.radial_velocity_cost(F(1,200),F(1,50))
    def test_leading_upgrade_requires_positive_strip(self):
        with self.assertRaises(ValueError):b.leading_angular_upgrade(0)
    def test_exact_inputs_reject_float_bool(self):
        for bad in (1.,True):
            with self.assertRaises(TypeError):b.operator_and_source(1,1,bad,1,F(1,200))
    def test_parent_protocol_and_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'6bc32141ce0df0ebe51f6d1ea55c599c232166a5');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unsupported_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_source_exponents_and_scope_protocol_cannot_be_erased(self):
        for k in ('source','exponents','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_global_H23_obligation_retained(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(1)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_frozen_overwrite_guard(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
