"""Third activated/core/source jets, affine coefficient bounds and strict scope."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import third_bounds as b

class ThirdCollarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def leading(self):return b.leading_third(2,F(1,8),4,F(1,100),F(1,10000))
    def test_actual_reference_third_envelope(self):self.assertEqual(self.leading()['reference_third_jet'],130375)
    def test_actual_collar_ordering_and_regular_lower_radius(self):
        v=self.leading();self.assertEqual(v['Xa'],1);self.assertEqual(v['xi_lower'],1);self.assertGreater(v['xi_upper']**2,v['X_upper'])
    def test_reference_third_core_majorant(self):self.assertTrue(self.r['symbolic']['checks']['Phi_third_majorant_below_S']);self.assertTrue(self.r['symbolic']['checks']['core_third_field_factor_below_two'])
    def test_native_shear_second_radial_derivative_bound(self):self.assertTrue(self.r['symbolic']['checks']['log_shear_second_bound'])
    def test_complete_third_swirl_exponential_product(self):self.assertTrue(self.r['symbolic']['checks']['complete_swirl_third_product'])
    def test_actual_axial_third_keeps_both_cutoff_derivatives(self):self.assertTrue(self.r['symbolic']['checks']['axial_activation_third'])
    def test_Jyyy_has_both_width_losses(self):self.assertEqual(self.leading()['J_yyy'],659600000)
    def test_original_squared_flat_shape_second_bound(self):self.assertTrue(self.r['symbolic']['checks']['actual_sigma_second_from_previous_bound'])
    def test_true_width_refinement_cost(self):self.assertTrue(self.r['checks']['true_activation_width_cost_retained'])
    def test_exponential_is_separated(self):self.assertEqual(self.leading()['exponential_power'],1);self.assertEqual(self.r['separated_operator_power'],2)
    def test_giant_exponential_not_materialized(self):
        r=b.leading_third(10**40,F(1,8),4,F(1,100),F(1,10000));self.assertGreater(r['exponential_power'],10**30);self.assertLess(len(str(r['leading_mixed_third_factor'].numerator)),1000)
    def test_pair_product_derivative_cross_terms(self):self.assertEqual(b.Pair(2,3)*b.Pair(5,7),b.Pair(10,29))
    def test_pair_sum_power_and_constant_division(self):self.assertEqual((b.Pair(2,3)+1)**2/3,b.Pair(3,6))
    def test_pair_no_signed_cancellation_allowed(self):
        for v,d in ((-1,0),(1,-1)):
            with self.assertRaises(ValueError):b.Pair(v,d)
        with self.assertRaises(ValueError):b.Pair(1)/0
        with self.assertRaises(ValueError):b.Pair(1)**True
    def test_all_six_rows_and_sources_have_gradient_bounds(self):
        r=b.operator_pairs(4,2,100,F(1,200));self.assertEqual(len(r['A0']),6);self.assertEqual(len(r['A1']),6);self.assertEqual(len(r['source']),6);self.assertGreater(r['Qxi'],0)
    def test_constant_kinematic_rows_have_zero_xi_derivative(self):
        r=b.operator_pairs(4,2,100,F(1,200));self.assertEqual([v['xi'] for v in r['A0'][:3]],['0','0','0'])
    def test_pressure_source_gradient_is_not_zero(self):
        r=b.operator_pairs(4,2,100,F(1,200));self.assertGreater(F(r['source'][3]['xi']),0)
    def test_full_shifted_axial_source_gradient_is_not_zero(self):
        r=b.operator_pairs(4,2,100,F(1,200));self.assertGreater(F(r['full_shifted_axial']['xi']),0);self.assertGreater(F(r['source'][5]['xi']),0)
    def test_source_bounds_enclosed_by_common_C_and_Cxi(self):
        r=b.operator_pairs(4,2,100,F(1,200))
        for row in r['A0']+r['A1']+r['source']:self.assertLessEqual(F(row['value']),r['Q']);self.assertLessEqual(F(row['xi']),r['Qxi'])
    def test_all_18_independent_retained_system_gradient_diagnostics(self):
        rows=self.r['retained_system_diagnostic']['rows'];self.assertEqual(len(rows),18);self.assertTrue(all(all(v['checks'].values()) for v in rows))
    def test_pressure_column_and_full_sources_preserved(self):self.assertTrue(self.r['retained_system_diagnostic']['checks']['both_pressure_column_couplings_retained']);self.assertTrue(self.r['retained_system_diagnostic']['checks']['pressure_and_both_shifted_sources_nonzero'])
    def test_finite_differences_explicitly_not_interval_certificates(self):self.assertIn('not interval',self.r['retained_system_diagnostic']['derivative_method'])
    def test_operator_factor_scaling_encloses_degree_two_exponential(self):
        lo=b.operator_pairs(4,2,7,F(1,200));hi=b.operator_pairs(4,2,35,F(1,200));self.assertLessEqual(hi['Q'],25*lo['Q']);self.assertLessEqual(hi['Qxi'],25*lo['Qxi'])
    def test_third_X_conversion_from_regular_P_and_Q(self):self.assertTrue(self.r['symbolic']['checks']['third_X_from_regular_P_Q'])
    def test_inverse_radius_derivative_omission_detected(self):self.assertTrue(self.r['symbolic']['checks']['missing_inverse_radius_terms_detected'])
    def test_singular_diagonal_and_angular_source_terms_in_affine_Wx(self):
        r=b.coefficient_third(2,3,1,F(1,4));self.assertEqual(r['W_xi'],(21,2));self.assertEqual(r['W_xi_eta'],(368,8))
    def test_affine_second_W_derivative_charges_all_sources(self):
        r=b.coefficient_third(2,3,1,F(1,4));self.assertEqual(r['W_xixi'],(871,29))
    def test_affine_J3_exact_coefficients(self):
        r=b.coefficient_third(2,3,1,F(1,4));self.assertEqual(r['scalar_radial_jets'][3],(F(937,8),F(35,8)))
    def test_no_operator_gradient_can_be_silently_discarded(self):self.assertTrue(self.r['checks']['missing_operator_derivative_underestimates_J3'])
    def test_eta_strip_losses_remain(self):self.assertTrue(self.r['checks']['eta_strip_loss_needed']);self.assertEqual(b.coefficient_third(2,3,1,F(1,4))['C_eta1_loss'],9)
    def test_Picard_toy_majorant_has_justified_exponent(self):self.assertTrue(self.r['symbolic']['checks']['explicit_Picard_diagnostic_bound']);self.assertIn('2^4000',self.r['auxiliary_M_upper'])
    def test_third_cutoff_input_not_erased_in_join(self):self.assertTrue(self.r['checks']['interface_keeps_positive_unknown_third_cutoff_cost'])
    def test_invalid_ranges_rejected(self):
        for values in ((1,1,1,F(1,100),F(1,10000)),(2,0,1,F(1,100),F(1,10000)),(2,1,0,F(1,100),F(1,10000)),(2,1,1,F(1,10),F(1,10000)),(2,1,1,F(1,100),F(1,2000))):
            with self.assertRaises(ValueError):b.leading_third(*values)
        for values in ((0,0,1,F(1,4)),(1,-1,1,F(1,4)),(1,1,0,F(1,4)),(1,1,1,2)):
            with self.assertRaises(ValueError):b.coefficient_third(*values)
        for values in ((4,1,100,F(1,200)),(4,2,0,F(1,200)),(4,2,1,F(1,10))):
            with self.assertRaises(ValueError):b.operator_pairs(*values)
    def test_float_and_bool_inputs_rejected(self):
        for v in (True,.1):
            with self.assertRaises(TypeError):b.Pair(v)
            with self.assertRaises(TypeError):b.coefficient_third(v,3,1,F(1,4))
    def test_parent_protocol_and_record(self):self.assertEqual(audit.protocol()['parent_commit'],'2fcf54abdc2b3faf661d73e0feae79197f96c565');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()));audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_source_bytes_preserved(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_all_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_false_formal_flag(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_protocol_cannot_be_erased(self):
        r=deepcopy(self.record);r['protocol']['pressure']='erased'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_actual_cutoff_obligation_not_dispatched_by_J3(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_frozen_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
