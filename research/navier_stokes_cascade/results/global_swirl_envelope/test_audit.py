"""Same-root C3 heat target, exact terminal taper and global swirl controls."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import swirl_bounds as w

class SwirlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def exact(self,k):self.assertTrue(self.r['exact']['checks'][k])
    def symbolic(self,k):self.assertTrue(self.r['symbolic']['checks'][k])
    def test_Gamma_first_three_derivative_caps(self):self.exact('gamma_derivatives_through_three_bounded_by_2_4_32_h')
    def test_Gamma_zero_derivative_cap_one(self):self.assertEqual(w.gamma_derivative_cap(F(1,100),0),1)
    def test_actual_positive_moment_formula(self):self.assertEqual(w.gamma_derivative_cap(F(1,100),2),F(1,100)*F(101,100)**2*F(201,100))
    def test_Gamma_C3_difference_cap(self):self.exact('angular_Gamma_difference_C3_below_512_h_over_X')
    def test_squared_heat_factor_norm(self):self.exact('square_difference_budget_valid_at_X100')
    def test_full_C3_inverse_normalizers(self):self.symbolic('full_C3_inverse_amplitude_normalizations')
    def test_pressure_and_energy_target_rows(self):self.exact('normalized_heat_pressure_and_energy_rows_fit_8192')
    def test_infinite_angular_target_row(self):self.exact('normalized_heat_angular_row_retains_h_integral');self.symbolic('infinite_angular_tail_retains_h_cancellation')
    def test_C3_target_value(self):self.assertEqual(w.heat_root_budget(F(1,100),100)['target_C3'],F(8192,10**10))
    def test_strict_C3_contraction(self):self.exact('strict_C3_heat_root_contraction')
    def test_new_C3_radius(self):
        r=w.heat_root_budget(F(1,100),100);self.assertEqual(r['root_C3'],30*r['target_C3']);self.assertEqual(r['common_pointwise_radius'],F(2**18,10**10))
    def test_same_root_common_pointwise_ball(self):self.exact('old_and_new_roots_in_common_pointwise_ball');self.exact('common_pointwise_uniqueness_lipschitz_below_one')
    def test_larger_C3_ball_not_inside_old_C1_ball(self):self.exact('new_C3_radius_not_claimed_inside_old_C1_radius')
    def test_heat_cutoff_keeps_inverse_width(self):self.exact('heat_cutoff_J3_below_two_million');self.assertGreater(F(self.r['exact']['heat_cutoff_J3_upper']),1000000)
    def test_interpolation_theta_J3(self):self.exact('interpolation_theta_J3_below_two')
    def test_original_terminal_normalization_and_width_two(self):self.symbolic('original_terminal_derivative_has_half_width');self.exact('exact_terminal_two_unit_normalized_factor_J3_below_two')
    def test_scalar_profile_three_jet_formula(self):self.exact('scalar_slope_and_two_derivatives_fit_4096');self.assertGreater(F(self.r['exact']['scalar_J3_coefficient']),3000)
    def test_positive_heat_argument_mixed_J3(self):self.exact('positive_heat_argument_J3_below_twenty_over_X')
    def test_Gamma_joint_not_only_angular_bound(self):self.exact('joint_Gamma_difference_J3_below_64_h_over_X')
    def test_actual_heat_multiplier_norm(self):self.exact('actual_heat_multiplier_J3_below_two')
    def test_both_relative_corrections_small_in_J3(self):self.exact('both_relative_edit_J3_budgets_below_one');self.assertEqual(self.r['exact']['actual_relative_edit_binary_exponents'],dict(heat=-1616,angular=-7235))
    def test_full_reference_product_constant(self):self.exact('reference_product_fits_two60_Ahi');self.exact('actual_H_dominates_reference_prefactor')
    def test_Gamma_third_chain(self):self.symbolic('Gamma_third_angular_chain_includes_Z_second')
    def test_omitted_second_argument_term_detected(self):self.symbolic('missing_Gamma_second_argument_term_detected')
    def test_F_first_inverse_radius_conversion(self):self.symbolic('F_radial_first_conversion')
    def test_F_second_inverse_radius_conversion(self):self.symbolic('F_radial_second_conversion')
    def test_geometric_log_factor_norm(self):self.symbolic('geometric_log_jet_factorial_sum_below_two')
    def test_all_ten_positive_Gamma_controls(self):self.assertEqual(len(self.r['Gamma_diagnostic']['rows']),10);self.assertTrue(self.r['Gamma_diagnostic']['checks']['all_ten_positive_integral_controls'])
    def test_Gamma_poles_have_nonzero_eta_derivatives(self):self.assertTrue(self.r['Gamma_diagnostic']['checks']['both_poles_have_nonzero_first_eta_heat_jet'])
    def test_no_negative_Gamma_argument_samples(self):self.assertTrue(all(r['checks']['nonnegative_heat_argument'] for r in self.r['Gamma_diagnostic']['rows']))
    def test_Global_F_value_geometric_cost(self):self.assertEqual(w.swirl_inputs(1,16,7,101)['F_C3'],404)
    def test_Global_F_X_cost(self):self.assertEqual(w.swirl_inputs(F(1,2),16,7,101)['F_X_C2'],1212)
    def test_Global_F_XX_keeps_lower_log_derivative(self):self.assertEqual(w.swirl_inputs(F(1,2),16,7,101)['F_XX_C1'],7272)
    def test_axis_common_bound_not_divided_by_zero(self):
        r=w.swirl_inputs(1,16,10000,1);self.assertEqual(r['F_C3'],10000);self.assertEqual(r['F_X_C2'],10000);self.assertEqual(r['F_XX_C1'],10000)
    def test_complete_shifted_swirl_operator_coefficients(self):
        r=w.swirl_inputs(1,16,7,101);self.assertEqual(r['shifted_F_C1'],64*16**2*r['F_XX_C1']+704*16*r['F_X_C2']+1664*r['F_C3'])
    def test_complete_first_constant_positive(self):self.assertTrue(self.r['composition']['checks']['full_first_constant_positive'])
    def test_shifted_swirl_viscosity_omission_detected(self):self.assertTrue(self.r['composition']['checks']['shifted_swirl_viscosity_enters_theta_source'])
    def test_pressure_source_and_shifted_U_omissions_detected(self):self.assertTrue(self.r['composition']['checks']['leading_pointwise_source_and_shifted_U_not_erased'])
    def test_same_first_operator_interface(self):self.assertTrue(self.r['composition']['checks']['full_swirl_operator_same_coefficient_ledger'])
    def test_invalid_heat_h_or_radius(self):
        for h,X in ((0,100),(F(1,50),100),(F(1,100),99)):
            with self.assertRaises(ValueError):w.heat_root_budget(h,X)
    def test_invalid_swirl_geometry(self):
        for v in ((0,16,7,101),(2,1,7,101),(1,16,-1,101)):
            with self.assertRaises(ValueError):w.swirl_inputs(*v)
    def test_gamma_order_validation(self):
        for j in (-1,True,F(2)):
            with self.assertRaises(ValueError):w.gamma_derivative_cap(F(1,100),j)
    def test_exact_rational_inputs_only(self):
        for h in (True,.01):
            with self.assertRaises(TypeError):w.heat_root_budget(h,100)
    def test_parent_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'ca0e6a3aec5b491bfd0069418261a092d8fdd03e');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_historical_bytes(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_scope_promotions_rejected(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_scope_type_mutation_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_cutoff_obligation_not_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_root_budget_mutation_rejected(self):
        r=deepcopy(self.record);r['reviews']['exact']['worst_parameter_heat_root']['root_C3']='0'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
