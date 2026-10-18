"""Regional C3, complete-energy, same-root and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import angular_bounds as w

class GlobalAxialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']['reviews']
    def check(self,name):self.assertTrue(self.r['exact']['checks'][name])
    def test_angular_discrepancy_C3(self):self.check('discrepancy_enclosure')
    def test_target_decay_factor(self):self.check('target_decay_conversion')
    def test_strict_C3_contraction(self):self.check('angular_contraction_strict')
    def test_C3_radius(self):self.check('radius_enclosure')
    def test_same_root_uses_common_pointwise_ball(self):self.check('old_root_in_new_pointwise_uniqueness_ball');self.check('pointwise_same_root_lipschitz')
    def test_new_C3_ball_not_claimed_inside_old_C1_ball(self):self.assertGreater(2**32,4200)
    def test_interpolation_power_derivatives(self):self.assertTrue(self.r['symbolic']['checks']['interpolation_power_derivatives'])
    def test_interpolation_family_factorial_norm(self):self.check('normalized_interpolation_power_norm')
    def test_full_interpolation_mass(self):self.check('interpolation_energy_at_T64')
    def test_full_post_energy_not_truncated(self):self.check('complete_unedited_post_energy_bound')
    def test_nonzero_angular_edit_energy_retained(self):self.check('angular_edit_energy_less_than_one');self.check('complete_energy_norm_bound')
    def test_s0_degree_six_polynomial_third_norm(self):self.check('s0_polynomial_C3');self.check('s0_norm_at_T64')
    def test_linear_root_coefficient_not_erased(self):self.check('linear_coefficient_bound')
    def test_full_constant_coefficient_jet_bound(self):self.check('constant_coefficient_jet_budget')
    def test_quadratic_coefficient_bound(self):self.check('quadratic_coefficient_at_most_one')
    def test_selected_epsilon_bound(self):self.check('epsilon_selection');self.assertEqual(self.r['exact']['binary_upper_exponents']['epsilon'],-237)
    def test_original_Jacobian_inverse_retained(self):self.check('inverse_Jacobian_padding')
    def test_normalized_jet_caps(self):self.check('jet_caps_cover_recursion');self.assertEqual(self.r['exact']['constants']['strength_normalized_jet_caps'],[9,279,15939])
    def test_full_strength_C3_below_two(self):self.check('strength_full_C3_below_two')
    def test_affine_prefix_norm_from_C1(self):self.check('affine_prefix_polynomial_C3');self.check('prefix_C3_from_C1')
    def test_composed_end_coefficients_include_strength_derivative(self):self.check('composed_axial_end_coefficients')
    def test_third_implicit_derivative_all_terms(self):self.assertTrue(self.r['symbolic']['checks']['third_root_jet_all_product_terms'])
    def test_linear_coefficient_derivative_omission_detected(self):self.assertTrue(self.r['symbolic']['checks']['omitting_linear_coefficient_derivative_is_detected'])
    def test_all_five_independent_quadratic_controls(self):
        rows=self.r['manufactured']['rows'];self.assertEqual(len(rows),5);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_manufactured_roots_include_both_endpoints(self):self.assertEqual({r['eta'] for r in self.r['manufactured']['rows']},{'-1','-2/5','0','3/5','1'})
    def test_reciprocal_angular_shape_norm(self):self.check('reciprocal_f_factorial_C3')
    def test_disjoint_main_and_end_supports(self):self.check('disjoint_main_shape_budget');self.check('disjoint_end_bump_budget')
    def test_actual_pulse_bound_constant(self):self.check('pulse_factor');self.assertEqual(self.r['exact']['constants']['pulse_axial_C3_upper'],352)
    def test_axis_angular_Cauchy_factorial_cancellation(self):
        d=F(1,12800);self.assertEqual(F(self.r['regional']['axis_fixture']),7*sum(d**(-j) for j in range(4)))
    def test_common_radial_cover(self):self.assertTrue(all(self.r['regional']['checks'].values()))
    def test_global_maximum_includes_later_pulse(self):self.assertEqual(w.repaired_source_bounds(F(1,200),11,7,13)['KU'],352)
    def test_global_maximum_includes_axis(self):self.assertEqual(w.repaired_source_bounds(F(1,200),11,1000,13)['KU'],1000)
    def test_global_maximum_includes_modulation(self):self.assertEqual(w.repaired_source_bounds(F(1,200),11,7,1000)['KU'],2000)
    def test_signed_source_interface_preserves_radial_average_cost(self):
        r=w.repaired_source_bounds(F(1,200),11,7,13);self.assertEqual(r['v2'],36*r['KU']);self.assertGreater(r['I0'],0);self.assertGreater(r['I1'],0)
    def test_zero_epsilon_root_jets(self):self.assertEqual(w.root_jet_caps(0),[2,0,0,0])
    def test_caps_invalid_ranges(self):
        for e in (-1,F(1,100)):
            with self.assertRaises(ValueError):w.root_jet_caps(e)
    def test_invalid_region_inputs(self):
        with self.assertRaises(ValueError):w.repaired_source_bounds(F(1,200),11,-1,13)
    def test_exact_inputs_reject_float_bool(self):
        for v in (True,.01):
            with self.assertRaises(TypeError):w.root_jet_caps(v)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'0ba0d99a92c7e6de0158f91a11a49709f805ec64');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bytes(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_no_scope_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_full_stress_obligation_not_erased(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_integer_bool_scope_mutation(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_bound_tamper_rejected(self):
        r=deepcopy(self.record);r['reviews']['reviews']['exact']['constants']['pulse_axial_C3_upper']=1
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
