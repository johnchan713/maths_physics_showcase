"""Exact forward pressure kernels, conditional majorants and scope guards."""
from copy import deepcopy
from fractions import Fraction as Q
import json,tempfile,unittest
from pathlib import Path
import sympy as s
import audit
import volterra as v

class VolterraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def test_twelve_both_forward_identities_exact(self):self.assertTrue(self.r['kernels']['checks']['twelve_exact_both_forward_kernel_identities'])
    def test_three_distinct_radii(self):self.assertEqual({r['X'] for r in self.r['kernels']['rows']},{'3/2','2','3'})
    def test_poles_and_interior_angles(self):self.assertEqual({r['eta'] for r in self.r['kernels']['rows']},{'-1','0','1/2','1'})
    def test_axis_and_incoming_fifth_memory_nonzero(self):self.assertTrue(self.r['kernels']['checks']['nonzero_incoming_fifth_memory_not_erased'])
    def test_endpoint_D_eta_M5_survives(self):self.assertTrue(self.r['kernels']['checks']['both_angular_endpoints_keep_D_eta_M5'])
    def test_Cp_cannot_be_dropped(self):self.assertTrue(self.r['kernels']['checks']['complete_Cp_high_jet_is_needed_for_swirl_kernel'])
    def test_swirl_kernel_zero_only_at_current_radius(self):self.assertTrue(self.r['kernels']['checks']['local_swirl_kernel_vanishes_only_at_current_radius'])
    def test_cases_not_actual_selected_profile(self):self.assertTrue(self.r['kernels']['manufactured_not_selected_profile'])
    def test_exact_logarithmic_antiderivative(self):
        x,X=s.symbols('x X',positive=True);self.assertEqual(v.radial_integral(1/x,x,X),s.log(X))
    def test_exact_half_power_antiderivative(self):
        x,X=s.symbols('x X',positive=True);self.assertEqual(s.simplify(v.radial_integral(s.sqrt(2*x),x,X)-(2*s.sqrt(2)/3)*(X**s.Rational(3,2)-1)),0)
    def test_nonpolynomial_radial_input_rejected(self):
        x,X=s.symbols('x X',positive=True)
        with self.assertRaises(ValueError):v.radial_integral(s.exp(x),x,X)
    def test_all_four_lower_product_partitions(self):
        e=s.Symbol('e');f=1+2*e+3*e**2+5*e**3+7*e**4+11*e**5
        self.assertEqual(v.lower_product(f,f,e,0),2*2*7+2*3*5)
    def test_fifth_normalization_is_120(self):
        e=s.Symbol('e');self.assertEqual(v.jet5(7*e**5,e,0),7)
    def fixture(self,**change):
        d=dict(A=3,B=10,left=1,X=2,eta=Q(1,2),h=Q(1,200),entry=(1,2,3,4,5,6),FE=7,FU=8,YE=9,YU=10);d.update(change);return v.certificate(**d)
    def test_quadratic_kernel_cap(self):self.assertEqual(self.fixture()['kappa_cap'],Q(1,2))
    def test_zero_length_keeps_entry(self):
        q=self.fixture(left=2);self.assertEqual(q['T2'],q['incoming_T2']);self.assertGreater(q['T2'],0)
    def test_halved_length_quarter_kernel_cap(self):self.assertEqual(self.fixture(X=Q(3,2))['kappa_cap'],self.fixture()['kappa_cap']/4)
    def test_pressure_axis_bound_has_same_coefficient_as_Cp(self):
        a=self.fixture(entry=(0,0,0,0,1,0));b=self.fixture(entry=(0,0,0,0,0,1));self.assertEqual(a['T2'],b['T2'])
    def test_T1_affine_in_both_high_field_jets(self):
        z=self.fixture(FE=0,FU=0);a=self.fixture(FE=7,FU=0);b=self.fixture(FE=0,FU=8);q=self.fixture();self.assertEqual(q['T1'],a['T1']+b['T1']-z['T1'])
    def test_T2_affine_in_both_high_field_jets(self):
        z=self.fixture(FE=0,FU=0);a=self.fixture(FE=7,FU=0);b=self.fixture(FE=0,FU=8);q=self.fixture();self.assertEqual(q['T2'],a['T2']+b['T2']-z['T2'])
    def test_pole_still_has_axial_fifth_cost(self):self.assertGreater(self.fixture(eta=1,FU=8)['T2'],self.fixture(eta=1,FU=0)['T2'])
    def test_swirl_fifth_cost_vanishes_at_pole(self):self.assertEqual(self.fixture(eta=1,FE=7)['T2'],self.fixture(eta=1,FE=0)['T2'])
    def test_pressure_fourth_factor_ten_A(self):
        q=self.fixture();self.assertEqual(q['pressure2_C4'],10+v.remainders(3,10)['common']+30*q['T2'])
    def test_common_remainder_bounds_both_rows(self):
        for A in (1,2,5):
            for B in (1,3,10):
                q=v.remainders(A,B);self.assertLess(q['R1'],q['common']);self.assertLess(q['R2'],q['common'])
    def test_explicit_remainder_constants(self):
        q=v.remainders(1,1);self.assertEqual((q['W'],q['Q'],q['N'],q['R1'],q['R2'],q['common']),(21,156,91,2496,7280,10000))
    def test_original_forward_lower_remainders(self):self.assertTrue(self.r['lower_shear']['checks']['actual_original_forward_lower_remainder_controls'])
    def test_L_inverse_geometric_cost_below_two(self):self.assertTrue(self.r['lower_shear']['checks']['C4_L_inverse_norm_geometric_series_25_over23_below_two'])
    def test_high_mixed_shear_coefficient_linear(self):self.assertTrue(self.r['lower_shear']['checks']['mixed_shear_highest_jet_is_linear_with_same_E_value'])
    def test_shear_native_E_inverse_cost(self):
        a=self.fixture(YE=0);b=self.fixture(YE=1);self.assertEqual(b['shear_a_high4']-a['shear_a_high4'],6)
    def test_axial_mixed_shear_cost_separate(self):
        a=self.fixture(YU=0);b=self.fixture(YU=1);self.assertEqual(b['shear_b_high4']-a['shear_b_high4'],6);self.assertEqual(a['shear_a_high4'],b['shear_a_high4'])
    def test_native_J4_not_relabelled_J5(self):self.assertTrue(self.r['lower_shear']['checks']['original_J4_only_bounds_slow_y_jets_through_three_not_four'])
    def test_negative_high_inputs_rejected(self):
        for key in ('FE','FU','YE','YU'):
            with self.assertRaises(ValueError):self.fixture(**{key:-1})
    def test_invalid_geometry_rejected(self):
        for change in (dict(A=0),dict(B=0),dict(left=0),dict(left=3),dict(X=4),dict(eta=2),dict(h=0)):
            with self.assertRaises(ValueError):self.fixture(**change)
    def test_all_six_entry_bounds_required(self):
        with self.assertRaises(ValueError):self.fixture(entry=(0,0,0,0,0))
        with self.assertRaises(ValueError):self.fixture(entry=(0,0,0,0,0,-1))
    def test_float_bool_inputs_rejected(self):
        for x in (3.,True):
            with self.assertRaises(TypeError):v.remainders(x,10)
    def test_no_gate_acceptance_from_manufactured_bounds(self):self.assertFalse(self.fixture()['actual_gate_accepted']);self.assertFalse(self.fixture()['actual_native_high_jets_supplied'])
    def test_parent_protocol_frozen_record(self):self.assertEqual(audit.protocol()['parent_commit'],'e85d339f19d696a915a047a923cb0303f4824acd');audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bindings(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_unproved_scope_promotions_rejected(self):
        for k,vv in audit.FLAGS.items():
            if not vv:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_numeric_zero_not_boolean_false(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_pressure_memory_kernel_shear_and_scope_not_erased(self):
        for k in ('moments','T1','T2','memory','kernel','pressure','shear','scope'):
            r=deepcopy(self.record);r['protocol'][k]='erased'
            with self.assertRaises(ValueError):audit.validate(r)
    def test_native_high_obligation_not_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_output_and_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
