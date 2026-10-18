"""Exact envelope, retained-system, provenance and scope controls."""
from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
import envelope as w

class CollarEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews']
    def fixture(self):return w.bounds(2,F(1,8),4,F(1,100),F(1,10000))
    def test_actual_cutoff_is_squared_flat_step(self):self.assertIn('-1/s^2',self.r['symbolic']['cutoff'])
    def test_cutoff_derivative_bound(self):self.assertTrue(self.r['symbolic']['checks']['flat_squared_cutoff_derivative_64'])
    def test_cutoff_extremum_has_correct_sign(self):self.assertTrue(self.r['symbolic']['checks']['maximum_x3_exp_minus_x2_on_x_ge_2_at_2'])
    def test_core_reference_envelope(self):self.assertEqual(self.fixture()['reference_jet'],26075)
    def test_original_axis_radius(self):self.assertEqual(self.fixture()['Xa'],1)
    def test_actual_activated_endpoint_not_core_endpoint(self):self.assertGreater(self.fixture()['X_upper'],self.fixture()['Xa'])
    def test_xi_upper_encloses_sqrt_radius(self):
        r=self.fixture();self.assertGreater(r['xi_upper']**2,r['X_upper'])
    def test_log_shear_envelope(self):self.assertEqual(self.fixture()['log_shear'],100)
    def test_log_shear_second_derivative_envelope(self):self.assertEqual(self.fixture()['log_shear_y'],8000)
    def test_J_second_retains_width_cost(self):self.assertEqual(self.fixture()['J_yy'],324000)
    def test_second_product_rule_retains_Jy_squared(self):self.assertTrue(self.r['symbolic']['checks']['activated_swirl_second_product_rule'])
    def test_log_radial_chain_rule(self):self.assertTrue(self.r['symbolic']['checks']['independent_log_chain_rule'])
    def test_omitted_lower_log_derivative_detected(self):self.assertTrue(self.r['symbolic']['checks']['dropping_lower_log_derivative_changes_result'])
    def test_width_refinement_increases_envelope(self):self.assertTrue(self.r['checks']['width_cost_visible'])
    def test_no_division_by_flat_ea_or_axis_coordinate(self):
        r=self.fixture();self.assertGreater(r['radial_factor'],r['reference_jet']);self.assertEqual(r['Xa'],1)
    def test_parameter_strip_shrink(self):self.assertEqual(self.fixture()['delta'],F(1,12800))
    def test_second_parameter_Cauchy_factor_retained(self):
        r=self.fixture();self.assertEqual(r['angular_loss'],1+12800+2*12800**2)
    def test_regular_radial_average_jet_bound(self):self.assertTrue(self.r['symbolic']['checks']['regular_v_radial_bound'])
    def test_regular_v_eta_retains_reciprocal_L_derivative(self):self.assertTrue(self.r['symbolic']['checks']['regular_v_eta_bound'])
    def test_complex_L_has_no_pole(self):self.assertTrue(self.r['symbolic']['checks']['reciprocal_L_no_pole'])
    def test_full_source_quadratic_cost(self):
        a=w.primitive_operator_bounds(1,1,2);b=w.primitive_operator_bounds(1,1,4);self.assertEqual(b['w0'],4*a['w0'])
    def test_full_shifted_axial_second_envelope(self):self.assertEqual(w.primitive_operator_bounds(1,1,2)['axial_second'],1024)
    def test_A0_polynomial_domination(self):self.assertTrue(self.r['symbolic']['checks']['A0'])
    def test_A1_polynomial_domination(self):self.assertTrue(self.r['symbolic']['checks']['A1'])
    def test_full_source_polynomial_domination(self):self.assertTrue(self.r['symbolic']['checks']['source'])
    def test_every_radial_source_term_enclosed(self):self.assertTrue(self.r['symbolic']['checks']['all_radial_source_terms_fit_W'])
    def test_every_shifted_axial_product_term_enclosed(self):self.assertTrue(self.r['symbolic']['checks']['all_shifted_axial_product_terms_fit_Z'])
    def test_A0_row_five_all_terms(self):self.assertTrue(self.r['symbolic']['checks']['A0_row_five_with_all_terms'])
    def test_A0_row_six_pressure_term(self):self.assertTrue(self.r['symbolic']['checks']['A0_row_six_with_pressure'])
    def test_A1_row_five_all_terms(self):self.assertTrue(self.r['symbolic']['checks']['A1_row_five_with_angular_couplings'])
    def test_A1_row_six_pressure_coupling(self):self.assertTrue(self.r['symbolic']['checks']['A1_row_six_with_pressure_coupling'])
    def test_all_18_independent_retained_system_controls(self):
        rows=self.r['retained_system_diagnostic']['rows'];self.assertEqual(len(rows),18);self.assertTrue(all(all(r['checks'].values()) for r in rows))
    def test_system_axis_and_angular_endpoints_included(self):
        rows=self.r['retained_system_diagnostic']['rows'];self.assertEqual({r['eta'] for r in rows},{'-1.0','0.0','1.0'});self.assertIn('0.0',{r['xi'] for r in rows})
    def test_nonzero_pressure_source_retained(self):self.assertTrue(self.r['retained_system_diagnostic']['checks']['nonzero_full_pressure_source'])
    def test_nonzero_axial_source_retained(self):self.assertTrue(self.r['retained_system_diagnostic']['checks']['nonzero_shifted_axial_source'])
    def test_pressure_stays_separate_system_variable(self):self.assertTrue(self.r['retained_system_diagnostic']['checks']['pressure_is_a_separate_system_variable'])
    def test_exp_separated_representation(self):
        r=self.fixture();self.assertEqual(r['exponential_power'],1);self.assertEqual(r['operator_power'],2)
    def test_large_exponent_does_not_materialize_giant_power(self):
        r=w.bounds(10**40,F(1,8),4,F(1,100),F(1,10000));self.assertGreater(r['exponential_power'],10**30);self.assertLess(len(str(r['operator_factor'].numerator)),1000)
    def test_picard_strip_condition(self):
        r=self.fixture();self.assertLess(r['Delta']+r['d'],r['sigma']);self.assertEqual(r['Delta'],r['delta']/4)
    def test_bad_parameter_ranges_rejected(self):
        for values in ((1,1,1,F(1,100),F(1,10000)),(2,0,1,F(1,100),F(1,10000)),(2,1,0,F(1,100),F(1,10000)),(2,1,1,F(1,10),F(1,10000)),(2,1,1,F(1,100),F(1,2000))):
            with self.assertRaises(ValueError):w.bounds(*values)
    def test_exact_arithmetic_refuses_float_and_bool(self):
        for v in (2.,True):
            with self.assertRaises(TypeError):w.bounds(v,1,1,F(1,100),F(1,10000))
    def test_operator_primitive_invalid_ranges_rejected(self):
        for values in ((-1,1,1),(1,-1,1),(1,1,F(1,2))):
            with self.assertRaises(ValueError):w.primitive_operator_bounds(*values)
    def test_parent_and_protocol(self):self.assertEqual(audit.protocol()['parent_commit'],'1e910ccc6aa4da81eefb1d05e76e010943d0e91b');self.assertEqual(audit.protocol(),json.loads((audit.HERE/'protocol.json').read_text()))
    def test_frozen_record(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
    def test_historical_bytes_unchanged(self):self.assertTrue(all(self.record['inherited_provenance']['bindings'].values()))
    def test_no_scope_promotion(self):
        for k,v in audit.FLAGS.items():
            if not v:
                r=deepcopy(self.record);r['flags'][k]=True
                with self.assertRaises(ValueError):audit.validate(r)
    def test_complete_pressure_obligation_cannot_be_removed(self):
        r=deepcopy(self.record);r['remaining_obligations'].pop(0)
        with self.assertRaises(ValueError):audit.validate(r)
    def test_integer_boolean_promotion_rejected(self):
        r=deepcopy(self.record);r['flags']['formal_certificate']=0
        with self.assertRaises(ValueError):audit.validate(r)
    def test_constant_tampering_rejected(self):
        r=deepcopy(self.record);r['reviews']['exact_fixture']['operator_factor']='1'
        with self.assertRaises(ValueError):audit.validate(r)
    def test_fresh_replay_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(p)]);audit.validate(json.loads(p.read_text()))
            with self.assertRaises(SystemExit):audit.main(['--output',str(p)])

if __name__=='__main__':unittest.main()
