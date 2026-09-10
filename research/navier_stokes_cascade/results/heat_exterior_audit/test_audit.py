"""Heat compensation, whole-tail integration and reference-scope controls."""
from copy import deepcopy
import json
import unittest
import mpmath as mp
from bounds import heat_bound,HERE
from construction import TailMoments,Compensation,Rule,heat_polynomial,heat_integral,pole_control,relative
from exterior import Exterior
from audit import validate,compare,collar_limit_gap,GATES,STATUS


class BoundsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 160
        cls.b = heat_bound()

    def test_all_outward_bounds(self):
        self.assertTrue(all(v is True for v in self.b['checks'].values()))

    def test_matrix_uniform_inverse(self):
        self.assertLess(max(self.b['inverse_row_upper']),15)

    def test_positive_tiny_budgets(self):
        for key in ('absolute_w_error_upper','relative_Q_error_upper','coefficient_universal_upper'):
            self.assertTrue(0<self.b[key]<mp.mpf('1e-6000'))

    def test_finite_cone(self):
        self.assertGreater(self.b['finite_normalized_cone_gap_lower'],mp.mpf('.8'))
        self.assertEqual(self.b['a_minus_two_floor'],'h')

    def test_moments_and_scope(self):
        self.assertEqual(len(self.b['exact_moments']),5)
        self.assertIn('regular axis',self.b['scope'])
        self.assertIn('global-blowup-unverified',STATUS)

    def test_invalid_parameters(self):
        for kw in ({'md':3},{'digits':20},{'radius':99}):
            with self.assertRaises(ValueError): heat_bound(**kw)


class ConstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 160
        cls.h = mp.mpf('.005')
        cls.tail = TailMoments(cls.h,mp.log(10000),32,8,24)
        cls.matrix = Compensation('.001',Rule(32,8))
        cls.log_rh = 4*mp.log(cls.h)-14
        cls.row = cls.tail.target('.5',mp.log(100),cls.log_rh)
        cls.root = cls.matrix.newton(cls.row['target'])
        cls.graded = cls.matrix.graded(cls.row['target'])

    def test_actual_three_moments(self):
        residual = self.matrix.B*self.root+self.matrix.quadratic(self.root,self.root)-mp.matrix(self.row['target'])
        self.assertLess(max(abs(v) for v in residual)/self.graded['scale'],mp.mpf('1e-100'))

    def test_nonlinear_grade_and_remainder(self):
        surrogate = [mp.fsum(g[j] for g in self.graded['components']) for j in range(3)]
        gap = max(abs(v-w) for v,w in zip(surrogate,self.root))
        self.assertTrue(0<gap<self.graded['remainder'])
        self.assertGreater(max(abs(v) for v in self.graded['components'][1]),0)

    def test_independent_source_ODE(self):
        result = self.matrix.independent_moments(self.root)
        self.assertLess(max(result['gaps']),mp.mpf('1e-9'))
        self.assertGreater(result['quadratic_omission_upper'],0)

    def test_target_and_root_eta_derivatives(self):
        delta = mp.mpf('1e-20')
        plus = self.tail.target(mp.mpf('.5')+delta,mp.log(100),self.log_rh)
        minus = self.tail.target(mp.mpf('.5')-delta,mp.log(100),self.log_rh)
        exact = self.matrix.derivative(self.root,self.row['target_eta'])
        for p,m,e in zip(plus['target'],minus['target'],self.row['target_eta']):
            self.assertLess(relative((p-m)/(2*delta),e),mp.mpf('1e-35'))
        for p,m,e in zip(self.matrix.newton(plus['target']),self.matrix.newton(minus['target']),exact):
            self.assertLess(relative((p-m)/(2*delta),e),mp.mpf('1e-35'))

    def test_pole_zero_value_nonzero_derivative(self):
        pole = self.tail.target(1,mp.log(100),self.log_rh)
        self.assertTrue(all(v==0 for v in pole['target']))
        self.assertTrue(any(v!=0 for v in pole['target_eta']))

    def test_tiny_nonlinearity_is_not_zero(self):
        g = self.matrix.graded([mp.mpf('1e-1200'),mp.mpf('1e-1000'),mp.mpf('1e-410')])
        self.assertGreater(max(abs(v) for v in g['components'][1]),0)
        self.assertTrue(all(a+b==a for a,b in zip(g['components'][0],g['components'][1])))

    def test_positive_Taylor_target_budget(self):
        self.assertTrue(all(v>0 for v in self.row['taylor_target_C1_error']))

    def test_independent_Gamma_integral(self):
        q = heat_polynomial(self.h,'.001',24)
        reference = heat_integral(self.h,'.001',48)
        self.assertLess(abs(q['normalized_delta']-reference),q['remainder']/(self.h*mp.mpf('.001')))

    def test_polynomial_is_not_exact_heat_solution(self):
        q = heat_polynomial(self.h,'.001',24)
        self.assertNotEqual(q['ode_residual'],0)
        self.assertLess(relative(q['ode_direct'],q['ode_residual']),mp.mpf('1e-50'))

    def test_finite_h_pole_stress(self):
        with mp.workdps(320):
            h = mp.mpf('1e-220')
            p = pole_control(h)
            self.assertLess(relative(p['ps1_minus_two'],2*h),mp.mpf('1e-80'))
            self.assertGreater(p['ps1_minus_two'],0)
            self.assertLess(p['omitted_derivative_stress_over_F'],-1)


class ExteriorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 160
        cls.e = Exterior(order=32)

    def test_positive_stress_and_direction(self):
        for y in ('.5','1.5','2.5'):
            r = self.e.direct(y,'.5')
            self.assertGreater(r['theta_over_rho'],0)
            self.assertTrue(0<abs(r['ratio'])<r['ratio_upper'])
            self.assertGreater(r['a_minus_two'],self.e.h)
            self.assertGreater(r['directional_gap'],mp.mpf('1.99'))

    def test_independent_flat_integral(self):
        a,b = self.e.direct('2.5','.5'),self.e.flat('.5','.5')
        self.assertLess(relative(a['ratio'],b['ratio']),mp.mpf('1e-20'))

    def test_smooth_edge_direction(self):
        b = self.e.flat('.001','.5')
        self.assertGreater(b['btheta'],0)
        self.assertLess(relative(b['ratio_over_delta6'],b['expected_ratio_over_delta6']),mp.mpf('.002'))

    def test_flat_factor_survives_binary64_underflow(self):
        b = self.e.flat('.01','.5')
        self.assertGreater(b['positive_flat_factor'],0)
        self.assertEqual(b['binary64_flat_factor'],0)

    def test_zero_edge_uses_extended_direction(self):
        b = self.e.endpoint('.5')
        self.assertEqual((b['theta'],b['axial']),(0,0))
        self.assertEqual(b['direction'],[1,0])
        self.assertEqual(b['directional_gap'],2)

    def test_endpoint_ratio_rejected(self):
        with self.assertRaises(ValueError): self.e.direct(3,0)
        with self.assertRaises(ValueError): self.e.flat(0,0)

    def test_collar_selection_survives_precision_change(self):
        with mp.workdps(80):
            rows = [dict(delta=mp.mpf('.001'),ratio_over_delta6=mp.mpf(1),expected_ratio_over_delta6=mp.mpf(1))]
        with mp.workdps(160):
            self.assertNotEqual(rows[0]['delta'],mp.mpf('.001'))
            self.assertEqual(collar_limit_gap(rows,['.5'],1),0)
            with self.assertRaises(ValueError): collar_limit_gap([],['.5'],1)


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads((HERE/'evidence.json').read_text())

    def test_record_valid(self):
        validate(self.record)
        self.assertEqual(set(self.record['gates']),GATES)

    def test_full_proof_promotion_rejected(self):
        bad = deepcopy(self.record);bad['status']='blowup-proved'
        with self.assertRaises(ValueError): validate(bad)

    def test_missing_gate(self):
        bad = deepcopy(self.record);bad['gates'].pop('pole_omission_failure')
        with self.assertRaises(ValueError): validate(bad)

    def test_zeroed_small_bound(self):
        bad = deepcopy(self.record);bad['certificate']['absolute_w_error_upper']='0'
        with self.assertRaises(ValueError): validate(bad)

    def test_erased_quadratic_grade(self):
        bad = deepcopy(self.record);bad['selected']['rows'][0]['nonlinear_grade_norm']='0'
        with self.assertRaises(ValueError): validate(bad)

    def test_erased_pole_derivative(self):
        bad = deepcopy(self.record);bad['selected']['rows'][-1]['root_eta']=['0','0','0']
        with self.assertRaises(ValueError): validate(bad)

    def test_erased_ODE_cadence(self):
        bad = deepcopy(self.record);bad['independent_ODE'].pop()
        with self.assertRaises(ValueError): validate(bad)

    def test_tiny_record_comparison(self):
        with mp.workdps(120):
            compare('1e-8000','1.00000000000000000001e-8000','/physical')
            with self.assertRaises(ValueError): compare('1e-8000','0','/physical')


if __name__=='__main__':
    unittest.main()
