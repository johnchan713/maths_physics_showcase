#!/usr/bin/env python3
"""Check continuum enclosures, source algebra, and evidence integrity."""
import copy
from fractions import Fraction
import json
import unittest
import mpmath as mp
from bounds import axial_bound, derivative_upper, elementary_checks, step_at, lo, hi
from diagnostics import Moments, independent_ode, HERE, step, step_prime
from audit import validate, compare_record, relative


class IntervalTests(unittest.TestCase):
    def setUp(self):
        mp.mp.dps = 80
        mp.iv.dps = 50

    def test_step_point_enclosures(self):
        for x in [Fraction(0),Fraction(1,1000),Fraction(1,4),Fraction(1,2),Fraction(3,4),Fraction(999,1000),Fraction(1)]:
            exact = step(mp.mpf(x.numerator)/x.denominator)
            enclosure = step_at(x)
            self.assertLessEqual(lo(enclosure),exact)
            self.assertGreaterEqual(hi(enclosure),exact)

    def test_derivative_box_enclosures(self):
        for a,b in [(Fraction(0),Fraction(1,100)),(Fraction(1,5),Fraction(3,10)),(Fraction(49,100),Fraction(51,100)),(Fraction(9,10),Fraction(1))]:
            bound = derivative_upper(a,b)
            for i in range(21):
                x = a+(b-a)*Fraction(i,20)
                self.assertLessEqual(step_prime(mp.mpf(x.numerator)/x.denominator),bound)

    def test_step_derivative_independent_differentiation(self):
        for x in map(mp.mpf,['.03','.24','.43','.5','.61','.97']):
            self.assertLess(abs(mp.diff(step,x)-step_prime(x)),mp.mpf('1e-65'))

    def test_global_derivative_polynomial_identity(self):
        # Exact rational arithmetic checks the identity used in the proof.
        for i in range(51):
            u = Fraction(i,50)
            lhs = (1-u)**4+64*u-(1+3*u)*(1-u)
            rhs = 58*u+9*u*u-4*u**3+u**4
            self.assertEqual(lhs,rhs)
            self.assertGreaterEqual(rhs,0)

    def test_finite_constants(self):
        self.assertTrue(all(elementary_checks().values()))

    def test_selected_bound_and_finite_remainders(self):
        result = axial_bound(64,256)
        self.assertLess(result['maximum_bsw_absolute_upper'],mp.mpf('.6'))
        self.assertGreater(result['second_A24_margin_lower'],mp.mpf('.8'))
        self.assertGreater(result['relative_Q_loss_upper'],0)
        self.assertGreater(result['bs_squared_upper'],0)

    def test_inconclusive_does_not_mean_rejected(self):
        result = axial_bound(32,64)
        self.assertEqual(result['status'],'bound-inconclusive')

    def test_invalid_domains(self):
        for args in [(3,32,50),(4,0,50),(4,32,20),(4.5,32,50)]:
            with self.assertRaises(ValueError): axial_bound(*args)
        with self.assertRaises(ValueError): derivative_upper(Fraction(1),Fraction(0))


class SourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 80
        cls.m = Moments()

    def setUp(self):
        mp.mp.dps = 80

    def test_moment_constant_bounds(self):
        self.assertGreaterEqual(1-self.m.rI1,3/(8*mp.e))
        self.assertLessEqual(1-self.m.rI1,mp.mpf(3)/8)
        self.assertLess(self.m.energy1,mp.mpf('1.5'))

    def test_positive_Q_identity_and_symmetry(self):
        for eta in ['0','1e-12','.5','1']:
            row = self.m.axial(mp.expm1(2),eta)
            opposite = self.m.axial(mp.expm1(2),-mp.mpf(eta))
            self.assertGreater(row['Q'],0)
            self.assertLess(relative(row['Q'],row['Q_from_positive_identity']),mp.mpf('1e-65'))
            self.assertEqual(row['Q'],opposite['Q'])
            self.assertEqual(row['n'],-opposite['n'])
            self.assertEqual(row['bsw'],opposite['bsw'])

    def test_axis_memory_is_positive_and_shear_zero(self):
        row = self.m.axial(self.m.s.td,0)
        self.assertGreater(row['Q'],0)
        self.assertEqual(row['bsw'],0)
        self.assertEqual(row['bs_squared'],0)

    def test_axial_failure_against_source_ODE(self):
        y = mp.expm1(2)
        row = self.m.axial(y,'.5')
        source = independent_ode(self.m,'.5',[float(y)])[0]
        self.assertLess(row['Pc_over_ps1'],-1)
        for k in ('Q','n'):
            self.assertLess(relative(row[k],source[k]),mp.mpf('1e-9'))

    def test_later_lambda_failure_is_retained(self):
        row = self.m.intermediate_start('.5')
        self.assertGreater(row['Q'],0)
        self.assertGreater(row['lambda_w_squared'],mp.mpf('1e27'))

    def test_double_precision_long_stage_is_rejected(self):
        with self.assertRaises(ValueError): independent_ode(self.m,'.5',[64.])


class RecordTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = json.loads((HERE/'evidence.json').read_text())

    def test_pinned_record(self):
        validate(self.record)

    def test_scope_promotion_rejected(self):
        record = copy.deepcopy(self.record)
        record['status']='complete-profile-verified'
        with self.assertRaises(ValueError): validate(record)

    def test_missing_gate_rejected(self):
        record = copy.deepcopy(self.record)
        record['gates'].pop(next(iter(record['gates'])))
        with self.assertRaises(ValueError): validate(record)

    def test_failure_erasure_rejected(self):
        record = copy.deepcopy(self.record)
        record['old_Md_failure']['Pc_over_ps1']='1'
        with self.assertRaises(ValueError): validate(record)

    def test_provenance_change_rejected(self):
        record = copy.deepcopy(self.record)
        record['provenance'][next(iter(record['provenance']))]='0'*64
        with self.assertRaises(ValueError): validate(record)

    def test_record_numeric_change_rejected(self):
        record = copy.deepcopy(self.record)
        record['selected']['maximum_bsw_absolute_upper']='.01'
        with self.assertRaises(ValueError): compare_record(self.record,record)


if __name__ == '__main__':
    unittest.main()
