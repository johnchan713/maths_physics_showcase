import unittest,json
from fractions import Fraction as Q
from math import factorial
import sympy as s
import haar_norms as h
import audit

class HaarNormTests(unittest.TestCase):
    def test_raw_is_not_through_norm(self):
        x=s.symbols('x0:5');f=s.exp(sum(x));at=dict.fromkeys(x,0)
        raw=[f.subs(at)]+[s.diff(f,v).subs(at) for v in x]
        self.assertEqual(max(raw),1);self.assertEqual(sum(raw),6)
        self.assertEqual(h.count(1),6)
        self.assertEqual(h.count(2),Q(37,2))

    def test_mixed_coordinate_majorant(self):
        R,Z=s.symbols('R Z',positive=True)
        chart=h.reverse.ReverseChart(1,2)
        for d in (s.Rational(1,2),s.Rational(2),s.Rational(12,5)):
            f=(R**d)**2*Z**2
            for n in range(4):
                raw=sum(abs(s.diff(f,R,i,Z,j).subs({R:1,Z:1})) for i in range(n+1) for j in range(n+1-i))
                self.assertLessEqual(raw,4*chart.transfer(n).evaluate(1))

    def test_distinct_covariance_and_background_conversions(self):
        p=h.common.normal.roots.fixture();stress=h.common.normal.stress.fixture
        for n in range(3):
            v=h.r_inputs(n,p,stress,Rmin=1,Rmax=2)
            raw=h.common.common_inputs(n,p,stress,dict(Rmin=1,Rmax=2))
            self.assertEqual(v['W'].c,(v['covariance_transfer']*raw['W']).c)
            self.assertEqual(v['higher'].c,(h.count(n)*raw['higher']).c)
            self.assertEqual(v['D'].c,(v['Dcurl']+v['higher']).c)
            if n:self.assertNotEqual(v['D'].c,raw['D'].c)
            self.assertEqual(v['pole'],12*n+18)

    def test_scalar_pressure_requests_next_covariance_jet(self):
        trace=[]
        def inputs(j):trace.append(j);return dict(W=h.P(j+1))
        v=h.initial_scalars(2,inputs,{-1:1,0:2,1:3,2:4},h.P(5))
        self.assertEqual(trace,[2,3])
        self.assertEqual(v['P'].evaluate(1),150)
        self.assertEqual(v['Jtheta'].evaluate(1),60)
        self.assertEqual(v['Jz'].evaluate(1),210)

    def test_pressure_cross_term_does_not_vanish(self):
        R,Z=s.symbols('R Z');wrz=R*Z
        term=-s.diff(s.integrate(wrz,(R,1,2)),Z)
        self.assertEqual(term,-s.Rational(3,2))

    def test_derivative_shift_count(self):
        # f=z^(n+1) at zero saturates the derivative-shift factor.
        for n in range(6):
            self.assertEqual(Q(factorial(n+1),factorial(n)),n+1)
        with self.assertRaises(ValueError):h.initial_scalars(0,None,{0:1},h.P(1))
        with self.assertRaises(TypeError):h.initial_scalars(0,None,{-1:1,0:1,1:1,2:1},1)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['full_source_residual_or_blowup'])

if __name__=='__main__':unittest.main()
