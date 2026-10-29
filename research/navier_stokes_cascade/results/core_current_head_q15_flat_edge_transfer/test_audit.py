from decimal import Decimal, localcontext
from fractions import Fraction as Q
import json
import unittest
import audit
import edge

def dec(q):
    q=Q(q);return Decimal(q.numerator)/Decimal(q.denominator)

class EdgeTests(unittest.TestCase):
    def test_ceiling_covers_fractional_poles(self):
        for p in (Q(0),Q(1,3),Q(7,2),Q(400009)):
            for s in (Q(0),Q(1,3),Q(9)):
                r=edge.request(p,s)
                self.assertTrue(0<=Q(r['power_margin'])<2)
    def test_rational_majorant_all_delta(self):
        for w in (Q(1),Q(1,100),Q(1,10**30)):
            for p,s in ((0,0),(3,4),(8,5)):
                for delta in (Q(1),Q(1,3),Q(1,10**30)):
                    self.assertLessEqual(edge.normalized_majorant(w,p,s,delta),edge.fixture_constant(w,p,s))
    def test_independent_exponential_diagnostics(self):
        with localcontext() as ctx:
            ctx.prec=80
            for w in (Q(1),Q(1,10)):
                for L,x in ((Q(2),Q(1,10)),(Q(2),Q(19,10)),(Q(100),Q(50))):
                    delta=min(Q(1),x,L-x)
                    zeta=dec(-w*w/(x*x)-4/(L-x)**2).exp()
                    for p,s in ((0,0),(3,4),(8,5)):
                        lhs=zeta*dec(delta)**(-p)
                        rhs=dec(edge.fixture_constant(w,p,s))*dec(delta)**s
                        self.assertLessEqual(lhs,rhs)
    def test_interior_case_cannot_use_edge_exponential(self):
        w,x,L=Q(1),Q(50),Q(100)
        self.assertEqual(min(1,x,L-x),1)
        self.assertLess(w*w/x**2+4/(L-x)**2,w*w)
        self.assertGreaterEqual(edge.fixture_constant(w,0),1)
    def test_fixed_width_is_essential(self):
        for ell in (1,2,10,1000):
            w=Q(1,ell)
            self.assertLessEqual(1+4/(2-w)**2,5)
            self.assertEqual(w**-2,ell**2)
    def test_huge_request_stays_symbolic(self):
        r=edge.request(10**100,12)
        self.assertEqual(r['n'],(10**100+12)//2)
        self.assertFalse(r['materialized_constant'])
        with self.assertRaises(ValueError):edge.fixture_constant(Q(1,10),10**100)
    def test_invalid_inputs(self):
        for p in (-1,0.5,True):
            with self.assertRaises(ValueError):edge.request(p)
        for w in (0,-1,2,0.1):
            with self.assertRaises(ValueError):edge.fixture_constant(w,1)
    def test_frozen_record_and_scope(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['scope']['complete_representative_bounds_proved'])
        self.assertFalse(r['scope']['q_exponent_reserve_spent'])

if __name__=='__main__':unittest.main()
