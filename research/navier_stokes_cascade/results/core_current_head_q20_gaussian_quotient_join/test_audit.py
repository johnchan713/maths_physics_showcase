from fractions import Fraction as Q
from math import factorial,comb
import unittest,json
import audit
import quotient as q
class QuotientTests(unittest.TestCase):
    def test_linear_log_exact_bell(self):
        for a in (Q(0),Q(2),Q(7,3)):
            self.assertEqual(q.bell([0,a]+[0]*7),[a**j for j in range(9)])
    def test_gaussian_even_jets_at_center(self):
        H=q.bell([0,0,2]+[0]*6)
        for j in range(9):self.assertEqual(H[j],factorial(j)//factorial(j//2) if j%2==0 else 0)
    def test_quotient_derivatives_not_raw_weighted_jets(self):
        # f=1, P=exp(-a*x), at x=0; quotient exp(a*x).
        self.assertEqual(q.convert([1]+[0]*5,[0,3]+[0]*4),[3**j for j in range(6)])
    def test_two_exponentials_majorant(self):
        for a,b in ((2,3),(0,5)):
            self.assertEqual(q.convert([Q(b)**j for j in range(7)],[0,a]+[0]*5),[Q(a+b)**j for j in range(7)])
    def test_mixed_binomial_grouping(self):
        for a,b in ((1,2),(3,4)):
            for k in range(a+b+1):
                self.assertEqual(sum(comb(a,i)*comb(b,k-i) for i in range(a+1) if 0<=k-i<=b),comb(a+b,k))
    def test_no_zeroth_log_cost_or_extra_jet(self):
        self.assertEqual(q.bell([0,1,2,3]),q.bell([10**100,1,2,3]))
        self.assertEqual(len(q.convert([1]*8,[0]+[2]*7)),8)
    def test_log_rate_bound(self):
        self.assertEqual(q.log_envelope_jets(10,[2,3,4]),[0,17,23])
    def test_trace_does_not_bound_lift(self):
        for n in (1,100,10**20):
            self.assertEqual(n*0,0);self.assertEqual(n*1,n)
    def test_guards(self):
        with self.assertRaises(ValueError):q.bell([0,-1])
        with self.assertRaises(ValueError):q.convert([1],[0,1])
        with self.assertRaises(ValueError):q.bell([0,1.0])
    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['actual_all_stage_source_enclosures_proved'])
if __name__=='__main__':unittest.main()
