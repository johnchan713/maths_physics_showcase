import json, unittest
from fractions import Fraction as Q
import sympy as s
import angular as a
import audit

class Tests(unittest.TestCase):
    def test_two_updates_nonconstant_graph(self):
        R,M=s.symbols('R M',positive=True)
        # A single Fourier mode: L and Nt have different nonzero symbols.
        l=3*s.I; inv=1/(5*s.I); j=R+1; W=R**3
        dr=lambda f:s.diff(f,R)+M*j*l*f
        lap=lambda f:dr(dr(f))+dr(f)/R-f/R**2
        f0=dr(W)+2*W/R
        u0=-inv*f0
        f1=-lap(u0)
        u1=-inv*f1
        final=-lap(u1)
        pol=s.Poly(final.expand(),M)
        self.assertEqual(pol.degree(),5)
        self.assertEqual(s.expand(pol.LC()-j**5*inv**2*l**5*W),0)
        self.assertNotEqual(pol.nth(4),0)
    def test_rate_cap(self):
        k=Q(1,100000);d=a.ledger(k)
        self.assertEqual(d['final_leading'],3-5*k)
        self.assertEqual(d['final_error']-d['final_leading'],k)
        self.assertEqual(a.necessary_order(3-5*k+k/2),1+k/2)
        self.assertEqual(a.necessary_order(3),1+k)
    def test_two_divisors(self):
        b=s.sqrt(2)-1
        for x,y in ((2,5),(5,12),(12,29),(-12,-29)):
            l=x-b*y;t=b*x+y;r=abs(x)+abs(y)
            self.assertEqual(abs(s.expand(l*t)),b*abs((x+y)**2-2*y*y))
            self.assertTrue(bool(abs(t**2/l**5)<27**2*r**7))
    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
