import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q82
class Tests(unittest.TestCase):
    def test_radial_adjoint(self):
        R=s.symbols('R',positive=True);V=R**3+2*R
        h=(R-1)**3*(2-R)**3;a=2*V/R
        direct=s.integrate(a*(s.diff(h,R)+2*h/R),(R,1,2))
        self.assertEqual(direct,s.integrate((6*V/R**2-2*s.diff(V,R)/R)*h,(R,1,2)))
        self.assertNotEqual(direct,s.integrate(6*V*h/R**2,(R,1,2)))
    def test_shells(self):
        for start,r in [(Q(1),0),(Q(1),1),(2-3*q82.K,0)]:
            b=q82.row(start,r)
            self.assertEqual(start+b['p']*q82.K,q82.TARGET)
            for (j,q),power in zip(b['radial_auxiliary_pairs'],(2,-2)):
                self.assertEqual(j-r-q,power);self.assertEqual(j+q,b['extra_source'])
                self.assertEqual(-2*j+power,-2*b['p'])
    def test_fallback_gap(self):
        g=q82.TARGET-1
        self.assertEqual(g,2-5*q82.K);self.assertLess(g,2)
        # At epsilon=1/N take B=N^2; no constant pays a gap >2 as N grows.
        for N in (4,16,64):self.assertGreater(Q(N**3,2*N*N),1)
    def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
