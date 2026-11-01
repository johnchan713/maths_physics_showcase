import json,unittest
from fractions import Fraction as Q
import sympy as s
import audit,q73

class Tests(unittest.TestCase):
    def test_varying_quadratic_saturates_positive_recurrence(self):
        x=s.symbols('x');c0=s.Rational(1,10)
        # c-(1+2x+3x²)c²=z. Independent formal Newton iteration.
        c=c0;z=c0-c0**2+sum(2*(8*x)**k for k in range(1,6));q=1+2*x+3*x*x
        for _ in range(3):
            c=(c-(c-q*c*c-z)/(1-2*q*c)).series(x,0,6).removeO().expand()
        self.assertEqual(q73.root_bounds(5,2,Q(1,8),Q(1,10),Q(5,4),(1,2,3)),
                         [c.coeff(x,k) for k in range(6)])

    def test_ledger_normalizer_and_axis_cancellation(self):
        z=s.symbols('z');axis=z**3-3*z
        cp1=1+z;cp2=2-z*z;f=1/(1+z*z)
        row=((axis+cp1)-(axis+cp2))/f**2
        self.assertEqual(s.expand(row),s.expand((cp1-cp2)*(1+z*z)**2))
        self.assertNotEqual(s.diff(row,z,2).subs(z,0),s.diff(cp1-cp2,z,2).subs(z,0))

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
