import json, unittest
from math import factorial
import sympy as s
import audit,q66

class Tests(unittest.TestCase):
    def test_polynomial_recurrence(self):
        u,w=s.symbols('u w');p=s.Integer(1)
        for n in range(9):
            poly=s.Poly(p,u,w)
            self.assertLessEqual(poly.total_degree(),2*n)
            self.assertLessEqual(sum(abs(c) for c in poly.coeffs()),2**n*factorial(n))
            bound=sum(abs(c)*factorial(i)*factorial(j) for (i,j),c in poly.terms())
            self.assertLessEqual(bound,8**n*factorial(n)**3)
            p=s.expand(-u*u*s.diff(p,u)+w*w*s.diff(p,w)+(u*u-w*w)*p)

    def test_plateau_amplitude_composition(self):
        x=s.symbols('x');a=1+x+x*x;b=s.Rational(1,4)+x*x
        v=a+b*b/a;delta=s.Rational(1,2)
        # v(0)=17/16 lies strictly in the K=1 plateau.
        actual=s.sqrt((2+delta/2-v)/a).series(x,0,5).removeO().expand()
        bound=q66.amplitude_majorant(x+x*x,x*x,s.Integer(1),s.Rational(1,4),delta,x,4)
        for j in range(5):
            self.assertTrue(abs(actual.coeff(x,j))<=bound.coeff(x,j))

    def test_normalized_cutoff_scaling(self):
        x=s.symbols('x');delta=s.Rational(1,3);v=x+x*x
        # Independent exact polynomial K surrogate with derivatives at the bound.
        y=s.symbols('y')
        K=sum(s.Rational(q66.cutoff_bound(j),factorial(j))*y**j for j in range(5))
        direct=K.subs(y,8*v/delta).expand()
        self.assertEqual(direct.coeff(x,2),q66.cutoff_bound(1)*24+q66.cutoff_bound(2)*24**2//2)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
