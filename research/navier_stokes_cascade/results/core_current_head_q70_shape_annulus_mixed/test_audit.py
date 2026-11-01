import json,unittest
from fractions import Fraction as Q
from math import factorial
import sympy as s
import audit,q70

class Tests(unittest.TestCase):
    def test_majorant_exponential(self):
        z=s.symbols('z');T=2000;J=3
        log=z/10+sum((J+6)*s.Rational(q70.step(k),factorial(k)*T**k)*z**k for k in range(1,5))
        p=s.exp(log).series(z,0,5).removeO()
        for k in range(5): self.assertEqual(p.coeff(z,k),q70.shape(k,T,J))

    def test_relative_edit_leibniz(self):
        y=s.symbols('y');b=Q(572,10**16)
        beta=sum(s.Rational(q70.bump(k))*y**k/factorial(k) for k in range(5))
        f=s.exp(y/10)*(1+3*b*beta)
        for k in range(5): self.assertEqual(s.diff(f,y,k).subs(y,0),q70.relative(k,b))

    def test_partial_pressure_reconstruction(self):
        X,z=s.symbols('X z', positive=True)
        g=X*(1+z*z)**2/2  # E=X(1+eta²)
        P=z**3-s.integrate(g,(X,X,2))
        self.assertEqual(P.subs(X,2),z**3)
        for i in range(1,4):
            for j in range(5):
                self.assertEqual(s.simplify(s.diff(P,X,i,z,j)-s.diff(g,X,i-1,z,j)),0)
        self.assertNotEqual(s.diff(P-z**3,z).subs({X:1,z:1}),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
