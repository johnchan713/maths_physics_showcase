import json,unittest
from fractions import Fraction as Q
import sympy as s
import q68,audit

class Tests(unittest.TestCase):
    def test_ideal_trace_densities(self):
        x,z,P,R=s.symbols('x z P R',positive=True); f=1/(1+z*z)
        E=P*f*x**s.Rational(1,10);U=4*z;X=R*x
        I=5*s.sqrt(2)*R**s.Rational(3,2)*P*f*x**s.Rational(8,5)/8
        primitives=[4*z*X,I,4*z*I,R*(16*z*z*x-5*P**2*f**2*x**s.Rational(6,5)/12),5*P**2*f**2*x**s.Rational(1,5)/2]
        for v,g in zip(primitives,[U,s.sqrt(2*X)*E,s.sqrt(2*X)*E*U,U**2-E**2/2,E**2/(2*X)]):
            self.assertEqual(s.simplify(s.diff(v,x)/R-g),0)

    def test_forward_integral_control(self):
        # Integrate the disk density majorants on [2,8], including the
        # logarithmic pressure row; these are not a sampled physical field.
        R,B,P,C=2,8,2,2
        m=q68.moments(R,B,P,C,2,4,2)
        start=q68.moments(R,R,P,C,0,2,2)
        X=s.symbols('X',positive=True)
        for u,v,g in zip(start,m,[6,4*s.sqrt(2*X),24*s.sqrt(2*X),44,8/X]):
            self.assertTrue(s.integrate(g,(X,R,B))<=v-u)

    def test_atom_and_pressure_derivative(self):
        z=s.symbols('z'); f=1/(1+z*z)
        axis=-s.Rational(5,2)*f**2-f
        merged=s.simplify(axis+s.Rational(5,2)*f**2)
        self.assertEqual(s.simplify(merged+f),0)
        self.assertNotEqual(s.diff(merged,z).subs(z,s.Rational(1,2)),0)
        self.assertGreater(1-2*Q(1,100)*Q(3,2)**2,Q(1,2))
        self.assertGreaterEqual(1-Q(1,2)**2,Q(3,4))

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
