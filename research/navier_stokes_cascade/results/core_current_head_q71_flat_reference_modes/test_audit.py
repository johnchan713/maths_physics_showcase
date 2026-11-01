import json, unittest
import sympy as s
from math import factorial
from fractions import Fraction as Q
import audit,q71

class Tests(unittest.TestCase):
    def test_reference_source_powers(self):
        X,B,z=s.symbols('X B z',positive=True)
        F=1+z*z;u=1-3*z;M,I,J,S,P=[z**j for j in range(1,6)]
        M+=u*(X-B);I+=F*(X**2-B**2);J+=F*u*(X**2-B**2)
        S+=u**2*(X-B)-F**2*(X**2-B**2)/2;P+=F**2*(X-B)
        h=s.Rational(1,100);D=s.Rational(1,2)-h;A=1-D;d=1-z*z;L=1-2*h*z*z
        W=1-(2*D*z*M+d*s.diff(M,z))/X
        N=(1-h)*I-D*z*s.diff(I,z)-d*s.diff(J,z)+2*(h-D)*z*J
        p=(-X*W+N/(2*X*F))/L
        n=(-W*u+(D*(M-z*s.diff(M,z))+4*h*z*S-d*s.diff(S,z))/X+4*A*z*P-d*s.diff(P,z))/L
        for f in (p,n):
            self.assertEqual(s.Poly(s.cancel(X*f),X).degree(),2)

    def test_mode_interpolation(self):
        for q,rows in [(s.Integer(2),[8,9,2]),(s.Rational(1,2),[2,9,8])]:
            inv=s.Matrix([[q**(j*t) for j in (-1,0,1)] for t in range(3)]).inv()
            self.assertEqual([sum(map(abs,inv.row(j))) for j in range(3)],rows)
            a=s.Matrix([3,-7,2])
            self.assertEqual(inv*s.Matrix([3*q**(-t)-7+2*q**t for t in range(3)]),a)

    def test_switch_product_derivatives(self):
        y=s.symbols('y');delta=Q(1,16)
        t=sum(q71.switch(i,delta)*y**i/factorial(i) for i in range(5))
        # Positive majorants, with derivatives of exp(2y) equal to 2**r.
        for i in range(5):
            a,b=q71.slopes(i,delta,110,1,1)
            self.assertEqual(s.diff(delta*10*t*s.exp(y)+s.Rational(4,5)*t,y,i).subs(y,0),a)
            self.assertEqual(s.diff(delta*19*110*t*s.exp(2*y)/2,y,i).subs(y,0),b)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
