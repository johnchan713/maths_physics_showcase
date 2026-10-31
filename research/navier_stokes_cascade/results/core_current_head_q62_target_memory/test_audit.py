from fractions import Fraction as F
from math import factorial
import json, unittest
import sympy as s
import audit, q62


class Tests(unittest.TestCase):
    def test_full_stress(self):
        E,U,t=s.symbols('E U t',nonzero=True)
        M,I,J,S,P,Me,Ie,Je,Se,Pe=s.symbols('M I J S P Me Ie Je Se Pe')
        X,r,h=s.Integer(2),s.Integer(2),s.Rational(1,100)
        D,A=s.Rational(1,2)-h,s.Rational(1,2)+h
        L,d=1-2*h*t*t,1-t*t
        W=1-(2*D*t*M+d*Me)/X
        Q=-W+((1-h)*I-D*t*Ie-d*Je+2*(h-D)*t*J)/(X*r*E)
        N=-W*U+(D*(M-t*Me)+4*h*t*S-d*Se)/X+4*A*t*P-d*Pe
        direct=(E/r*(X*Q/L-(1-2*X*3/E)),E/r*(X*N/(L*E)+2*X*7/E))
        alg=q62.stress(X,r,h,t,E,U,3,7,M,I,J,S,P,Me,Ie,Je,Se,Pe)
        self.assertEqual([s.simplify(a-b) for a,b in zip(alg,direct)],[0,0])
        self.assertEqual(alg[1].diff(Pe).subs(t,0),-X/r)

    def test_real_endpoint_heat_jet(self):
        z=s.symbols('z'); h=s.Rational(1,100)
        for B in (1,2,5):
            # At eta=1, zeta=0; exact normalized Gamma moments, one-sided.
            arg=2*(1-(1+z)**2)/B
            p=sum((-1)**j*s.rf(h,j)*s.rf(1+h,j)*arg**j/factorial(j)
                  for j in range(1,6)).expand()
            C,G=q62.heat(5,h,B)
            self.assertLessEqual(sum(abs(p.coeff(z,j)) for j in range(1,6)),C/B)
            self.assertEqual(G,1+C/B)

    def test_mixed_primitive(self):
        x,t=s.symbols('x t');g=(1+x+t)**4;P=t**3-s.integrate(g,(x,x,3))
        for i in range(4):
            for j in range(4):
                def norm(f,a,b,X):
                    return abs(f.diff(x,a).diff(t,b).subs({x:X,t:1}))/factorial(a)/factorial(b)
                val=q62.primitive(i,j,3,2,lambda a,b:norm(g,a,b,3 if i==0 else 1))
                self.assertLessEqual(norm(P,i,j,1),val)
                if i:self.assertEqual(norm(P,i,j,1),val)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__=='__main__':unittest.main()
