from fractions import Fraction as Q
from math import factorial,comb
import json,unittest
import sympy as s
import audit,q64


class Tests(unittest.TestCase):
    def test_contraction_disk(self):
        for m in (Q(1,8),Q(1),Q(2)):
            for D in (Q(1,10),Q(3),Q(10)):
                rho,delta=q64.disk(m,D)
                self.assertLessEqual(D*rho/m,Q(1,4))
                self.assertLessEqual(delta/m+D*rho*rho/(2*m),rho/2)

    def test_frequency_factorials(self):
        for r in range(9):
            k=2; N=7;delta=Q(1,13)
            H=lambda a,b:Q(1,1+a+b)
            raw=sum(comb(r,j)*Q(N**j,N)*factorial(r-j)*factorial(k)*
                    H(r-j,k)*factorial(j)*(2/delta)**j for j in range(r+1))
            self.assertEqual(q64.pullback(r,k,N,delta,H),raw)

    def test_moving_inverse_has_slow_term(self):
        x,z=s.symbols('x z');a=1+x
        # Explicit local nonlinear phase phi=a*theta+theta^2/2.
        theta=s.sqrt(a*a+2*z)-a
        phi=a*theta+theta**2/2
        self.assertEqual(s.simplify(phi-z),0)
        self.assertEqual(s.simplify(theta.diff(x)+theta/(a+theta)),0)
        self.assertEqual(s.simplify(theta.diff(z)-1/(a+theta)),0)
        self.assertNotEqual(theta.diff(x).subs({x:0,z:1}),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__=='__main__':unittest.main()
