import unittest,json
from math import factorial
import sympy as s
import audit,q78

class Tests(unittest.TestCase):
    def test_mixed_table_independent_polynomial(self):
        p,q=s.symbols('p q');exp=sum(p**i/s.factorial(i) for i in range(4))
        J=2*((exp-1)+exp*(2*q+q*q))
        W=s.Poly(4+sum(s.factorial(a+1)*J**a for a in range(1,6)),p,q)
        table=q78.heat_table(3,2)
        for i in range(4):
            for j in range(3): self.assertEqual(table[i,j],W.coeff_monomial(p**i*q**j))

    def test_real_pole_mixed_derivative(self):
        # At eta=1, H_eta=-4 H'(0)/X for every radius:
        # taking y derivatives preserves its magnitude, not zero.
        h=s.Rational(1,8);X=7
        for i in range(5):
            exact=4*h*(1+h)/(X*factorial(i))
            self.assertLessEqual(exact,h*q78.heat_table(i,1)[i,1]/X)

    def test_terminal_product_cross_term(self):
        y=s.symbols('y');r=s.Rational(1,8)
        a=(1-r*(1-y/2))/(1-r);b=1+y*y/3
        self.assertEqual(s.expand((a*b).diff(y,2)-a*b.diff(y,2)-a.diff(y,2)*b),s.expand(2*a.diff(y)*b.diff(y)))
        self.assertNotEqual(a.diff(y)*b.diff(y),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
