from fractions import Fraction as F
import json, unittest
import sympy as s
import audit, q61


class Tests(unittest.TestCase):
    def test_uniform_scale_split(self):
        for k in range(1,20):
            e=F(1,k)
            for B in (1,2,8,32,200):
                case,factor=q61.branch(e,1,B,1)
                if case=='compare':
                    self.assertEqual(factor*(1-e*e*B),1)
                    self.assertLessEqual(factor,2)
                else:
                    self.assertLessEqual(1/e**2,factor)

    def test_gaussian_weighted_propagator(self):
        # V=P(v)/P(w); constant perturbation delta, L=7, nonconstant P.
        v,w=s.symbols('v w')
        P=s.exp(-(v-3)**2); delta=s.Rational(1,16)
        V=P/P.subs(v,w)*s.exp(delta*(v-w))
        self.assertEqual(s.simplify(V.diff(v)-(P.diff(v)/P+delta)*V),0)
        case,factor=q61.branch(F(1,4),1,7,1)
        self.assertEqual(case,'compare')
        self.assertLess(s.exp(delta*7),factor)

    def test_column_difference_with_sign_change(self):
        x,e=s.symbols('x e', positive=True)
        H=s.Matrix([[1,x],[0,1]])
        H0=H+s.Matrix([[0,e**2*(1+x)],[0,0]])
        T=H*s.Matrix([e**2/2,1])
        z=H.inv()*T; z0=H0.inv()*T
        self.assertEqual(s.simplify(z0-z-H0.inv()*(H-H0)*z),s.zeros(2,1))
        self.assertEqual(s.simplify((z0-z)[0].diff(x)), -e**2)
        self.assertLess(s.simplify(z0[0]),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__=='__main__':unittest.main()
