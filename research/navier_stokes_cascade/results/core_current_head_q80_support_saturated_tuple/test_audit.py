import unittest,json
from fractions import Fraction as Q
import sympy as s
import audit,q80
class Tests(unittest.TestCase):
    def test_clipped_hull_and_exact_exterior(self):
        args=(Q(1,200),1,1);f=(7,11);u=(5,17)
        cut=q80.clipped(*args,1,10,3,f,u,zero_total=True)
        full=q80.radial.reconstruct(*args,200,Q(1,2),20,f,u)['through'][1]
        self.assertLess(cut['bound'],full)
        self.assertEqual(cut['native_upper'],3)
        self.assertEqual(q80.clipped(*args,8,10,3,f,u,zero_total=True)['bound'],0)
        with self.assertRaises(ValueError):q80.clipped(*args,1,10,3,f,u,zero_total=False)
        with self.assertRaises(ValueError):q80.clipped(Q(1,200),0,1,1,10,3,f,u,zero_total=True)
    def test_primitive_cancellation_not_local_velocity(self):
        x,e=s.symbols('x e');M=e*x**2*(1-x)**3;U=s.diff(M,x)
        self.assertEqual(s.integrate(U,(x,0,1)),0)
        for k in range(3):self.assertEqual(s.diff(M,e,k).subs(x,1),0)
        self.assertEqual(U.subs(x,1),0);self.assertEqual(s.diff(U,x).subs(x,1),0)
        # A compact velocity with nonzero total leaves a primitive and radial tail.
        bad=s.integrate(e*x**2*(1-x)**2,(x,0,1))
        h=s.Rational(1,200);D=s.Rational(1,2)-h
        V=(-2*(D+2*h)*e*bad-(1-e**2)*s.diff(bad,e))/(1-2*h*e**2)
        self.assertNotEqual(V.subs(e,0),0)
    def test_pressure_survives_velocity_zero(self):
        a,b,mu,c,q,h,K1,K2=s.symbols('a b mu c q h K1 K2',positive=True)
        P=lambda X:c*q**(2*h)*(K1-mu**2/(4*X))+q**(4*h)*K2
        self.assertEqual(s.simplify(P(b)-P(a)-c*q**(2*h)*mu**2*(1/a-1/b)/4),0)
    def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
