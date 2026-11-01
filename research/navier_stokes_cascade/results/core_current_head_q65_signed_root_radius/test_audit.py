from fractions import Fraction as Q
import json,unittest
import sympy as s
import audit,q65


class Tests(unittest.TestCase):
    def test_uniform_disks(self):
        for G,E in ((2,9),(3,27),(4,81)):
            r=q65.radii(G,E,Q(1,32))
            A,C,b,R,tau,M,m,a,P,sig=(r[k] for k in ('A','C','b','R','tau','M','m','a','P','sigma'))
            self.assertLess(C*C*((1-A*R)**(-2)-1),b/2)
            self.assertLessEqual(2*G*tau+tau*tau,R/2)
            self.assertLessEqual(12*M*a/tau**2,m/4)
            self.assertLessEqual((1+P)*sig/m,a/4)

    def test_removable_variance(self):
        z=s.symbols('z')
        I=sum(z**(2*k)/(4**k*s.factorial(k)**2) for k in range(7))
        V=((I.subs(z,2*z)/I**2-1)/z**2).series(z,0,8).removeO()
        self.assertEqual(V.subs(z,0),s.Rational(1,2))
        self.assertEqual(s.simplify(V.subs(z,-z)-V),0)
        mu,p=s.symbols('mu p')
        F=mu*s.sqrt(V.subs(z,mu*p))
        self.assertEqual(s.simplify(F.subs(p,0)),mu/s.sqrt(2))

    def test_finite_two_input_composition(self):
        x,y=s.symbols('x y')
        p=x+x*y;v=y+x*x
        # psi(p,s)=s/(1-p); coefficients <=1 on its germ.
        scale=s.symbols('t')
        exact=(v/(1-p)).subs({x:scale*x,y:scale*y}).series(scale,0,7).removeO().subs(scale,1).expand()
        bound=sum((p+v)**k for k in range(1,7)).expand()
        for i in range(7):
            for j in range(7-i):
                self.assertLessEqual(abs(exact.coeff(x,i).coeff(y,j)),bound.coeff(x,i).coeff(y,j))

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))


if __name__=='__main__':unittest.main()
