import json,unittest
import sympy as s
import audit,q76

class Tests(unittest.TestCase):
    def test_gamma_endpoint_and_composition(self):
        z=s.symbols('z');h=s.Rational(1,8);X=3
        for k in range(1,7):
            # At eta=1 the exact Taylor polynomial in Z gives the true jets.
            series=sum((-1)**i*s.rf(h,i)*s.rf(1+h,i)*(-4*z/X-2*z*z/X)**i/s.factorial(i) for i in range(1,k+1)).expand()
            self.assertLessEqual(sum(abs(series.coeff(z,j)) for j in range(1,k+1)),h*q76.kernel(k)/X)

    def test_tail_integrals(self):
        x,h=s.symbols('x h',positive=True)
        for p,answer in [(3+2*h,1/(2+2*h)),(2+2*h,1/(1+2*h)),(1+h,1/h)]:
            self.assertEqual(s.simplify(s.integrate(x**(-p),(x,1,s.oo))-answer),0)

    def test_relative_bump_changes_pressure_before_exit(self):
        t,z=s.symbols('t z');E=(1+z*z)**-1*s.exp(-t);b=z*t*(1-t)
        inc=(E**2*((1+b)**2-1)/2).expand()
        self.assertEqual(s.simplify(inc-E**2*(b+b*b/2)),0)
        self.assertNotEqual(s.diff(inc,z).subs({z:0,t:s.Rational(1,2)}),0)
        self.assertEqual(inc.subs(t,1),0)

    def test_evidence(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
