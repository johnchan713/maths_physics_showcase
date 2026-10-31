import json, unittest
import sympy as s
import audit
import q58_weights as weights


class Tests(unittest.TestCase):
    def test_nonconstant_scale_pressure_moment(self):
        x,z=s.symbols('x z');scale=1+z**2
        chi=(x-2)**2*(3-2*(x-2))  # smoothstep on [2,3]
        mu=3-s.integrate(chi,(x,2,3))
        val=s.integrate(scale*x*s.diff(chi,x),(x,2,3))
        self.assertEqual(s.simplify(val-scale*mu),0)
        self.assertNotEqual(s.diff(val,z),0)

    def test_native_swirl_weight(self):
        U,C=s.symbols('U C',positive=True);R=U**s.Rational(2,3)
        F=U**3-U
        w=C-R**2/2; a=2*w*F
        val=s.diff(a,U)-2*a*s.diff(R,U)/R
        reduced=(2*C-R**2)*s.diff(F,U)-4*C*s.diff(R,U)*F/R
        self.assertEqual(s.simplify(val-reduced),0)
        self.assertNotEqual(s.simplify(val-reduced+4*C*s.diff(R,U)*F/R),0)

    def test_radial_bound_on_fractional_chart(self):
        U=s.symbols('U',positive=True);R=U**s.Rational(2,3)
        coeff,base=weights.geometry(1,2,3)
        for m in range(6):
            for i,expr in enumerate((R*s.diff(R,U),s.diff(R,U)*(3-R**2/2)/R)):
                self.assertLess(abs(s.diff(expr,U,m).subs(U,1)),coeff[i]*s.factorial(m)*base**m)

    def test_record(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
