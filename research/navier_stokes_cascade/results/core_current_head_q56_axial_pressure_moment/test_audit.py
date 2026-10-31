import json
import unittest
import sympy as s
import audit
import moment


class Tests(unittest.TestCase):
    def test_moving_two_branch_pressure_on_characteristic(self):
        # All quantities are pulled back to the SAME path; allow nonzero mass.
        r,x,z=s.symbols('r x z')
        chi=(r-1)+z*(r-1)*(2-r)
        g=x**2+z*x
        integ=lambda f:s.integrate(f,(r,1,2))
        mass=s.integrate(g,(x,1,2))
        pressure=s.integrate(g,(x,1,r))-chi*mass
        c=2-integ(r*chi)
        original=integ(r*pressure)
        reduced=integ((c-r**2/2)*g.subs(x,r))
        self.assertEqual(s.expand(original-reduced),0)
        frozen=integ((c-r**2/2)*s.diff(g,z).subs(x,r))
        self.assertEqual(s.simplify(s.diff(original,z)-frozen-s.diff(c,z)*mass),0)
        self.assertNotEqual(s.diff(c,z)*mass,0)

    def test_adjoint_on_resonant_bump(self):
        r,c=s.symbols('r c',positive=True)
        # Envelope after phase cancellation on the characteristic.
        F=r*(r-1)**3*(2-r)**3
        V=r**2+1;G=r**3;w=c-r**2/2
        d=moment.weights(r,c,V,G)
        integ=lambda f:s.integrate(s.expand(f),(r,1,2))
        self.assertEqual(integ(w*s.diff(s.diff(F,r)+F/r,r)),0)
        self.assertEqual(s.simplify(integ(-2*V*w/r*(s.diff(F,r)+2*F/r)-d['rt']*F)),0)
        self.assertEqual(s.simplify(integ(-w*(s.diff(F,r)+F/r)-d['rr']*F)),0)
        self.assertEqual(s.simplify(sum(d[i] for i in ('rr','tt','zz'))),0)
        self.assertNotEqual(integ(d['rt']*F),0)

    def test_frozen(self):
        self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
