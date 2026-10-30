import unittest,json
from fractions import Fraction as Q
import sympy as s
import fixed_maps as f
import audit

class FixedTests(unittest.TestCase):
    def test_gap_unchanged_by_scale(self):
        eta=f.gap(1,16,2,8)
        self.assertEqual(eta,Q(2,3))
        for scale in (Q(1,4),1,4):
            self.assertEqual(f.gap(scale,16*scale,2*scale,8*scale),eta)

    def test_band_dependent_collar_has_no_fixed_floor(self):
        values=[f.gap(1,16,1+Q(1,n),8) for n in (2,10,100)]
        self.assertGreater(values[0],values[1]);self.assertGreater(values[1],values[2])
        self.assertLess(values[2],Q(1,50))

    def test_log_lower_bound_identity(self):
        r=s.symbols('r',positive=True)
        self.assertEqual(s.factor(s.diff(s.log(r)-2*(r-1)/(r+1),r)),(r-1)**2/(r*(r+1)**2))

    def test_full_angular_matrix_and_inverse_norm(self):
        g=s.Rational(2);v=s.Rational(1,3);nodes=[g*g,v/(g*g),v]
        V=s.Matrix([[z**j for j in range(3)] for z in nodes]);B=s.diag(2,6,-5)*V
        D=f.determinants(2,Q(1,3))
        self.assertEqual(abs(V.det()),D['angular'])
        self.assertEqual(abs(B.det()),60*D['angular'])
        inv=B.inv();bound=max(sum(abs(B.cofactor(j,i)) for j in range(3)) for i in range(3))/abs(B.det())
        self.assertEqual(bound,max(sum(abs(inv[i,j]) for j in range(3)) for i in range(3)))
        for target in (s.Matrix([0,-1,0]),s.Matrix([0,0,-1])):
            self.assertEqual(B*(inv*target),target)

    def test_axial_conditioning(self):
        for v in (Q(1,2),Q(9,10),Q(99,100)):
            V=s.Matrix([[1,2],[1,2*v]])
            self.assertEqual(abs(V.det()),f.determinants(2,v)['axial'])
        self.assertEqual(s.Matrix([[1,2],[1,2]]).det(),0)
        self.assertGreater(1/f.determinants(2,Q(99,100))['axial'],1/f.determinants(2,Q(1,2))['axial'])

    def test_hull_contains_bumps_and_gaps(self):
        lo,hi=f.hull(2,3,2);self.assertEqual((lo,hi),(2,12))
        self.assertTrue(lo<7<hi)
        self.assertFalse(any(a<=7<=b for a,b in ((2,3),(4,6),(8,12))))
        for root_s in (Q(1,2),1,2):
            self.assertEqual((root_s*hi)**2/(2*root_s**2),hi**2/2)

    def test_zero_moment_compact_potential(self):
        x=s.symbols('x');F=1-s.Rational(9,14)*x
        self.assertEqual(s.integrate(x*F,(x,1,2)),0)
        y=s.symbols('y');Psi=s.integrate(y*F.subs(x,y),(y,1,x))/x
        self.assertEqual(s.simplify(s.diff(Psi,x)+Psi/x-F),0)
        self.assertEqual(Psi.subs(x,1),0);self.assertEqual(Psi.subs(x,2),0)
        # Polynomial interval fixture checks moments, not smooth extension.

    def test_guards(self):
        with self.assertRaises(ValueError):f.gap(1,4,1,3)
        with self.assertRaises(ValueError):f.determinants(2,1)
        with self.assertRaises(ValueError):f.hull(2,3,1)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['full_pressure_identified_with_leading'])
if __name__=='__main__':unittest.main()
