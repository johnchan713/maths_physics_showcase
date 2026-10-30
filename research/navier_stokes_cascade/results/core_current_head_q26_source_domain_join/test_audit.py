from fractions import Fraction as Q
import unittest,json
import sympy as s
import domain as d
import audit

class DomainTests(unittest.TestCase):
    def test_dyadic_closed_boundaries(self):
        for star in (4,16,30):
            q=Q(1,2**(star+2))
            for ell in range(star,star+5):
                r=d.coverage(star,1,q,ell)
                self.assertTrue(r['accepted']);self.assertTrue(r['background_inside'])
                self.assertEqual(r['q_unchanged'],q)
            with self.assertRaises(ValueError):d.coverage(star,1,q,star-1)

    def test_background_cap_controls_full_enlargement(self):
        qstar=Q(1,2**40);q=qstar/64
        r=d.coverage(16,qstar,q,44)
        self.assertEqual(r['background_upper'],qstar/4)
        self.assertTrue(r['background_inside'])

    def test_unrestricted_exterior_slow_box_can_add_earlier_label(self):
        star=16;Q0=Q(1,2**star);q=Q0/4
        self.assertTrue(d.coverage(star,1,q,star)['accepted'])
        qoutside=4*Q0;earlier=2*Q0
        self.assertLessEqual(qoutside/earlier,4)
        self.assertGreater(qoutside,q)
        with self.assertRaises(ValueError):d.coverage(star,1,qoutside,star-1)

    def test_shifted_radial_integral_keeps_q_and_shift(self):
        r,H,q,M,a=s.symbols('r H q M a')
        F=q+(r+a)*(H+M*a)
        value=s.integrate(F,(a,0,1))
        self.assertEqual(s.diff(value,q),1)
        self.assertNotEqual(s.expand(value-value.subs(M,0)),0)
        self.assertEqual(s.expand(value),q+r*H+H/2+M*r/2+M/3)

    def test_repeated_fiber_maps_do_not_import_other_q(self):
        # Discrete exact analogue: nonlocal r/H sums, with q fixed.
        data={(q,r,h):Q(q+r+h) for q in (1,2,3) for r in (0,1) for h in (0,1)}
        other=dict(data)
        for key in other:
            if key[0]>1:other[key]+=10**10
        def op(f):return {(q,r,h):sum(f[q,R,H] for R in (0,1) for H in (0,1))+f[q,r,h]**2 for q,r,h in f}
        for _ in range(4):
            data,other=op(data),op(other)
            for r in (0,1):
                for h in (0,1):self.assertEqual(data[1,r,h],other[1,r,h])

    def test_radial_primitive_fills_gap_at_same_q(self):
        q=s.symbols('q');source=[q,0,-q]
        primitive=[sum(source[:i+1]) for i in range(3)]
        self.assertEqual(source[1],0);self.assertEqual(primitive[1],q)
        self.assertEqual(primitive[-1],0)

    def test_physical_cutoff_has_nonzero_commutator(self):
        t=s.symbols('t');q=1-t;f=1+t*t;chi=q*q
        comm=s.expand(s.diff(chi*f,t)-chi*s.diff(f,t))
        self.assertEqual(comm,s.expand(-2*q*f));self.assertNotEqual(comm,0)

    def test_whole_moving_pressure_product(self):
        z=s.symbols('z');rho=1+z*z;P=1+z
        self.assertEqual(s.diff(rho*P,z)-rho*s.diff(P,z),2*z*P)

    def test_guards(self):
        with self.assertRaises(TypeError):d.coverage(16,1,0.001,16)
        with self.assertRaises(ValueError):d.coverage(16,1,Q(1,2**10),16)

    def test_evidence(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['source_norm_proved'])
        self.assertFalse(r['symbolic']['historical_localization_identified'])
if __name__=='__main__':unittest.main()
