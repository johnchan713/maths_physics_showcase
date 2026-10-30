from fractions import Fraction as Q
from math import ceil, floor, comb
import json
import unittest
import audit
import support as s

class SupportTests(unittest.TestCase):
    def test_projected_count_majorant(self):
        for S in (1,16**2,10**12):
            for length in (0,Q(1,7),1,Q(19,3)):
                self.assertLessEqual(s.projected_count(S,length),s.count_constant(length)*S**3)
        self.assertEqual(s.projected_count(1,0),54*2744*3)
    def test_mesh_count_including_closed_boundaries(self):
        # A mesh support [k-1,k+1] meets [a,a+l] iff k in [a-1,a+l+1].
        for a in (Q(0),Q(1,7),Q(-4,3)):
            for length in (Q(0),Q(1,7),Q(2),Q(19,3)):
                actual=floor(a+length+1)-ceil(a-1)+1
                self.assertLessEqual(actual,ceil(length)+3)
        self.assertEqual(floor(3)-ceil(-1)+1,5) # three-point endpoints matter
    def test_aggregate_polynomial(self):
        native=[2,Q(1,3),5]
        aggregate=s.aggregate_polynomial(native,Q(7,2))
        for S in (1,256,10**8):
            self.assertLessEqual(s.projected_count(S,Q(7,2))*s.evaluate(native,S),s.evaluate(aggregate,S))
        self.assertEqual(len(aggregate),len(native)+3)
    def test_filled_labels_can_add(self):
        # Disjoint input radial bumps of unit integral produce equal scalar moments.
        # A common reserved profile times those moments has N-fold size.
        for N in (2,17,100):
            moments=[Q(1)]*N
            reserved_value=Q(3,7)
            self.assertEqual(sum(m*reserved_value for m in moments),N*reserved_value)
            self.assertGreater(sum(m*reserved_value for m in moments),reserved_value)
    def test_leibniz_is_exact_for_exponentials(self):
        for m in range(8):
            a,b=Q(2),Q(3)
            self.assertEqual(s.product_jet_bound([a**k for k in range(m+1)],[b**k for k in range(m+1)]),(a+b)**m)
        self.assertEqual(s.product_jet_bound([1,2],[3,5]),11)
    def test_compact_pressure_and_moment(self):
        # g=R(1-R)(2+R), rho=6R(1-R); p and p' vanish at both ends.
        g=[0,2,-1,-1];rho=[0,6,-6]
        P,p=s.pressure_fixture(g,rho)
        self.assertEqual(p[0],0);self.assertEqual(sum(p),0)
        self.assertNotEqual(P,0)
        self.assertEqual(s.integral(p,1),-s.integral(g,2)/2+s.integral(rho,2)*P/2)
        # The uncompensated primitive has a nonzero outgoing constant.
        self.assertNotEqual(s.integral(g),0)
    def test_gap_fill_by_zero_total_primitive(self):
        # Unit-mass 6u(1-u) on left quarter; its negative on right quarter.
        def bump_primitive(u): return 3*u*u-2*u*u*u
        left_mass=bump_primitive(Q(1))-bump_primitive(Q(0))
        right_mass=-left_mass
        self.assertEqual(left_mass+right_mass,0)
        # At the midpoint both sources vanish; the primitive equals the left mass.
        self.assertEqual(left_mass,1)
    def test_moving_pressure_product(self):
        # c(Z)=mu*(1+Z), P(Z)=1+Z^2, Z=2; frozen-c derivative misses mu*P.
        mu=Q(3,5);Z=Q(2)
        whole=mu*(1+2*Z+3*Z*Z)
        frozen=mu*(1+Z)*2*Z
        self.assertEqual(whole-frozen,mu*(1+Z*Z))
        self.assertNotEqual(whole,frozen)
    def test_log_degree_and_guards(self):
        for x in (1,2,100):
            self.assertLessEqual((16*x*x)**3,4096*(1+x)**6)
        for S in (0,-1,1.0,True):
            with self.assertRaises(ValueError):s.projected_count(S,1)
        with self.assertRaises(ValueError):s.aggregate_polynomial([],1)
        with self.assertRaises(ValueError):s.pressure_fixture([1],[2])
    def test_provenance_and_scope(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertTrue(r['scope']['projected_mean_pressure_count_derived'])
        self.assertFalse(r['scope']['actual_normalized_coefficient_estimates_proved'])
        self.assertFalse(r['scope']['unconditional_residual_or_blowup'])

if __name__=='__main__':unittest.main()
