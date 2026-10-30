from fractions import Fraction as Q
import json
import unittest
import audit
import transfer as t

class TransferTests(unittest.TestCase):
    def test_corrected_rates_all_components(self):
        k=Q(1,100000);delta=Q(1,2)-3*k
        expected=[1,Q(3,2),Q(3,2)-k,Q(3,2)-k,Q(3,2)-k,Q(5,2)-k]
        for j in (1,2,10**30):
            r=t.rows(Q(1,513),j)
            self.assertEqual(list(r),list(t.KINDS))
            self.assertEqual([v['normalized_order']-delta*j for v in r.values()],expected)
    def test_pressure_is_smallest_margin(self):
        k=Q(1,100000)
        for h in (Q(1,513),Q(1,10**100),Q(99,10000)):
            for j,m in ((1,0),(2,3),(10**30,9)):
                rows=t.rows(h,j,m)
                margin=h*(Q(3,2)-k)+m*(1-5*h/2)
                self.assertEqual(min(r['margin'] for r in rows.values()),margin)
                self.assertEqual(rows['complete_fixed_gauge_pressure_difference']['margin'],margin)
                self.assertGreater(margin,0)
    def test_normalization_dictionary(self):
        h=Q(1,513);r=t.rows(h);A=Q(1,2)+h
        self.assertEqual(r['wave_potential']['normalization_loss'],h)
        self.assertEqual(r['direct_axisymmetric_swirl']['normalization_loss'],A)
        self.assertEqual(r['complete_fixed_gauge_pressure_difference']['normalization_loss'],2*A)
    def test_conservative_derivative_cost(self):
        for h in (Q(0),Q(1,513),Q(1,100)):
            nu=1+5*h/2
            self.assertGreaterEqual(nu,1+3*h/2)
            self.assertGreaterEqual(nu,Q(1,2)+h/2)
            self.assertLess(nu,2)
    def test_polynomial_log_bound(self):
        coeffs=[Q(3),Q(1,7),Q(5)];K=10**30
        C,P=t.polynomial_log_bound(coeffs,K)
        for x in (1,2,100):
            for S in (0,x*x,16*x*x):
                lhs=sum(c*(S+K)**i for i,c in enumerate(coeffs))
                self.assertLessEqual(lhs,C*(1+x)**P)
    def test_band_power_conversion(self):
        for E in (-19,-1,0,1,27):
            for ratio in (Q(1,4),Q(1),Q(4)):
                q=Q(1,100);band=q/ratio
                self.assertLessEqual(band**E,4**abs(E)*q**E)
    def test_guards(self):
        for h in (0,1,0.005,True):
            with self.assertRaises(ValueError):t.rows(h)
        with self.assertRaises(ValueError):t.rows(Q(1,513),0)
        with self.assertRaises(ValueError):t.polynomial_log_bound([1,-1],0)
    def test_provenance_scope(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['scope']['actual_complete_normalized_inputs_proved'])
        self.assertFalse(r['scope']['mean_pressure_overlap_verified'])

if __name__=='__main__':unittest.main()
