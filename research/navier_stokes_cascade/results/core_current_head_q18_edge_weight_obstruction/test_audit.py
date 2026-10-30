from fractions import Fraction as Q
import json
import unittest
import audit
import weights as w

class WeightTests(unittest.TestCase):
    def test_source_sum_retains_half_weight(self):
        self.assertEqual(w.sum_weight(Q(1,2),w.product(Q(1,2),Q(1,2)),w.product(1,Q(1,2))),Q(1,2))
        self.assertFalse(w.embeds(Q(1,2),1))
    def test_covariance_has_full_weight(self):
        self.assertEqual(w.product(Q(1,2),Q(1,2)),1)
        self.assertTrue(w.embeds(1,Q(1,2)))
    def test_signed_target_threshold(self):
        self.assertEqual(w.signed_quotient(1),Q(1,2))
        self.assertEqual(w.signed_quotient(Q(1,2)),0)
        self.assertLess(w.signed_quotient(Q(1,4)),0)
    def test_divergent_lower_bound_all_finite_poles(self):
        for p in (0,1,100,400009,10**30):
            for n in (max(1,4*p),max(1,4*p)+100):
                self.assertGreaterEqual(w.radial_log_lower(n,p),Q(n*n,4))
    def test_finite_band_constant_cannot_repair_edge(self):
        for constant in (1,10**5,10**40):
            n=4*constant+4
            self.assertGreater(w.radial_log_lower(n,constant)-constant,n)
    def test_zero_harmonic_mean_not_zero_source(self):
        # Laurent coefficients of cos(theta): mean zero, value at theta=0 is one.
        f={-1:Q(1,2),1:Q(1,2)}
        self.assertEqual(f.get(0,0),0);self.assertEqual(sum(f.values()),1)
        square={}
        for k,a in f.items():
            for l,b in f.items(): square[k+l]=square.get(k+l,0)+a*b
        self.assertEqual(square[0],Q(1,2))
    def test_guards(self):
        for x in (-1,0.5,True):
            with self.assertRaises(ValueError):w.weight(x)
        with self.assertRaises(ValueError):w.sum_weight()
    def test_provenance_scope(self):
        r=audit.generate()
        self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertFalse(r['symbolic']['actual_source_violates_full_weight_bound_proved'])
        self.assertFalse(r['symbolic']['typed_all_stage_induction_completed'])

if __name__=='__main__':unittest.main()
