from fractions import Fraction as Q
import unittest,json
import audit
import typed as t
class TypedTests(unittest.TestCase):
    def test_cycle_invariant(self):
        wave,mean=Q(1,2),Q(1)
        for _ in range(100):
            r=t.cycle(wave,mean);wave,mean=r['wave'],r['mean']
            self.assertEqual(wave,Q(1,2));self.assertEqual(mean,1)
            self.assertEqual(r['signed_target'],1)
            self.assertEqual(r['harmonic_source'],Q(1,2))
    def test_half_target_not_accepted(self):
        with self.assertRaises(ValueError):t.cycle(target_floor=Q(1,2))
        with self.assertRaises(ValueError):t.cycle(mean=Q(1,2))
        with self.assertRaises(ValueError):t.cycle(stress=0)
    def test_angular_linear_projection(self):
        f={-2:Q(3),1:Q(4),3:Q(-2)}
        for _ in range(5):
            self.assertEqual(f.get(0,0),0);f=t.derivative(f)
    def test_opposite_harmonics_survive(self):
        f={-1:Q(1,2),1:Q(1,2)}
        self.assertEqual(t.multiply(f,f)[0],Q(1,2))
        self.assertNotEqual(t.multiply(f,f)[0],f.get(0,0)**2)
    def test_complete_covariance_both_crosses(self):
        u={-1:Q(1),1:Q(1)};v={-1:Q(2),1:Q(2)}
        total={k:u[k]+v[k] for k in u}
        delta=t.multiply(total,total)[0]-t.multiply(u,u)[0]
        self.assertEqual(delta,2*t.multiply(u,v)[0]+t.multiply(v,v)[0])
        self.assertNotEqual(delta,t.multiply(u,v)[0]+t.multiply(v,v)[0])
    def test_angular_and_auxiliary_projections_differ(self):
        # mode index is (angular,auxiliary); one survives Haar, the other angular mean.
        f={(1,0):Q(2),(0,1):Q(3)}
        self.assertEqual(sum(v for (a,b),v in f.items() if b==0),2)
        self.assertEqual(sum(v for (a,b),v in f.items() if a==0),3)
    def test_wave_pressure_not_promoted(self):
        r=t.cycle();self.assertEqual(r['mean_pressure'],1)
        self.assertEqual(r['harmonic_pressure'],Q(1,2))
    def test_evidence_scope(self):
        r=audit.generate();self.assertEqual(r,json.loads((audit.HERE/'evidence.json').read_text()))
        self.assertTrue(r['symbolic']['conditional_radial_cycle_closed'])
        self.assertFalse(r['symbolic']['actual_inverse_root_and_rate_contracts_certified'])
if __name__=='__main__':unittest.main()
