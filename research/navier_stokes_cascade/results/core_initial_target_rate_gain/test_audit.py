from copy import deepcopy
from pathlib import Path
from fractions import Fraction as Q
import json,unittest
import audit
from rate_gain import divisor_data,exponents,final_row,shell_inverse_bound,target_gain


class InitialTargetRateGainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=audit.generate()
    def control(self,key):
        found=[v for name,v in self.r['review']['checks'].items() if key in name]
        self.assertTrue(found,key);self.assertTrue(all(found),key)
    def test_profile_gate(self):self.control('finite_derivative');self.control('zero_collars')
    def test_target_gain(self):
        self.control('Haar_rows');self.control('moment_bumps');self.control('compact_target');self.control('principal_and_curl');self.control('pair_remainder')
    def test_source_and_final_margin(self):
        self.control('graph_loss');self.control('fast_source');self.control('radial_viscosity');self.control('A1_bound');self.control('strict_kappa')
    def test_guards(self):
        for mode in ((0,0),(1.0,0)):
            with self.assertRaises(ValueError):divisor_data(*mode)
        for bad in (0,-1,1):
            with self.assertRaises(ValueError):exponents(bad)
        with self.assertRaises(ValueError):shell_inverse_bound(6)
    def test_exact_rate_values(self):
        k=Q(1,100000);g=target_gain(k);f=final_row(k)
        self.assertEqual(g['target_gain'],k);self.assertEqual(f['strict_margin'],k)
    def test_frozen_record(self):
        audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
    def test_tampering(self):
        for key in ('flags','review','source_hashes'):
            r=deepcopy(self.r);r[key]['tamper']=True
            with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
