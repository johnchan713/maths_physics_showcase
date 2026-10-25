from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from profile_gate import actual_profile_implication,exponents,fourier_injectivity,mode_divisors


class InitialActualProfileGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=audit.generate()
    def control(self,key):
        found=[v for name,v in self.r['review']['checks'].items() if key in name]
        self.assertTrue(found,key);self.assertTrue(all(found),key)
    def test_irrational_mode_gate(self):
        self.control('radial_divisor');self.control('temporal_divisor');self.control('injective');self.control('Pell_norm')
    def test_actual_profile_distinction(self):
        self.control('signed_multiplier');self.control('strict_transverse_support');self.control('internal_graph');self.control('zero_target')
    def test_rate_scope(self):
        self.control('orders_match');self.control('not_promoted');self.control('remaining_strict_gate')
    def test_guards(self):
        for k in ((0,0),(1.0,0)):
            with self.assertRaises(ValueError):mode_divisors(*k)
        with self.assertRaises(ValueError):fourier_injectivity({(0,0):1})
        for k in (0,-1,1):
            with self.assertRaises(ValueError):exponents(k)
    def test_nonzero_target_is_a_required_branch_hypothesis(self):
        self.assertFalse(actual_profile_implication(False)['exact_profile_cancellation_excluded'])
    def test_frozen_record(self):
        audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
    def test_tampering(self):
        for key in ('flags','review','source_hashes'):
            r=deepcopy(self.r);r[key]['tamper']=True
            with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
