from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from witness import Mode,exponents,local_witness


class InitialViscositySharpnessWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=audit.generate()
    def control(self,key):
        found=[v for name,v in self.r['review']['checks'].items() if key in name]
        self.assertTrue(found,key);self.assertTrue(all(found),key)
    def test_covariance_to_source(self):self.control('principal_cross');self.control('source_order_D');self.control('connection_part')
    def test_fast_inverse_and_viscosity(self):self.control('temporal_inverse');self.control('viscosity_derivatives');self.control('laplacian_parts')
    def test_obstruction_scope(self):self.control('both_missing');self.control('not_a_consequence');self.control('special_identity');self.control('zero_Haar')
    def test_fourier_algebra(self):
        m=Mode(3,5);self.assertEqual(m.derivative(),Mode(5,-3));self.assertEqual(m.inverse_derivative().derivative(),m)
    def test_parameter_guard(self):
        for k in (0,-1,1):
            with self.assertRaises(ValueError):exponents(k)
    def test_frozen_record(self):
        audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
    def test_tampering(self):
        for key in ('flags','review','source_hashes'):
            r=deepcopy(self.r);r[key]['tamper']=True
            with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
