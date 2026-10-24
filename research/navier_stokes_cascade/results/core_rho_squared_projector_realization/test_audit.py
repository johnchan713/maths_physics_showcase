from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import RhoSquaredProjectors
from controls import fixture
class ProjectorTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_signed_profile_and_normalization(self):self.control('not_assumed_positive');self.control('normalization_inverse');self.control('mass_one')
 def test_nonnegative_square_and_actual_support(self):self.control('comes_from_square');self.control('support_endpoints')
 def test_projector_algebra(self):self.control('idempotence')
 def test_actual_implicit_mixed_jets(self):self.control('implicit_native_mixed')
 def test_all_product_jets(self):self.control('product_jet_majorant')
 def test_floor_and_logarithmic_chart(self):self.control('single_moving_support');self.control('single_logarithmic')
 def test_historical_identity_and_scope(self):self.control('not_silently');self.control('no_complete')
 def test_domains(self):
  v=RhoSquaredProjectors(fixture())
  with self.assertRaises(ValueError):v.native(5,1)
  with self.assertRaises(ValueError):v.native(0,3)
  with self.assertRaises(TypeError):RhoSquaredProjectors(object())
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
