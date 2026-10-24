from pathlib import Path
from copy import deepcopy
import json,unittest
import audit
from bounds import SignedTarget
from controls import fixture
class SignedTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_shell_and_symbolic_primitive(self):self.control('same_shell_ratio');self.control('three_branch')
 def test_two_edge_derivative_orders(self):self.control('two_edge_derivative')
 def test_full_reciprocal_weight_jets(self):self.control('reciprocal_weight_jet')
 def test_nonlinear_native_composition(self):self.control('native_finite');self.control('native_log')
 def test_unchanged_signed_root(self):self.control('original_signed_root');self.control('root_object')
 def test_complete_curl_next_jet(self):self.control('signed_curl')
 def test_symbolic_chain_and_scope(self):self.control('finite_coefficients');self.control('scope_still_open')
 def test_derivative_and_type_contracts(self):
  b=fixture()
  with self.assertRaises(ValueError):b.native(3,1)
  with self.assertRaises(ValueError):b.native(0,3)
  with self.assertRaises(ValueError):b.signed_curl(2)
  with self.assertRaises(ValueError):SignedTarget(b.source,4)
  with self.assertRaises(TypeError):SignedTarget(object(),0)
 def test_frozen_sources(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
