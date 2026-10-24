import json,unittest
from copy import deepcopy
from pathlib import Path
import audit
from controls import fixture,data
from bounds import P,Prefast
class SourceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v);self.assertTrue(all(v))
 def test_all_complete_historical_rows(self):self.control('historical_zero_mean')
 def test_whole_pressure(self):self.control('whole_pressure')
 def test_no_Haar_rate_promotion(self):self.control('zero_Haar');self.control('graph_derivative')
 def test_derived_graph_constants(self):self.control('graph_and_radius')
 def test_full_source_jet_routing(self):self.control('both_covariance_orders');self.control('raw_tangential');self.control('raw_radial')
 def test_prefast_scope(self):self.control('postfast_claim')
 def test_missing_inputs_and_small_pole_rejected(self):
  b=fixture()
  with self.assertRaises(ValueError):b.rows(0,{},54)
  with self.assertRaises(ValueError):b.rows(0,data(3),53)
  with self.assertRaises(TypeError):Prefast(object())
 def test_record_and_hashes(self):
  audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering_rejected(self):
  for k in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[k]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
