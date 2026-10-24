from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import ExpSum,weighted_maps,fast_coefficients,Q
class MapTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v);self.assertTrue(all(v))
 def test_historical_constants(self):self.control('historical_')
 def test_actual_plateau_gap(self):self.control('same_collar')
 def test_exact_arithmetic(self):self.control('exact_exponential')
 def test_giant_width_not_materialized(self):self.control('tiny_width');self.control('giant_diagnostic')
 def test_actual_shift(self):self.control('fixed_shift');self.control('shifted_graph')
 def test_axial_A1_retained(self):self.control('actual_axial')
 def test_radial_input_order(self):self.control('radial_uses');self.control('source_orders')
 def test_invalid_inputs(self):
  for kw in (dict(gap=0),dict(width=2),dict(p=-1)):
   args=dict(p=3,L=4,width=Q(1,2),gap=1,jacobian=2,KT_profile=3,KA_profile=4);args.update(kw)
   with self.assertRaises(ValueError):weighted_maps(**args)
  with self.assertRaises(ValueError):ExpSum([(-1,0)])
  with self.assertRaises(TypeError):ExpSum([(0.5,0)])
 def test_reject_tampered_record(self):
  for k in ('flags','source_hashes','review'):
   r=deepcopy(self.r);r[k]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
 def test_scope_and_five_sources(self):
  self.assertEqual(len(self.r['dependency_evidence_sha256']),5);self.assertTrue(all(self.r['gates'].values()));self.assertFalse(self.r['flags']['actual_source_slow_mean_full_residual_or_blowup_certified'])
if __name__=='__main__':unittest.main()
