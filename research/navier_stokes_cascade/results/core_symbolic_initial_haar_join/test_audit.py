import json,unittest
from pathlib import Path
from copy import deepcopy
import audit
from controls import fixture
from bounds import EP,InitialHaar,row_budget
class HaarTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_all_original_complete_rows(self):self.control('original_row_polynomials')
 def test_all_slots_and_sharp_weighted_fields(self):self.control('input_slots');self.control('sharp_slow');self.control('complete_fast')
 def test_complete_covariance_and_higher_stress(self):self.control('inserted_once')
 def test_pressure_input_and_covariance_orders(self):self.control('order_n_plus_one');self.control('n_plus_thirteen')
 def test_whole_pressure_changes_only_axial_bound(self):self.control('whole_pressure_defect');self.control('spuriously')
 def test_distinct_projector_channels(self):self.control('axial_projector');self.control('angular_projector')
 def test_projected_source_and_fixed_pullback_gate(self):self.control('row_plus_bump');self.control('not_silently_unit')
 def test_scope_and_rate(self):self.control('no_signed_root');self.control('original_order')
 def test_invalid_inputs(self):
  b=fixture()
  with self.assertRaises(ValueError):b.rows(1)
  with self.assertRaises(ValueError):row_budget(0,{})
  with self.assertRaises(TypeError):InitialHaar(object(),0,lambda n,p:1,lambda n:1)
  with self.assertRaises(ValueError):b.projected_source(0,3,pullback=lambda n,e:1,bump=lambda n,e:1)
 def test_frozen_sources(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
