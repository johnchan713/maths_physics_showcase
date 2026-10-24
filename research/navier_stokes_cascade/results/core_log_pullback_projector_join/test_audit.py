from pathlib import Path
from copy import deepcopy
import json,unittest
import audit
from bounds import LogPullback,ProjectedSource
from controls import fixture
class PullbackTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_formal_exponential_and_composition(self):self.control('exponential_recurrence');self.control('composition_coefficients')
 def test_mixed_actual_coordinate_form(self):self.control('mixed_actual_form')
 def test_multiplier_and_nonlinear_jets(self):self.control('distinct_axial');self.control('second_jet')
 def test_normalized_positive_support_obstruction(self):self.control('positive_family');self.control('mass_one');self.control('grows_exactly');self.control('not_identically_zero')
 def test_actual_join_substitution(self):self.control('substituted');self.control('inserted_once')
 def test_full_pressure_and_scope(self):self.control('full_pressure')
 def test_invalid_hull_width_and_orders(self):
  v=fixture();p=v.haar.initial.budget.primary;h=dict(l=2,u=3,Xa=1,Xb=8,t1=p.primary.root.width)
  for change in (dict(t1=2*h['t1']),dict(l=1),dict(u=5,Xb=20)):
   with self.assertRaises(ValueError):ProjectedSource(v.haar,v.sigma_norm,dict(h,**change))
  with self.assertRaises(ValueError):LogPullback(0)
  with self.assertRaises(ValueError):v.pull.pullback(-1,1)
  with self.assertRaises(ValueError):v.pull.pullback(0,3)
  with self.assertRaises(TypeError):ProjectedSource(object(),v.sigma_norm,h)
 def test_frozen_sources(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
