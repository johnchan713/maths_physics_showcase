import json,unittest
from copy import deepcopy
from pathlib import Path
import audit
from controls import fixture
from bounds import EP,FreshMeans,ExpSum
class FreshTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_original_rows_and_coefficient_ring(self):self.control('historical_complete');self.control('exact_')
 def test_constructed_fast_state_and_fresh_scalars(self):self.control('constructed_fast');self.control('fresh_pressure_and_scalars')
 def test_pressure_memory_and_A1(self):self.control('moving_pressure');self.control('axial_A1')
 def test_fresh_slow_response(self):self.control('response_uses_fresh');self.control('next_fresh');self.control('stale_angular')
 def test_old_mean_not_reset(self):self.control('not_reset')
 def test_pole_and_jet_ledger(self):self.control('doubled_pressure');self.control('covariance_order_twelve')
 def test_moving_hull_symbolic_floor(self):self.control('hull_floor')
 def test_scope(self):self.control('no_later_signed');self.control('coarse_slow');self.assertTrue(all(self.r['gates'].values()))
 def test_declared_order_rejected(self):
  b,_=fixture()
  with self.assertRaises(ValueError):b.slow_increment(1)
  with self.assertRaises(ValueError):b.rows(2)
 def test_wrong_hull_or_changed_width_rejected(self):
  b,_=fixture();kw=dict(l=2,u=3,Xa=1,Xb=8,t1=b.primary.primary.root.width)
  for change in (dict(Xa=2),dict(t1=2*kw['t1'])):
   with self.assertRaises(ValueError):FreshMeans(b.primary,b.fixed,b.compact,b.slow,dict(kw,**change),0,162)
 def test_invalid_positive_ring(self):
  for v in (-1,{0:-1},{-1:2}):
   with self.assertRaises(ValueError):EP(v)
  with self.assertRaises(ValueError):EP(1).rational_upper(-1)
 def test_tampering_and_frozen_sources(self):
  audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()))
  for key in ('review','source_hashes','flags'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
