from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import rates,requests,current_postfast,EP,Q
class RefreshTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_real_rates_and_zero_strict_margin(self):self.control('fine_error');self.control('zero_strict');self.control('positive_margin')
 def test_fresh_scalar_and_whole_pressure_rates(self):self.control('postfast_scalar');self.control('whole_pressure_orders')
 def test_compact_counts(self):self.control('least_strict');self.control('differ_from')
 def test_whole_request_graph_and_poles(self):self.control('requests_traced');self.control('high_jet');self.control('pressure_domain');self.control('graph_covered')
 def test_current_scalar_no_reset_or_duplicate(self):self.control('exactly_once');self.control('duplicate_wave');self.control('old_fast_correlations')
 def test_fresh_slow_response(self):self.control('current_pressure');self.control('next_scalar_order')
 def test_full_differences_and_actual_fast(self):self.control('viscosity_and_pressure');self.control('old_plus_actual');self.control('no_extra_Qcirc')
 def test_limits(self):self.control('feedback_does_not')
 def test_domains(self):
  for k in (0,-1,Q(1,2)):
   with self.assertRaises(ValueError):rates(k)
  with self.assertRaises(ValueError):requests(0,old_pole=-1,target_allowance=1,retained_wave_pole=1)
  with self.assertRaises(ValueError):current_postfast(0,dict(updated_scalars={}),{}, {},EP(1))
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
