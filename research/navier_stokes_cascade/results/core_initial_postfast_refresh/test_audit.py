from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import InitialPostfastRefresh,complete_pressure_envelope,sharp_slow_fields,a1_time_error_request,EP
from controls import fixture

class InitialPostfastRefreshTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_complete_pressure_memory(self):self.control('five_sequential');self.control('not_reset')
 def test_scalar_refresh_without_double_count(self):self.control('exactly_once');self.control('double_counting');self.control('distinct_actual_axial')
 def test_sharp_slow_response(self):self.control('angular_and_axial');self.control('radial_field')
 def test_new_A1_request(self):self.control('least_strict');self.control('sharper_but');self.control('source_state');self.control('dependencies');self.control('does_not_reach')
 def test_rates_and_limits(self):self.control('first_real_rates');self.control('zero_strict_margin');self.control('no_complete')
 def test_domain_and_type_guards(self):
  with self.assertRaises(TypeError):InitialPostfastRefresh(object())
  with self.assertRaises(ValueError):sharp_slow_fields(0,dict(P=1),dict(P=1),{}, {})
  with self.assertRaises(ValueError):a1_time_error_request(0,0)
  with self.assertRaises(ValueError):complete_pressure_envelope(1,2,3,4,5)['pm'].rational_upper(-1)
 def test_short_fixture_rejects_uncovered_refresh(self):
  v=InitialPostfastRefresh(fixture())
  with self.assertRaises(ValueError):v.pressure_envelope(0)
  with self.assertRaises(ValueError):v.scalars(0)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
