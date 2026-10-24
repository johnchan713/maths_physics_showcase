from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from controls import fixture,operators
from bounds import RealCovariance,EP,signed
class CovarianceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_all_six_entries_and_signed_squares(self):self.control('all_six');self.control('not_zero')
 def test_real_phase_and_rate_limits(self):self.control('cosine_cross');self.control('rate_limiter');self.control('cannot_be_promoted')
 def test_common_atlas_and_one_reverse_transfer(self):self.control('atlas_physical');self.control('one_chart_transfer')
 def test_complete_change_and_pressure(self):self.control('contains_principal');self.control('complete_pressure_target')
 def test_new_operator_domain(self):self.control('domain_is_requested')
 def test_current_scalar_and_whole_pressure_additions(self):self.control('added_not_reset');self.control('moving_pressure');self.control('includes_old')
 def test_scope(self):self.control('no_later_complex')
 def test_stale_operator_or_missing_input_rejected(self):
  b=fixture()
  with self.assertRaises(ValueError):b.pressure(0,lambda n,p:dict(operators(n,p),pole=b.target.pole))
  with self.assertRaises(ValueError):b.pressure(0,lambda n,p:dict(operators(n,p),order=0))
  with self.assertRaises(ValueError):b.pressure(0,lambda n,p:{})
  with self.assertRaises(ValueError):b.pressure(1,operators)
  with self.assertRaises(ValueError):b.coefficients(3)
  with self.assertRaises(TypeError):RealCovariance(object(),0)
 def test_giant_new_pressure_constant_remains_symbolic(self):
  def op(n,p):return dict(operators(n,p),Kp=signed.maps.ExpSum([(2,10**12)]))
  v=fixture().pressure(0,op)['changes']['pm_next'];self.assertTrue(v.record())
  with self.assertRaises(ValueError):v.rational_upper(1)
 def test_frozen_sources(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
