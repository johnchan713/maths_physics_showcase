import json,unittest
from pathlib import Path
from copy import deepcopy
import audit
from bounds import SharpInitial
from controls import fixture
class SharpTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_initial_plus_fast_squares(self):self.control('fresh_scalar');self.control('quadratic_not_zero')
 def test_exact_fast_and_slow_formulas(self):self.control('exact_historical')
 def test_scalar_and_moment_input_orders(self):self.control('covariance_order')
 def test_moving_whole_pressure(self):self.control('whole_axial_moving');self.control('moment_not_deletable')
 def test_weighted_reserved_support(self):self.control('weight_floor_inserted')
 def test_stronger_rates_use_squares_not_linear_promotion(self):self.control('squares_supply');self.control('order_three')
 def test_all_inherited_exact_pressure_identities(self):self.control('pressure:')
 def test_initial_zero_Haar_reduction_and_later_failure(self):self.control('prefix:')
 def test_later_or_residual_claim_rejected(self):self.control('no_later_prefix')
 def test_order_and_hull_contracts(self):
  a=fixture()
  with self.assertRaises(ValueError):a.output(2)
  with self.assertRaises(ValueError):a.axial_moment(1,1)
  with self.assertRaises(ValueError):a.preslow(5)
  with self.assertRaises(ValueError):a.slow_fields(4)
  with self.assertRaises(ValueError):SharpInitial(a.budget,0,(0,2))
  with self.assertRaises(ValueError):SharpInitial(a.budget,0,(1,3))
  with self.assertRaises(TypeError):SharpInitial(object(),0,(1,2))
 def test_frozen_sources(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering_rejected(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
