from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import InitialFinalRowObstruction,a1_phi_recipe,complete_final_pressure,strict_fast_closure
from controls import fixture

class InitialFinalRowObstructionTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_exact_final_row_and_pressure(self):self.control('adds_Haar');self.control('old_plus');self.control('pressure_and_viscosity');self.control('extra_Qcirc');self.control('all_five')
 def test_A1_raw_recipe(self):self.control('least_strict');self.control('raw_profiles');self.control('symbolic_not')
 def test_sharp_rate_obstruction(self):self.control('source_rate_bottleneck');self.control('radial_viscosity');self.control('sharp_slow');self.control('independent_frozen');self.control('requires_positive');self.control('either_positive')
 def test_scope(self):self.control('not_complete_residual')
 def test_domain_and_type_guards(self):
  with self.assertRaises(TypeError):InitialFinalRowObstruction(object())
  with self.assertRaises(ValueError):a1_phi_recipe(0,geometry=dict(L=0,width='1/2',gap='1'),pole=1)
  with self.assertRaises(ValueError):complete_final_pressure({},1)
  with self.assertRaises(ValueError):strict_fast_closure(-1,0)
 def test_short_fixture_rejects_uncovered_final_row(self):
  v=InitialFinalRowObstruction(fixture())
  with self.assertRaises(ValueError):v.final_rows(0,dict(Etheta=1,Ez=1),1)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
