from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from controls import fixture,remaining
from bounds import MovingCompensator
class PressureTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_dilation_and_composition(self):self.control('dilation_has');self.control('full_mixed_compensator')
 def test_actual_implicit_coordinate_local_jets(self):self.control('actual_implicit_coordinate');self.control('Z_derivative')
 def test_signed_mass_and_moving_second_moment(self):self.control('signed_fixed');self.control('moment_bound_uses');self.control('exact_moving_pressure')
 def test_new_pole_moment_weights(self):self.control('radial_weight');self.control('integrated_slow_jets')
 def test_derived_inputs_routed_to_complete_pressure(self):self.control('three_independent');self.control('complete_pressure')
 def test_whole_current_pressure_memory(self):self.control('whole_current_pressure')
 def test_moving_support_and_scope(self):self.control('moving_support');self.control('gaps_remain')
 def test_missing_next_profile_jet_and_impossible_unit_mass_rejected(self):
  v=fixture()
  with self.assertRaises(ValueError):MovingCompensator(v.update,v.hull,[1])
  with self.assertRaises(ValueError):MovingCompensator(v.update,v.hull,[0,0])
  with self.assertRaises(ValueError):v.rho(5)
  with self.assertRaises(ValueError):v.moment(-1)
  with self.assertRaises(TypeError):MovingCompensator(object(),v.hull,[1,1])
 def test_same_width_and_whole_support_required(self):
  v=fixture()
  with self.assertRaises(ValueError):MovingCompensator(v.update,dict(v.hull,t1=v.hull['t1']/2),[1,1])
  with self.assertRaises(ValueError):MovingCompensator(v.update,dict(v.hull,l=1),[1,1])
 def test_remaining_map_coverage_and_keys_rejected(self):
  v=fixture()
  for f in (lambda n,p:dict(remaining(n,p),order=0),lambda n,p:dict(remaining(n,p),pole=0),lambda n,p:dict(remaining(n,p),Kmom=1)):
   with self.assertRaises(ValueError):v.pressure(0,f)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
