from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import FullCompactPressure,Q
from controls import fixture,bound_fixture,profiles
class PressureResponseTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_angular_and_auxiliary_averages_differ(self):self.control('not_auxiliary_Haar')
 def test_actual_shift_and_compact_A0_identity(self):self.control('coefficient_derivative');self.control('coefficient_derivatives');self.control('exact_A0');self.control('not_identically_zero');self.control('zero_Haar')
 def test_polynomial_diagnostic_endpoint_scope(self):self.control('extension_jets')
 def test_original_raw_profile_constants_and_five_variable_conversion(self):self.control('original_raw');self.control('next_cutoff');self.control('five_variable')
 def test_recomputed_full_response_and_compensation(self):self.control('complete_row');self.control('uncompensated');self.control('order_are_requested');self.control('changed_radial_row')
 def test_remaining_old_row_and_A0_obligations(self):self.control('remain_explicit')
 def test_Haar_only_pressure_and_missing_profile_rejected(self):
  v=bound_fixture()
  with self.assertRaises(ValueError):FullCompactPressure(v.state,1,v.pole,Q(1,10),v.state.maps(1,v.pole))
  with self.assertRaises(TypeError):FullCompactPressure(object(),1,v.pole,Q(1,10),profiles())
  with self.assertRaises(ValueError):FullCompactPressure(v.state,1,v.pole,Q(1,10),{})
 def test_full_pressure_order_pole_and_actual_gap_required(self):
  with self.assertRaises(ValueError):bound_fixture(order=0).covariance_response(0)
  with self.assertRaises(ValueError):bound_fixture(pole=0).covariance_response(0)
  state=fixture()
  with self.assertRaises(ValueError):FullCompactPressure(state,1,56,0,profiles())
  with self.assertRaises(ValueError):FullCompactPressure(state,1,56,100,profiles())
  with self.assertRaises(ValueError):FullCompactPressure(state,1,56,Q(1,10),dict(profiles(),f0=-1))
 def test_same_choice_identity_and_covariance_ledger(self):
  v=bound_fixture();v.state.update.target=v.state.family.compensator.update.target
  with self.assertRaises(ValueError):v.covariance_response(0)
  with self.assertRaises(ValueError):bound_fixture().covariance_response(1)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
