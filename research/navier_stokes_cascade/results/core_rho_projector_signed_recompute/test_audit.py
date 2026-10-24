from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import RecomputedSigned,family_module
from controls import fixture
class RecomputeTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_chosen_source_and_single_log_bump(self):self.control('log_bump_inserted');self.control('row_and_moment');self.control('not_arbitrary')
 def test_new_signed_amplitude_and_original_root(self):self.control('not_reused');self.control('amplitude_recipe')
 def test_all_covariance_terms(self):self.control('squares_and_curl')
 def test_native_pressure_projectors_and_profile(self):self.control('not_logarithmic');self.control('same_choice_fixed')
 def test_full_pressure_scalars_and_memory(self):self.control('all_scalar_changes');self.control('added_not_reset');self.control('original_and_new')
 def test_scope_and_historical_objects(self):self.control('still_unchanged');self.control('not_continuum')
 def test_changed_profile_and_pressure_data_rejected(self):
  v=fixture();v.covariance(0);v.family.compensator.profile=tuple(x+1 for x in v.family.compensator.profile)
  with self.assertRaises(ValueError):v.pressure(0)
  v=fixture();v.compensator.profile=tuple(x+1 for x in v.compensator.profile)
  with self.assertRaises(ValueError):v.pressure(0)
 def test_choice_and_primary_swaps_rejected(self):
  v=fixture();v.update.target=v.family.compensator.update.target
  with self.assertRaises(ValueError):v.covariance(0)
  v=fixture();v.source.sigma_norm=lambda n,e:1
  with self.assertRaises(ValueError):v.covariance(0)
  v=fixture();v.target.root=object()
  with self.assertRaises(ValueError):v.covariance(0)
 def test_whole_finite_ledger_and_types(self):
  v=fixture();v.haar.max_output=2
  with self.assertRaises(ValueError):RecomputedSigned(v.haar,v.family,0)
  v=fixture();v.family.compensator.max_order=3
  with self.assertRaises(ValueError):RecomputedSigned(v.haar,v.family,0)
  v=fixture()
  with self.assertRaises(ValueError):v.pressure(1)
  with self.assertRaises(TypeError):RecomputedSigned(v.haar,object(),0)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
