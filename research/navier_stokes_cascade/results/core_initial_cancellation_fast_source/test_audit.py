from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import CurrentFastSource,prefast_zero_haar_rows
from controls import fixture
class SourceTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_exact_cancellation_and_projection(self):self.control('cancellation_identity');self.control('annihilated');self.control('tail_is_axial');self.control('deleting_original')
 def test_actual_old_fields_and_sharp_slow_path(self):self.control('not_silently');self.control('sharp_weighted')
 def test_full_differences_and_current_covariance(self):self.control('both_differences');self.control('viscosity_differences')
 def test_original_prefast_and_A1_source_program(self):self.control('source_program_is_explicit');self.control('generic_later_stage')
 def test_raw_A1_obligation_and_actual_axial_gap(self):self.control('strict_gain');self.control('not_truncated');self.control('existing_pole');self.control('gap_not_copied')
 def test_no_unjustified_fine_source_or_inverse(self):self.control('prevents_unconditional');self.control('falsely_certified');
 def test_new_and_old_request_ledgers(self):self.control('next_fast');self.control('nested_primary')
 def test_materialized_inverse_is_blocked(self):
  with self.assertRaises(NotImplementedError):fixture().inverse_source()
 def test_prefix_and_domain_changes_rejected(self):
  v=fixture();v.budget.pole+=1
  with self.assertRaises(ValueError):v.rows(0)
  v=fixture();v.pressure.order=0
  with self.assertRaises(ValueError):v.rows(0)
  v=fixture();v.initial.max_output=-1
  with self.assertRaises(ValueError):v.rows(0)
  with self.assertRaises(TypeError):CurrentFastSource(object(),1)
  v=fixture()
  with self.assertRaises(ValueError):CurrentFastSource(v.pressure,0)
  with self.assertRaises(ValueError):v.rows(1)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
