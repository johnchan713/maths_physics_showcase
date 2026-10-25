from copy import deepcopy
from pathlib import Path
import json,unittest
import audit
from bounds import InitialQuotientFastInverse,fast_majorants,EP
from controls import fixture

class QuotientFastTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.r=audit.generate()
 def control(self,k):
  v=[v for n,v in self.r['review']['checks'].items() if k in n];self.assertTrue(v,k);self.assertTrue(all(v),k)
 def test_exact_radial_quotient_and_pressure_memory(self):
  self.control('gr_equals');self.control('Wrr_pressure');self.control('rhoP_disappears');self.control('full_pressure_complement')
 def test_positive_direct_tail_program(self):
  self.control('program_and_orders');self.control('compensated_pressure_memory');self.control('profile_recipe')
 def test_correction_of_high_jet_route(self):
  self.control('strictly_beats');self.control('200001_jet')
 def test_whole_source_preserves_regular_part(self):self.control('regular_source_is_unchanged')
 def test_actual_fast_map_formula(self):
  self.control('without_second_projection');self.control('gd_minus');self.control('output_orders')
 def test_limits_and_next_obligation(self):self.control('does_not_claim');self.control('remain_next')
 def test_domain_and_type_guards(self):
  with self.assertRaises(TypeError):InitialQuotientFastInverse(object())
  v=fixture();v.budget.compact=lambda n,p:dict(KT=EP(1),KA=EP(2))
  with self.assertRaises(ValueError):v.compact(0)
  v=fixture();v.budget.pole+=1
  with self.assertRaises(ValueError):v.source(0)
  with self.assertRaises(ValueError):fast_majorants(0,1,1,1,dict(KT=1,KA=1),dict(KT=1,KA=1,I=1))
 def test_short_fixture_rejects_uncovered_fast_output(self):
  with self.assertRaises(ValueError):fixture().fast_fields(0)
 def test_frozen_record(self):audit.validate(json.loads((Path(__file__).parent/'evidence.json').read_text()));self.assertTrue(all(self.r['gates'].values()))
 def test_tampering(self):
  for key in ('flags','review','source_hashes'):
   r=deepcopy(self.r);r[key]['tamper']=True
   with self.assertRaises(ValueError):audit.validate(r)
if __name__=='__main__':unittest.main()
