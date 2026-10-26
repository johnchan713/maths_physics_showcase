from copy import deepcopy
from fractions import Fraction as F
import json,tempfile,unittest
from pathlib import Path
import audit
from borel_join import EDGES,NODES,geometric_tail,normalized_increment_bound,stage_margins,support_count,topological_order

class LeadingBackgroundBorelJoinTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.record=audit.generate()
 def test_route_is_acyclic(self):self.assertEqual(set(topological_order()),set(NODES))
 def test_cycle_is_rejected(self):
  with self.assertRaises(ValueError):topological_order(edges=EDGES+(("summed_background","core"),))
 def test_rule_margins(self):
  m=stage_margins(3,F(1,5),F(2,3),4,7,5,20);self.assertTrue(all(v>0 for v in m.values()))
 def test_normalized_increment(self):self.assertTrue(normalized_increment_bound(3,F(1,5),F(2,3),4,20)['bounded_by_2_minus_n'])
 def test_geometric_tail(self):
  self.assertEqual(geometric_tail(4,1),F(1,8))
  with self.assertRaises(ValueError):geometric_tail(0,1)
 def test_local_finiteness(self):
  self.assertEqual(support_count(F(1,20),(1,2,4,8,16,32)),5)
  with self.assertRaises(ValueError):support_count(F(1,20),(1,3,4))
 def test_scope_correction(self):
  f=self.record['flags'];self.assertTrue(f['leading_complex_analyticity_restricted_to_protected_collar']);self.assertFalse(f['global_heat_endpoint_claimed_complex_analytic'])
  self.assertTrue(f['rule_defined_borel_sum_identified_with_actual_background']);self.assertFalse(f['same_band_phase_normal_value_cone_contract_independently_closed_here']);self.assertFalse(f['unconditional_pulse_residual_force_or_blowup_proved'])
 def test_dependencies(self):
  self.assertEqual(len(self.record['dependency_evidence_sha256']),12);self.assertTrue(self.record['inherited_source_checks']);self.assertTrue(all(self.record['inherited_source_checks'].values()))
 def test_frozen_record_and_tampering(self):
  frozen=json.loads((audit.HERE/'evidence.json').read_text());audit.validate(frozen);bad=deepcopy(frozen);bad['flags']['unconditional_pulse_residual_force_or_blowup_proved']=True
  with self.assertRaises(ValueError):audit.validate(bad)
 def test_no_overwrite(self):
  with tempfile.TemporaryDirectory(dir=audit.HERE) as d:
   p=Path(d)/'kept.json';p.write_text('keep')
   with self.assertRaises(SystemExit):audit.main(['--output',str(p)])
   self.assertEqual(p.read_text(),'keep')
if __name__=='__main__':unittest.main()
