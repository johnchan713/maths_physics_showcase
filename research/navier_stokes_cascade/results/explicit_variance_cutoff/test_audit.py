#!/usr/bin/env python3
"""Exact tests for the explicit source-admissible variance cutoff."""
from copy import deepcopy
import json,tempfile,unittest
from pathlib import Path
import audit
import cutoff
class ExplicitCutoffTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.record=audit.generate();cls.r=cls.record['reviews'];cls.s=cls.r['selection']
 def test_exact_derivative_majorants(self):self.assertEqual([x['integer_majorant'] for x in self.s['derivative_rows']],[1,4,80,3600,292320])
 def test_exact_term_counts(self):self.assertEqual([x['term_count'] for x in self.s['derivative_rows']],[1,2,5,10,18])
 def test_cutoff_definition_and_scale(self):self.assertIn('integral_0^t',self.s['cutoff']);self.assertEqual(self.s['scaled'],'zetaL(v)=K(8*(v-2)/deltaL-1)')
 def test_integral_lower_bound(self):self.assertEqual(self.s['integral_lower'],'exp(-9/2)/3');self.assertEqual(self.s['integral_reciprocal_upper'],729)
 def test_raw_derivative_bounds(self):self.assertEqual(self.s['raw_derivative_strict_upper_bounds'],[729,2916,58320,2624400,213101280])
 def test_Z4_Z5(self):self.assertEqual(self.s['Z4_strict_upper'],2624400);self.assertEqual(self.s['Z5_strict_upper'],213101280)
 def test_all_exact_checks(self):self.assertTrue(all(self.s['checks'].values()))
 def test_universal_minimum(self):self.assertTrue(all(cutoff.exact_review(2**16)['checks'].values()))
 def test_larger_exact_A(self):self.assertTrue(all(cutoff.exact_review(2**17)['checks'].values()))
 def test_invalid_A(self):
  for a in (True,2**16-1,2**16.):
   with self.assertRaises(ValueError):cutoff.exact_review(a)
 def test_correction_is_explicit(self):self.assertIs(self.record['flags']['inherited_claim_that_a_quantitative_actual_cutoff_was_already_fixed_corrected'],True)
 def test_selected_not_silent_replacement(self):self.assertIs(self.r['consequences']['source_previously_specified_only_an_admissible_class'],True);self.assertIs(self.r['consequences']['explicit_member_now_fixed_without_changing_deltaL_or_plateaus'],True)
 def test_both_gates_discharged(self):self.assertIs(self.r['consequences']['Z4_gate_and_actual_target_S_Ceta4_below_G'],True);self.assertIs(self.r['consequences']['Z5_gate_with_native_R5_A32'],True)
 def test_same_repair_gate(self):self.assertIs(self.record['flags']['same_frequency_compact_modulation_and_exact_five_moment_repair_gate_discharged'],True)
 def test_no_FINAL_or_global_promotion(self):
  for k in ('complete_FINAL_joint_leading_and_global_physical_scale_joins_proved','uniform_actual_pulses_complete_covariance_axial_pressure_means_residual_or_blowup_verified'):
   self.assertIs(self.record['flags'][k],False)
 def test_fixed137_history(self):self.assertEqual(len(self.record['inherited_provenance']['frozen_records']),137);self.assertEqual(self.record['inherited_provenance']['fixed_manifest_sha256'],audit.MANIFEST_SHA)
 def test_historical_bindings(self):self.assertTrue(self.record['inherited_provenance']['all_bindings_verified']);self.assertGreater(self.record['inherited_provenance']['binding_count'],2464)
 def test_exact_frozen_evidence(self):audit.validate(json.loads((audit.HERE/'evidence.json').read_text()))
 def test_false_scope_promotions_rejected(self):
  for k,v in audit.FLAGS.items():
   if not v:
    r=deepcopy(self.record);r['flags'][k]=True
    with self.assertRaises(ValueError):audit.validate(r)
 def test_zero_cannot_replace_false(self):
  r=deepcopy(self.record);r['flags']['complete_FINAL_joint_leading_and_global_physical_scale_joins_proved']=0
  with self.assertRaises(ValueError):audit.validate(r)
 def test_all_protocol_bound(self):
  for k in audit.protocol():
   r=deepcopy(self.record);r['protocol'][k]='erased'
   with self.assertRaises(ValueError):audit.validate(r)
 def test_remaining_obligations_bound(self):
  r=deepcopy(self.record);r['remaining_obligations'].pop(0)
  with self.assertRaises(ValueError):audit.validate(r)
 def test_fresh_replay_and_overwrite_refusal(self):
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'fresh.json';audit.main(['--verify-record',str(audit.HERE/'evidence.json'),'--output',str(out)]);audit.validate(json.loads(out.read_text()))
   with self.assertRaises(SystemExit):audit.main(['--output',str(out)])
if __name__=='__main__':unittest.main()
