import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import q146 as q
import audit

class Tests(unittest.TestCase):
 def post(self):
  p=q.previous.previous.previous
  c=p.previous.selected_caps(2,4,16,Q(1,100),Q(1,2**20),10,100,p.previous.IDENTITY)
  short=p.previous.build(c,Q(1,100),Q(1,2**20),p.previous.IDENTITY)
  args=(Q(1,100),Q(1,2**20),110,4,10,10,p.IDENTITY)
  return p.post_taper(short,p.reference_sources(short,*args),*args)
 def bound(self):
  post=self.post();s=post['state']
  caps=q.cap_tuple(Q(1,8),2,1,Q(1,100),1,1,10,200,16,2010,s['E']['P'],s['invE']['P'],q.IDENTITY)
  return q.bind(post,caps,q.IDENTITY),caps
 def test_exact_cap_specialization(self):
  b,c=self.bound();self.assertGreaterEqual(c['Ecap'],c['endpoint_E']);self.assertGreaterEqual(c['invEcap'],c['endpoint_invE'])
  self.assertGreaterEqual(c['Ecap'],4*c['P']);self.assertGreaterEqual(c['invEcap'],20/c['P'])
  self.assertTrue(b['same_profile_complex_cap_formula_certified']);self.assertFalse(b['exact_repair_identity_and_cone_join_certified'])
 def test_symbolic_selected_cap_and_guards(self):
  c=q.cap_tuple(Q(1,8),2,1,Q(1,100),10**6,1,10,200,16,2010,1,1,q.IDENTITY,False)
  self.assertFalse(c['selected_huge_integer_materialized']);self.assertIn('3^n',c['formulas']['Ecap'])
  with self.assertRaises(ValueError):q.cap_tuple(Q(1,8),2,1,Q(1,50),1,1,10,200,16,2010,1,1,q.IDENTITY)
 def test_repair_request_retains_all_memory(self):
  b,_=self.bound();r=q.repair_request(b,q.REPAIR_IDENTITY)
  self.assertEqual(r['rows'],q.ROWS);self.assertIn('Pi',r['baseline_required']);self.assertTrue(r['inverse_lambda_retained'])
  self.assertFalse(r['same_modulation_memory_certified']);self.assertFalse(r['exact_repair_and_inner_outer_cone_join_certified'])
 def test_stage_shape_cannot_certify_actual_baseline(self):
  s=q.validate_baseline_shape({k:(1,2) for k in q.ROWS});self.assertTrue(s['pressure_row_retained']);self.assertFalse(s['actual_function_inequalities_certified'])
  with self.assertRaises(ValueError):q.validate_baseline_shape({k:(1,2) for k in q.ROWS if k!='Pi'})
 def test_guards(self):
  b,_=self.bound()
  with self.assertRaises(ValueError):q.repair_request(b,'wrong')
  bad=dict(b);bad['recipe']=dict(b['recipe']);bad['recipe']['full_p2_Pi_and_Pi_eta_retained']=False
  with self.assertRaises(ValueError):q.repair_request(bad,q.REPAIR_IDENTITY)
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
