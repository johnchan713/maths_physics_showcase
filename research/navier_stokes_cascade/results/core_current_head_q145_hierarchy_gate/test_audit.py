import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import q145 as q
import audit

class Tests(unittest.TestCase):
 def recipe(self):
  p=q.previous.previous
  c=p.previous.selected_caps(2,4,16,Q(1,100),Q(1,2**20),10,100,p.previous.IDENTITY)
  short=p.previous.build(c,Q(1,100),Q(1,2**20),p.previous.IDENTITY)
  args=(Q(1,100),Q(1,2**20),110,4,10,10,p.IDENTITY)
  post=p.post_taper(short,p.reference_sources(short,*args),*args)
  s=post['state'];E=max(Q(100),2*s['E']['P']);I=max(Q(100),2*s['invE']['P'])
  return q.previous.shape_annulus(post,Q(1,8),2,1,Q(1,100),10,E,I,200,16,2010,q.previous.IDENTITY)
 def test_exact_symbolic_hierarchy(self):
  h=q.hierarchy_inequalities();self.assertTrue(all(h['exact_checks'].values()))
  self.assertEqual(h['definitions']['N'],'1+floor(H^32)');self.assertFalse(h['hierarchy_materialized'])
  with self.assertRaises(ValueError):q.hierarchy_inequalities(q.A_MIN-1)
 def test_scalar_gate_reduction_keeps_pressure(self):
  r=q.reduce_gate(self.recipe(),q.IDENTITY)
  self.assertTrue(r['scalar_gate']['G_cap']);self.assertEqual(r['scalar_gate']['Eexp_choice'],'original H')
  self.assertTrue(all(r['pressure_memory'].values()));self.assertGreater(r['a_floor'],0)
  self.assertFalse(r['same_profile_complex_caps_certified']);self.assertFalse(r['actual_q140_q138_cone_join_certified'])
 def test_remaining_obstruction_is_profile_binding(self):
  r=q.reduce_gate(self.recipe(),q.IDENTITY)
  self.assertIn('same-profile Q73/Q143/Q127',r['remaining_inequality'])
  self.assertNotIn('Eexp>=exp(G)',r['remaining_inequality'])
 def test_guards_reject_changed_identity_and_pressure_reset(self):
  recipe=self.recipe()
  with self.assertRaises(ValueError):q.reduce_gate(recipe,'wrong')
  recipe=dict(recipe);recipe['full_p2_Pi_and_Pi_eta_retained']=False
  with self.assertRaises(ValueError):q.reduce_gate(recipe,q.IDENTITY)
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()

