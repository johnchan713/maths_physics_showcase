import json,tempfile,unittest
from fractions import Fraction as Q
from pathlib import Path
import q143 as q
import audit

class Tests(unittest.TestCase):
 def short(self,P=16):
  c=q.previous.selected_caps(2,4,P,Q(1,100),Q(1,2**20),10,100,q.previous.IDENTITY)
  return q.previous.build(c,Q(1,100),Q(1,2**20),q.previous.IDENTITY)
 def recipe(self,P=16):
  short=self.short(P);args=(Q(1,100),Q(1,2**20),110,4,10,10,q.IDENTITY)
  src=q.reference_sources(short,*args)
  return q.post_taper(short,src,*args)
 def test_exact_mode_sources_and_fixed_recurrence(self):
  r=self.recipe();self.assertEqual(r['reference_sources']['p1']['A'],1)
  self.assertEqual(r['reference_sources']['Xns']['A'],2)
  self.assertTrue(r['reference_sources']['all_five_reference_rows_reintegrated'])
  self.assertTrue(r['reference_sources']['full_reference_I_J_Pi_and_Pi_eta_memory_retained'])
  self.assertFalse(r['field']['order_dependent_norm_exponentiated'])
  self.assertFalse(r['inherited_exp_Tm_squared_cap_used'])
 def test_actual_memory_and_complete_pressure(self):
  a=self.recipe(16);b=self.recipe(32)
  self.assertTrue(a['all_five_actual_incoming_rows_retained'])
  self.assertTrue(a['full_p2_Pi_and_Pi_eta_retained'])
  self.assertGreater(b['state']['p']['P'],a['state']['p']['P'])
 def test_floor_and_order_zero_cap(self):
  r=self.recipe();self.assertEqual(r['a_floor'],Q(1,10*2**20))
  self.assertEqual(r['field']['exp_value_cap'],3)
  self.assertTrue(r['selected_order_zero_floor_and_exp_value_cap_certified'])
 def test_old_cap_obstruction_and_new_binary_budget(self):
  r=self.recipe();old=q.old_q73_log_lower(Q(1,100),8)
  new=q.bound(r,'p',8,8,True)['upper']
  self.assertGreater(old,new)
  self.assertIsInstance(q.bound(r,'p',100005,100005,True)['upper'],int)
 def test_guards(self):
  short=self.short();args=(Q(1,100),Q(1,2**20),110,4,10,10,q.IDENTITY)
  with self.assertRaises(ValueError):q.reference_sources(short,*args[:-1],'wrong')
  src=q.reference_sources(short,*args)
  with self.assertRaises(ValueError):q.post_taper(short,src,Q(1,100),Q(1,2**10),110,4,10,10,q.IDENTITY)
  with self.assertRaises(ValueError):q.reference_sources(short,Q(1,100),Q(1,2**20),110,4,4,10,q.IDENTITY)
  with self.assertRaises(ValueError):q.post_taper(self.short(32),src,*args)
  with self.assertRaises(ValueError):q.bound({},'p',1,1)
 def test_00_record(self):
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
