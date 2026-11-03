import unittest,json
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
import audit,q120

class Tests(unittest.TestCase):
 def test_closed_bounds_dominate_independent_old_recurrences(self):
  old=q120.borel.step_bounds(10)
  for j in range(11):
   self.assertLessEqual(q120.borel.gaussian_raw(j),12**j*factorial(j)**3)
   self.assertLessEqual(old['reciprocal'][j],81*2048**j*factorial(j)**3)
   self.assertLessEqual(old['coarse'][j],82*2048**j*factorial(j)**3)
  for m in range(1,11):
   self.assertLessEqual(q120.borel.cutoff_bounds(m)['K'][m],82*4097**m*factorial(m)**3)
   self.assertLessEqual(q120.implemented_bound(m),2**q120.cutoff_binary_upper(m))

 def test_reciprocal_induction_and_ordered_partition_bounds(self):
  self.assertLess(Q(162)*Q(12,2048)/(1-Q(12,2048)),1)
  self.assertLess(Q(81)/(1-Q(12,2048)),82)
  for m in range(1,13):
   for j in range(1,m+1):
    S=q120.borel.stirling(m)[j]
    self.assertLessEqual(S*factorial(j),factorial(m)*comb(m-1,j-1))
   self.assertEqual(sum(comb(m-1,j-1)*4096**j for j in range(1,m+1)),4096*4097**(m-1))

 def test_fixed_function_low_certificates_and_claim_scope(self):
  self.assertEqual([q120.implemented_bound(j) for j in range(3)],[1,64,4928])
  self.assertTrue(q120.claim_bracket(2,3)['impossible_by_midpoint'])
  self.assertTrue(q120.claim_bracket(2,100)['undecided_by_these_bounds'])
  self.assertTrue(q120.claim_bracket(2,4928)['certified_by_upper'])
  self.assertFalse(q120.claim_bracket(2,4928)['actual_field_identity_certified'])
  with self.assertRaises(ValueError):q120.cutoff_provider('compact-X')
  with self.assertRaises(ValueError):q120.cutoff_provider('variance')

 def test_retained_constant_slack_transfer_without_scale_selection(self):
  for h,n,m,M,oldK in ((Q(1,200),4,3,Q(7),Q(100)),(Q(1,200),2,1,Q(0),Q(6000)),(Q(1,200),1,0,Q(2),Q(10000))):
   first=q120.retained_row(h,n,m,M,oldK,0)
   self.assertLessEqual(first['new_constant'],first['inflation']*first['old_constant'])
   self.assertFalse(first['sufficient_direct_new_row_gate'])
   bits=q120.borel.ceil_log2(first['old_constant'])+first['extra_binary_slack']
   supplied=Q(n+bits,n*h)
   r=q120.retained_row(h,n,m,M,oldK,supplied)
   self.assertTrue(r['sufficient_constant_transfer_gate']);self.assertTrue(r['sufficient_direct_new_row_gate'])
   self.assertEqual(r['supplied_log2_lower'],supplied)
   self.assertFalse(r['scales_selected_or_changed']);self.assertFalse(r['infinite_schedule_certified'])

 def test_deferred_large_row_uses_only_binary_budgets(self):
  r=q120.deferred_row(Q(1,200),100010,100010,123456,0)
  self.assertFalse(r['sufficient_new_row_gate']);self.assertFalse(r['factorials_materialized'])
  self.assertGreater(r['normalized_constant_binary_upper'],123456)
  self.assertEqual(q120.cutoff_recipe(100011)['factorial_argument'],100011)
  with self.assertRaises(ValueError):q120.materialize(q120.cutoff_recipe(100011))
  for n,m,Mbits in ((2,1,10),(5,4,20)):
   r=q120.deferred_row(Q(1,200),n,m,Mbits,0)
   C=max(Q(1),2*q120.implemented_bound(m+1)*(1+2*n*Q(1,200))**m*2**Mbits)
   self.assertLessEqual(C,2**r['normalized_constant_binary_upper'])

 def test_current_background_accepts_identified_cutoff_norm_provider(self):
  current=q120.current;trace=set()
  b=current.ClosedBackground(Q(1,200),lambda m:Q(m+2),lambda n,m:Q(n+m+1),
       q120.cutoff_provider(q120.IDENTITY),1,8,'euler',trace)
  old=current.base.BackgroundBudget(b.h,b.leading,b.profiles,b.cutoff,1,8,'euler')
  for m in range(3):self.assertEqual(b.bound(m),old.bound(m))
  v=b.bound(4)
  expected=current.base.increment_constant(b.h,1,4,Q(6),q120.implemented_bound(5))
  self.assertEqual(v['head']['head'],sum(current.base.increment_constant(b.h,n,4,Q(n+5),q120.implemented_bound(5)) for n in range(1,4)))
  self.assertGreater(expected,1);self.assertIn(('cutoff',5),trace)
  self.assertEqual(v['chart']['kind'],'closed-real-slow-chart-majorant')

 def test_invalid_rows_and_inputs(self):
  for args in ((Q(0),1,0,1,1,0),(Q(1,200),1,2,1,1,0),(Q(1,200),1,0,-1,1,0)):
   with self.assertRaises(ValueError):q120.retained_row(*args)
  with self.assertRaises(ValueError):q120.cutoff_recipe(True)
  with self.assertRaises(ValueError):q120.deferred_row(Q(1,200),2,1,True,0)

 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
