import unittest,importlib.util
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
import q125 as q
import audit
class Tests(unittest.TestCase):
 def test_uniform_reciprocal_at_all_checked_orders(self):
  import sympy as s
  z=s.symbols('z')
  for h in (Q(1,200),Q(99,10000)):
   g=1/(1-2*s.Rational(h.numerator,h.denominator)*z*z)
   for m in range(13):
    norm=sum(abs(s.diff(g,z,j).subs(z,1))/s.factorial(j) for j in range(m+1))
    self.assertLessEqual(norm,1/(1-8*s.Rational(h.numerator,h.denominator)))
 def test_independent_full_pressure_and_moments_with_tighter_caps(self):
  p=q.HERE.parent/'core_current_head_q124_full_pressure_discrepancies/test_audit.py'
  import sys
  previous=sys.modules.get('q124');sys.modules['q124']=q.current
  try:
   sp=importlib.util.spec_from_file_location('q125_full_source_control',p);t=importlib.util.module_from_spec(sp);sp.loader.exec_module(t)
  finally:
   if previous is None:sys.modules.pop('q124',None)
   else:sys.modules['q124']=previous
  t.q=q
  t.Tests('test_independent_full_order_two_pressure_and_five_moments').test_independent_full_order_two_pressure_and_five_moments()
 def test_tighter_than_historical_and_binary_bounds(self):
  for n,m in ((1,0),(2,3),(4,7)):
   U=[Q(i+2) for i in range(n)];F=[Q(i+3) for i in range(n)];prev=None if n==1 else Q(5)
   new=q.bounds(Q(1,200),n,m,1,1,2,3,U,F,prev,7,11)
   old=q.current.bounds(Q(1,200),n,m,1,1,2,3,U,F,prev,7,11)
   for k in ('I0','I1','axis'):self.assertLessEqual(new[k],old[k])
   for k in new['discrepancies']:self.assertLessEqual(new['discrepancies'][k],old['discrepancies'][k])
   bits=q.bounds(Q(1,200),n,m,1,1,2,3,[2]*n,[3]*n,None if n==1 else 3,3,4,binary=True)
   exact=q.bounds(Q(1,200),n,m,1,1,2,3,[4]*n,[8]*n,None if n==1 else 8,8,16)
   for k,v in q.raw_discrepancies(exact).items():self.assertLessEqual(v,2**bits['raw_discrepancy_binary_upper'][k])
  huge=q.bounds(Q(1,200),2,100005,1,1,2,3,[2,3],[5,7],11,13,17)
  self.assertFalse(huge['factorials_materialized'])
  with self.assertRaises(ValueError):q.raw_discrepancies(huge)
 def test_first_seed_recovers_c3_and_bounds_exact_analytic_profile(self):
  delta=Q(1,4)
  r=q.first_seed(3,delta,7,11,q.SEED_IDENTITY)
  self.assertEqual(r['upper'],18*(1+2/delta+4/delta**2+8/delta**3))
  import sympy as s
  z=s.symbols('z');f=(1+z/2)**4
  M=(1+(1+Q(3,2)*delta)/2)**4
  for m in range(8):
   norm=sum(s.diff(f,z,j).subs(z,1)/s.factorial(j) for j in range(m+1))
   cap=q.first_seed(m,delta,M,0,q.SEED_IDENTITY)['upper']
   self.assertLessEqual(norm,s.Rational(cap.numerator,cap.denominator))
  for m in range(12):
   bits=q.first_seed(m,delta,3,4,q.SEED_IDENTITY,binary=True)['upper']
   self.assertLessEqual(q.first_seed(m,delta,8,16,q.SEED_IDENTITY)['upper'],2**bits)
 def test_first_full_discrepancy_join_and_high_order_counts(self):
  r=q.first_discrepancies(Q(1,200),5,1,1,2,3,7,11,Q(1,4),13,17,q.SEED_IDENTITY)
  J=q.first_seed(5,Q(1,4),13,17,q.SEED_IDENTITY)['upper']
  manual=q.bounds(Q(1,200),1,5,1,1,2,3,[7],[11],None,J,J)
  self.assertEqual(r['moments'],manual);self.assertEqual(r['leading_U_input_order'],7)
  high=q.first_discrepancies(Q(1,200),100005,1,1,2,3,7,11,Q(1,4),13,17,q.SEED_IDENTITY,binary=True)
  self.assertEqual(high['leading_U_input_order'],100007)
  self.assertFalse(high['actual_leading_norms_and_scale_certified'])
  self.assertFalse(high['seed']['actual_original_strip_or_solution_bounds_certified'])
 def test_scope_guards(self):
  with self.assertRaises(ValueError):q.first_seed(2,Q(1,4),1,1,'new unrelated collar')
  with self.assertRaises(ValueError):q.first_seed(100005,Q(1,4),1,1,q.SEED_IDENTITY)
  with self.assertRaises(ValueError):q.first_seed(2,0,1,1,q.SEED_IDENTITY)
  with self.assertRaises(ValueError):q.first_seed(True,Q(1,4),1,1,q.SEED_IDENTITY)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
