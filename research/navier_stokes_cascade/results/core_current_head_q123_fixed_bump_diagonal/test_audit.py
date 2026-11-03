import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q123 as q
import audit
class Tests(unittest.TestCase):
 def test_shifted_polynomial_full_coefficients(self):
  for m in range(13):
   for alpha in (Q(0),Q(1,2)):
    coeff=q.old.log_coefficients(m,alpha)
    self.assertEqual(sum(abs(c) for c in coeff),__import__('functools').reduce(lambda a,i:a*(1+alpha+i),range(m),Q(1)))
    self.assertLessEqual(sum(abs(c) for c in coeff),(1+alpha)**m*factorial(m))
  self.assertEqual(q.old.log_coefficients(3,Q(1,2)),(Q(-15,8),Q(23,4),Q(-9,2),Q(1)))
 def test_independent_exact_logarithmic_pullback(self):
  import sympy as s
  x=s.symbols('x',positive=True);y=s.symbols('y');b=1+y+y**3+y**5
  for alpha in (s.Rational(0),s.Rational(1,2)):
   for m in range(7):
    c=q.old.log_coefficients(m,Q(alpha));rhs=sum(s.Rational(v.numerator,v.denominator)*s.diff(b,y,j).subs(y,s.log(x)) for j,v in enumerate(c))*x**(-m-alpha)
    self.assertEqual(s.simplify(s.diff(x**(-alpha)*b.subs(y,s.log(x)),x,m)-rhs),0)
 def test_low_actual_bump_certificates_preserved(self):
  self.assertEqual(tuple(2**q.bump_bits(j)>=q.LOW[j] for j in range(4)),(True,)*4)
  B=q.old.repair_bump_bounds(3)
  self.assertEqual(tuple(B),q.LOW)
  for p,rp in ((Q(1,4),Q(1,2)),(Q(2),Q(2))):
   for m in range(4):
    bound=q.moderate_radial(m,p,rp,q.IDENTITY)
    oldU=q.old.bump_radial(B[:m+1],p,rp,'U')
    oldF=q.old.bump_radial(B[:m+1],p,rp,'F')
    self.assertLessEqual(max(oldU),bound['U']);self.assertLessEqual(max(oldF),bound['E_over_R'])
 def test_binary_bounds_and_large_order(self):
  for m in range(12):
   r=q.radial_bits(m,Q(1,4),Q(1,2),q.IDENTITY);v=q.moderate_radial(m,Q(1,4),Q(1,2),q.IDENTITY)
   for k in ('U','E_over_R'):self.assertLessEqual(v[k],2**r[k])
  r=q.radial_bits(100006,Q(1,4),Q(1,2),q.IDENTITY)
  self.assertFalse(r['factorials_materialized'])
  with self.assertRaises(ValueError):q.moderate_radial(100006,1,1,q.IDENTITY)
 def test_current_full_join_dominates_moderate_exact_recipe(self):
  h=Q(1,200);D={k:Q(7) for k in q.repair.KEYS}
  for n,m in ((1,0),(3,2),(4,4)):
   r=q.current_row(h,n,m,8,1,4,2,Q(1,100),3,D,3,4,Q(1,4),Q(1,2),q.repair.IDENTITY,q.IDENTITY,0)
   BF=q.moderate_radial(m,Q(1,4),Q(1,2),q.IDENTITY)['E_over_R']
   BU=q.moderate_radial(m+1,Q(1,4),Q(1,2),q.IDENTITY)['U']
   exact=q.repair.finalized_tuple(h,n,m,8,1,4,2,Q(1,100),3,D,8,16,BF,BU,q.repair.IDENTITY)
   self.assertLessEqual(q.q121.materialize(exact['tuple_recipe']),2**r['tuple_binary_upper'])
   self.assertFalse(r['row']['sufficient_new_row_gate'])
   supplied=Q(n+r['row']['normalized_constant_binary_upper'],n*h)
   passed=q.current_row(h,n,m,8,1,4,2,Q(1,100),3,D,3,4,Q(1,4),Q(1,2),q.repair.IDENTITY,q.IDENTITY,supplied)
   self.assertTrue(passed['row']['sufficient_new_row_gate']);self.assertFalse(passed['row']['scales_selected_or_changed'])
  huge=q.current_row(h,100005,100005,8,1,4,2,Q(1,100),3,D,3,4,Q(1,4),Q(1,2),q.repair.IDENTITY,q.IDENTITY,0)
  self.assertEqual(huge['discrepancy_order'],100006);self.assertFalse(huge['factorials_materialized'])
  self.assertFalse(huge['actual_seed_discrepancy_identity_support_or_scale_inputs_certified'])
 def test_input_guards(self):
  for args in ((True,1,1,q.IDENTITY),(1,0,1,q.IDENTITY),(1,1,2,q.IDENTITY),(1,1,1,'Borel cutoff')):
   with self.assertRaises(ValueError):q.radial_bits(*args)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
