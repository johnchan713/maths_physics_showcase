import unittest
from fractions import Fraction as Q
from pathlib import Path
import q122 as q
import audit
class Tests(unittest.TestCase):
 def test_exact_normalizer_jets(self):
  import sympy as s
  z=s.symbols('z');f=sum((j+1)*z**j for j in range(8));e=s.Rational(3)
  for k in range(10):
   rhs=2*s.diff(f,z,k)+2*k*s.diff(f,z,k-1) if k else 2*f
   if k>=2:rhs+=k*(k-1)*s.diff(f,z,k-2)
   # Pointwise identity includes eta-dependent factors, unlike its supremum bound.
   exact=(1+z*z)*s.diff(f,z,k)+(2*k*z*s.diff(f,z,k-1) if k else 0)+(k*(k-1)*s.diff(f,z,k-2) if k>=2 else 0)
   self.assertEqual(s.expand(s.diff((1+z*z)*f/e,z,k)-exact/e),0)
   jets=[Q(sum(abs(c) for c in s.Poly(s.diff(f,z,j),z).all_coeffs())) for j in range(k+1)]
   bound=q.normalizer_raw(k,jets,3)
   for eta in (-1,Q(1,3),1):self.assertLessEqual(abs(s.diff((1+z*z)*f/e,z,k).subs(z,s.Rational(eta.numerator,eta.denominator))),s.Rational(bound.numerator,bound.denominator))
 def test_five_row_linear_solve_mixed_targets(self):
  import sympy as s
  z=s.symbols('z');d=[(j+1)*(1+z+z**(j+2)) for j in range(5)]
  P=(1+z*z)/3;rho=2;lam=s.Rational(1,100);Kg=s.Rational(1,4)
  target=s.Matrix([d[0]/rho,(d[0]/rho-Kg*P*d[3])/lam,Kg*d[1],-P*d[4]/rho,P*d[2]])
  B=s.diag(1,2,3,4,5);c=-B.inv()*target
  self.assertEqual(B*c+target,s.zeros(5,1))
  for m in range(9):
   D={f'm{i+1}':Q(max(sum(abs(a) for a in s.Poly(s.diff(f,z,k),z).all_coeffs()) for k in range(m+1))) for i,f in enumerate(d)}
   r=q.coefficient_bound(2,m,2,Q(1,100),3,D,q.IDENTITY);C=r['coefficient_raw_upper']
   for f in c:
    for k in range(m+1):
     for eta in (-1,0,1):self.assertLessEqual(abs(s.diff(f,z,k).subs(z,eta)),s.Rational(C.numerator,C.denominator))
 def test_raw_through_order_cost_and_factorial_convention(self):
  for m in range(12):
   self.assertEqual(q.normalizer_raw(m,[Q(1)]*(m+1),2),Q(m*m+m+2,2))
  # The normalized polynomial coefficients are (2,2,1)/e at worst;
  # the inherited complete factorial C_eta^3 normalizer norm is 5/e.
  self.assertEqual(sum([Q(2,3),Q(2,3),Q(1,3)]),Q(5,3))
 def test_current_tuple_join_and_large_order_recipe(self):
  D={k:Q(7) for k in q.KEYS}
  r=q.finalized_tuple(Q(1,200),1,4,8,1,4,2,Q(1,100),3,D,5,6,7,8,q.IDENTITY)
  C=r['repair']['coefficient_raw_upper']
  expected=q.current.recipe(Q(1,200),1,4,8,1,4,5+7*C,6+8*C)
  self.assertEqual(r['tuple_recipe'],expected)
  large=q.finalized_tuple(Q(1,200),1,100005,8,1,4,2,Q(1,100),3,D,5,6,7,8,q.IDENTITY)
  self.assertEqual(large['discrepancy_order'],100006)
  with self.assertRaises(ValueError):q.current.materialize(large['tuple_recipe'])
  self.assertFalse(large['actual_seed_bump_moment_or_support_inputs_certified'])
 def test_scope_and_validation(self):
  D={k:Q(1) for k in q.KEYS}
  for args in ((0,2,2,Q(1,100),3,D,q.IDENTITY),(1,2,2,Q(1,100),3,D,'generic smooth e'),(1,2,2,Q(1,100),3,{'m1':1},q.IDENTITY),(1,2,2,0,3,D,q.IDENTITY)):
   with self.assertRaises(ValueError):q.coefficient_bound(*args)
  r=q.coefficient_bound(1,2,2,Q(1,100),3,D,q.IDENTITY)
  self.assertFalse(r['actual_full_pressure_discrepancy_bounds_certified']);self.assertFalse(r['leading_root_or_support_certified'])
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
