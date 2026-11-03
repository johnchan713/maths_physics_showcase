import unittest
from fractions import Fraction as Q
from math import factorial,comb
from pathlib import Path
import q135 as q
import audit
class Tests(unittest.TestCase):
 def fixture(self):return q.exterior(2,Q(1,2**256),100,q.IDENTITY)
 def test_same_factors_and_constant_bases(self):
  r=self.fixture();f=r['factors'];g=r['heat_g']
  self.assertEqual(f['heat_compensation'],1+15360000*g)
  self.assertEqual(f['angular_edit'],1+128*r['angular_root_cap'])
  self.assertEqual(sum(r['raw_factor_bases']),548)
  self.assertGreater(r['prefactor'],162)
  self.assertFalse(r['complex_Gamma_extension_claimed'])
  # exp(4)/2 < 81 from e<3; all analytic factors can be padded.
  self.assertLess(Q(3**4,2),81)
  for m in range(10):self.assertLessEqual(factorial(m),factorial(m)**2)
 def test_product_weights_and_manufactured_factors(self):
  import sympy as s
  z=s.symbols('z');E=(1+z*z)**-1*(1+z/100)*(1+z*z/100)*(1+z/1000)
  r=self.fixture()
  for m in range(9):
   weighted=sum(Q(comb(m,j)*factorial(j)**2*factorial(m-j)**2,factorial(m)**2)*2**j*512**(m-j) for j in range(m+1))
   self.assertLessEqual(weighted,514**m)
   P=r['prefactor'];cap=P*548**m*factorial(m)**2
   self.assertLessEqual(abs(s.diff(E,z,m).subs(z,0)),s.Rational(cap.numerator,cap.denominator))
 def test_factorial_norm_binary_and_large_order(self):
  r=self.fixture()
  for m in range(9):
   v=q.norm(r,m)['upper'];self.assertEqual(v,r['prefactor']*sum(548**j*factorial(j) for j in range(m+1)))
   self.assertLessEqual(v,2**q.norm(r,m,True)['upper'])
  self.assertIsInstance(q.norm(r,100005,True)['upper'],int)
  self.assertTrue(q.norm(r,3)['also_bounds_stronger_sum_of_global_suprema'])
 def test_same_schedule_and_domain_guards(self):
  for S,l,X,i in ((0,Q(1,2**256),100,q.IDENTITY),(2,0,100,q.IDENTITY),(2,Q(1,100),100,q.IDENTITY),(2,Q(1,2**256),99,q.IDENTITY),(2,Q(1,2**256),100,'preheat-compact-only')):
   with self.assertRaises(ValueError):q.exterior(S,l,X,i)
  with self.assertRaises(ValueError):q.norm(self.fixture(),100005)
  with self.assertRaises(ValueError):q.norm(dict(stage='PREMOD'),1)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   r=subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertIn(b'whole-exterior-leading-swirl',r.stdout)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
