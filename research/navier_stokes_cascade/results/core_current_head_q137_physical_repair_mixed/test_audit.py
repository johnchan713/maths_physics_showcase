import unittest
from fractions import Fraction as Q
from math import factorial,comb
from pathlib import Path
import q137 as q
import audit
class Tests(unittest.TestCase):
 def fixture(self):
  r=q.previous.axial.repair.root(dict(stage='same-leading-full-target-Gevrey3',Z=Q(1,1000),L=Q(2)),Q(1,100),beta=1,qnorm=1)
  return q.mixed(q.previous.axial.repair.field(r,2),Q(1,4),q.IDENTITY)
 def test_stirling_cycle_formula_and_all_order_cap(self):
  from sympy.functions.combinatorial.numbers import stirling
  for D in (1,3,32600):
   for r in range(13):
    total=sum(int(stirling(r,k,kind=1))*D**k*factorial(k+1)**3 for k in range(r+1))
    self.assertLessEqual(total,(8*(1+D))**r*factorial(r)**3)
    if r:
     for k in range(1,r+1):self.assertLessEqual(int(stirling(r,k,kind=1))*factorial(k),factorial(r)*comb(r-1,k-1))
 def test_signed_physical_chain_and_multiplier(self):
  import sympy as s
  from sympy.functions.combinatorial.numbers import stirling
  X,y,z=s.symbols('X y z',positive=True);f=y**4+2*y*y
  for r in range(7):
   converted=X**(-r)*sum(stirling(r,k,kind=1,signed=True)*s.diff(f,y,k).subs(y,s.log(X)) for k in range(r+1))
   self.assertEqual(s.simplify(s.diff(f.subs(y,s.log(X)),X,r)-converted),0)
  recipe=self.fixture();F=s.log(X)**2*(z+z*z)/(1000*s.sqrt(2*X))
  for r in range(4):
   for k in range(4):
    value=abs(s.diff(F,X,r,z,k).subs({X:1,z:0}));cap=q.bound(recipe,'F',r,k)['upper']
    self.assertLessEqual(value,s.Rational(cap.numerator,cap.denominator))
 def test_total_raw_binary_and_tuple_contribution(self):
  recipe=self.fixture()
  for kind in ('E','F','U'):
   for m in range(7):
    cap=q.raw_total(recipe,kind,m)['upper']
    for r in range(m+1):
     for k in range(m-r+1):self.assertLessEqual(q.bound(recipe,kind,r,k)['upper'],cap)
    self.assertLessEqual(cap,2**q.raw_total(recipe,kind,m,True)['upper'])
   self.assertIsInstance(q.raw_total(recipe,kind,100005,True)['upper'],int)
   self.assertIsInstance(q.bound(recipe,kind,100005,100005,True)['upper'],int)
  t=q.repair_tuple(recipe,Q(1,1000),2,2,1,2)
  self.assertEqual(t['tuple_recipe']['U_order'],3)
  self.assertFalse(t['actual_complete_leading_tuple_cone_or_scale_certified'])
 def test_angular_radial_obstruction_and_memory(self):
  import sympy as s
  X,M=s.symbols('X M',positive=True);f=s.sin(M*X)
  self.assertEqual(s.diff(f,X).subs(X,0),M)
  # Value bound |Acal|<=N alone does not control the radial modulation derivative.
  self.assertEqual(s.diff(s.exp(s.sin(M*X)),X).subs(X,0),M)
  t=s.symbols('t');b=t*t*(1-t)**2
  self.assertEqual(b.subs(t,1),0);self.assertEqual(s.integrate(b,(t,0,1)),s.Rational(1,30))
 def test_scope_and_order_guards(self):
  r=self.fixture()
  self.assertFalse(r['complete_modulated_profile_or_pressure_support_certified'])
  with self.assertRaises(ValueError):q.mixed(dict(stage='PREMOD'),1,q.IDENTITY)
  with self.assertRaises(ValueError):q.mixed(dict(stage='same-fixed-bump-repair-field-Gevrey3'),0,q.IDENTITY)
  with self.assertRaises(ValueError):q.bound(r,'F',100005,0)
  with self.assertRaises(ValueError):q.bound(r,'Pi',1,1)
  with self.assertRaises(ValueError):q.repair_tuple(r,Q(1,1000),12,2,1,2)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   p=subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertIn(b'physical-leading-repair-mixed',p.stdout);self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
