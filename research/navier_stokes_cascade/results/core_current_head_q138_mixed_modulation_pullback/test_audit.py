import unittest
from fractions import Fraction as Q
from math import factorial,comb
from pathlib import Path
import q138 as q
import audit
class Tests(unittest.TestCase):
 def fixture(self,N=8):
  p=q.pullback(2,2,2,2,1,N,q.IDENTITY)
  return p,q.physical(q.modulation(p,2,2,2,2,2,q.IDENTITY),Q(1,4),q.IDENTITY)
 def test_exact_carrier_chain_and_radial_countercontrol(self):
  import sympy as s
  y,z,t=s.symbols('y z t');P=(1+y+y*y)*(1+z*z)*s.sin(t);N=8
  for r in range(5):
   for k in range(4):
    direct=s.diff(P.subs(t,N*y)/N,y,r,z,k)
    chain=sum(comb(r,j)*N**j*s.diff(P,y,r-j,z,k,t,j).subs(t,N*y)/N for j in range(r+1))
    self.assertEqual(s.simplify(direct-chain),0)
  for N in (1,8,64):
   self.assertEqual(s.diff(s.sin(N*y)/N,y).subs(y,0),1)
   recipe=q.pullback(2,2,2,2,1,N,q.IDENTITY)
   for r in range(6):
    raw=abs(s.diff(s.sin(N*y)/N,y,r).subs(y,0))
    self.assertLessEqual(raw,recipe['primitive_prefactor']*recipe['radial_base']**r*factorial(r)**3)
  self.assertEqual(s.diff(s.sin(8*y)/8,y,3).subs(y,0),-64)
 def test_exponential_positive_majorant_and_difference(self):
  import sympy as s
  x,z=s.symbols('x z');p,_=self.fixture();T=p['primitive_prefactor'];u=1/(8*(1+T));S=T*((1-u)**-2-1)
  self.assertLessEqual(S,Q(16,49));self.assertLessEqual(S,T)
  # Independent finite bivariate exponential of a slow jet with zero real base.
  D=(x+z+x*z)/8
  for r in range(5):
   for k in range(5):
    raw=abs(s.diff(s.exp(D)-1,x,r,z,k).subs({x:0,z:0}))
    cap=p['difference_prefactor']*p['exp_radial_base']**r*p['exp_angular_base']**k*(factorial(r)*factorial(k))**3
    self.assertLessEqual(raw,s.Rational(cap.numerator,cap.denominator))
 def test_stirling_physical_conversion(self):
  from sympy.functions.combinatorial.numbers import stirling
  for R in (Q(1),Q(3,2),Q(100)):
   for r in range(12):
    total=sum(int(stirling(r,j,kind=1))*R**j*factorial(j)**3 for j in range(r+1))
    self.assertLessEqual(total,(1+R)**r*factorial(r)**3)
  p,m=self.fixture();self.assertEqual(m['radial_bases']['F']-m['radial_bases']['E'],4)
 def test_final_repair_join_and_binary(self):
  p,m=self.fixture()
  root=q.previous.previous.axial.repair.root(dict(stage='same-leading-full-target-Gevrey3',Z=Q(1,1000),L=Q(2)),Q(1,100),beta=1,qnorm=1)
  rep=q.previous.mixed(q.previous.previous.axial.repair.field(root,2),Q(1,4),q.previous.IDENTITY)
  joined=q.join_repair(m,rep,q.IDENTITY)
  for kind in ('F','U'):
   self.assertEqual(joined['prefactors'][kind],m['prefactors'][kind]+rep['prefactors'][kind])
   for r in range(5):
    for k in range(5):
     self.assertGreaterEqual(q.bound(joined,kind,r,k)['upper'],q.bound(m,kind,r,k)['upper'])
     self.assertLessEqual(q.bound(joined,kind,r,k)['upper'],2**q.bound(joined,kind,r,k,True)['upper'])
   self.assertIsInstance(q.bound(joined,kind,100005,100005,True)['upper'],int)
  self.assertFalse(joined['actual_same_profile_cover_complete_pressure_target_or_cone_certified'])
 def test_identity_value_and_stage_guards(self):
  for vals in ((2,2,2,2,9,8,q.IDENTITY),(2,2,2,2,1,Q(3,2),q.IDENTITY),(2,2,2,2,1,8,'angular-only')):
   with self.assertRaises(ValueError):q.pullback(*vals)
  p,m=self.fixture()
  with self.assertRaises(ValueError):q.physical(dict(stage='PREMOD'),1,q.IDENTITY)
  with self.assertRaises(ValueError):q.bound(m,'F',100005,0)
  with self.assertRaises(ValueError):q.join_repair(m,dict(stage='four-row-repair'),q.IDENTITY)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   r=subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertIn(b'mixed-Ny-pullback',r.stdout);self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
