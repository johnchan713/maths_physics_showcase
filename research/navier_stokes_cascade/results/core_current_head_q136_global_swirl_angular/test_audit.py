import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q136 as q
import audit
class Tests(unittest.TestCase):
 def fixture(self):
  p=dict(stage='real-zero-mean-primitive',prefactor=Q(4),base=Q(2))
  mod=q.axial.repair.previous.modulation(4,3,4,Q(1,4),p,p,1,8,q.axial.repair.previous.IDENTITY)
  root=q.axial.repair.root(dict(stage='same-leading-full-target-Gevrey3',Z=Q(1,1000),L=Q(2)),Q(1,100),beta=1,qnorm=1)
  return mod,q.axial.repair.field(root,2)
 def inputs(self,m,binary=False):
  c=q.compact(*self.fixture(),Q(1,4),100,q.IDENTITY)
  return q.axis(3 if binary else 8,Q(1,8),m,binary),q.compact_norm(c,m,binary),q.exterior_norm(2,Q(1,2**256),100,m,binary)
 def test_regular_axis_and_singular_quotient(self):
  import sympy as s
  X,z=s.symbols('X z',positive=True);F=1/(1+z*z)
  E=s.sqrt(2*X)*F
  self.assertEqual(s.limit(E,X,0,dir='+'),0);self.assertEqual(s.limit(E/s.sqrt(2*X),X,0,dir='+'),F)
  for m in range(8):
   cap=q.axis(2,Q(1,2),m)['upper']
   self.assertLessEqual(sum(abs(s.diff(F,z,j).subs(z,0))/factorial(j) for j in range(m+1)),cap)
 def test_final_compact_and_E_to_F_loss(self):
  mod,r=self.fixture();c=q.compact(mod,r,Q(1,4),100,q.IDENTITY)
  self.assertEqual(c['prefactor'],4*(mod['prefactors']['E']+q.axial.BUMP0*r['prefactor']))
  for x in (Q(1,100),Q(1,4),Q(1),Q(100)):
   D=max(Q(1),1/x);self.assertGreaterEqual(D*D,1/(2*x))
  for m in range(8):self.assertEqual(q.compact_norm(c,m)['upper'],c['prefactor']*sum(c['base']**j*factorial(j)**2 for j in range(m+1)))
 def test_cover_and_intended_maximum(self):
  for m in range(7):
   items=self.inputs(m);joined=q.join(*items,Q(1,2),100,q.IDENTITY)
   self.assertEqual(joined['upper'],max(v['upper'] for v in items))
   self.assertFalse(joined['stronger_sum_of_global_suprema_certified'])
  items=self.inputs(2)
  with self.assertRaises(ValueError):q.join(items[0],dict(items[1],xmax=200),q.exterior_norm(2,Q(1,2**256),200,2),1,100,q.IDENTITY)
  with self.assertRaises(ValueError):q.join(*items,Q(1,8),100,q.IDENTITY)
  with self.assertRaises(ValueError):q.join(*items,Q(1,2),101,q.IDENTITY)
 def test_binary_and_high_orders(self):
  for m in range(7):
   a=q.join(*self.inputs(m),Q(1,2),100,q.IDENTITY);b=q.join(*self.inputs(m,True),Q(1,2),100,q.IDENTITY)
   self.assertLessEqual(a['upper'],2**b['upper'])
  self.assertIsInstance(q.join(*self.inputs(100005,True),Q(1,2),100,q.IDENTITY)['upper'],int)
 def test_stages_orders_and_first_pair(self):
  m=3;F=q.join(*self.inputs(m),Q(1,2),100,q.IDENTITY)
  U=dict(stage='whole-leading-pointwise-X-angular-algebra-norm',order=m+2,upper=123,binary=False)
  self.assertEqual(q.first_inputs(F,U,m)['U_order'],5)
  for bad in (dict(U,order=m),dict(U,stage='PREMOD'),dict(U,binary=True)):
   with self.assertRaises(ValueError):q.first_inputs(F,bad,m)
  items=self.inputs(m)
  with self.assertRaises(ValueError):q.join(items[0],dict(items[1],order=2),items[2],1,100,q.IDENTITY)
  mod,r=self.fixture()
  with self.assertRaises(ValueError):q.compact(dict(mod,stage='PREMOD'),r,1,100,q.IDENTITY)
  with self.assertRaises(ValueError):q.compact(mod,r,0,100,q.IDENTITY)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   r=subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertIn(b'regular-axis-final-compact',r.stdout);self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
