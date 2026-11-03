import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q133 as q
import audit
def fixture():
 p=dict(stage='real-zero-mean-primitive',prefactor=Q(4),base=Q(2))
 mod=q.repair.previous.modulation(4,3,4,Q(1,4),p,p,1,8,q.repair.previous.IDENTITY)
 root=q.repair.root(dict(stage='same-leading-full-target-Gevrey3',Z=Q(1,1000),L=Q(2)),Q(1,100),beta=1,qnorm=1)
 return mod,q.repair.field(root,2)
class Tests(unittest.TestCase):
 def test_exact_smooth_norm_order_countercontrol(self):
  import sympy as s
  X,z=s.symbols('X z',real=True);U=100*(1-X+X*s.sin(4*z)/4);d=s.diff(U,z)
  self.assertEqual(U.subs({X:0,z:0}),100)
  self.assertEqual(U.subs({X:s.Rational(1,2),z:s.pi/8}),s.Rational(125,2))
  self.assertEqual(d.subs({X:s.Rational(1,2),z:0}),50)
  self.assertEqual(d.subs({X:1,z:0}),100)
  value=100-75*X;deriv=100*X
  self.assertEqual(s.expand(value+deriv),100+25*X)
  c=q.counterexample();self.assertEqual(c['weak_global'],125)
  self.assertEqual(max(c['weak_axis'],c['weak_compact']),c['weak_global'])
  self.assertGreater(c['strong_global'],max(c['strong_axis'],c['strong_compact']))
  old=q.regional.current.global_U(1,c['strong_axis'],c['strong_compact'],q.regional.current.IDENTITY)
  self.assertEqual(old['upper'],c['invalid_strong_max']);self.assertLess(old['upper'],c['strong_global'])
  # The SAME historical maximum is correct on its intended pointwise-X input bounds.
  valid=q.regional.current.global_U(1,c['weak_axis'],c['weak_compact'],q.regional.current.IDENTITY)
  self.assertEqual(valid['upper'],c['weak_global'])
 def test_final_field_and_normalized_jets(self):
  import sympy as s
  mod,field=fixture();r=q.compact(mod,field,q.IDENTITY)
  self.assertEqual(r['prefactor'],mod['prefactors']['U']+q.BUMP0*field['prefactor'])
  z=s.symbols('z');c=(s.sqrt(1+s.Rational(1,250)*s.sin(z))-1)/2
  U=1+z+s.sin(z)/8+2*c/(1+z*z) # fixed bump value 1 <= actual B0
  for m in range(7):
   cap=r['prefactor']*r['base']**m*factorial(m)**3
   self.assertLessEqual(abs(s.diff(U,z,m).subs(z,0)),s.Rational(cap.numerator,cap.denominator))
   norm=sum(abs(s.diff(U,z,j).subs(z,0))/s.factorial(j) for j in range(m+1))
   b=q.compact_norm(r,m)['upper'];self.assertLessEqual(norm,s.Rational(b.numerator,b.denominator))
 def test_regional_weak_strong_and_orderwise_gluing(self):
  mod,field=fixture();r=q.compact(mod,field,q.IDENTITY)
  for m in range(7):
   a=q.regional.axis_U(m,7,Q(1,8),q.regional.AXIS_ID);c=q.compact_norm(r,m)
   weak=q.regional_max(m,a,c,q.regional.current.IDENTITY)
   strong=q.regional_sum(m,a,c,q.regional.current.IDENTITY)
   jet=q.jetwise(m,7,Q(1,8),r,q.regional.current.IDENTITY)
   ext=q.regional.current.axial_angular(m,q.regional.current.IDENTITY)['upper']
   self.assertEqual(weak['upper'],max(a['upper'],c['upper'],8,ext))
   old_join=q.regional.leading_U(m,7,Q(1,8),c,q.regional.AXIS_ID,q.regional.current.IDENTITY)
   self.assertEqual(old_join['global_U']['upper'],weak['upper'])
   self.assertEqual(strong['upper'],a['upper']+c['upper']+8+ext)
   self.assertLessEqual(weak['upper'],strong['upper']);self.assertLessEqual(jet['upper'],strong['upper'])
   self.assertFalse(weak['stronger_sum_of_global_suprema_certified'])
 def test_binary_budgets_and_high_order(self):
  r=q.compact(*fixture(),q.IDENTITY)
  for m in range(7):
   c=q.compact_norm(r,m);cb=q.compact_norm(r,m,True)
   self.assertLessEqual(c['upper'],2**cb['upper'])
   a=q.regional.axis_U(m,8,Q(1,8),q.regional.AXIS_ID);ab=q.regional.axis_U(m,3,Q(1,8),q.regional.AXIS_ID,True)
   for method in (q.regional_max,q.regional_sum):
    self.assertLessEqual(method(m,a,c,q.regional.current.IDENTITY)['upper'],2**method(m,ab,cb,q.regional.current.IDENTITY,True)['upper'])
  m=100005;c=q.compact_norm(r,m,True);a=q.regional.axis_U(m,3,Q(1,8),q.regional.AXIS_ID,True)
  self.assertEqual(q.regional_sum(m,a,c,q.regional.current.IDENTITY,True)['upper'],q.regional_max(m,a,c,q.regional.current.IDENTITY,True)['upper']+2)
 def test_final_stage_source_order_and_scope_guards(self):
  r=q.compact(*fixture(),q.IDENTITY);m=5;c=q.compact_norm(r,m);a=q.regional.axis_U(m,7,Q(1,8),q.regional.AXIS_ID)
  g=q.regional_max(m,a,c,q.regional.current.IDENTITY);self.assertEqual(q.first_input(g,3)['input_order'],5)
  with self.assertRaises(ValueError):q.first_input(g,5)
  with self.assertRaises(ValueError):q.first_input(q.regional_sum(m,a,c,q.regional.current.IDENTITY),3)
  with self.assertRaises(ValueError):q.compact_norm(dict(stage='PREMOD'),m)
  with self.assertRaises(ValueError):q.regional_max(m,a,dict(c,stage='PREMOD'),q.regional.current.IDENTITY)
  with self.assertRaises(ValueError):q.regional_max(m,a,dict(c,order=4),q.regional.current.IDENTITY)
  with self.assertRaises(ValueError):q.regional_max(m,a,dict(c,binary=True),q.regional.current.IDENTITY)
  with self.assertRaises(ValueError):q.compact_norm(r,100005)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
