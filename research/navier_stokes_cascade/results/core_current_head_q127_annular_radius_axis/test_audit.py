import unittest
from fractions import Fraction as Q
from pathlib import Path
import q127 as q
import audit
class Tests(unittest.TestCase):
 def test_uniform_complex_contraction_gates(self):
  for H,Z in ((Q(1,8),Q(2)),(Q(1,1024),Q(10**20)),(Q(1,4),Q(0))):
   r=q.radius(H,Z,q.PREMOD_ID)
   self.assertTrue(all(r['checks'].values()));self.assertLessEqual(r['sigma'],H/2)
   self.assertEqual(r['real_value_cap'],q.D)
   self.assertFalse(r['actual_target_and_same_root_contracts_certified'])
 def test_actual_varying_quadratic_jets_and_frozen_countercontrol(self):
  d=q.D;v=d/100;qcoef=[Q(1,100),0,Q(1,50),0,Q(1,100)]
  # Scalar positive-coefficient control with the same eta-dependent quadratic shape.
  a=[Q(0)];frozen=[Q(0)]
  for k in range(1,9):
   def lower(coeff,root):
    return sum(coeff[i]*root[j]*root[k-i-j] for i in range(min(k,4)+1) for j in range(k-i+1) if j<len(root) and k-i-j<len(root))
   a.append((v if k==1 else 0)-lower(qcoef,a))
   frozen.append((v if k==1 else 0)-lower([Q(1,100),0,0,0,0],frozen))
  self.assertNotEqual(a[4],frozen[4])
  R=q.radius(Q(1,8),2*v,q.PREMOD_ID)
  for k in range(1,9):self.assertLessEqual(abs(a[k]),q.A*R['sigma']**(-k))
  # Nonzero base: Q_eta derivatives acting on c0^2 cannot be discarded.
  c=[d/10,v]+[Q(0)]*7
  def coeff(poly,j):return poly[j] if j<len(poly) else 0
  target=[coeff(c,k)+sum(coeff(qcoef,i)*coeff(c,j)*coeff(c,k-i-j) for i in range(min(k,4)+1) for j in range(k-i+1)) for k in range(9)]
  restored=[c[0]]
  for k in range(1,9):
   cost=sum(coeff(qcoef,i)*coeff(restored,j)*coeff(restored,k-i-j) for i in range(min(k,4)+1) for j in range(k-i+1))
   restored.append((target[k]-cost)/(1+2*qcoef[0]*c[0]))
  self.assertEqual(restored,c)
  self.assertNotEqual(target[2]-(qcoef[0]*v*v),0)
 def test_inverse_radius_and_premod_binary_budget(self):
  for hb,zb in ((3,1),(10,20),(1,0)):
   H=Q(1,2**hb);Z=Q(2**zb);r=q.radius(H,Z,q.PREMOD_ID)
   self.assertLessEqual(1/r['sigma'],2**q.inverse_radius_bits(hb,zb))
   for m in range(8):
    exact=q.premod_U(m,8,H,Z,q.PREMOD_ID)
    bits=q.premod_U(m,3,hb,zb,q.PREMOD_ID,binary=True)
    self.assertLessEqual(exact['upper'],2**bits['upper'])
  high=q.premod_U(100005,3,10,20,q.PREMOD_ID,binary=True)
  self.assertEqual(high['stage'],'PREMOD');self.assertFalse(high['actual_prefix_target_or_final_modulated_norm_certified'])
 def test_same_axis_cauchy_cover_and_bits(self):
  delta=Q(1,8);K=7
  self.assertEqual(q.axis_U(3,K,delta,q.AXIS_ID)['upper'],K*(1+1/delta+1/delta**2+1/delta**3))
  for m in range(10):
   self.assertLessEqual(q.axis_U(m,8,delta,q.AXIS_ID)['upper'],2**q.axis_U(m,3,delta,q.AXIS_ID,binary=True)['upper'])
 def test_final_compact_input_cannot_be_replaced_by_premod(self):
  m=3;premod=q.premod_U(m,7,Q(1,8),2,q.PREMOD_ID)
  with self.assertRaises(ValueError):q.leading_U(m,7,Q(1,8),premod,q.AXIS_ID,q.current.IDENTITY)
  final=dict(stage=q.FINAL_STAGE,order=m,upper=11,binary=False)
  joined=q.leading_U(m,7,Q(1,8),final,q.AXIS_ID,q.current.IDENTITY)
  self.assertEqual(joined['global_U']['upper'],max(joined['axis']['upper'],11,8,40*(2**(m+1)-1)))
  self.assertFalse(joined['actual_final_compact_norm_and_scale_certified'])
 def test_identity_and_scope_guards(self):
  for H,Z,identity in ((0,1,q.PREMOD_ID),(Q(1,8),-1,q.PREMOD_ID),(Q(1,8),1,'later constant compact map')):
   with self.assertRaises(ValueError):q.radius(H,Z,identity)
  with self.assertRaises(ValueError):q.premod_U(100005,1,Q(1,8),1,q.PREMOD_ID)
  with self.assertRaises(ValueError):q.axis_U(True,1,Q(1,8),q.AXIS_ID)
  with self.assertRaises(ValueError):q.axis_U(3,1,Q(1,8),'first-coefficient seed')
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
