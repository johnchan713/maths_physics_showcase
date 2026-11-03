import unittest
from fractions import Fraction as Q
from pathlib import Path
import q126 as q
import audit
class Tests(unittest.TestCase):
 def test_all_uniform_value_inequalities(self):
  checks=q.value_checks();self.assertTrue(all(checks.values()));self.assertGreaterEqual(len(checks),14)
 def test_exact_strength_jets_retain_odd_affine_end_energy(self):
  eps=Q(1,2**256);k=Q(1,5);p=[Q(0),Q(1),Q(0),Q(1)];a=[Q(1)];n=8
  def sq(v,j):return sum(v[i]*v[j-i] for i in range(j+1) if i<len(v) and j-i<len(v))
  for j in range(1,n+1):
   z=[(a[i] if i<j else 0)+(p[i] if i<len(p) else 0) for i in range(j+1)]
   cost=k*sq(a,j)+eps*sq(z,j)
   a.append(-cost/(2*(k+eps)))
  major=q.pulse.energy_root(n,Q(1),1/(2*(k+eps)),k,[Q(0)]*(n+1),[(eps,1,p)])
  self.assertEqual(a[1],-eps/(k+eps));self.assertNotEqual(a[1],0)
  for j in range(1,n+1):
   z=[a[i]+(p[i] if i<len(p) else 0) for i in range(n+1)]
   self.assertEqual(k*sq(a,j)+eps*sq(z,j),0)
   self.assertLessEqual(abs(a[j]),major[j]);self.assertLessEqual(major[j],Q(3,2)*2**j)
 def test_positive_future_mass_control_and_angular_cap(self):
  import sympy as s
  z=s.symbols('z');tiny=s.Rational(1,2**256)**30
  # A finite positive measure with mass 512T, real radial weights independent of eta.
  # Polynomial powers are a subset of the exact pre-heat interpolation powers.
  T=64;future=512*T*((1+z*z)**2+(1+z*z)+1)/3*(1+64*tiny*z)**2
  for m in range(7):
   for eta in (-1,0,1):
    norm=sum(abs(s.diff(future,z,j).subs(z,eta))/s.factorial(j) for j in range(m+1))
    self.assertLessEqual(norm,2**14*T*(2**(m+1)-1))
 def test_closed_axial_and_global_binary_bounds(self):
  for m in range(16):
   a=q.axial_angular(m,q.IDENTITY);b=q.axial_angular(m,q.IDENTITY,binary=True)
   self.assertLessEqual(a['upper'],2**b['upper'])
   g=q.global_U(m,7,11,q.IDENTITY)
   self.assertEqual(g['upper'],max(7,11,8,a['upper']))
  high=q.global_U(100007,13,17,q.IDENTITY,binary=True)
  self.assertEqual(high['upper'],100014)
  self.assertFalse(high['exterior']['heat_tail_complex_extension_claimed'])
  self.assertFalse(high['actual_axis_compact_matching_or_support_inputs_certified'])
 def test_current_first_discrepancy_join(self):
  h=Q(1,200);m=3
  r=q.current_first(h,m,1,1,2,3,7,11,13,Q(1,4),17,19,q.IDENTITY,q.current.SEED_IDENTITY)
  globalU=q.global_U(m+2,7,11,q.IDENTITY)['upper']
  expected=q.current.first_discrepancies(h,m,1,1,2,3,globalU,13,Q(1,4),17,19,q.current.SEED_IDENTITY)
  self.assertEqual(r['first'],expected)
  high=q.current_first(h,100005,1,1,2,3,7,11,13,Q(1,4),17,19,q.IDENTITY,q.current.SEED_IDENTITY,binary=True)
  self.assertEqual(high['global_U']['order'],100007)
  self.assertFalse(high['actual_leading_profile_or_scale_certified'])
 def test_scope_and_identity_guards(self):
  with self.assertRaises(ValueError):q.axial_angular(4,'complex Gamma heat continuation')
  with self.assertRaises(ValueError):q.axial_angular(True,q.IDENTITY)
  with self.assertRaises(ValueError):q.axial_angular(100005,q.IDENTITY)
  with self.assertRaises(ValueError):q.global_U(3,-1,1,q.IDENTITY)
  self.assertFalse(q.axial_angular(3,q.IDENTITY)['physical_wave_bound_certified'])
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
