import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q139 as q
import audit
class Tests(unittest.TestCase):
 def fixture(self):
  disk=q.phase_disk(q.scalar.parameter_disk(2,8,q.scalar.IDENTITY),q.IDENTITY)
  return disk,q.bind(3,2,2,disk,q.IDENTITY)
 def test_joint_inverse_gates(self):
  for r,W,F,G in ((Q(1,4),1,1,2),(Q(1,100),10**6,10**4,100)):
   d=q.phase_disk(dict(stage='analytic-scalar-parameter-disk',r=r,W=W,F=F,G=G),q.IDENTITY)
   self.assertLess(d['lipschitz'],Q(1,16));self.assertLessEqual(d['forcing'],d['theta_ball']/4)
   self.assertLess(d['forcing']+d['lipschitz']*d['theta_ball'],d['theta_ball']/2)
   self.assertLess(d['phase_radius'],1);self.assertFalse(d['off_constraint_complex_periodicity_assumed'])
  self.fixture()
 def test_fixed_origin_inverse_motion(self):
  import sympy as s
  a,t=s.symbols('a t',real=True)
  phi=(t+a*s.sin(t))/(2*s.pi);Theta_a=-s.diff(phi,a)/s.diff(phi,t)
  g_a=s.cos(t)*Theta_a
  self.assertEqual(s.simplify(g_a.subs({a:0,t:s.pi/4})),-s.Rational(1,2))
  self.assertNotEqual(g_a.subs({a:0,t:s.pi/4}),0)
  self.assertEqual(phi.subs(t,0),0)
 def test_real_mean_constraint_not_off_constraint_identity(self):
  import sympy as s
  a,ts,mu,z=s.symbols('a ts mu z',positive=True);t=s.symbols('t',real=True)
  v=a*(1+ts*ts+mu*mu/2);loop=ts+mu*s.sin(t)
  mass=s.integrate(a*(1+loop*loop)/(2*s.pi*(v+z)),(t,0,2*s.pi))
  self.assertEqual(s.simplify(mass-v/(v+z)),0)
  self.assertEqual(s.simplify(mass.subs(z,0)),1)
  weighted_source=-a*(1-mass)/2
  self.assertEqual(s.simplify(weighted_source.subs(z,0)),0)
  self.assertNotEqual(s.simplify(weighted_source),0)
 def test_independent_finite_mixed_composition(self):
  import sympy as s
  x,y=s.symbols('x y');disk=dict(stage='same-joint-scalar-phase-inverse-disk',scalar_radius=Q(1,8),phase_radius=Q(1,8),primitive_complex_cap=6)
  recipe=q.bind(3,2,2,disk,q.IDENTITY)
  # Analytic scalar/phase model cap <=4<6 on the disks; actual slow input is polynomial.
  poly=s.Poly(sum((4*(x+y+x*y))**n for n in range(7)),x,y)
  for i in range(4):
   for k in range(4):
    for j in range(4):
     raw=abs(poly.coeff_monomial(x**i*y**k))*factorial(i)*factorial(k)*4**j*factorial(j)
     cap=q.bound(recipe,i,k,j)['upper'];self.assertLessEqual(raw,s.Rational(cap.numerator,cap.denominator))
 def test_binding_binary_and_guards(self):
  disk,r=self.fixture();p=q.pullback(r,1,8,q.IDENTITY)
  self.assertEqual(p['radial_base'],r['A']+8*r['C'])
  self.assertFalse(r['actual_mixed_scalar_state_root_augmented_variance_and_pressure_binding_certified'])
  for i in range(5):
   for k in range(5):
    for j in range(4):self.assertLessEqual(q.bound(r,i,k,j)['upper'],2**q.bound(r,i,k,j,True)['upper'])
  self.assertIsInstance(q.bound(r,100005,100005,100005,True)['upper'],int)
  with self.assertRaises(ValueError):q.bind(3,2,2,dict(stage='angular-only'),q.IDENTITY)
  with self.assertRaises(ValueError):q.phase_disk(dict(stage='PREMOD'),q.IDENTITY)
  with self.assertRaises(ValueError):q.pullback(r,1,8,'unweighted-mean')
  with self.assertRaises(ValueError):q.bound(r,100005,0,0)
 def test_00_record(self):
  import json,tempfile
  with tempfile.TemporaryDirectory() as t:
   out=Path(t)/'replay.json';frozen=audit.HERE/'evidence.json'
   out.write_text(json.dumps(audit.generate(),indent=2,sort_keys=True)+'\n')
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
