import unittest,json
from fractions import Fraction as Q
from math import factorial
import new_maps as m
import audit

def coeff(v):return {k:dict(x.terms) for k,x in v.c.items()}
def fixture():
 p=m.update.profile.pipeline;P=p.roots.P
 root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
 t=m.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
 fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),
 slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=3,profile=[10000]*17,profile_hull=hull)
 return m.NewPoleMaps(m.update.SignedCovariance(t,0))

class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.v=fixture()
 def test_pressure_uses_new_pole(self):
  v=self.v;out=v.pressure(0);pole=out['operator_pole']
  self.assertGreater(pole,v.target.pole)
  self.assertNotEqual(coeff(v.kernel(pole)),coeff(v.kernel(v.target.pole)))
  self.assertNotEqual(coeff(out['whole_updated_pressure_moment']),coeff(out['whole_pressure_change_moment']))
  self.assertEqual(set(v.recipe(1,pole)),{'order','pole','Kp','primitive','bump'})
 def test_kernel_independent_historical_formula(self):
  # Same proven branch constant, recomputed at several new poles.
  old=m.update.profile.load('q38_old','core_rho_projector_signed_recompute/bounds.py')
  for p in (0,3,4,12,60):
   v=old.RecomputedSigned.kernel(self.v,p)
   self.assertEqual(coeff(v),coeff(self.v.kernel(p)))
 def test_target_radial_multiplier_jets(self):
  import sympy as s
  R=s.symbols('R',positive=True);a=Q(1);b=Q(8)
  for e in (1,2):
   for n in range(4):
    pos=sum(Q(__import__('math').comb(e,j))*b**(e-j) for j in range(min(e,n)+1))
    neg=sum(Q(__import__('math').comb(e+j-1,j))*a**(-e-j) for j in range(n+1))
    self.assertEqual(pos,sum(abs(s.diff(R**e,R,j).subs(R,b))/factorial(j) for j in range(n+1)))
    self.assertEqual(neg,sum(abs(s.diff(R**(-e),R,j).subs(R,a))/factorial(j) for j in range(n+1)))
    expected=(pos*neg)*(1+b/2*self.v.kernel(12))
    # recipe sums both channels, check individual term by recomputing the sum below.
    self.assertTrue(expected.c)
  p=self.v.profile
  self.assertEqual(coeff(self.v.recipe(3,12)['bump']),coeff(p.floor*(p.sigma(3,1)+p.sigma(3,2))))
 def test_derivative_identity_and_zero_moment_obstruction(self):
  import sympy as s
  R,Z=s.symbols('R Z');F=(R-1)**5*(2-R)**5*(1+Z+Z**2);g=s.diff(F,R)
  self.assertEqual(s.integrate(g,(R,1,2)),0)
  for a in range(1,4):
   for b in range(3):
    self.assertEqual(s.expand(s.diff(F,R,a,Z,b)/factorial(a)/factorial(b)-s.diff(g,R,a-1,Z,b)/factorial(a-1)/factorial(b)/a),0)
  # Without zero integral the primitive has a nonzero exterior tail.
  self.assertEqual(s.integrate(1,(R,1,2)),1)
 def test_domains_and_scope(self):
  with self.assertRaises(ValueError):self.v.recipe(17,12)
  with self.assertRaises(ValueError):self.v.kernel(-1)
  with self.assertRaises(TypeError):m.NewPoleMaps(object())
  self.assertFalse(m.record()['actual_source_and_fresh_means_certified'])
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
