import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q129 as q
import audit
class Tests(unittest.TestCase):
 def test_contraction_and_composition_disks(self):
  for W,B,C,d in [(1,1,1,Q(1,8)),(100,3,5,Q(1,1000)),(2**400,2**200,2**300,Q(1,2**500))]:
   r=q.inverse(W,B,C,d,q.IDENTITY);k,t,e=r['k'],r['t'],r['eps']
   self.assertLessEqual(k*e/(1-e),t/14)
   self.assertLessEqual(k*(1/((1-e)*(1-t)**2)-1),Q(1,4))
   self.assertLessEqual(8*k*(e+t),Q(17,128))
   p=q.compose(r,7,3,9,q.IDENTITY);v=2*t*p['radius']/(C*r['radius'])
   self.assertLessEqual(3*p['radius'],Q(1,4));self.assertLessEqual(9*v,Q(1,4))
 def test_independent_implicit_root(self):
  import sympy as s
  h=s.symbols('h');u=-2+2*s.sqrt(1-h) # u+h+u^2/4=0
  self.assertEqual(s.simplify(u+h+u*u/4),0)
  r=q.inverse(1,1,1,Q(1,2),q.IDENTITY)
  for m in range(1,9):
   b=q.bound(r,m)['upper'];self.assertLessEqual(abs(s.diff(u,h,m).subs(h,0)),s.Rational(b.numerator,b.denominator))
 def test_finite_gevrey_majorant(self):
  def mul(a,b,n):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(n+1)]
  # Actual smooth jets phi=u+sum_i (i!)^2 h^i+u^2; all finite identities.
  # The sign-free inverse recurrence is exact for phi=u-sum_i (i!)^2 h^i-u^2.
  n=8;raw=[Q(0)]*(n+1);major=[Q(0)]*(n+1)
  for m in range(1,n+1):
   raw[m]=factorial(m)**2+mul(raw,raw,n)[m]
   major[m]=1+mul(major,major,n)[m]
   self.assertLessEqual(raw[m]/factorial(m)**2,major[m])
   if m>=3:self.assertGreater(raw[m],major[m]) # missing factorial fails
  for i in range(7):
   for j in range(7):self.assertLessEqual(factorial(i)**2*factorial(j)**2,factorial(i+j)**2)
 def test_full_majorant_against_independent_algebraic_root(self):
  import sympy as s
  x=s.symbols('x');r=q.inverse(1,1,1,1,q.IDENTITY);k=r['k']
  kk=s.Rational(k.numerator,k.denominator)
  T=(1-s.sqrt(1-4*kk*(1+kk)*x/(1-x)))/(2*(1+kk))
  self.assertEqual(s.simplify(T-(1+kk)*T*T-kk*x/(1-x)),0)
  series=s.series(T,x,0,8).removeO().expand()
  for m in range(1,8):
   c=series.coeff(x,m);b=r['t']*r['eps']**(-m)
   self.assertGreater(c,0);self.assertLessEqual(c,s.Rational(b.numerator,b.denominator))
 def test_fixed_origin_and_weighted_mean(self):
  import sympy as s
  theta,eta=s.symbols('theta eta',real=True);e=s.Rational(1,4)
  phi=(theta+e*eta*(1-s.cos(theta)))/(2*s.pi);w=s.diff(phi,theta)
  theta_eta=-s.diff(phi,eta)/w
  self.assertEqual(theta_eta.subs({eta:0,theta:s.pi}),-2*e)
  moving=(s.diff(s.sin(theta),theta)*theta_eta).subs({eta:0,theta:s.pi})
  self.assertEqual(moving,2*e);self.assertNotEqual(moving,0)
  weighted=(1+e*s.cos(theta))/(2*s.pi)
  self.assertEqual(s.integrate(s.cos(theta),(theta,0,2*s.pi)),0)
  self.assertEqual(s.integrate(s.cos(theta)*weighted,(theta,0,2*s.pi)),e/2)
  ts=s.Rational(1,3);a=s.Rational(3,2);tloop=ts+e*s.sin(theta)
  vs=a*(1+ts**2);rho=a*e**2/2;v=vs+rho;b=-a*ts
  density=a*(1+tloop**2)/(2*s.pi*v)
  self.assertEqual(s.simplify(s.integrate(density,(theta,0,2*s.pi))),1)
  aL=v/(1+tloop**2);bL=-v*tloop/(1+tloop**2)
  for f in (-(aL-a)/2,(bL-b)/2):
   self.assertEqual(s.simplify(s.integrate(s.cancel(f*density),(theta,0,2*s.pi))),0)
  self.assertNotEqual(s.simplify(s.integrate(density*v/vs,(theta,0,2*s.pi))),1)
  z=s.symbols('z',real=True);g=s.sin(2*s.pi*z)
  v=s.symbols('v');P=s.integrate(g,(z,0,z))-s.integrate(s.integrate(g,(z,0,v)),(v,0,1))
  self.assertEqual(s.simplify(s.diff(P,z)-g),0);self.assertEqual(s.integrate(P,(z,0,1)),0)
  self.assertEqual(s.simplify(P.subs(z,1)-P.subs(z,0)),0)
 def test_budgets_join_and_guards(self):
  r=q.inverse(2,3,4,Q(1,16),q.IDENTITY);p=q.compose(r,5,2,3,q.IDENTITY)
  for recipe in (r,p):
   for m in range(1,7):self.assertLessEqual(q.bound(recipe,m)['upper'],2**q.bound(recipe,m,True)['upper'])
   self.assertEqual(q.bound(recipe,100005,True)['factorial_power'],3)
  self.assertEqual(q.axial_join(8,p,7,3)['raw_upper'],8+q.bound(p,3)['upper']/7)
  self.assertFalse(q.axial_join(8,p,7,3)['actual_final_repair_or_cone_certified'])
  with self.assertRaises(ValueError):q.axial_join(1,r,7,1)
  with self.assertRaises(ValueError):q.inverse(1,1,1,2,q.IDENTITY)
  with self.assertRaises(ValueError):q.bound(r,0)
  with self.assertRaises(ValueError):q.bound(p,100005)
  with self.assertRaises(ValueError):q.inverse(1,1,1,0,q.IDENTITY)
  with self.assertRaises(ValueError):q.compose(r,1,1,1,'unweighted theta mean')
  with self.assertRaises(ValueError):q.axial_join(1,p,Q(3,2),1)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
