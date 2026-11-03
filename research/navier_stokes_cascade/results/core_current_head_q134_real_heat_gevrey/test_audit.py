import unittest
from fractions import Fraction as Q
from math import factorial, comb
from pathlib import Path
import q134 as q
import audit

def poch(h,m):
 v=Q(1)
 for j in range(m):v*=h+j
 return v
class Tests(unittest.TestCase):
 def test_real_pole_chain_and_zero_Gamma_radius(self):
  # Exact pole coefficients: (-1)^a from H and DeltaZ cancel.
  for h in (Q(1,100),Q(1,10000)):
   for X in (100,1000):
    recipe=q.heat(h,X)
    for m in range(1,13):
     coeff=Q(0)
     for a in range((m+1)//2,m+1):
      k=m-a
      coeff+=poch(h,a)*poch(1+h,a)/factorial(a)*comb(a,k)*4**(a-k)*2**k/X**a
     raw=coeff*factorial(m)
     self.assertLessEqual(raw,q.bound(recipe,m)['upper'])
    ratios=[(m+h)*(m+1+h)/(m+1) for m in (10,100,1000)]
    self.assertGreater(ratios[-1],900);self.assertTrue(ratios[0]<ratios[1]<ratios[2])
  self.assertEqual(4*Q(1,16)+2*Q(1,16)**2,Q(33,128))
 def test_square_multiplier_and_polynomial_normalization(self):
  import sympy as s
  v=s.symbols('v');eta=s.symbols('eta',real=True)
  for degree,cost in ((1,5),(2,25)):
   poly=s.Poly(s.expand((1+(eta+v)**2)**degree),v)
   coeffs=[s.Abs(c.subs(eta,1)) for c in reversed(poly.all_coeffs())]
   self.assertEqual(sum(coeffs),cost)
   for m in range(13):
    normalized=sum(coeffs[j]*64**(m-j)*factorial(m-j) for j in range(min(m,2*degree)+1))
    self.assertLessEqual(normalized,cost*64**m*factorial(m))
  h=Q(1,100);X=100;A=4*h/X
  self.assertLessEqual(2*A+A*A,16*h/X)
  for m in range(13):
   square=sum(Q(comb(m,j)*factorial(j)**2*factorial(m-j)**2,factorial(m)**2)*32**m for j in range(m+1))
   self.assertLessEqual(square,64**m)
 def test_infinite_tail_integrals_and_nonuniform_truncation(self):
  import sympy as s
  x=s.symbols('x',positive=True);h=s.symbols('h',positive=True)
  for exponent,expect in ((-3-2*h,1/(2+2*h)),(-2-2*h,1/(1+2*h)),(-1-h,1/h)):
   self.assertEqual(s.simplify(s.integrate(x**exponent,(x,1,s.oo))-expect),0)
  self.assertEqual(s.simplify(h*s.integrate(x**(-1-h),(x,1,s.oo))),1)
  # At h=1/100, even radius 2^100 leaves HALF the marginal I mass.
  self.assertEqual((s.Integer(2)**100)**(-s.Rational(1,100)),s.Rational(1,2))
  t=q.tail(Q(1,100),100,q.IDENTITY)
  self.assertEqual(t['row_prefactors']['I'],40*t['g'])
  self.assertEqual(t['Z'],800*t['g']);self.assertTrue(all(v<=t['Z'] for v in t['row_prefactors'].values()))
 def test_same_root_majorant_and_weighted_recurrence(self):
  r=q.root(q.tail(Q(1,100),100,q.IDENTITY),q.IDENTITY);a,b=r['a'],r['b']
  self.assertLessEqual(8*a*b,1);self.assertLess(a/7+b*a*a,a);self.assertLessEqual(2*b*a,Q(1,4))
  self.assertEqual(r['base'],512)
  Z=a/r['gamma'];d=[Q(0)]
  for m in range(1,21):
   d.append(r['gamma']*(Z*64**m*factorial(m)+r['qnorm']*sum(d[j]*d[m-j] for j in range(1,m))))
   self.assertLessEqual(d[m]*factorial(m),a*512**m*factorial(m)**2)
  # Exact scalar B=1/50,Q=100,z=g*eta, on the small accepted branch.
  # c=(-1+sqrt(1+10^6*g*eta))/10000; its endpoint value is safely small.
  import sympy as s
  z=s.symbols('z');g=s.Rational(1,10**10);c=(-1+s.sqrt(1+10**6*g*z))/10000
  self.assertEqual(s.simplify(c/50+100*c*c-g*z),0)
  for m in range(1,9):
   cap=q.bound(r,m)['upper'];self.assertLessEqual(abs(s.diff(c,z,m).subs(z,0)),s.Rational(cap.numerator,cap.denominator))
  self.assertLess(1-s.sqrt(1-10**6*g),s.Rational(6400,10**10)*10000)
 def test_relative_field_and_binary_budget(self):
  t=q.tail(Q(1,100),100,q.IDENTITY);r=q.root(t,q.IDENTITY);f=q.field(r,2)
  self.assertEqual(f['prefactor'],30720000*2*t['g']);self.assertEqual(f['base'],514)
  for recipe in (q.heat(Q(1,100),100),r,f):
   for m in range(9):self.assertLessEqual(q.bound(recipe,m)['upper'],2**q.bound(recipe,m,True)['upper'])
   self.assertIsInstance(q.bound(recipe,100005,True)['upper'],int)
  self.assertFalse(r['complex_slow_angle_extension_claimed'])
  self.assertFalse(t['actual_reference_attenuation_complete_tail_and_target_binding_certified'])
 def test_identity_domain_and_order_guards(self):
  t=q.tail(Q(1,100),100,q.IDENTITY)
  for h,X in ((0,100),(Q(1,10),100),(Q(1,100),99)):
   with self.assertRaises(ValueError):q.heat(h,X)
  with self.assertRaises(ValueError):q.tail(Q(1,100),100,'five-row-compact')
  with self.assertRaises(ValueError):q.root(dict(t,stage='same-leading-full-target-Gevrey3'),q.IDENTITY)
  with self.assertRaises(ValueError):q.root(dict(t,Z=0),q.IDENTITY)
  with self.assertRaises(ValueError):q.bound(q.heat(Q(1,100),100),100005)
  with self.assertRaises(ValueError):q.bound(q.heat(Q(1,100),100),True)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   child=subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertIn(b'real-heat-endpoint-Gevrey2',child.stdout)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
