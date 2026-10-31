import unittest,json
from fractions import Fraction as Q
from math import factorial
import chart as c
import audit

def coeff(v):return {j:dict(x.terms) for j,x in v.c.items()}
class Tests(unittest.TestCase):
 def test_series_oracles(self):
  import sympy as s
  t=s.symbols('t')
  for a in (Q(1,2),Q(2),Q(3,2)):
   for n in range(4):
    A=4*t+sum(3*abs(s.binomial(s.Rational(a.numerator,a.denominator),j))*t**j for j in range(1,n+1))
    for normalized,f in ((True,1/(1-A)),(False,s.exp(A))):
     series=s.series(f,t,0,n+1).removeO().expand()
     expected=factorial(n)*sum(series.coeff(t,j) for j in range(n+1))
     self.assertEqual(c.transfer(n,a,3,1,normalized),expected)
 def test_actual_mixed_jets(self):
  import sympy as s
  from itertools import product
  R,U,Z,T,y,z=s.symbols('R U Z T y z',positive=True)
  F=R**3+R*Z+T**2+y*z
  G=F.subs(R,s.sqrt(U));xs=(R,Z,T,y,z);us=(U,Z,T,y,z)
  for radius in (1,2):
   point={R:radius,U:radius**2,Z:1,T:1,y:1,z:1}
   for n in range(4):
    indices=[i for i in product(range(n+1),repeat=5) if sum(i)<=n]
    nr=0;mu=0
    for I in indices:
     fr=F;gu=G;fac=1
     for x,u,k in zip(xs,us,I):fr=s.diff(fr,x,k);gu=s.diff(gu,u,k);fac*=factorial(k)
     nr+=abs(fr.subs(point))/fac;mu=max(mu,abs(gu.subs(point)))
    forward=c.transfer(n,Q(1,2),2,1,True)
    reverse=c.transfer(n,2,4,1,False)
    self.assertLessEqual(mu,forward*nr)
    self.assertLessEqual(nr,reverse*mu)
 def test_roundtrip_provider(self):
  p=c.provider
  raw=p.CompactProvider(lambda n,v:dict(order=n,pole=v,cutoff_order=n+1,bounds={k:1 for k in p.KEYS}),lambda n,v:3,4,9)
  out=c.ChartedCompact(raw,d=2,Rmin=1,Rmax=2,Umin=1,Umax=4)
  for n in (0,1,3):
   a,b=out.factors(n);x=raw(n,7);y=out(n,7)
   for k in ('KT','KA'):self.assertEqual(coeff(y[k]),coeff(a*b*x[k]))
   self.assertEqual(coeff(y['I']),coeff(x['I']));self.assertEqual((y['order'],y['pole']),(n,7))
  self.assertEqual(out.factors(0),(1,1))
  with self.assertRaises(ValueError):out(5,7)
 def test_guards(self):
  with self.assertRaises(ValueError):c.transfer(2,1,2,0,True)
  with self.assertRaises(TypeError):c.transfer(2,1,2,1,1)
  with self.assertRaises(TypeError):c.ChartedCompact(lambda n,p:None,d=2,Rmin=1,Rmax=2,Umin=1,Umax=4)
 def test_norm_distinction(self):
  # Already for the identity chart, normalized input does not imply raw <= input.
  self.assertGreater(c.transfer(3,1,1,1,True),c.transfer(3,1,1,1,False))
  self.assertEqual(c.record()['derivative_order_loss'],0)
  self.assertFalse(c.record()['complete_source_residual_or_blowup_certified'])
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((c.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
