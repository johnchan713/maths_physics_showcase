import unittest,json
from fractions import Fraction as Q
import radial as r
import audit

def coeff(v):return {j:dict(x.terms) for j,x in v.c.items()}
class Tests(unittest.TestCase):
 def test_hulls_both_sides_of_one(self):
  import sympy as s
  for d in (Q(1,2),Q(1),Q(2),Q(7,3)):
   for a,b in ((Q(1,4),Q(1,2)),(Q(1,2),Q(3)),(Q(2),Q(4))):
    g=r.RadialGeometry(d,a,b)
    for x in (a,b):
     value=s.Rational(x.numerator,x.denominator)**s.Rational(d.numerator,d.denominator)
     self.assertLessEqual(g.Umin,value);self.assertLessEqual(value,g.Umax)
 def test_exact_fractional_multiplier_jets(self):
  import sympy as s
  U=s.symbols('U',positive=True)
  for d,lo,hi in ((Q(1,2),Q(1),Q(4)),(Q(2),Q(1,2),Q(2))):
   g=r.RadialGeometry(d,lo,hi);D=s.Rational(d.numerator,d.denominator)
   for e in (0,1,2):
    fs=dict(f=U**((e-D+1)/D)/D,g=U**(-e/D),k=D*U**((D-1-e)/D))
    for n in range(5):
     bounds=g.multipliers(n,e)
     for key,F in fs.items():
      for j in range(n+1):
       for R in (lo,hi):
        point=s.Rational(R.numerator,R.denominator)**D
        self.assertLessEqual(abs(s.diff(F,U,j).subs(U,point)),bounds[key])
 def test_raw_falling_factorial(self):
  g=r.RadialGeometry(1,1,1)
  self.assertEqual(g.power(4,-1),24)
  self.assertEqual(g.power(4,0),1)
  self.assertGreater(g.power(4,-1),g.power(0,-1))
 def test_bound_provider_same_hull(self):
  def remaining(n,p):return dict(order=n,pole=p,cutoff_order=n+1,bounds={k:1 for k in r.REMAINING})
  out=r.build(Q(1,2),Q(1,4),2,remaining,lambda n,p:3,4,10)
  self.assertEqual((out.Umin,out.Umax),(Q(1,64),8))
  g=r.RadialGeometry(Q(1,2),Q(1,4),2)
  for n in (0,2):
   expected=r.provider.constants(n,dict(remaining(n,7)['bounds'],**g.multipliers(n)))
   a,b=out.factors(n);actual=out(n,7)
   for k in ('KT','KA'):self.assertEqual(coeff(actual[k]),coeff(a*b*expected[k]))
   self.assertEqual((actual['order'],actual['pole']),(n,7))
 def test_guards_and_no_override(self):
  with self.assertRaises(ValueError):r.RadialGeometry(3,1,2)
  with self.assertRaises(ValueError):r.RadialGeometry(1,0,2)
  with self.assertRaises(ValueError):r.RadialGeometry(1,1,2).multipliers(0,3)
  def bad(n,p):return dict(order=n,pole=p,cutoff_order=n+1,bounds={k:1 for k in r.provider.KEYS})
  with self.assertRaises(ValueError):r.build(1,1,2,bad,lambda n,p:1,1,1)(0,1)
 def test_evidence(self):
  self.assertEqual(audit.generate(),json.loads((r.HERE/'evidence.json').read_text()))
  self.assertFalse(r.record()['initial_A1_tail_full_source_residual_or_blowup_certified'])
if __name__=='__main__':unittest.main()
