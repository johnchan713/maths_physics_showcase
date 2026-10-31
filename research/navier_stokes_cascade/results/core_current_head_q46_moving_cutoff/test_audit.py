import unittest,json
from fractions import Fraction as Q
import moving as m
import audit
class Tests(unittest.TestCase):
 def test_exp_and_shape_coefficients(self):
  import sympy as s
  c=m.Moving(2,1,2,3,Q(1,4),[1,2,3,4])
  t=s.symbols('t');L=sum((Q(1,j)+sum(Q(32**i*64**(j-i),3) for i in range(j+1)))*t**j for j in (1,2,3))
  V=3*(s.series(s.exp(L),t,0,4).removeO()-1)
  outer=s.Poly(1+2*V+3*V**2+4*V**3,t)
  self.assertEqual(c.raw(3),tuple(outer.nth(i)*s.factorial(i) for i in range(4)))
 def test_mixed_time_radial_slice(self):
  import sympy as s
  U,T=s.symbols('U T',positive=True)
  # On the Z=0 slice s=T, independently of the physical h.
  X=U**2/(2*T);F=X**2
  c=m.Moving(1,1,2,8,1,[64,16,1,0,0]);raw=c.raw(4)
  for n in range(5):
   for k in range(n+1):
    self.assertLessEqual(abs(s.diff(F,U,k,T,n-k).subs({U:2,T:Q(1,2)})),raw[n])
  self.assertNotEqual(s.diff(F,U,T).subs({U:2,T:1}),0)
 def test_extra_jet_and_constant(self):
  c=m.Moving(1,1,2,8,1,[2,0,0]);self.assertEqual(c(1,7)['bounds'],dict(chi=2,one_minus_chi=3,chi_U=0))
  with self.assertRaises(ValueError):c(2,7)
 def test_full_binding(self):
  kw=dict(d=1,Rmin=1,Rmax=8,Xa=1,Xb=8,c0=2,c1=4,width=Q(1,4),time_inverse=lambda n,p:2,max_order=2,max_pole=3)
  with self.assertRaises(ValueError):m.build(**kw,shape=[1,2,3])
  v=m.build(**kw,shape=[1,2,3,4])(2,3);self.assertEqual((v['order'],v['pole']),(2,3))
  self.assertTrue(v['KA'].c[0].terms);self.assertTrue(v['KT'].c[0].terms)
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((m.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
