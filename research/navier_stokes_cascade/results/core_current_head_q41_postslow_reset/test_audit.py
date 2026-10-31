import unittest,json
from fractions import Fraction as Q
import postslow as p
import audit

def coeff(v):return {k:dict(x.terms) for k,x in v.c.items()}
def fixture():
 r=p.refresh;z=r.maps.update.profile.pipeline;P=z.roots.P
 root=z.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
 t=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
 fixed=lambda n,pole:{k:1 for k in z.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),
 slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=14,profile=[10000]*28,profile_hull=hull)
 return r.maps.NewPoleMaps(r.maps.update.SignedCovariance(t,11))

class Tests(unittest.TestCase):
 def test_independent_reset_remainders(self):
  # Small manufactured fields test the inherited complete scalar polynomial.
  old=dict(x=p.EP(2),u=p.EP(3),z=p.EP(5));inc=dict(x=p.EP(7),u=p.EP(11),z=p.EP(13))
  n=2;b=p.EP(17);M=19;v=p.reset(n,old,inc,b,M)
  S=2*(17+2)*7+7**2;H=(17+2)*13+5*7+7*13;T=2*3*11+11**2
  radial=(n+1)*(n+3)*7+S+T+(n+1)*H
  expected=dict(P=M*radial,Jtheta=M*(5*11+3*13+11*13),Jz=M*(2*5*13+13**2)+Q(M*radial,2))
  for k in expected:self.assertEqual(coeff(v[k]),coeff(p.EP(expected[k])))
 def test_only_slow_crosses_can_use_Haar(self):
  import sympy as s
  t=s.symbols('t',real=True);old=3+s.cos(t);slow=2;fast=2*s.cos(t)
  avg=lambda f:s.integrate(f,(t,0,2*s.pi))/(2*s.pi)
  self.assertEqual(avg(old*slow),avg(old)*slow)
  self.assertNotEqual(avg(old*fast),avg(old)*avg(fast))
  self.assertEqual(avg(fast),0)
 def test_exact_whole_moving_pressure_product(self):
  import sympy as s
  Z=s.symbols('Z');c=1+Z**2;P=2+Z;J=3-Z
  self.assertEqual(s.expand(s.diff(J+c*P,Z)-(s.diff(J,Z)+s.diff(c,Z)*P+c*s.diff(P,Z))),0)
  self.assertNotEqual(s.expand(s.diff(J+c*P,Z)-(s.diff(J,Z)+c*s.diff(P,Z))),0)
 def test_factory_full_ledger_and_identity(self):
  current=fixture();tail=lambda n,v:dict(order=n,pole=v,Ez=1)
  compact=lambda n,v:dict(order=n,pole=v,KT=1,KA=1,I=1)
  out=p.build(current,tail,compact,0);c=out.current
  self.assertIs(c.fast.rows.current,c.current)
  self.assertEqual((c.fast.rows.max_row,c.fast.max_output,c.max_output),(10,5,3))
  self.assertEqual(current.target.max_order,14)
  self.assertEqual(current.profile.max_order,27)
  self.assertEqual(out.moment,(6-1)*(1+1+6+36))
  # Guards are checked without evaluating the enormous coefficient programs.
  with self.assertRaises(ValueError):p.build(current,tail,compact,1)
  with self.assertRaises(ValueError):out.scalars(2)
  with self.assertRaises(ValueError):out.axial(1)
  independent=p.refresh.CurrentScalars(current,lambda n:{},3)
  with self.assertRaises(TypeError):p.Postslow(independent,0)
 def test_scope_and_rate_boundary(self):
  a=2-Q(3,100000);r=p.record()
  self.assertEqual(r['postslow_order'],str(a+1))
  self.assertEqual(r['whole_axial_order'],str(a+2))
  self.assertLess(a+1,2*a)
  self.assertFalse(r['actual_high_jets_full_rows_source_residual_or_blowup_certified'])
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
