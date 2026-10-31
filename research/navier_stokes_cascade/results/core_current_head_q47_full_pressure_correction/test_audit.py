import unittest,json
from fractions import Fraction as Q
import full_rows as q
import audit
def coefficients(v):return {j:dict(c.terms) for j,c in v.c.items()}
def fixture(M):
 r=q.refresh;p=r.maps.update.profile.pipeline;P=p.roots.P
 root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
 t=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
 fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),
 slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=M+3,profile=[10000]*(M+17),profile_hull=hull)
 return r.maps.NewPoleMaps(r.maps.update.SignedCovariance(t,M))


class Tests(unittest.TestCase):
 def test_three_full_pressure_contributions(self):
  c=fixture(1);r=q.FullRows(c,lambda n,p:dict(order=n,pole=p,Ez=7),0,d=2,c0=2,c1=4,shape=[1]*5,time_inverse=lambda n,p:1)
  o=r.rows(0);K=r.pressure(3,r.pole)['KT'];rho=c.profile.weighted_rho(3)
  self.assertNotEqual(coefficients(K),coefficients(c.recipe(3,r.pole)['Kp']))
  for key in ('old_fast','old_slow','signed_wave'):
   v=o[key];self.assertEqual(coefficients(v['pm']),coefficients(K*(v['gr']+rho*v['P'])))
  self.assertEqual(coefficients(o['Ez']),coefficients(2*(7+sum((o[k]['Ez'] for k in ('old_fast','old_slow','signed_wave')),q.EP()))))
  self.assertTrue(o['old_fast']['viscosity'].c)
 def test_full_factory_and_ceiling(self):
  c=fixture(11);args=dict(d=1,c0=2,c1=4,time_inverse=lambda n,p:1)
  tail=lambda n,p:dict(order=n,pole=p,Ez=7)
  with self.assertRaises(ValueError):q.build(c,tail,0,shape=[1]*14,**args)
  a=q.build(c,tail,0,shape=[1]*15,**args)
  rows=a.current.fast.rows;self.assertIsInstance(rows,q.FullRows)
  self.assertIs(rows.current,c)
  self.assertNotEqual(coefficients(rows.pressure(0,rows.pole)['KT']),coefficients(a.current.fast.compact(0,rows.pole)['KT']))
 def test_haar_compatibility_is_not_modewise(self):
  import sympy as s
  R,y=s.symbols('R y');v=30*R**2*(1-R)**2*s.cos(y)
  # Finite integral control; the proof uses an interior smooth bump.
  self.assertEqual(s.integrate(v,(y,0,2*s.pi)),0)
  self.assertEqual(s.integrate(v,(R,0,1)),s.cos(y))
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((q.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
