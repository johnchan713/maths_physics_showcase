import unittest,json
from fractions import Fraction as Q
from math import comb,factorial
import current_fast as m
import audit

def coeff(v):return {k:dict(x.terms) for k,x in v.c.items()}
def fixture():
 r=m.refresh;p=r.maps.update.profile.pipeline;P=p.roots.P
 root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
 t=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
 fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),
 slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=4,profile=[10000]*18,profile_hull=hull)
 return r.maps.NewPoleMaps(r.maps.update.SignedCovariance(t,1))

class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.current=fixture();cls.v=m.CurrentRows(cls.current,lambda n,p:dict(order=n,pole=p,Ez=7),0)
  cls.out=cls.v.rows(0)
 def test_all_three_current_contributions_and_pressure(self):
  o=self.out
  for k in ('old_fast','old_slow','signed_wave'):
   self.assertTrue(o[k]['Ez'].c)
   self.assertTrue(o[k]['pm'].c)
  self.assertTrue(o['old_fast']['viscosity'].c)
  self.assertEqual(coeff(o['Etheta']),coeff(2*sum((o[k]['Etheta'] for k in ('old_fast','old_slow','signed_wave')),m.EP())))
 def test_initial_A1_survives_projection(self):
  zero=m.CurrentRows(self.current,lambda n,p:dict(order=n,pole=p,Ez=0),0).rows(0)
  self.assertEqual(coeff(self.out['Ez']),coeff(zero['Ez']+14))
  self.assertEqual(coeff(self.out['Etheta']),coeff(zero['Etheta']))
 def test_pole_and_fixed_root(self):
  v=self.v
  self.assertEqual(v.q,max(6*v.initial.budget.pole,27,24+v.current.target.pole))
  self.assertGreaterEqual(v.pole,v.current.state.covariance(3)['pole'])
  self.assertIs(v.current.target.root,v.current.state.primary.primary.root)
 def test_independent_fast_factor_and_derivative_calls(self):
  calls=[];ops=[]
  def rows(n):calls.append(n);return dict(Etheta=m.EP(2),Ez=m.EP(3))
  def compact(n,p):ops.append((n,p));return dict(order=n,pole=p,KT=n+2,KA=5,I=7)
  n=2;v=m.fast_recipe(n,rows,compact,100)
  t=comb(n+5,5)*24*2**(n+4)*factorial(n+4)*7
  x=comb(n+5,5)*24*2**(n+5)*factorial(n+5)*7*(n+3)
  self.assertEqual(coeff(v['u']),coeff(m.EP({1:2*t})))
  self.assertEqual(coeff(v['z']),coeff(m.EP({1:18*t})))
  self.assertEqual(coeff(v['x']),coeff(m.EP({1:3*x})))
  self.assertEqual(calls,[n+4,n+5]);self.assertEqual(ops,[(n,100),(n+1,100)])
 def test_domain_guards_and_scope(self):
  with self.assertRaises(ValueError):m.CurrentRows(self.current,lambda *a:{},1)
  with self.assertRaises(ValueError):self.v.rows(1)
  with self.assertRaises(ValueError):m.FastFields(self.v,lambda *a:{},0)
  with self.assertRaises(ValueError):m.CurrentRows(self.current,lambda n,p:dict(order=n,pole=p+1,Ez=1),0).rows(0)
  with self.assertRaises(ValueError):m.fast_recipe(0,lambda n:{},lambda n,p:dict(order=n,pole=p-1,KT=1,KA=1,I=1),100)
  self.assertFalse(m.record()['actual_fast_source_postslow_or_full_residual_certified'])
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
