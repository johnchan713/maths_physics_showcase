import unittest,json
from fractions import Fraction as Q
import refresh_scalars as r
import audit

def coeff(v):return {k:dict(x.terms) for k,x in v.c.items()}
def fixture():
 p=r.maps.update.profile.pipeline;P=p.roots.P
 root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
 t=r.maps.update.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),Rmin=1,Rmax=8,
 fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},compact=lambda n,pole:dict(KT=1,KA=1,I=6),
 slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},hull=hull,fast_hull=(1,8),max_output=4,profile=[10000]*18,profile_hull=hull)
 return r.maps.NewPoleMaps(r.maps.update.SignedCovariance(t,1))

class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.calls=[]
  def fast(n):
   cls.calls.append(n);return dict(order=n,pole=6,x=1,u=2,z=3)
  cls.v=r.CurrentScalars(fixture(),fast,0)
 def test_current_scalar_addition_once(self):
  v=self.v;n=0;o=v.old_full(1);f=v.fast_unweighted(1);M=v.moment
  x,u,z=(f[k] for k in ('x','u','z'));ox,ou,oz=(o[k] for k in ('x','u','z'))
  S=2*ox*x+x*x;T=2*ou*u+u*u;H=ox*z+oz*x+x*z
  change=dict(P=M*(S+T+H),Jtheta=M*(oz*u+ou*z+u*z),Jz=M*(2*oz*z+z*z)+M*(S+T+H)*Q(1,2))
  current=v.current.pressure(0)['updated_scalars']
  for k in current:self.assertEqual(coeff(v.postfast(0)[k]),coeff(current[k]+change[k]))
  self.assertNotEqual(coeff(change['P']),coeff(M*(x*x+u*u+x*z)))
 def test_full_old_mean_and_separate_poles(self):
  v=self.v;b=v.initial.budget;old=v.old_full(1)
  edge=r.maps.update.signed.maps.edge_absorption(b.pole,v.current.target.width)
  f=b._fast(1);s=v.initial.slow_fields(1,False)
  for k,j in [('x','X'),('u','U'),('z','Z')]:self.assertEqual(coeff(old[k]),coeff(edge*r.convert(f[k])+r.convert(s[j])))
  self.assertNotEqual(b.pole,6)
  E=r.maps.update.signed.maps.edge_absorption(6,v.current.target.width)
  self.assertEqual(coeff(v.fast_unweighted(1)['z']),coeff(r.EP(3*E)))
 def test_zero_haar_does_not_remove_old_correlations(self):
  import sympy as s
  t=s.symbols('t',real=True);old=s.cos(t);inc=2*s.cos(t)
  avg=lambda f:s.integrate(f,(t,0,2*s.pi))/(2*s.pi)
  self.assertEqual(avg(inc),0)
  self.assertEqual(avg((old+inc)**2-old**2),4)
  self.assertEqual(avg(inc**2),2)
 def test_slow_uses_fresh_next_scalar_and_same_hull(self):
  v=self.v;a=v.postfast(0);an=v.postfast(1);out=v.slow(0)
  self.assertEqual(coeff(out['u']),coeff(a['P']+a['Jz']))
  self.assertEqual(coeff(out['z']),coeff(a['Jtheta']))
  self.assertEqual(coeff(out['x']),coeff(an['Jtheta']))
  self.assertEqual(set(self.calls),{1,2})
  weighted=v.slow(0,True);floor=r.convert(v.initial.budget.floor)
  for k in out:self.assertEqual(coeff(weighted[k]),coeff(floor*out[k]))
 def test_guards_and_rates(self):
  v=self.v
  with self.assertRaises(ValueError):r.CurrentScalars(v.current,lambda n:{},1)
  with self.assertRaises(ValueError):v.old_full(-1)
  with self.assertRaises(ValueError):v.fast_unweighted(-1)
  with self.assertRaises(ValueError):v.postfast(2)
  with self.assertRaises(ValueError):v.slow(1)
  bad=r.CurrentScalars(v.current,lambda n:dict(order=0,pole=6,x=1,u=1,z=1),0)
  with self.assertRaises(ValueError):bad.fast_unweighted(1)
  self.assertEqual(r.record()['postfast_order'],str(2-Q(3,100000)))
  self.assertFalse(r.record()['postslow_reset_full_source_residual_or_blowup_certified'])
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
