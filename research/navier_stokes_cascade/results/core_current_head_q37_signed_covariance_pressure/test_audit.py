import unittest,json
from fractions import Fraction as Q
import signed_update as c
import audit

def coefficients(v):return {k:dict(x.terms) for k,x in v.c.items()}

def fixture():
 p=c.profile.pipeline;P=p.roots.P
 root=p.roots.RootBudget(H=lambda j:P(int(j==0)),target=lambda j:P(int(j==0)),wave=lambda j:P(int(j==0)),curl=lambda j:P(int(j==0)),loghalf=lambda j:P(),K0=1,beta=1,R=1,mass=2,width=Q(1,4))
 hull=dict(l=2,u=3,Xa=1,Xb=8,t1=root.width)
 return c.profile.build(root,lambda n:dict(target=Q(1),correction=Q(1),pole=3*n+3),
  Rmin=1,Rmax=8,fixed=lambda n,pole:{k:1 for k in p.fresh.FIXED if k!='rho'},
  compact=lambda n,pole:dict(KT=1,KA=1,I=6),slow=lambda n:{k:1 for k in ('H','Kp','Kz','Kg','Kpsi')},
  hull=hull,fast_hull=(1,8),max_output=3,profile=[10000]*17,profile_hull=hull)

class Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.v=c.SignedCovariance(fixture(),0)
  cls.calls=[]
  def maps(n,p):
   cls.calls.append((n,p));return dict(order=n,pole=p,Kp=1,primitive=1,bump=1)
  cls.out=cls.v.pressure(0,maps)
 def test_exact_six_entry_phase_identity(self):
  import sympy as s
  A=s.Matrix([1,2,3]);B=s.Matrix([2,-1,4]);L=s.Matrix([-3,1,2]);S=s.Matrix([1,4,-2])
  old=(A*A.T+B*B.T)/2
  new=((A+L)*(A+L).T+(B+S)*(B+S).T)/2
  principal=(A*L.T+L*A.T)/2
  error=(B*S.T+S*B.T+L*L.T+S*S.T)/2
  for i in range(3):
   for j in range(i,3):self.assertEqual((new-old-principal-error)[i,j],0)
  self.assertNotEqual(error,(B*S.T+S*B.T+S*S.T)/2)
  # A later cosine error has an additional surviving cross with L.
  E=s.Matrix([1,0,0]);self.assertNotEqual((E*L.T+L*E.T)/2,s.zeros(3))
 def test_order_and_pressure_pole(self):
  k=Q(1,100000);alpha=2-3*k
  self.assertEqual(min(alpha+1-2*k,2*alpha-1,2*alpha-2*k),3-6*k)
  v=self.v;pole=2*max(21,18+v.target.pole)
  self.assertEqual(self.calls,[(1,pole)])
  self.assertGreater(pole,v.target.pole)
  self.assertEqual(v.covariance(2)['pole'],pole)
  with self.assertRaises(ValueError):v.pressure(0,lambda n,p:dict(order=n,pole=v.target.pole,Kp=1,primitive=1,bump=1))
 def test_full_pressure_memory(self):
  v=self.v;o=self.out;old=v.target.source.haar.initial
  for k,key in [('P','moment_P'),('Jtheta','Jtheta'),('Jz','Jz')]:
   a=c.profile.rebind(old.output(0)[k],c.signed.projector.haar.sharp.fresh)
   self.assertEqual(coefficients(o['updated_scalars'][k]),coefficients(a+o['changes'][key]))
  self.assertNotEqual(coefficients(o['whole_updated_pressure_moment']),coefficients(o['whole_pressure_change_moment']))
 def test_all_covariance_terms_and_transfer(self):
  v=self.v;a=v.coefficients(0);t=c.EP(v.primary.chart.transfer(0))
  expected=t*(2*a['B']*a['S']+a['L']*a['L']+a['S']*a['S'])
  self.assertEqual(coefficients(v.covariance(0)['error']),coefficients(expected))
  self.assertTrue(a['S'].c)
  self.assertIs(v.target.root,v.primary.primary.root)
 def test_guards_and_scope(self):
  with self.assertRaises(ValueError):c.SignedCovariance(self.v.target,1)
  with self.assertRaises(ValueError):self.v.coefficients(3)
  with self.assertRaises(ValueError):self.v.pressure(1,lambda *a:{})
  self.assertFalse(c.record()['complete_residual_or_blowup'])
 def test_evidence(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))

if __name__=='__main__':unittest.main()
