import unittest,json,math
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
import audit,q116

def constant_coefficients(N=3,flat=2):
 G=[1,3,6,6][:N+1]
 return q116.coefficients(G,[1]+[0]*N,[0]*(N+1),[0]*(N+1),[0]*(N+1),[1]+[0]*N,
                         q0=1,T=1,L=1,dmin=1,entrance_flat_through=flat)

class Tests(unittest.TestCase):
 def test_exact_second_and_mixed_pulse_projection_identity(self):
  a,t=s.symbols('a t',real=True);n=s.Matrix([1,a+t+a*t,0]);np=n.diff(t)
  N=(n.T*n)[0];Pi=s.eye(3)-n*n.T/N
  K=s.Matrix([[1,a,0],[t,2,1],[0,a,3]]);A=-Pi*K-n*np.T/N
  d=2+a*a+t*t;p=a+t;m=2;B=A-(m*m*d+p)*s.eye(3)
  u=t**4*Pi*s.Matrix([1,a,t]);g=B*u-u.diff(t)+t**4*n*(a+1)
  Qlift=Pi*(K*n+np)/N
  @lru_cache(None)
  def D(name,alpha):
   v=dict(n=n,u=u,A=A,d=d,p=p,Pi=Pi,g=g)[name]
   return s.diff(v,a,alpha[0],t,alpha[1])
  for alpha in ((2,0),(1,1),(0,2),(0,3)):
   indices=set(beta for beta,_,_ in q116.partitions(alpha))
   tables={name:{b:D(name,b) for b in indices} for name in ('n','u','A','d','p','Pi','g')}
   out=q116.projected_forcing(alpha,**tables,Qlift=Qlift,m=m)
   raw=D('u',alpha);w=Pi*raw
   residual=w.diff(t)-B*w-out['forcing']
   reconstruction=raw-w+n/N*out['constraint']
   for point in ({a:0,t:0},{a:Q(1,3),t:Q(1,2)}):
    self.assertEqual(residual.subs(point).applyfunc(s.simplify),s.zeros(3,1))
    self.assertEqual(reconstruction.subs(point).applyfunc(s.simplify),s.zeros(3,1))
   self.assertEqual(raw.subs(t,0),s.zeros(3,1))

 def test_scalar_recurrence_matches_damped_zero_entry_model(self):
  out=constant_coefficients();t=s.symbols('t',nonnegative=True)
  self.assertEqual(out['velocity_over_residence'],list(map(Q,[1,3,6,6])))
  for m in (1,-1,2,-2,5,20):
   v=q116.specialize(out,m);lam=m*m
   u=-t**3/lam+3*t*t/lam**2-6*t/lam**3+6*(1-s.exp(-lam*t))/lam**4
   self.assertEqual(s.simplify(s.diff(u,t)+lam*u+t**3),0)
   for j in range(4):
    for x in (0,Q(1,4),1):
     actual=abs(float(s.diff(u,t,j).subs(t,x)))
     self.assertLessEqual(actual,float(v['velocity'][j])+1e-12)
   self.assertFalse(v['actual_source_or_high_jet_contracts_certified'])

 def test_entrance_gate_and_boundary_layer_obstruction(self):
  with self.assertRaises(ValueError):constant_coefficients(flat=1)
  with self.assertRaises(ValueError):constant_coefficients(flat=True)
  out=constant_coefficients()
  self.assertFalse(out['pressure_zero_extension_flatness_sufficient'])
  self.assertTrue(constant_coefficients(flat=3)['pressure_zero_extension_flatness_sufficient'])
  t=s.symbols('t',nonnegative=True)
  for m in (2,10,100):
   u=-(1-s.exp(-m*m*t))/(m*m)
   self.assertEqual(s.diff(u,t).subs(t,0),-1)
   self.assertLess(q116.first.residence(m,1,1),Q(1))

 def test_full_normal_pressure_has_no_residence_factor(self):
  # n=e1, K=0, g=a^2 e1, u=0: second normalized pressure derivative is 2i/(km).
  for m in (1,2,10):
   bounds=q116.pressure([1,2,2],[1,0,0],[0,0,0],[0,0,0],[1,0,0],[0,0,0],m)
   self.assertEqual(bounds,[Q(1,m),Q(2,m),Q(2,m)])
   self.assertGreater(bounds[2],2*q116.first.residence(m,1,1)/m if m>1 else 0)

 def test_nonconstant_geometry_and_all_coefficient_leaves(self):
  N=4
  out=q116.coefficients([1]*(N+1),[2]*(N+1),[3]*(N+1),[2]*(N+1),[1]*(N+1),[4]*(N+1),
    q0=2,T=3,L=5,dmin=Q(1,4),entrance_flat_through=4)
  self.assertEqual(len(out['velocity_over_residence']),N+1)
  self.assertTrue(all(x>0 for x in out['velocity_over_residence']))
  for m in (1,2,100):
   r=q116.specialize(out,m)
   self.assertLessEqual(m*m*r['residence'],out['m2_residence_upper'])
   self.assertEqual(r['velocity'],[r['residence']*v for v in out['velocity_over_residence']])
   self.assertTrue(all(x>0 for x in r['pressure_times_carrier']))
  self.assertFalse(out['actual_native_or_entrance_contracts_certified'])

 def test_multiindex_combinatorics_and_invalid_inputs(self):
  for alpha in ((2,3),(1,1,2)):
   self.assertEqual(sum(c for b,g,c in q116.partitions(alpha)),2**sum(alpha))
  for alpha in ((),(True,),(-1,)):
   with self.assertRaises(ValueError):list(q116.partitions(alpha))
  with self.assertRaises(ValueError):q116.jets([1.0])
  with self.assertRaises(ValueError):q116.convolution([1],[1,2])
  with self.assertRaises(ValueError):q116.specialize(constant_coefficients(),0)

 def test_record(self):self.assertEqual(audit.generate(),json.loads((audit.HERE/'evidence.json').read_text()))
if __name__=='__main__':unittest.main()
