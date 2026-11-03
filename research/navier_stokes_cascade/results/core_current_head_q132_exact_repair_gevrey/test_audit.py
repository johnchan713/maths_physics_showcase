import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q132 as q
import audit
def memory():
 p=dict(stage='real-zero-mean-primitive',prefactor=Q(4),base=Q(2))
 r=q.previous.modulation(4,3,4,Q(1,4),p,p,1,8,q.previous.IDENTITY)
 return q.previous.memory(r,1,3,{k:(Q(1),Q(2)) for k in q.ROWS})
def target():return q.target(memory(),{k:(Q(0),Q(1)) for k in q.ROWS},2,Q(1,2),Q(1,100),q.IDENTITY)
class Tests(unittest.TestCase):
 def test_normalizers_and_five_physical_rows(self):
  import sympy as s
  z=s.symbols('z');estar=2;rho=s.Rational(1,2);lam=s.Rational(1,100);e=estar/(1+z*z)
  for p,caps in [(1/e,[2,2,1]),(1/e**2,[4,8,8,4,1])]:
   scale=estar if p==1/e else estar**2
   for m in range(7):
    for eta in (-1,0,1):
     value=abs(s.diff(p,z,m).subs(z,eta))/s.factorial(m)
     cap=s.Rational(caps[m],scale) if m<len(caps) else 0
     self.assertLessEqual(value,cap)
  d=[(i+1)*(1+z+z**2)/100 for i in range(5)]
  scales=[rho*e,s.sqrt(2)*rho**s.Rational(3,2)*e,s.sqrt(2)*rho**s.Rational(3,2)*e**2,rho*e**2,e**2]
  m,i,j,ss,cp=[f/w for f,w in zip(d,scales)]
  mapped=[m,(m-j)/lam,i,-ss,cp]
  restored=[mapped[0]*scales[0],mapped[2]*scales[1],(mapped[0]-lam*mapped[1])*scales[2],-mapped[3]*scales[3],mapped[4]*scales[4]]
  for f,g in zip(d,restored):self.assertEqual(s.simplify(f-g),0)
  self.assertNotEqual(s.simplify(j-d[2]/(scales[2]*e**(-1))),0)
 def test_full_target_pressure_and_lambda_cost(self):
  base={k:(Q(0),Q(1)) for k in q.ROWS};m=memory();a=q.target(m,base,2,Q(1,2),Q(1,100),q.IDENTITY)
  base['Pi']=(Q(7),Q(3));b=q.target(m,base,2,Q(1,2),Q(1,100),q.IDENTITY)
  self.assertEqual(b['transformed_row_prefactors'][4]-a['transformed_row_prefactors'][4],175)
  other=q.target(m,base,2,Q(1,2),Q(1,200),q.IDENTITY)
  self.assertEqual(other['transformed_row_prefactors'][1],2*b['transformed_row_prefactors'][1])
  with self.assertRaises(ValueError):q.target(m,dict(M=(0,1),I=(0,1),J=(0,1),S=(0,1)),2,1,Q(1,100),q.IDENTITY)
 def test_majorant_contraction_zero_and_linear_limits(self):
  for Z in (Q(0),Q(1,2**100),Q(1),Q(2**100)):
   for qq in (0,1,100000):
    r=q.root(dict(stage='same-leading-full-target-Gevrey3',Z=Z,L=Q(3)),Q(1,10**10),beta=1,qnorm=qq)
    a,b,T,R=r['a'],r['b'],r['majorant_ball'],r['majorant_radius']
    if a:
     self.assertLessEqual(a*3*R/(1-3*R),T/6);self.assertLessEqual(b*T*T,T/8);self.assertLessEqual(2*b*T,Q(1,4))
    else:self.assertEqual(q.bound(r,5)['upper'],0)
  import sympy as s
  x=s.symbols('x');a=s.Rational(1,10);b=2;L=3
  V=(1-s.sqrt(1-4*b*a*L*x/(1-L*x)))/(2*b)
  self.assertEqual(s.simplify(V-a*L*x/(1-L*x)-b*V*V),0)
 def test_same_small_scalar_root_large_high_jets(self):
  import sympy as s
  z=s.symbols('z');k=16;eps=s.Rational(1,k**3);c=(s.sqrt(1+4*eps*s.sin(k*z))-1)/2
  r=q.root(dict(stage='same-leading-full-target-Gevrey3',Z=Q(1,k**3),L=Q(k)),Q(2,k*k),beta=1,qnorm=1)
  self.assertEqual(s.simplify(c+c*c-eps*s.sin(k*z)),0)
  self.assertGreater(abs(s.diff(c,z,3).subs(z,0)),s.Rational(2,k*k))
  for m in range(1,8):
   cap=q.bound(r,m)['upper'];self.assertLessEqual(abs(s.diff(c,z,m).subs(z,0)),s.Rational(cap.numerator,cap.denominator))
  f=q.field(r,2);ec=2*c/(1+z*z)
  for m in range(7):
   cap=f['prefactor']*f['base']**m*factorial(m)**3
   self.assertLessEqual(abs(s.diff(ec,z,m).subs(z,0)),s.Rational(cap.numerator,cap.denominator))
 def test_independent_vector_quadratic_jets_and_square_controls(self):
  def quad(a,b):return [0,0,0,-a[0]*b[0]-a[1]*b[1]+sum(a[j]*b[j]/2 for j in range(2,5)),sum(a[j]*b[j]/2 for j in range(2,5))]
  n=6;c=[[Q(0)]*5 for _ in range(n+1)];z=[[Q(0)]*5 for _ in range(n+1)]
  z[1]=[Q(j+1,10**6) for j in range(5)];z[2]=[Q(j+1,10**7) for j in range(5)]
  r=q.root(dict(stage='same-leading-full-target-Gevrey3',Z=Q(1,100),L=Q(2)),Q(1,100),beta=1,qnorm=5)
  for m in range(1,n+1):
   nonlinear=[sum(quad(c[j],c[m-j])[i] for j in range(1,m)) for i in range(5)]
   c[m]=[z[m][i]-nonlinear[i] for i in range(5)]
   self.assertEqual([c[m][i]+nonlinear[i]-z[m][i] for i in range(5)],[0]*5)
   for v in c[m]:self.assertLessEqual(abs(v)*factorial(m),q.bound(r,m)['upper'])
  self.assertNotEqual(quad(c[1],c[1])[4],0)
  self.assertNotEqual(-c[1][0]**2-c[1][1]**2,0)
 def test_guards_and_binary_budgets(self):
  r=q.root(target(),Q(1,10**12))
  self.assertEqual(q.bound(r,0)['upper'],Q(1,10**12));self.assertFalse(r['slow_analytic_extension_claimed'])
  for m in range(7):self.assertLessEqual(q.bound(r,m)['upper'],2**q.bound(r,m,True)['upper'])
  self.assertEqual(q.bound(r,100005,True)['factorial_power'],3)
  with self.assertRaises(ValueError):q.bound(r,100005)
  with self.assertRaises(ValueError):q.root(target(),1)
  with self.assertRaises(ValueError):q.root(target(),0,identity='earlier eta-varying quadratic')
  with self.assertRaises(ValueError):q.target(memory(),{k:(0,1) for k in q.ROWS},2,1,Q(1,100),'positive-order-linear')
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
