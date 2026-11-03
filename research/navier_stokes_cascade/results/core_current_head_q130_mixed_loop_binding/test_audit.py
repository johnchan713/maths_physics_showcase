import unittest
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb
import q130 as q
import audit
def fixture():
 state=q.cutoff.state_caps(Q(1,200),1,4,Q(1,8),dict(E=4,U=8,EX=2,UX=2,invE=3),dict(M=3,I=2,J=3,S=3,Pi=7))
 recipe=q.cutoff.parameters(state,Q(1,2),Q(1,4),2,9,Q(1,32),q.cutoff.IDENTITY)
 return state,recipe,q.scalar_jets(state,recipe,4,Q(1,2))
class Tests(unittest.TestCase):
 def test_parameter_disk_margins(self):
  for G,E in [(2,9),(10,3**10),(2**100,2**300)]: # third is a manufactured algebra control, not an exp bound
   d=q.parameter_disk(G,E,q.IDENTITY);r0,T,r=d['r0'],d['T'],d['r']
   dz=2*G*r0+r0*r0
   self.assertLessEqual(dz,3*G*r0);self.assertLessEqual(2*E*dz,Q(3,32))
   dt=14*T*r/r0;self.assertLessEqual(dt,1/(8*Q(G+1)))
   self.assertLess(dt*(2*G+dt),Q(1,4))
   self.assertLessEqual(4*G*G*T+4*G**3,d['F']);self.assertLessEqual(3*G,d['F'])
 def test_scalar_jet_binding_and_augmented_variance(self):
  import sympy as s
  state,recipe,jets=fixture();z=s.symbols('z');a=2+z*z;b=z/4;am=z*z;mu=z
  functions=dict(a=a,ts=-b/a,p=z,mu=mu,E=(1+z*z)/2,v=a+b*b/a+a*am*am)
  for name,f in functions.items():
   for m in range(7):
    raw=abs(s.diff(f,z,m).subs(z,0));bound=jets['prefactors'][name]*jets['L']**m*factorial(m)**3
    self.assertLessEqual(raw,s.Rational(bound.numerator,bound.denominator))
  self.assertEqual(jets['N'],sum(jets['prefactors'].values()))
  self.assertNotEqual(s.diff(functions['v']-(a+b*b/a),z,4).subs(z,0),0)
 def test_finite_mixed_composition(self):
  def mul(a,b,n):return [sum(a[l]*b[k-l] for l in range(k+1)) for k in range(n+1)]
  for i in range(1,9):
   g=[Q(0)]+[Q(2)*3**n*factorial(n)**2 for n in range(1,i+1)]
   power=[Q(1)]+[Q(0)]*i;value=Q(0)
   for k in range(1,i+1):
    power=mul(power,g,i);value+=power[i]
   for j in (0,1,4):
    self.assertLessEqual(value*5**j,9**i*5**j*factorial(i)**2)
    if i>=4:self.assertGreater(value*5**j,9**i*5**j)
 def test_actual_zero_pressure_loop_mixed_jets(self):
  import sympy as s
  z,theta=s.symbols('z theta',real=True);d0=s.Rational(1,4)
  # At p=0, q_z(0,theta)=sin(theta), including mu=0.
  a=2+z*z;ts=z/4;mu=z;E=1+z;v=a*(1+ts**2+(d0*mu)**2/2)
  t=ts+d0*mu*s.sin(theta);w=a*(1+t*t)/(2*s.pi*v)
  aL=v/(1+t*t);bL=-v*t/(1+t*t)
  sources=[w,-(aL-a)/2,E*(bL+a*ts)/2]
  disk=q.parameter_disk(2,9,q.IDENTITY)
  # This local control's six parameter raw jets at eta=0 admit prefactor sum 100, base 4.
  for i in range(8):
   norm=sum(abs(s.diff(f,z,i).subs(z,0))/s.factorial(i) for f in (a,ts,s.Integer(0),mu,E,v))
   self.assertLessEqual(norm,100*4**i*factorial(i)**2)
  jets=dict(stage='six-real-scalar-Gevrey3',N=Q(100),L=Q(4));bound=q.bind(jets,disk)
  for f,cap in zip(sources,[bound['W'],bound['F'],bound['F']]):
   for i,j in [(0,0),(1,0),(2,0),(1,1),(0,3)]:
    for th in (0,s.pi/2,s.pi):
     value=abs(s.diff(f,z,i,theta,j).subs({z:0,theta:th}))/(s.factorial(i)*s.factorial(j))
     upper=cap*bound['B']**i*bound['C']**j*factorial(i)**2
     self.assertLessEqual(value,s.Rational(upper.numerator,upper.denominator))
 def test_joined_inverse_primitive_and_guards(self):
  state,recipe,jets=fixture();disk=q.parameter_disk(2,9,q.IDENTITY);r=q.bind(jets,disk)
  self.assertEqual(r['density_floor'],Q(1,32));
  import sympy as s
  self.assertLess(1/(4*s.pi),s.Rational(1,8)) # generic G=2 pointwise caps do not give G^-3
  self.assertGreater(1/(4*s.pi),s.Rational(1,32))
  self.assertEqual(r['C'],1/disk['r']);self.assertEqual(r['B'],jets['L']*(1+jets['N']/disk['r']))
  for kind in ('inverse','primitive'):
   for m in range(1,5):self.assertLessEqual(q.previous.bound(r[kind],m)['upper'],2**q.previous.bound(r[kind],m,True)['upper'])
   self.assertEqual(q.previous.bound(r[kind],100005,True)['factorial_power'],3)
  self.assertFalse(r['actual_final_swirl_repair_cone_or_scale_certified'])
  with self.assertRaises(ValueError):q.parameter_disk(2,9,'new pressure reset')
  with self.assertRaises(ValueError):q.scalar_jets(state,recipe,4,1)
  with self.assertRaises(ValueError):q.bind(dict(stage='wrong'),disk)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
