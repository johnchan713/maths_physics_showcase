import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q131 as q
import audit
def primitive(M=4,L=2):return dict(stage='real-zero-mean-primitive',prefactor=Q(M),base=Q(L))
def fixture():return q.modulation(4,3,4,Q(1,4),primitive(),primitive(3),1,8,q.IDENTITY)
class Tests(unittest.TestCase):
 def test_exponential_finite_series_and_value_gate(self):
  def mul(a,b,m):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(m+1)]
  for m in range(1,9):
   t=Q(1,2);g=[0]+[t*2**k*factorial(k)**2 for k in range(1,m+1)];power=[Q(1)]+[Q(0)]*m;value=0
   for j in range(1,m+1):power=mul(power,g,m);value+=power[m]/factorial(j)
   self.assertLessEqual(value,3*t*(2*(1+t))**m*factorial(m)**2)
  # The value gate uses V0, even when the derivative prefactor is enormous.
  r=q.modulation(1,1,1,Q(1,2),primitive(2**500),primitive(),1,8,q.IDENTITY)
  self.assertEqual(r['exponential_value_factor'],3)
  self.assertEqual(r['real_swirl_floor_factor'],Q(1,3))
  with self.assertRaises(ValueError):q.modulation(1,1,1,Q(1,2),primitive(),primitive(),9,8,q.IDENTITY)
 def test_actual_signed_exponential_field_jets(self):
  import sympy as s
  z=s.symbols('z');E=(2+z*z);U=1+z;f=s.sin(z);g=(z+z*z)/2;N=8
  r=fixture();EN=E*s.exp(f/N);UN=U+g/N
  cases=dict(deltaE=EN-E,deltaInvE=1/EN-1/E,deltaU=UN-U,E=EN,invE=1/EN,U=UN)
  for kind,expr in cases.items():
   for m in range(7):
    value=abs(s.diff(expr,z,m).subs(z,0));cap=q.bound(r,kind,m)['upper']
    self.assertLessEqual(value,s.Rational(cap.numerator,cap.denominator))
  self.assertEqual(q.bound(r,'deltaE',0)['upper'],Q(6)) # e*3*max(1,4)/8
 def test_independent_five_density_expansions(self):
  import sympy as s
  z,x=s.symbols('z x',positive=True);E=2+z*z;U=1+z;dE=z/8;dU=(1+z)/16
  en=E+dE;un=U+dU
  actual=[un-U,s.sqrt(2*x)*(en-E),s.sqrt(2*x)*(en*un-E*U),un**2-en**2/2-U**2+E**2/2,(en**2-E**2)/(2*x)]
  expected=[dU,s.sqrt(2*x)*dE,s.sqrt(2*x)*(E*dU+U*dE+dE*dU),2*U*dU+dU*dU-E*dE-dE*dE/2,(E*dE+dE*dE/2)/x]
  for f,g in zip(actual,expected):self.assertEqual(s.expand(f-g),0)
  r=fixture();mem=q.memory(r,1,3,{k:(Q(1),Q(2)) for k in q.ROWS})
  # Compare actual normalized derivatives at fixed real X to the derived density majorants.
  for k,f in zip(('M','I','J','S','Pi'),actual):
   for m in range(5):
    value=abs(s.diff(f,z,m).subs({z:0,x:2}));cap=mem['density_prefactors'][k]*(2*r['base'])**m*factorial(m)**3
    self.assertLessEqual(value,s.Rational(cap.numerator,cap.denominator))
 def test_incoming_pressure_and_shifted_jets(self):
  r=fixture();traces={k:(Q(0),Q(1)) for k in q.ROWS};a=q.memory(r,1,3,traces)
  traces['Pi']=(Q(7),Q(3));b=q.memory(r,1,3,traces)
  self.assertEqual(b['prefactors']['Pi']-a['prefactors']['Pi'],7)
  self.assertEqual(set(b['prefactors']),q.ROWS)
  for m in range(13):
   self.assertLessEqual(b['prefactors']['Pi']*b['base']**(m+1)*factorial(m+1)**3,b['shifted_derivative_prefactors']['Pi']*b['shifted_derivative_base']**m*factorial(m)**3)
  import sympy as s
  z,x,e=s.symbols('z x e');h=s.Rational(1,200);c=z
  change=x*(4*(s.Rational(1,2)+h)*z*c-(1-z*z)*s.diff(c,z))/((1-2*h*z*z)*e)
  self.assertEqual(change.subs(z,0),-x/e)
  # Integrate the actual exponential pressure-density difference from x0=1.
  density=2*(s.exp(z/4)-1)/x
  P=z+z*z+2*s.log(x)*(s.exp(z/4)-1)
  self.assertEqual(s.simplify(s.diff(P,x)-density),0)
  traces={k:(Q(0),Q(1)) for k in q.ROWS};traces['Pi']=(Q(3),Q(2));mem=q.memory(r,1,3,traces)
  for m in range(7):
   value=abs(s.diff(P,z,m).subs({z:0,x:3}).subs(s.log(3),2)) # |log 3|<2
   cap=mem['prefactors']['Pi']*mem['base']**m*factorial(m)**3
   self.assertLessEqual(value,s.Rational(cap.numerator,cap.denominator))

  with self.assertRaises(ValueError):q.memory(r,1,3,dict(M=(0,1),I=(0,1),J=(0,1),S=(0,1)))
 def test_radial_countercontrol_and_binary_budgets(self):
  import sympy as s
  y=s.symbols('y');N=16;f=s.sin(2*s.pi*N*y)/N
  self.assertEqual(s.simplify(f.subs(y,y+s.Rational(1,N))-f),0)
  self.assertEqual(abs(s.diff(f,y,3).subs(y,0)),(2*s.pi)**3*N**2)
  self.assertGreater(N**2,Q(1,N))
  r=fixture()
  for kind in r['prefactors']:
   for m in range(6):self.assertLessEqual(q.bound(r,kind,m)['upper'],2**q.bound(r,kind,m,True)['upper'])
   self.assertEqual(q.bound(r,kind,100005,True)['factorial_power'],3)
  with self.assertRaises(ValueError):q.bound(r,'deltaE',100005)
  with self.assertRaises(ValueError):q.modulation(1,1,1,Q(1,2),dict(stage='real-moving-inverse'),primitive(),1,8,q.IDENTITY)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
