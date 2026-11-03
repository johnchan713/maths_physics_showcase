import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q128 as q
import audit
FIELDS=dict(E=4,U=8,EX=2,UX=2,invE=3)
TRACES=dict(M=3,I=2,J=3,S=3,Pi=7)
def fixture():return q.state_caps(Q(1,200),1,4,Q(1,8),FIELDS,TRACES)
class Tests(unittest.TestCase):
 def test_full_forward_memory_caps_and_pressure_derivative(self):
  r=fixture();zero=dict(TRACES,Pi=0)
  z=q.state_caps(Q(1,200),1,4,Q(1,8),FIELDS,zero)
  A=Q(1,2)+Q(1,200)
  self.assertEqual(r['p']-z['p'],3*4*(8*A+10/Q(1,8))*7/(1-8*Q(1,200)))
  self.assertEqual(set(r['primitives']),q.TRACES)
  self.assertGreater(r['primitives']['Pi'],TRACES['Pi'])
  self.assertTrue(r['full_incoming_pressure_required'])
  import sympy as s
  x,e=s.symbols('x e',positive=True);eta=s.symbols('eta');h=s.Rational(1,200);L=1-2*h*eta**2;d=1-eta**2
  c=eta
  delta=x*(4*(s.Rational(1,2)+h)*eta*c-d*s.diff(c,eta))/(L*e)
  self.assertEqual(delta.subs(eta,0),-x/e)

 def test_independent_original_state_jets_on_polynomial_memory(self):
  import sympy as s
  x,z,lx=s.symbols('x z lx');h=s.Rational(1,200);D=s.Rational(1,2)-h;A=s.Rational(1,2)+h
  E=(1+z*z)/s.sqrt(x);U=x*z+z**3;L=1-2*h*z*z;d=1-z*z
  M=1+z+(x*x-1)*z/2+(x-1)*z**3
  S=1+z*z+(x**3-1)*z*z/3+(x*x-1)*z**4+(x-1)*z**6-(1+z*z)**2*lx/2
  Pi=3+z+z**4+(1+z*z)**2*(1-1/x)/2
  w=x-2*D*z*M-d*s.diff(M,z)
  a=1-2*x*s.diff(E,x)/E;b=2*x*s.diff(U,x)/E
  p=(-w*U+D*(M-z*s.diff(M,z))+4*h*z*S-d*s.diff(S,z)+x*(4*A*z*Pi-d*s.diff(Pi,z)))/(L*E)
  r=fixture()
  for name,f in [('a',a),('b',b),('p',p)]:
   for j in range(4):
    for xx in (1,4):
     for eta in (-1,0,1):
      v=s.diff(f,z,j).subs({x:xx,z:eta})/s.factorial(j)
      poly=s.Poly(s.expand(v),lx)
      cap=sum(abs(c)*3**ij[0] for ij,c in poly.terms()) # |log X|<=3 on [1,4]
      bound=r[name]*r['rho']**(-j)
      self.assertLessEqual(cap,s.Rational(bound.numerator,bound.denominator))

 def test_closed_Gevrey_composition_retains_factorials(self):
  def mul(a,b,m):return [sum(a[i]*b[k-i] for i in range(k+1)) for k in range(m+1)]
  for m in range(1,9):
   g=[0]+[Q(2)*3**j for j in range(1,m+1)];power=[Q(1)]+[Q(0)]*m;value=Q(0)
   for j in range(1,m+1):
    power=mul(power,g,m);value+=729*8**j*factorial(j)**2*power[m]
   self.assertLessEqual(value,729*(3*(1+8*2))**m*factorial(m)**2)
   if m>=4:self.assertGreater(value,729*(3*(1+8*2))**m)
  for m in range(1,8):
   g=[0]+[Q(2)*5**j*factorial(j)**2 for j in range(1,m+1)];power=[Q(1)]+[Q(0)]*m;value=Q(0)
   for j in range(1,m+1):
    power=mul(power,g,m);value+=Q(1,10)*3**j*power[m]
   self.assertLessEqual(value,Q(1,10)*(5*(1+3*2))**m*factorial(m)**2)

 def test_actual_plateau_radicand_derivatives(self):
  import sympy as s
  z=s.symbols('z');a=2+z*z;b=z/4;v=a+b*b/a
  t=s.sqrt((s.Rational(5,2)-v)/a)
  state=dict(stage='PREMOD-state',rho=Q(1,4),a=Q(4),b=Q(1,2),p=Q(4))
  recipe=q.parameters(state,1,1,2,9,Q(1,32),q.IDENTITY)
  # v(0)=2 lies strictly in the real K=1 plateau; no complex K is used.
  for m in range(7):
   raw=abs(s.diff(t,z,m).subs(z,0));cap=q.derivative_bound(recipe,m,'amplitude')['upper']
   self.assertLessEqual(raw,s.Rational(cap.numerator,cap.denominator))

 def test_closed_parameters_bits_and_zero_pressure_scope(self):
  recipe=q.parameters(fixture(),Q(1,2),Q(1,4),2,9,Q(1,32),q.IDENTITY)
  old=q.signed.radii(2,9,Q(1,32));self.assertEqual(recipe['scalar_inverse_radius'],old['sigma'])
  self.assertEqual(recipe['mu_prefactor'],old['a']);self.assertFalse(recipe['slow_analytic_extension_claimed'])
  for kind in ('mu','amplitude'):
   for m in range(6):
    exact=q.derivative_bound(recipe,m,kind)['upper'];bits=q.derivative_bound(recipe,m,kind,binary=True)['upper']
    self.assertLessEqual(exact,2**bits)
   high=q.derivative_bound(recipe,100005,kind,binary=True)
   self.assertEqual(high['factorial_power'],3);self.assertFalse(high['slow_analytic_extension_claimed'])
  self.assertEqual(q.derivative_bound(recipe,0,'mu')['upper'],2)
  state=dict(fixture(),p=Q(0))
  self.assertGreater(q.parameters(state,Q(1,2),Q(1,4),2,9,Q(1,32),q.IDENTITY)['mu_base'],0)

 def test_identity_trace_and_materialization_guards(self):
  with self.assertRaises(ValueError):q.state_caps(Q(1,200),1,4,Q(1,8),FIELDS,dict(M=1,I=1,J=1,S=1))
  with self.assertRaises(ValueError):q.parameters(fixture(),0,1,2,9,Q(1,32),q.IDENTITY)
  with self.assertRaises(ValueError):q.parameters(fixture(),1,1,2,9,Q(1,32),'analytic slow cutoff')
  recipe=q.parameters(fixture(),1,1,2,9,Q(1,32),q.IDENTITY)
  with self.assertRaises(ValueError):q.derivative_bound(recipe,100005,'mu')
  with self.assertRaises(ValueError):q.derivative_bound(recipe,True,'amplitude')

 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
