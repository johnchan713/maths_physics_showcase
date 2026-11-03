import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q121 as q
import audit
class Tests(unittest.TestCase):
 def test_finite_positive_majorants_all_components(self):
  r=q.radial
  for m in range(9):
   h=Q(1,200);n=7;Xhi=Q(3);lo=Q(1,2);hi=Q(3);H=Q(5)
   p=q.recipe(h,n,m,Xhi,lo,hi,H,H);tau=p['tau']
   keys=[(a,b) for a in range(m+1) for b in range(m+1-a)]
   U={ab:H/Q(factorial(ab[0])*factorial(ab[1])) for ab in keys}
   M={ab:(Xhi if ab[0]==0 else 1)*U[ab] for ab in keys}
   # Raw derivative bounds also apply to M_eta, using U through m+1.
   eta={(0,0):Q(1),(0,1):Q(1)};X={(0,0):Xhi,(1,0):Q(1)};d={(0,0):Q(1),(0,1):Q(2),(0,2):Q(1)}
   L={(0,0):1/(1-2*h)}
   for k in range(1,m+1):L[0,k]=(4*h*L[0,k-1]+(2*h*L[0,k-2] if k>=2 else 0))/(1-2*h)
   mul=lambda a,b:r.mul(a,b,m)
   c=Q(1,2)-h+2*n*h
   V=mul(L,r.add(r.scale(2,mul(eta,mul(X,U))),r.scale(2*c,mul(eta,M)),mul(d,M)))
   I=r.radius_series(m,lo,hi,True)
   polys=[U,mul(r.radius_series(m,lo,hi),U),mul(I,V),r.scale(2,mul(I,mul(L,mul(eta,M))))]
   for poly in polys:
    self.assertLessEqual(sum(v*tau**sum(ab) for ab,v in poly.items()),p['prefactor'])
    for ab,v in poly.items():self.assertLessEqual(v*factorial(ab[0])*factorial(ab[1]),q.materialize(p))

 def test_independent_exact_profile_jets(self):
  import sympy as s
  x,z=s.symbols('x z');h=s.Rational(1,200);F=x*x+z**3;U=x*z*z+z**3+1
  M=s.integrate(U,(x,0,x));L=1-2*h*z*z;R=s.sqrt(2*x)
  for m in range(5):
   def rawmax(f,k):
    return max(abs(s.diff(f,x,a,z,b).subs({x:8,z:1})) for a in range(k+1) for b in range(k+1-a))
   for n in (0,2):
    c=s.Rational(1,2)-h+2*n*h
    fields=[R*F,U,(2*z*x*U-2*c*z*M-(1-z*z)*s.diff(M,z))/(L*R),2*z*M/(L*R)]
    recipe=q.recipe(Q(1,200),n,m,8,1,4,Q(rawmax(F,m)),Q(rawmax(U,m+1)))
    bound=q.materialize(recipe)
    for field in fields:
     for a in range(m+1):
      for b in range(m+1-a):
       jet=s.diff(field,x,a,z,b)
       for eta in (-1,0,1):
        value=jet.subs({x:2,z:eta})
        self.assertLessEqual(abs(value),s.Rational(bound.numerator,bound.denominator))

 def test_extra_angular_order_and_axis_memory(self):
  r=q.recipe(Q(1,200),3,9,8,1,4,1,1)
  self.assertEqual((r['F_order'],r['U_order']),(9,10));self.assertTrue(r['whole_axis_interval_required'])
  import sympy as s
  z=s.symbols('z');mu=s.Rational(7);R=s.Rational(2);h=s.Rational(1,200);c=s.Rational(1,2)-h+6*h
  for k in range(1,6):
   M=mu*z**k/s.factorial(k)
   W=-(2*c*z*M+(1-z*z)*s.diff(M,z))/((1-2*h*z*z)*R)
   self.assertEqual(s.diff(W,z,k-1).subs(z,0),-mu/R)
  self.assertFalse(r['pressure_support_certified'])

 def test_binary_budget_and_current_diagonal(self):
  for m in range(12):
   r=q.recipe(Q(1,200),m,m,8,1,4,3,5)
   self.assertLessEqual(q.materialize(r),2**q.tuple_binary_upper(r))
  r=q.recipe(Q(1,200),100010,100010,8,1,4,3,5)
  self.assertGreater(q.tuple_binary_upper(r),100010)
  with self.assertRaises(ValueError):q.materialize(r)
  joined=q.diagonal_from_profiles(Q(1,200),100010,100010,8,1,4,3,5,0)
  self.assertFalse(joined['row']['sufficient_new_row_gate']);self.assertFalse(joined['row']['factorials_materialized'])
  self.assertFalse(joined['actual_profile_and_scale_provenance_certified'])

 def test_hull_and_input_guards(self):
  for args in ((0,0,1,8,1,4,1,1),(Q(1,200),-1,1,8,1,4,1,1),(Q(1,200),0,True,8,1,4,1,1),(Q(1,200),0,1,8,1,3,1,1),(Q(1,200),0,1,8,1,4,-1,1)):
   with self.assertRaises(ValueError):q.recipe(*args)
 def test_zero_input(self):
  r=q.recipe(Q(1,200),0,5,8,1,4,0,0)
  self.assertEqual(q.materialize(r),0)
  self.assertGreaterEqual(q.tuple_binary_upper(r),0)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
