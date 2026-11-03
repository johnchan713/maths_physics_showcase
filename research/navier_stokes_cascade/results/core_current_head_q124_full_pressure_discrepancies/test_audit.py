import unittest
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import q124 as q
import audit
class Tests(unittest.TestCase):
 def test_independent_full_order_two_pressure_and_five_moments(self):
  import sympy as s
  x,z=s.symbols('x z');h=s.Rational(1,200);D=s.Rational(1,2)-h;L=1-2*h*z*z;d=1-z*z
  primitive=x**3*(1-x)**5;rad=s.diff(primitive,x)
  U=[(i+1)*rad*(z+z**3) for i in range(2)];M=[(i+1)*primitive*(z+z**3) for i in range(2)]
  v=[s.cancel((2*z*U[i]-2*(D+2*h*i)*z*M[i]/x-d*s.diff(M[i],z)/x)/L) for i in range(2)]
  def Z(a,f):return (2*a*z*f-2*z*x*s.diff(f,x)+d*s.diff(f,z))/L
  k=1
  w=((1-2*h*k)*v[k]+x*s.diff(v[k],x)+D*z*s.diff(v[k],z))/L-2*(2*s.diff(v[k],x)+x*s.diff(v[k],x,2))
  w+=sum(v[i]*(v[k-i]/2+x*s.diff(v[k-i],x))+U[i]*Z(2*h*(k-i)-1,v[k-i]) for i in range(k+1))
  w-=Z(-1-D,Z(-1,v[0]))
  def integ(f):
   p=s.Poly(s.cancel(f),x)
   return s.factor(sum(c/s.Rational(j[0]+1) for j,c in p.terms()))
  I=[integ(w),integ(x*w)];F=[(i+1)*(1-x)**6*(1+z) for i in range(2)]
  Fn=3*(1-x)**6*(1+z*z);Un=3*rad*(1+z*z);px=2*F[0]*Fn+F[1]*F[1]-w/2
  densities=[Un,2*x*Fn,px,2*x*(U[0]*Fn+Un*F[0]+U[1]*F[1]),2*U[0]*Un+U[1]**2-x*px]
  moments=[integ(f) for f in densities]
  def norm(f,r):return Q(sum(sum(abs(c) for c in s.Poly(s.diff(f,z,j),x,z).coeffs())/s.factorial(j) for j in range(r+1)))
  for m in (0,1,3):
   r=q.bounds(Q(1,200),2,m,1,1,1,1,[norm(f,m+2) for f in U],[norm(f,m) for f in F],norm(U[0],m+3),norm(Fn,m),norm(Un,m))
   for expr,cap in list(zip(I,[r['I0'],r['I1']]))+list(zip(moments,r['discrepancies'].values())):
    for eta in (-1,0,1):
     value=sum(abs(s.diff(expr,z,j).subs(z,eta))/s.factorial(j) for j in range(m+1))
     self.assertLessEqual(value,s.Rational(cap.numerator,cap.denominator))
   self.assertGreater(r['axial'][1],0)
 def test_general_second_discrepancies_recover_inherited_formula(self):
  import importlib.util
  p=q.HERE.parent/'higher_pressure_source_moments/moments.py';sp=importlib.util.spec_from_file_location('q124_old_moments',p);old=importlib.util.module_from_spec(sp);sp.loader.exec_module(old)
  r=q.bounds(Q(1,200),2,1,1,2,3,7,[11,17],[7,13],19,5,5)
  self.assertEqual(tuple(r['discrepancies'].values()),old.second_repair_bounds(2,3,5,7,11,13,17,r['I0'],r['I1']))
 def test_binary_majorizes_exact_and_raw_conversion(self):
  for n,m in ((1,0),(2,1),(4,4)):
   ub=list(range(1,n+1));fb=list(range(2,n+2));pb=None if n==1 else 3
   exact=q.bounds(Q(1,200),n,m,1,2,3,7,[2**b for b in ub],[2**b for b in fb],None if pb is None else 2**pb,8,16)
   r=q.bounds(Q(1,200),n,m,1,2,3,7,ub,fb,pb,3,4,binary=True)
   for k in ('I0','I1','axis'):self.assertLessEqual(exact[k],2**r[k])
   raw=q.raw_discrepancies(exact)
   for k,v in raw.items():self.assertLessEqual(v,2**r['raw_discrepancy_binary_upper'][k])
   self.assertEqual(raw['m1'],factorial(m)*exact['discrepancies']['m1'])
 def test_shifted_and_axis_terms_not_dropped(self):
  r=q.bounds(Q(1,200),2,5,1,1,2,3,[2,3],[5,7],11,13,17)
  z=q.bounds(Q(1,200),2,5,1,1,2,3,[2,3],[5,7],0,13,17)
  self.assertEqual(r['I1']-z['I1'],r['axial'][1]);self.assertGreater(r['axis'],0)
  self.assertEqual(r['previous_U_input_order'],8)
  first=q.bounds(Q(1,200),1,5,1,1,2,3,[2],[5],None,13,17)
  self.assertEqual(first['axial'],[0,0]);self.assertIsNone(first['previous_U_input_order'])
 def test_large_order_and_scope_guards(self):
  r=q.bounds(Q(1,200),2,100005,1,1,2,3,[2,3],[5,7],11,13,17,binary=True)
  self.assertEqual(r['previous_U_input_order'],100008);self.assertFalse(r['factorials_materialized'])
  self.assertEqual(set(r['raw_discrepancy_binary_upper']),{'m1','m2','m3','m4','m5'})
  self.assertFalse(r['absolute_source_integral_certified']);self.assertFalse(r['actual_profile_source_support_or_scale_inputs_certified'])
  with self.assertRaises(ValueError):q.raw_discrepancies(r)
  with self.assertRaises(ValueError):q.bounds(Q(1,200),2,1,1,1,2,3,[2,3],[5,7],None,13,17)
  with self.assertRaises(ValueError):q.bounds(Q(1,200),2,1,1,3,2,3,[2,3],[5,7],11,13,17)
  with self.assertRaises(ValueError):q.bounds(Q(1,200),2,100005,1,1,2,3,[2,3],[5,7],11,13,17)
 def test_closed_reciprocal_and_double_angular_operator(self):
  import sympy as s
  z=s.symbols('z');h=s.Rational(1,200);L=1-2*h*z*z;d=1-z*z;D=s.Rational(1,2)-h;f=1+z+z**6
  for m in range(7):
   reciprocal=sum(abs(s.diff(1/L,z,j).subs(z,1))/s.factorial(j) for j in range(m+1))
   self.assertLessEqual(reciprocal,4**(m+1))
  for m,b in ((0,s.Rational(0)),(2,s.Rational(7)),(4,s.Rational(-2))):
   def T(a,g):return (d*s.diff(g,z)+2*a*z*g)/L
   expr=T(b-D,T(b,f))
   fn=sum(sum(abs(c) for c in s.Poly(s.diff(f,z,j),z).all_coeffs())/s.factorial(j) for j in range(m+3))
   cap=16*4**(m+1)*4**(m+2)*(m+1+abs(b-D))*(m+2+abs(b))*fn
   for eta in (-1,0,1):
    self.assertLessEqual(sum(abs(s.diff(expr,z,j).subs(z,eta))/s.factorial(j) for j in range(m+1)),cap)
 def test_record(self):
  import subprocess,sys,tempfile
  with tempfile.TemporaryDirectory() as tmp:
   out=Path(tmp)/'replay.json';frozen=audit.HERE/'evidence.json'
   subprocess.run([sys.executable,str(audit.HERE/'audit.py'),'--verify-record',str(frozen),'--output',str(out)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
   self.assertEqual(out.read_bytes(),frozen.read_bytes())
if __name__=='__main__':unittest.main()
