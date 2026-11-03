"""Full signed pressure moments at every fixed angular order; actual norms are inputs."""
from fractions import Fraction as Q
from math import factorial,prod
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q124_current',HERE.parent/'core_current_head_q123_fixed_bump_diagonal/q123.py');current=importlib.util.module_from_spec(s);s.loader.exec_module(current)
order=current.order;exact=current.exact;lg=current.lg

def validate(h,n,m,nu,Xcut,Xvel,Xv,U,F,previous,Jf,Ju,binary):
 order(n);order(m)
 if n<1:raise ValueError('Positive coefficient order required')
 h,nu,Xcut,Xvel,Xv=map(exact,(h,nu,Xcut,Xvel,Xv))
 if not 0<h<Q(1,100) or nu<0 or not 0<Xcut<=Xvel<=Xv:raise ValueError('Original h, viscosity and complete ordered supports required')
 U,F=tuple(U),tuple(F)
 if len(U)!=n or len(F)!=n:raise ValueError('All finalized lower orders 0..n-1 required')
 if n>=2 and previous is None:raise ValueError('Previous radial-velocity angular loss requires U_(n-2) through m+3')
 if n==1 and previous is not None:raise ValueError('No shifted predecessor at pressure order one')
 values=U+F+(Jf,Ju)+(() if previous is None else (previous,))
 if binary:
  for v in values:order(v)
 else:
  values=tuple(map(exact,values))
  if min(values)<0:raise ValueError('Nonnegative whole-support algebra norms required')
  U=tuple(map(exact,U));F=tuple(map(exact,F));Jf,Ju=map(exact,(Jf,Ju));previous=None if previous is None else exact(previous)
 return h,nu,Xcut,Xvel,Xv,U,F,previous,Jf,Ju

def _calculate(h,n,m,nu,xc,xvel,xv,U,F,previous,Jf,Ju,binary):
 if binary:
  const=lambda a:None if a==0 else lg(max(Q(1),Q(a)))
  def add(*v):
   v=[x for x in v if x is not None]
   return None if not v else max(v)+lg(Q(len(v)))
  def mul(*v):return None if any(x is None for x in v) else sum(v)
  R=lambda r:2*r+2
 else:
  const=lambda a:Q(a)
  add=lambda *v:sum(v,Q(0));mul=lambda *v:prod(v)
  R=lambda r:Q(4)**(r+1)
 D=Q(1,2)-h;A=Q(1,2)+h;k=n-1;lam=2*h*k
 V=[mul(R(m+1),const(4*(m+3+D+2*h*i)),U[i]) for i in range(n)]
 B=add(*(mul(U[i],V[k-i]) for i in range(n)))
 C=add(*(mul(V[i],V[k-i]) for i in range(n)))
 axial=[]
 for j in (0,1):
  if previous is None:axial.append(const(0));continue
  beta=2*h*(k-1)+j
  vp=mul(R(m+2),const(4*(m+4+D+2*h*(k-1))),previous)
  axial.append(mul(const(16*nu*(m+1+abs(beta-D))*(m+2+abs(beta))*xv**(j+1)/Q(j+1)),R(m),R(m+1),vp))
 I0=add(mul(R(m),const(xv),add(mul(const(lam+2*D*(m+1)),V[k]),mul(const(4*(m+1)+4*abs(lam-A)),B))),mul(const(xv/2),C),mul(const(2*nu),V[k]),axial[0])
 I1=add(mul(R(m),const(xv*xv/2),add(mul(const(1+lam+2*D*(m+1)),V[k]),mul(const(4*(m+1)+4*abs(lam+D)),B))),mul(const(xv*xv/4),C),axial[1])
 FF=add(*(mul(F[i],F[n-i]) for i in range(1,n)))
 UF=add(*(mul(U[i],F[n-i]) for i in range(1,n)))
 UU=add(*(mul(U[i],U[n-i]) for i in range(1,n)))
 mu=dict(m1=mul(const(xc),Ju),m2=mul(const(xc*xc),Jf),
  m3=add(mul(const(2*xc),F[0],Jf),mul(const(xvel),FF),mul(const(Q(1,2)),I0)),
  m4=add(mul(const(xc*xc),add(mul(U[0],Jf),mul(F[0],Ju))),mul(const(xvel*xvel),UF)),
  m5=add(mul(const(2*xc),U[0],Ju),mul(const(xc*xc),F[0],Jf),mul(const(xvel),UU),mul(const(xvel*xvel/2),FF),mul(const(Q(1,2)),I1)))
 return dict(I0=I0,I1=I1,axial=axial,axis=mul(const(2*nu),V[k]),discrepancies=mu)

def bounds(h,n,m,nu,Xcut,Xvel,Xv,U,F,previous_U,Jf,Ju,binary=False):
 args=validate(h,n,m,nu,Xcut,Xvel,Xv,U,F,previous_U,Jf,Ju,binary)
 h,nu,xc,xvel,xv,U,F,previous,Jf,Ju=args
 if not binary and m>40:raise ValueError('Use binary bounds at large derivative orders')
 r=_calculate(h,n,m,nu,xc,xvel,xv,U,F,previous,Jf,Ju,binary)
 r.update(angular_order=m,U_input_order=m+2,F_input_order=m,
          previous_U_input_order=None if n==1 else m+3,
          norm='factorial-angular-product-algebra',binary=binary,
          actual_profile_source_support_or_scale_inputs_certified=False,
          absolute_source_integral_certified=False)
 if binary:
  fact=m*lg(Q(m)) if m else 0
  r['raw_discrepancy_binary_upper']={k:max(0,0 if v is None else v)+fact for k,v in r['discrepancies'].items()}
  r['factorials_materialized']=False
 return r

def raw_discrepancies(r):
 if r['binary']:raise ValueError('Binary bounds are recipes, not materialized norms')
 return {k:factorial(r['angular_order'])*v for k,v in r['discrepancies'].items()}

def record():
 return dict(status='all-order-full-pressure-repair-discrepancy-bounds',
  reciprocal_angular_algebra='R_m=4^(m+1)',
  shifted_viscosity='16*R_m*R_(m+1)*(m+1+abs(beta-D))*(m+2+abs(beta))*|vprev|_(m+2)',
  all_five_full_pressure_discrepancies=True,ordinary_to_raw='D_i <= m!*mu_i',
  actual_profile_support_scale_cone_pulse_residual_or_blowup_certified=False)
