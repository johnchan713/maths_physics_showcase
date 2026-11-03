"""Uniform angular product-algebra cost and the same first-coefficient cut seed."""
from fractions import Fraction as Q
from math import factorial,prod
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q125_current',HERE.parent/'core_current_head_q124_full_pressure_discrepancies/q124.py');current=importlib.util.module_from_spec(s);s.loader.exec_module(current)
order=current.order;exact=current.exact;lg=current.lg
SEED_IDENTITY='same-first-core-collar-Cauchy-cover'
def _calculate(h,n,m,nu,xc,xvel,xv,U,F,previous,Jf,Ju,binary):
 reciprocal=1/(1-8*h)
 if binary:
  const=lambda a:None if a==0 else lg(max(Q(1),Q(a)))
  def add(*v):
   v=[x for x in v if x is not None]
   return None if not v else max(v)+lg(Q(len(v)))
  def mul(*v):return None if any(x is None for x in v) else sum(v)
  R=lambda r:const(reciprocal)
 else:
  const=lambda a:Q(a)
  add=lambda *v:sum(v,Q(0));mul=lambda *v:prod(v)
  R=lambda r:reciprocal
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
 h,nu,xc,xvel,xv,U,F,previous,Jf,Ju=current.validate(h,n,m,nu,Xcut,Xvel,Xv,U,F,previous_U,Jf,Ju,binary)
 r=_calculate(h,n,m,nu,xc,xvel,xv,U,F,previous,Jf,Ju,binary)
 r.update(angular_order=m,U_input_order=m+2,F_input_order=m,
  previous_U_input_order=None if n==1 else m+3,norm='factorial-angular-product-algebra',binary=binary,
  reciprocal_algebra_upper=1/(1-8*h),factorials_materialized=False,
  actual_profile_source_support_or_scale_inputs_certified=False,absolute_source_integral_certified=False)
 if binary:
  fact=m*lg(Q(m)) if m else 0
  r['raw_discrepancy_binary_upper']={k:max(0,0 if v is None else v)+fact for k,v in r['discrepancies'].items()}
 return r

def raw_discrepancies(r):
 if r['binary'] or r['angular_order']>40:raise ValueError('Use binary budgets for high raw derivatives')
 return current.raw_discrepancies(r)

def first_seed(m,delta,core,collar,identity,binary=False):
 order(m);delta=exact(delta)
 if delta<=0 or identity!=SEED_IDENTITY:raise ValueError('Same positive first-coefficient Cauchy radius and strip cover required')
 if binary:
  order(core);order(collar)
  value=lg(Q(m+1))+m*lg(max(Q(1),2/delta))+1+max(core,collar)
 else:
  core,collar=map(exact,(core,collar))
  if min(core,collar)<0 or m>40:raise ValueError('Nonnegative moderate first core/collar bounds; use bits at high order')
  value=(core+collar)*sum((2/delta)**j for j in range(m+1))
 return dict(coefficient=1,angular_order=m,upper=value,binary=binary,
             actual_original_strip_or_solution_bounds_certified=False)

def first_discrepancies(h,m,nu,Xcut,Xvel,Xv,U0,F0,delta,core,collar,identity,binary=False):
 seed=first_seed(m,delta,core,collar,identity,binary)
 r=bounds(h,1,m,nu,Xcut,Xvel,Xv,[U0],[F0],None,seed['upper'],seed['upper'],binary)
 return dict(seed=seed,moments=r,leading_U_input_order=m+2,leading_F_input_order=m,
             actual_leading_norms_and_scale_certified=False)

def record():
 return dict(status='uniform-angular-reciprocal-and-first-cut-seed',
  reciprocal_algebra='R=1/(1-8h) for every fixed m',
  pressure_operator_cost='polynomial of degree at most three in m with native input budgets held fixed',
  same_first_seed='(Mcore+Mcollar)*sum_(j=0)^m (2/delta)^j',
  leading_U_order='m+2 for the first incoming discrepancy, m+3 shifted predecessor at later orders',
  actual_global_profile_scale_cone_pulse_residual_or_blowup_certified=False)
