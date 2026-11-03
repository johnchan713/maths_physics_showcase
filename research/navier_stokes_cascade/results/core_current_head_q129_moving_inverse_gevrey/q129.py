"""Finite real angular Gevrey bounds for the same moving phase inverse."""
from fractions import Fraction as Q
from math import factorial
import importlib.util
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q129_previous',HERE.parent/'core_current_head_q128_real_cutoff_signed_root/q128.py')
previous=importlib.util.module_from_spec(s);s.loader.exec_module(previous)
exact=previous.exact;order=previous.order;lg=previous.lg
IDENTITY='original-fixed-origin-phase-and-density-weighted-source-mean'
def inverse(W,B,C,floor,identity):
 W,B,C,floor=map(exact,(W,B,C,floor))
 if identity!=IDENTITY or min(W,B,C)<1 or not 0<floor<=W:raise ValueError('Original phase identity, mixed density caps >=1 and positive real density floor required')
 A=8*W;k=A*C/floor;K=max(Q(1),k)
 t=min(Q(1,8),1/(64*K));eps=min(Q(1,8),t/(16*K))
 return dict(stage='real-moving-inverse',k=k,t=t,eps=eps,C=C,B=B,
  radius=eps/B,prefactor=t/C,base=B/eps,
  slow_analytic_extension_claimed=False,actual_mixed_density_inputs_certified=False)
def compose(inv,F,Bf,Cf,weighted_mean_identity):
 F,Bf,Cf=map(exact,(F,Bf,Cf))
 if inv.get('stage')!='real-moving-inverse' or min(F,Bf,Cf)<1 or weighted_mean_identity!=IDENTITY:raise ValueError('Same inverse, mixed source caps >=1 and original weighted mean identity required')
 R=inv['radius'];t=inv['t'];C=inv['C']
 r=min(R/2,R*C/(8*Cf*t),1/(4*Bf))
 return dict(stage='real-zero-mean-primitive',radius=r,prefactor=4*F,base=1/r,
  source_prefactor=2*F,zero_mean_requires_density_weight=True,
  actual_mixed_source_and_mean_inputs_certified=False,slow_analytic_extension_claimed=False)
def bound(recipe,m,binary=False):
 order(m)
 if recipe.get('stage') not in ('real-moving-inverse','real-zero-mean-primitive'):raise ValueError('Bound recipe required')
 if recipe['stage']=='real-moving-inverse' and m==0:raise ValueError('The inverse lift has no periodic value bound')
 pref,base=recipe['prefactor'],recipe['base']
 if binary:v=lg(max(Q(1),pref))+m*lg(max(Q(1),base))+(3*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary budget at high order')
  v=pref*base**m*factorial(m)**3
 return dict(order=m,upper=v,binary=binary,factorial_power=3,actual_profile_or_scale_certified=False)
def axial_join(Uraw,primitive,N,m):
 order(m);Uraw,N=exact(Uraw),exact(N)
 if primitive.get('stage')!='real-zero-mean-primitive' or Uraw<0 or N<1 or N.denominator!=1:raise ValueError('Raw incoming U norm and same fixed positive integer frequency required')
 return dict(raw_upper=Uraw+bound(primitive,m)['upper']/N,order=m,
  actual_final_repair_or_cone_certified=False)
def record():
 return dict(status='real-moving-phase-inverse-and-zero-mean-primitive-Gevrey3-bound',
  inverse_positive_raw='(t/C)*(B/eps)^m*(m!)^3, m>=1',
  primitive_raw='4*F*r^(-m)*(m!)^3, m>=0',
  pressure_input='Q128 complete incoming Pi and Pi_eta retained in the signed-root inputs',
  origin='phi(eta,theta)=integral_0^theta w(eta,s) ds; derivatives at fixed zeta',
  mean='integral_0^(2*pi) source(eta,theta)*w(eta,theta) dtheta = 0',
  actual_mixed_inputs_final_repair_scale_pulse_residual_or_blowup_certified=False,
  slow_analytic_extension_claimed=False)
