"""Real angular compact modulation bounds, before the later exact repair."""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q131_previous',HERE.parent/'core_current_head_q130_mixed_loop_binding/q130.py')
previous=importlib.util.module_from_spec(s);s.loader.exec_module(previous)
exact=previous.exact;order=previous.previous.order;lg=previous.previous.lg
IDENTITY='same-original-fixed-integer-frequency-and-real-primitive-value-cap'
ROWS={'M','I','J','S','Pi'}
def modulation(e,invE,u,rho,primitiveA,primitiveB,valuecap,N,identity):
 e,invE,u,rho,valuecap,N=map(exact,(e,invE,u,rho,valuecap,N))
 if identity!=IDENTITY or min(e,invE,rho)<=0 or u<0 or valuecap<0 or N<1 or N.denominator!=1 or valuecap>N or rho>1:raise ValueError('Original fixed integer frequency, real value cap <=N and actual PREMOD complex field caps required')
 for p in (primitiveA,primitiveB):
  if p.get('stage')!='real-zero-mean-primitive':raise ValueError('Same Q129/Q130 zero-mean primitive bounds required')
 MA,LA=map(exact,(primitiveA['prefactor'],primitiveA['base']));MB,LB=map(exact,(primitiveB['prefactor'],primitiveB['base']))
 if min(MA,MB)<0 or min(LA,LB)<1:raise ValueError('Nonnegative primitive prefactors and bases >=1 required')
 R=1/rho;t=MA/N;Lexp=LA*(1+t);B=max(R+Lexp,LB)
 D=3*max(valuecap,MA)/N;DE=e*D;DI=invE*D;DU=MB/N
 return dict(stage='angular-modulation-before-exact-repair',N=N,base=B,
  prefactors=dict(deltaE=DE,deltaInvE=DI,deltaU=DU,E=3*e*max(Q(1),t),invE=3*invE*max(Q(1),t),U=u+DU),
  original_field_caps=dict(E=e,invE=invE,U=u),
  real_value_caps=dict(E=3*e,invE=3*invE,U=u+valuecap/N),
  real_swirl_floor_factor=Q(1,3),exponential_value_factor=Q(3),
  actual_original_fields_primitives_value_gate_or_frequency_certified=False,
  actual_final_repair_cone_pulse_residual_or_blowup_certified=False,slow_analytic_extension_claimed=False)
def bound(recipe,kind,m,binary=False):
 order(m)
 if recipe.get('stage')!='angular-modulation-before-exact-repair' or kind not in recipe['prefactors']:raise ValueError('Same pre-repair modulation component required')
 P,B=recipe['prefactors'][kind],recipe['base']
 if binary:value=lg(max(Q(1),P))+m*lg(max(Q(1),B))+(3*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary factorial budget at high order')
  value=P*B**m*factorial(m)**3
 return dict(order=m,upper=value,factorial_power=3,actual_profile_or_scale_certified=False)
def memory(recipe,x0,x1,incoming):
 x0,x1=exact(x0),exact(x1)
 if recipe.get('stage')!='angular-modulation-before-exact-repair' or not 0<x0<=x1 or set(incoming)!=ROWS:raise ValueError('Same compact interval and all five actual incoming discrepancy caps required')
 traces={k:tuple(map(exact,v)) for k,v in incoming.items()}
 if any(len(v)!=2 or v[0]<0 or v[1]<1 for v in traces.values()):raise ValueError('Incoming raw Gevrey prefactor/base pairs required')
 e,u=recipe['original_field_caps']['E'],recipe['original_field_caps']['U'];DE,DU=recipe['prefactors']['deltaE'],recipe['prefactors']['deltaU']
 r=max(Q(1),2*x1)
 density=dict(M=DU,I=r*DE,J=r*(e*DU+u*DE+DE*DU),S=2*u*DU+DU*DU+e*DE+DE*DE/2,Pi=(e*DE+DE*DE/2)/x0)
 base=max(2*recipe['base'],*(v[1] for v in traces.values()))
 pref={k:traces[k][0]+(x1-x0)*density[k] for k in ROWS}
 return dict(stage='full-five-incoming-angular-modulation-memory',density_prefactors=density,prefactors=pref,base=base,
  shifted_derivative_prefactors={k:v*base for k,v in pref.items()},shifted_derivative_base=8*base,
  actual_incoming_memory_interval_caps_and_final_repair_certified=False,slow_analytic_extension_claimed=False)
def record():
 return dict(status='closed-angular-swirl-and-five-row-modulation-memory-bound',
  swirl='E_N=E*exp(Acal/N), inverse E_N=E^-1*exp(-Acal/N)',
  angular_difference='D=3*max(V0,MA)/N; delta E and inverse E use e*D and invE*D',
  value_gate='|Acal|,|Bcal|<=V0<=N implies exp(|Acal|/N)<3; N is the unchanged selected integer',
  memory='All five incoming differences retained, including Pi and its angular derivative',
  shifted_jets='D*L^(m+1)*((m+1)!)^3 <= D*L*(8*L)^m*(m!)^3',
  radial_warning='Pure angular fixed-y estimates only; no radial O(1/N) claim',
  actual_profile_frequency_gate_final_repair_cone_pulse_residual_or_blowup_certified=False,
  slow_analytic_extension_claimed=False)
