"""Closed mixed modulation at the original Ny carrier; actual primitive jets remain inputs."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q138_previous',HERE.parent/'core_current_head_q137_physical_repair_mixed/q137.py')
previous=importlib.util.module_from_spec(s);s.loader.exec_module(previous)
exact=previous.exact;order=previous.order;lg=previous.lg
IDENTITY='same-fixed-origin-weighted-zero-mean-primitives-at-original-Ny-and-full-pressure-state'
def pullback(M,A,B,C,V0,N,identity):
 M,A,B,C,V0,N=map(exact,(M,A,B,C,V0,N))
 if identity!=IDENTITY or min(M,A,B,C,N)<1 or N.denominator!=1 or not 0<=V0<=N:raise ValueError('Actual mixed primitive caps, real value gate and original integer Ny frequency required')
 T=M/N;R=A+C*N;L=8*(1+T)
 return dict(stage='same-Ny-mixed-primitive-pullback',primitive_prefactor=T,radial_base=R,angular_base=B,
  exp_prefactor=Q(9),difference_prefactor=9*max(V0,M)/N,exp_radial_base=L*R,exp_angular_base=L*B,
  N=N,M=M,A=A,B=B,C=C,V0=V0,
  actual_mixed_fixed_zeta_inputs_mean_and_value_gate_certified=False)
def modulation(p,e,invE,u,Ay,Be,identity):
 e,invE,u,Ay,Be=map(exact,(e,invE,u,Ay,Be))
 if identity!=IDENTITY or p.get('stage')!='same-Ny-mixed-primitive-pullback' or min(e,invE)<=0 or u<0 or min(Ay,Be)<1:raise ValueError('Same mixed pullback and actual PREMOD E/inverse/U mixed caps required')
 D=p['difference_prefactor'];T=p['primitive_prefactor'];R=Ay+p['exp_radial_base'];K=Be+p['exp_angular_base']
 return dict(stage='same-pre-repair-log-mixed-modulation-Gevrey3',radial_base=R,angular_base=K,
  prefactors=dict(E=9*e,invE=9*invE,U=u+T,deltaE=e*D,deltaInvE=invE*D,deltaU=T),
  N=p['N'],actual_premod_primitives_pressure_and_frequency_inputs_certified=False)
def physical(mod,xmin,identity):
 xmin=exact(xmin)
 if identity!=IDENTITY or mod.get('stage')!='same-pre-repair-log-mixed-modulation-Gevrey3' or xmin<=0:raise ValueError('Same modulated cover and positive lower radius required')
 pref={k:mod['prefactors'][k] for k in ('E','U','deltaE','deltaU')};pref['F']=pref['E']*max(Q(1),1/xmin);pref['deltaF']=pref['deltaE']*max(Q(1),1/xmin)
 R=(1+mod['radial_base'])/xmin
 return dict(stage='same-pre-repair-ordinary-X-mixed-modulation-Gevrey3',prefactors=pref,
  radial_bases={k:R+(1/xmin if k in ('F','deltaF') else 0) for k in pref},angular_base=mod['angular_base'],xmin=xmin,
  actual_whole_compact_mixed_inputs_supports_and_low_order_gates_certified=False)
def join_repair(mod,repair,identity):
 if identity!=IDENTITY or mod.get('stage')!='same-pre-repair-ordinary-X-mixed-modulation-Gevrey3' or repair.get('stage')!='same-leading-physical-repair-mixed-Gevrey3':raise ValueError('Same ordinary-X modulation and original full five-row repair required')
 return dict(stage='conditional-final-compact-raw-mixed-Gevrey3',prefactors={k:mod['prefactors'][k]+repair['prefactors'][k] for k in ('F','U')},
  radial_bases={k:max(mod['radial_bases'][k],repair['radial_bases'][k]) for k in ('F','U')},angular_base=max(mod['angular_base'],repair['angular_base']),
  actual_same_profile_cover_complete_pressure_target_or_cone_certified=False)
def bound(recipe,kind,r,k,binary=False):
 order(r);order(k)
 if recipe.get('stage') not in ('same-pre-repair-ordinary-X-mixed-modulation-Gevrey3','conditional-final-compact-raw-mixed-Gevrey3') or kind not in recipe['prefactors']:raise ValueError('Ordinary-X mixed component recipe required')
 P,R,K=recipe['prefactors'][kind],recipe['radial_bases'][kind],recipe['angular_base']
 if binary:v=lg(max(Q(1),P))+r*lg(max(Q(1),R))+k*lg(max(Q(1),K))+(3*r*lg(Q(r)) if r else 0)+(3*k*lg(Q(k)) if k else 0)
 else:
  if max(r,k)>12:raise ValueError('Use binary mixed budget at high order')
  v=P*R**r*K**k*(factorial(r)*factorial(k))**3
 return dict(component=kind,X_order=r,eta_order=k,upper=v,binary=binary,actual_complete_profile_or_scale_certified=False)
def record():
 return dict(status='mixed-Ny-pullback-exponential-and-physical-modulation-envelope',
  input='fixed-zeta primitive raw (i,k,j) <= M A^i B^k(i!k!)^3 j! C^j; actual mixed data, not Q129 angular data',
  pullback='P(y,eta,Ny)/N raw <= (M/N)(A+C*N)^r B^k(r!k!)^3',
  exponential='value gate |P|<=V0<=N; exp pref9, difference pref9max(V0,M)/N; bases 8(1+M/N)*(A+C*N,B)',
  ordinary_X='raw log pref P R^r K^k(r!k!)^3 becomes P[(1+R)/xmin]^r K^k(r!k!)^3',
  repair='add Q137 ordinary mixed prefactors and take max bases, retaining full five-row pressure target',
  remaining='bind actual mixed slow primitive and PREMOD jets, then join all regions/axis primitive and verify original gates/scales',
  radial_O_inverse_N_claimed=False,actual_complete_profile_cone_pulse_covariance_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
