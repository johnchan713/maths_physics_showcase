"""Closed real Gevrey jets of the SAME leading quadratic five-row repair."""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q132_previous',HERE.parent/'core_current_head_q131_swirl_memory_gevrey/q131.py')
previous=importlib.util.module_from_spec(s);s.loader.exec_module(previous)
exact=previous.exact;order=previous.order;lg=previous.lg
ROWS=previous.ROWS
IDENTITY='same-leading-constant-BQ-five-row-repair-with-e_star-over-one-plus-eta_squared'
def target(memory,baseline,e_star,rho,lam,identity):
 e_star,rho,lam=map(exact,(e_star,rho,lam))
 if memory.get('stage')!='full-five-incoming-angular-modulation-memory' or identity!=IDENTITY or set(baseline)!=ROWS or min(e_star,rho,lam)<=0 or lam>Q(1,100):raise ValueError('Complete five-row memory, all baseline target defects and SAME leading patch required')
 pairs={k:tuple(map(exact,v)) for k,v in baseline.items()}
 if any(len(v)!=2 or v[0]<0 or v[1]<1 for v in pairs.values()):raise ValueError('Baseline raw Gevrey prefactor/base pairs required')
 D={k:memory['prefactors'][k]+pairs[k][0] for k in ROWS};L=max(memory['base'],*(v[1] for v in pairs.values()))
 geom=max(Q(1),rho**(-2));P1=5/e_star;P2=25/e_star**2
 A=dict(M=geom*P1*D['M'],I=geom*P1*D['I'],J=geom*P2*D['J'],S=geom*P2*D['S'],Pi=geom*P2*D['Pi'])
 transformed=(A['M'],(A['M']+A['J'])/lam,A['I'],A['S'],A['Pi'])
 return dict(stage='same-leading-full-target-Gevrey3',Z=max(transformed),L=L,physical_defect_prefactors=D,
  transformed_row_prefactors=transformed,actual_target_binding_normalizer_and_baseline_certified=False)
def root(tgt,cstar,beta=1000,qnorm=100000,identity=IDENTITY):
 cstar,beta,qnorm=map(exact,(cstar,beta,qnorm))
 if tgt.get('stage')!='same-leading-full-target-Gevrey3' or identity!=IDENTITY or cstar<0 or beta<=0 or qnorm<0 or 2*beta*qnorm*cstar>=Q(1,2):raise ValueError('SAME already accepted small leading root and retained Jacobian gate required')
 Z,L=exact(tgt['Z']),exact(tgt['L'])
 if Z<0 or L<1:raise ValueError('Nonnegative full-target prefactor and base >=1 required')
 gamma=2*beta;a=gamma*Z;b=gamma*qnorm
 T=min(a,1/(8*b)) if b else a
 B=8*L*max(Q(1),8*a*b)
 R=1/B
 return dict(stage='same-accepted-leading-quadratic-root-Gevrey3',positive_prefactor=a,base=B,value_cap=cstar,
  gamma=gamma,a=a,b=b,majorant_radius=R,majorant_ball=T,
  actual_original_acceptance_and_positive_field_gates_certified=False,slow_analytic_extension_claimed=False)
def bound(recipe,m,binary=False):
 order(m)
 if recipe.get('stage')!='same-accepted-leading-quadratic-root-Gevrey3':raise ValueError('Same leading quadratic root recipe required')
 P=recipe['value_cap'] if m==0 else recipe['positive_prefactor'];B=recipe['base']
 if binary:v=lg(max(Q(1),P))+m*lg(max(Q(1),B))+(3*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary budget at high order')
  v=P*B**m*factorial(m)**3
 return dict(order=m,upper=v,factorial_power=3,actual_profile_or_scale_certified=False)
def field(recipe,e_star):
 e_star=exact(e_star)
 if recipe.get('stage')!='same-accepted-leading-quadratic-root-Gevrey3' or e_star<=0:raise ValueError('Original leading root and same positive amplitude required')
 return dict(stage='same-fixed-bump-repair-field-Gevrey3',prefactor=2*e_star*max(recipe['value_cap'],recipe['positive_prefactor']),base=4+recipe['base'],
  radial_factor='actual Q63 B_r; ordinary X derivatives also pay x_min^(-r)*r!',
  actual_patch_supports_low_order_positive_field_gates_and_final_cone_certified=False)
def record():
 return dict(status='closed-Gevrey3-normalization-and-same-leading-quadratic-repair',
  normalization='M/(rho*e), I/(sqrt(2)*rho^(3/2)*e), J/(sqrt(2)*rho^(3/2)*e^2), S/(rho*e^2), Pi/e^2; transformed (m,(m-j)/lambda,i,-s,cp)',
  polynomial_costs='5/e_star for 1/e, 25/e_star^2 for 1/e^2; inverse lambda retained',
  root='positive jets <= (2*beta*Z)*[8*L*max(1,8*(2*beta)^2*q*Z)]^m*(m!)^3',
  memory='all five baseline target defects AND all five incoming modulation differences retained',
  scope='SAME accepted leading constant-B,Q root; not the earlier eta-dependent root or positive-order linear repair',
  slow_analytic_extension_claimed=False,actual_target_gates_final_field_cone_pulse_residual_or_blowup_certified=False)
