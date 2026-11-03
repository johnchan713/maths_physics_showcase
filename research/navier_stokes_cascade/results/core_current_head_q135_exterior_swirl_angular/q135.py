"""Conditional angular envelope of the SAME complete exterior leading swirl."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('q135_heat',HERE.parent/'core_current_head_q134_real_heat_gevrey/q134.py')
heat=importlib.util.module_from_spec(s);s.loader.exec_module(heat)
exact=heat.exact;order=heat.order;lg=heat.lg
IDENTITY='same-whole-exterior-reference-factorization-and-accepted-heat-angular-roots'
def exterior(S0,lam,Xstar,identity):
 S0,lam,Xstar=map(exact,(S0,lam,Xstar))
 if identity!=IDENTITY or S0<=0 or not 0<lam<=Q(1,2**256) or Xstar<100:raise ValueError('Same full exterior schedule, positive scalar cap, selected lambda range and Xstar required')
 h=lam**2;g=h**4/Xstar;r=2**19*lam**30
 root=heat.root(heat.tail(h,Xstar,heat.IDENTITY),heat.IDENTITY)
 factors=dict(gamma=Q(81),heat_compensation=1+192*root['positive_prefactor'],angular_edit=1+128*r,terminal_taper=1/(1-h/1000),terminal_heat=1+4*h/Xstar)
 P=S0
 for v in factors.values():P*=v
 return dict(stage='same-whole-exterior-leading-E-and-F0-raw-angular-Gevrey2',prefactor=P,base=Q(548),factors=factors,
  raw_factor_bases=(2,512,2,0,32),heat_g=g,angular_root_cap=r,
  actual_reference_factorization_scalar_cap_roots_supports_and_matching_certified=False,
  complex_Gamma_extension_claimed=False,actual_mixed_radial_or_physical_pulse_certified=False)
def norm(recipe,m,binary=False):
 order(m)
 if recipe.get('stage')!='same-whole-exterior-leading-E-and-F0-raw-angular-Gevrey2':raise ValueError('Same complete exterior recipe required')
 P,B=recipe['prefactor'],recipe['base']
 if binary:v=lg(max(Q(1),P))+lg(Q(m+1))+m*lg(B)+(m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary budget at high angular order')
  v=P*sum(B**j*factorial(j) for j in range(m+1))
 return dict(stage='whole-exterior-leading-F0-factorial-angular-norm',order=m,upper=v,binary=binary,
  norm='sup_X>=Xstar sum_j sup_eta |d_eta^j F0|/j!',also_bounds_stronger_sum_of_global_suprema=True,
  actual_axis_final_compact_cover_mixed_radial_or_scale_certified=False)
def record():
 return dict(status='whole-exterior-leading-swirl-angular-Gevrey2-envelope',
  raw='E0 and F0 <= P*548^m*(m!)^2 for X>=Xstar',
  P='81*S0*(1+15360000*h^4/Xstar)*(1+128*2^19*lambda^30)*(1+4*h/Xstar)/(1-h/1000)',
  factorial_norm='P*sum_(j=0)^m 548^j*j!; binary ceilLog(max(1,P))+ceilLog(m+1)+10m+m*ceilLog(m)',
  schedule='same nonincreasing eta-independent scalar S<=S0, gamma, three heat bumps, two angular bumps, exact two-unit terminal taper and Gamma multiplier',
  pulse='axial profile strength changes U, not the prescribed scalar swirl on that stage; no physical wave estimate follows',
  pressure='same Q134 infinite heat target and Q126 pre-heat future-mass root; no incoming pressure or other moment reset',
  remaining_inequality='whole F0 angular norm requires actual regular-axis and finalized compact bounds on the SAME cover; Q121 also needs mixed X,eta derivatives',
  actual_whole_F0_cone_pulse_residual_or_blowup_certified=False,historical_frozen_records_modified=False)
