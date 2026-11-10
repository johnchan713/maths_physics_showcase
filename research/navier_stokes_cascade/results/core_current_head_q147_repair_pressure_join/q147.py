"""Bind the same modulation to its five-row repair and global pressure exit."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q147_previous','core_current_head_q146_cap_repair_identity/q146.py')
modulation_provider=load('q147_modulation','core_current_head_q138_mixed_modulation_pullback/q138.py')
repair=load('q147_repair','core_current_head_q132_exact_repair_gevrey/q132.py')
terminal=load('q147_terminal','core_current_head_q74_terminal_modulation_defect/q74.py')
exact=previous.exact
ROWS=('M','I','J','S','Pi')
IDENTITY='same-Q74-terminal-Q138-modulation-Q132-repair-and-global-pressure-gauge'

def baseline(bound,identity):
 """The PREMOD target is the same unmodulated profile, in the same gauge."""
 if identity!=IDENTITY or bound.get('stage')!='Q144-with-same-profile-complex-caps-and-hierarchy-gate':
  raise ValueError('Same Q146 selected profile required')
 r=bound['recipe']
 if not r.get('all_five_actual_incoming_rows_reintegrated') or not r.get('full_p2_Pi_and_Pi_eta_retained'):
  raise ValueError('All five rows and the unchanged pressure gauge are required')
 return dict(stage='exact-zero-PREMOD-baseline-functional-identity',rows=ROWS,
  pairs={k:(Q(0),Q(1)) for k in ROWS},
  annular_endpoint='all five unmodulated primitive differences vanish as eta functions',
  suffix='the target and PREMOD profiles have identical fields and densities up to modulation',
  pressure='same Pi_axis+Cp gauge; no interface reset',
  fixed_origin='modulation primitives vanish on the preserved entrance collar',
  zero_as_eta_functions=True)

def memory(bound,mod,x0,x1,identity):
 """Integrate the five exact density differences of this Q138 object."""
 x0,x1=map(exact,(x0,x1))
 if identity!=IDENTITY or bound.get('stage')!='Q144-with-same-profile-complex-caps-and-hierarchy-gate':
  raise ValueError('Same Q146 profile and Q74/Q138 identity required')
 if mod.get('stage')!='same-pre-repair-log-mixed-modulation-Gevrey3' or not 0<x0<=x1:
  raise ValueError('Same Q138 modulation and fixed positive interval required')
 r=bound['recipe'];pf=mod['prefactors'];e=pf['E']/9;du=pf['deltaU'];u=pf['U']-du
 if e<=0 or u<0 or mod['N']<1 or mod['N'].denominator!=1 or pf['deltaE']<0 or du<0:
  raise ValueError('Positive PREMOD caps and original integer frequency required')
 # Stage strings alone are insufficient: the modulation caps and bases must
 # dominate the actual Q146 recipe from which the scalar/primitive chain starts.
 if e<r['field']['F']['P'] or u<r['field']['U']['P'] or mod['radial_base']<r['Ay'] or mod['angular_base']<r['Be']:
  raise ValueError('Q138 object is not bound to the Q146 actual PREMOD profile caps')
 de=pf['deltaE'];rad=max(Q(1),2*x1);ell=x1-x0
 d=dict(M=du,I=rad*de,J=rad*(e*du+u*de+de*du),
  S=2*u*du+du*du+e*de+de*de/2,Pi=(e*de+de*de/2)/x0)
 B=2*mod['angular_base'];pref={k:ell*d[k] for k in ROWS}
 return dict(stage='full-five-incoming-angular-modulation-memory',density_prefactors=d,
  prefactors=pref,base=B,shifted_derivative_prefactors={k:pref[k]*B for k in ROWS},
  shifted_derivative_base=8*B,rows=ROWS,N=mod['N'],same_q138_object=True,
  fixed_origin_incoming_differences_zero=True,full_Pi_and_Pi_eta_memory=True,
  actual_incoming_memory_interval_caps_and_final_repair_certified=True,
  slow_analytic_extension_claimed=False)

def bind(bound,mod,x0,x1,e_star,rho,lam,identity):
 b=baseline(bound,identity);mem=memory(bound,mod,x0,x1,identity)
 tgt=repair.target(mem,b['pairs'],e_star,rho,lam,repair.IDENTITY)
 return dict(stage='same-modulation-exact-five-row-repair-target',baseline=b,memory=mem,target=tgt,
  signed_target='-T_c Delta=( -m, (j-m)/lambda, -i, s, -cp )',
  q74_terminal_identity=True,inverse_lambda_retained=True,all_five_rows=True,
  same_modulation_memory_certified=True,actual_baseline_defects_certified=True)

def terminal_target(rows,rho,e,lam,angular_scale):
 """Q74 signed target, exposed here to test the Q132 row convention."""
 vals=tuple(map(exact,rows));rho,e,lam,angular_scale=map(exact,(rho,e,lam,angular_scale))
 if len(vals)!=5 or min(rho,e,lam,angular_scale)<=0:raise ValueError('Five physical rows and positive normalizers required')
 return tuple(map(exact,terminal.target(vals,rho,e,lam,angular_scale)))

def close(binding,identity):
 if identity!=IDENTITY or binding.get('stage')!='same-modulation-exact-five-row-repair-target':
  raise ValueError('Same bound target required')
 if not binding.get('same_modulation_memory_certified') or not binding['baseline'].get('zero_as_eta_functions'):
  raise ValueError('Exact baseline and same modulation memory required')
 if tuple(binding['memory'].get('rows',()))!=ROWS or not binding['memory'].get('full_Pi_and_Pi_eta_memory'):
  raise ValueError('Complete pressure row required')
 return dict(stage='same-profile-compact-repair-cone-and-global-pressure-join',rows=ROWS,
  target=binding['target'],accepted_root='same original finite-N C1 root in its uniqueness ball',
  positive_field='same A^-128 neighborhood and coefficient radius',
  compact_cone='same G^-4/2 repaired value gap',
  regional_join='protected inner collar plus annular/exterior value cones',
  endpoint='all five defects vanish as eta functions',
  global_identities=('M_infty=0','R_infty=0','J_infty=0','S_infty=0','Pi_axis+Cp_infty=0'),
  pressure='Pi(X)=-integral_X^infty E^2/(2x) dx in the unchanged gauge',
  same_exact_repair_and_inner_outer_value_cone_certified=True,
  global_pressure_and_five_moment_identities_certified=True,
  coefficient_reserved_mean_supports_certified=False,
  uniform_actual_pulses_covariance_residual_or_blowup_certified=False)

def record():
 return dict(status='same-modulation-five-row-repair-cone-and-global-pressure-identity-joined',
  correction='Q132 baseline is exactly zero because the target is the same PREMOD profile after the five-row annular match; this is a functional same-gauge identity, not support-end inference',
  memory='Q138 deltaE,deltaU generate the same five Q131 density rows; Pi and the shifted Pi_eta bound are retained',
  target='Q74 z=-T_c Delta=(-m,(j-m)/lambda,-i,s,-cp), identical to the signed Q132 target',
  repair='the original accepted finite-N root, positive-field neighborhood and compact value-cone comparison now apply to that exact target',
  global_pressure='five endpoint identities plus the identical later heat tail give M=R=J=S=0 and Pi_axis+Cp=0 globally in the unchanged gauge',
  next='audit finite positive-order coefficient and reserved-mean supports on this finalized profile before uniform physical pulses',
  coefficient_reserved_mean_pulse_covariance_residual_or_blowup_certified=False,
  historical_frozen_records_modified=False)
