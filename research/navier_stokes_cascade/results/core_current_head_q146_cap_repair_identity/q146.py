"""Specialize Q144 caps and retain the unresolved five-row repair identity."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('q146_previous','core_current_head_q145_hierarchy_gate/q145.py')
repair=load('q146_repair','core_current_head_q132_exact_repair_gevrey/q132.py')
exact=previous.previous.exact
IDENTITY='same-Q73-strip-Q143-endpoint-Q127-root-and-full-pressure-profile'
REPAIR_IDENTITY='same-Q138-modulation-Q131-five-memory-and-Q132-leading-repair'
ROWS=('M','I','J','S','Pi')

def cap_tuple(H,Z,B0,J6,e,v,C,R,P,length,endpoint_E,endpoint_invE,identity,materialize=True):
 vals=tuple(map(exact,(H,Z,B0,J6,e,v,C,R,P,length,endpoint_E,endpoint_invE)))
 H,Z,B0,J6,e,v,C,R,P,length,endpoint_E,endpoint_invE=vals
 if identity!=IDENTITY or not 0<H<=1 or Z<0 or B0<1 or not 0<=J6<=B0/100:
  raise ValueError('Same Q73 strip, selected shape hierarchy and Q127 root required')
 if min(e,v,C,R,P,length,endpoint_E,endpoint_invE)<=0 or C<1 or P<16 or R<110:
  raise ValueError('Positive actual entry, geometry, annular and endpoint caps required')
 if length.denominator!=1 or length<2000*B0+10 or R/110>2**length.numerator:
  raise ValueError('Same certified shape/annular logarithmic length required')
 L0=J6+e*v;x=200*B0+L0+4
 n=(x.numerator+x.denominator-1)//x.denominator
 out=dict(stage='same-profile-Q73-Q127-cap-specialization',H=H,Z=Z,B0=B0,J6=J6,L0=L0,
  exponent=n,C=C,R=R,P=P,length=length,endpoint_E=endpoint_E,endpoint_invE=endpoint_invE,
  formulas=dict(Ecap='max(endpoint_E,3^n/C,4P)',invEcap='max(endpoint_invE,C*3^n,20/P)'),
  same_profile_complex_cap_formula_certified=True,selected_huge_integer_materialized=False)
 if materialize:
  if n>1000000:raise ValueError('Keep the selected enormous cap symbolic')
  p3=Q(3)**n;out.update(Ecap=max(endpoint_E,p3/C,4*P),invEcap=max(endpoint_invE,C*p3,20/P),
   selected_huge_integer_materialized=True)
 return out

def bind(post,caps,identity):
 if identity!=IDENTITY or caps.get('stage')!='same-profile-Q73-Q127-cap-specialization' or not caps.get('selected_huge_integer_materialized'):
  raise ValueError('Moderate exact specialization required for executable Q144 control')
 r=previous.previous.shape_annulus(post,caps['H'],caps['Z'],caps['B0'],caps['J6'],caps['L0'],
  caps['Ecap'],caps['invEcap'],caps['R'],caps['P'],caps['length'],previous.previous.IDENTITY)
 reduced=previous.reduce_gate(r,previous.IDENTITY)
 return dict(stage='Q144-with-same-profile-complex-caps-and-hierarchy-gate',recipe=r,reduction=reduced,
  same_profile_complex_cap_formula_certified=True,selected_profile_huge_cap_materialized=False,
  exact_repair_identity_and_cone_join_certified=False)

def repair_request(bound,identity):
 if identity!=REPAIR_IDENTITY or bound.get('stage')!='Q144-with-same-profile-complex-caps-and-hierarchy-gate':
  raise ValueError('Same cap-bound modulation region and full-pressure repair identity required')
 r=bound['recipe']
 if not r.get('all_five_actual_incoming_rows_reintegrated') or not r.get('full_p2_Pi_and_Pi_eta_retained'):
  raise ValueError('All five rows and Pi/Pi_eta pressure memory required')
 return dict(stage='actual-five-row-exact-repair-binding-request',rows=ROWS,
  baseline_required={k:'T_%s,L_%s with |d_%s^(m)|<=T_%s L_%s^m (m!)^3'%(k,k,k,k,k) for k in ROWS},
  normalized_target='(m,(m-j)/lambda,i,-s,cp)',inverse_lambda_retained=True,
  same_modulation_memory_certified=False,actual_baseline_defects_certified=False,
  accepted_root_certified=False,positive_field_neighborhood_certified=False,
  exact_repair_and_inner_outer_cone_join_certified=False)

def validate_baseline_shape(baseline):
 if set(baseline)!=set(ROWS):raise ValueError('Every M,I,J,S,Pi baseline row is required')
 for k,v in baseline.items():
  if not isinstance(v,(tuple,list)) or len(v)!=2:raise ValueError('Each row needs a prefactor/base pair')
  T,L=map(exact,v)
  if T<0 or L<1:raise ValueError('Nonnegative prefactor and Gevrey base at least one required')
 return dict(rows=ROWS,pressure_row_retained='Pi' in baseline,complete_shape_only=True,
  actual_function_inequalities_certified=False)

def record():
 return dict(status='same-profile-complex-cap-specialization-and-exact-repair-identity-obstruction',
  cap='L0=J6+e*v; n=ceil(200B0+L0+4); Ecap=max(endpoint,3^n/C,4P); invEcap=max(endpoint inverse,C*3^n,20/P)',
  correction='Q73 and Q127 already analytically discharge Q145 same-profile cap inequality; selected huge caps need not be materialized',
  repair='stage compatibility does not bind Q138 modulation to Q131 five-row memory and Q132 actual baseline target',
  blocking_inequality='for every r in M,I,J,S,Pi prove |partial_eta^m d_r^base|<=T_r L_r^m(m!)^3 for the same modulation, then retain (m,(m-j)/lambda,i,-s,cp) and the original root/positive-field gates',
  exact_repair_cone_join_pulse_covariance_residual_or_blowup_certified=False,
  historical_frozen_records_modified=False)
