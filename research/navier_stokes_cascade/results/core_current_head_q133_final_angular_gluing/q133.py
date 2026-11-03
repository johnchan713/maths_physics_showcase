"""Final compact angular assembly with explicit order of spatial suprema."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import importlib.util
HERE=Path(__file__).resolve().parent
def load(name,path):
 s=importlib.util.spec_from_file_location(name,HERE.parent/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
repair=load('q133_repair','core_current_head_q132_exact_repair_gevrey/q132.py')
regional=load('q133_regional','core_current_head_q127_annular_radius_axis/q127.py')
exact=repair.exact;order=repair.order;lg=repair.lg
IDENTITY='same-final-leading-compact-cover-including-modulation-repair-and-gaps'
# Q63's conservative exact all-order B_r at r=0; no numerical bump maximum.
BUMP0=Q(162*32600)
def compact(modulation,repair_field,identity):
 if identity!=IDENTITY or modulation.get('stage')!='angular-modulation-before-exact-repair' or repair_field.get('stage')!='same-fixed-bump-repair-field-Gevrey3':raise ValueError('Same final leading compact cover, Q131 modulation and Q132 repair field required')
 P=modulation['prefactors']['U']+BUMP0*repair_field['prefactor']
 B=max(modulation['base'],repair_field['base'])
 return dict(stage='final-leading-compact-raw-angular-Gevrey3',prefactor=P,base=B,
  actual_whole_compact_PREMOD_caps_flat_extensions_target_and_acceptance_gates_certified=False,
  slow_analytic_extension_claimed=False)
def compact_norm(recipe,m,binary=False):
 order(m)
 if recipe.get('stage')!='final-leading-compact-raw-angular-Gevrey3':raise ValueError('FINAL compact raw angular recipe required; PREMOD alone is not sufficient')
 P,B=recipe['prefactor'],recipe['base']
 if binary:v=lg(max(Q(1),P))+lg(Q(m+1))+m*lg(max(Q(1),B))+(2*m*lg(Q(m)) if m else 0)
 else:
  if m>12:raise ValueError('Use binary budget at high angular order')
  v=P*sum(B**j*factorial(j)**2 for j in range(m+1))
 return dict(stage='final-leading-compact-modulation-and-repair',order=m,upper=v,binary=binary,
  norm='sum_j sup_compact_X,eta |d_eta^j U|/j!',controls_pointwise_X_angular_norm=True,actual_profile_and_scale_certified=False)
def regional_sum(m,axis,final_compact,pulse_identity,binary=False):
 order(m)
 for x,stage in ((axis,'leading-axis'),(final_compact,'final-leading-compact-modulation-and-repair')):
  if not isinstance(x,dict) or x.get('stage')!=stage or x.get('order')!=m or x.get('binary') is not binary:raise ValueError('Same-order, same-convention axis and FINAL compact norm recipes required')
 exterior=regional.current.axial_angular(m,pulse_identity,binary)
 if binary:
  order(axis['upper']);order(final_compact['upper']);v=max(axis['upper'],final_compact['upper'],3,exterior['upper'])+2
 else:
  a,c=exact(axis['upper']),exact(final_compact['upper'])
  if min(a,c)<0:raise ValueError('Nonnegative regional norms required')
  v=a+c+8+exterior['upper']
 return dict(stage='strong-whole-leading-factorial-angular-norm',order=m,upper=v,binary=binary,
  norm='sum_j sup_whole_profile |d_eta^j U0|/j!',gluing='regional sum for stronger sum-of-global-suprema norm',
  actual_original_axis_exterior_compact_cover_and_gates_certified=False,actual_mixed_radial_tuple_or_cone_certified=False)
def jetwise(m,axis_value,delta,final_compact,pulse_identity):
 order(m);axis_value,delta=exact(axis_value),exact(delta)
 if m>12 or axis_value<0 or not 0<delta<=1 or final_compact.get('stage')!='final-leading-compact-raw-angular-Gevrey3' or pulse_identity!=regional.current.IDENTITY:raise ValueError('Same original regional derivative caps and moderate FINAL compact order required')
 P,B=final_compact['prefactor'],final_compact['base']
 value=sum(max(axis_value*delta**(-j),P*B**j*factorial(j)**2,Q(4) if j<=1 else Q(0),Q(40)*2**j) for j in range(m+1))
 return dict(stage='strong-orderwise-whole-leading-angular-norm',order=m,upper=value,
  actual_regional_caps_and_cover_certified=False)
def regional_max(m,axis,final_compact,pulse_identity,binary=False):
 # The regional sum call validates stage/order/convention; the intended Q125 norm permits max.
 regional_sum(m,axis,final_compact,pulse_identity,binary)
 exterior=regional.current.axial_angular(m,pulse_identity,binary)
 value=max(axis['upper'],final_compact['upper'],3 if binary else Q(8),exterior['upper'])
 return dict(stage='whole-leading-pointwise-X-angular-algebra-norm',order=m,upper=value,binary=binary,
  norm='sup_X sum_j sup_eta |d_eta^j U0(X,eta)|/j!',
  actual_original_caps_cover_target_and_acceptance_certified=False,stronger_sum_of_global_suprema_certified=False)
def first_input(recipe,discrepancy_order):
 order(discrepancy_order)
 if recipe.get('stage')!='whole-leading-pointwise-X-angular-algebra-norm' or recipe.get('order')!=discrepancy_order+2:raise ValueError('Q125 requires the FINAL whole-support pointwise-X angular norm through m+2')
 return dict(U0_upper=recipe['upper'],input_order=recipe['order'],binary=recipe['binary'],
  actual_F0_core_collar_support_scale_or_discrepancies_certified=False)
def counterexample():
 # Smooth rectangle control, not the physical Navier-Stokes profile.
 return dict(weak_axis=Q(225,2),weak_compact=Q(125),weak_global=Q(125),
  strong_axis=Q(150),strong_compact=Q(325,2),strong_global=Q(200),
  invalid_strong_max=Q(325,2),safe_strong_two_region_sum=Q(625,2))
def record():
 return dict(status='final-leading-angular-assembly-and-norm-order-boundary',
  final_compact='P=Q131_U_pref+B0*Q132_repair_field_pref; B=max(the two bases); raw <= P*B^m*(m!)^3',
  factorial_norm='P*sum_(j=0)^m B^j*(j!)^2; not the raw m-th derivative bound',
  intended_Q125_norm='sup_X sum_j sup_eta |d_eta^j U0(X,eta)|/j!; regional maximum is valid',
  stronger_norm='sum_j sup_X,eta |d_eta^j U0|/j!; regional sum or orderwise maximum required',
  correction_scope='The smooth counterexample blocks a stronger reinterpretation, not the historical intended pointwise-X maximum rule',
  binary='stronger sum gluing costs at most two extra bits for four regions; weak intended norm uses max',
  pressure='Q132 full baseline and incoming five-row target premises remain; no pressure reset',
  actual_profile_gates_mixed_radial_support_scale_pulse_residual_or_blowup_certified=False,
  historical_frozen_records_modified=False)
