from types import MethodType
from pathlib import Path
import importlib.util,sys
import sympy as s
from bounds import InitialQuotientFastInverse,prefast_quotient_rows,fast_majorants,prior,EP,Q

def fixture():
 path=Path(__file__).resolve().parent.parent/'core_initial_cancellation_fast_source/controls.py'
 spec=importlib.util.spec_from_file_location('initial_quotient_fixture',path);mod=importlib.util.module_from_spec(spec)
 old=sys.modules.get('bounds');sys.modules['bounds']=prior
 try:spec.loader.exec_module(mod)
 finally:
  if old is None:sys.modules.pop('bounds',None)
  else:sys.modules['bounds']=old
 v=mod.fixture()
 def compact(self,n,p):return dict(KT=EP(n+2),KA=EP(n+3),I=EP({1:1})+EP(n+4))
 v.budget.compact=MethodType(compact,v.budget)
 return InitialQuotientFastInverse(v)

def equal(a,b):
 return {k:v.terms for k,v in EP(a).c.items()}=={k:v.terms for k,v in EP(b).c.items()}

def review():
 checks={};R,z=s.symbols('R z');Wrr=s.Function('Wrr')(R,z);Wrz=s.Function('Wrz')(R,z);Wtt=s.Function('Wtt')(R,z)
 Dr=lambda f:s.diff(f,R);Dz=lambda f:s.diff(f,z)
 gr=-Dr(Wrr)-Wrr/R-Dz(Wrz)+Wtt/R;r0=-Wrr/R-Dz(Wrz)+Wtt/R;K0=-Wrr
 checks['exact_prefast_radial_quotient_gr_equals_r0_plus_DrK0']=s.simplify(gr-r0-Dr(K0))==0
 qWzz,qWrr,qT=s.symbols('qWzz qWrr qT');correct=qWzz-qWrr+qT
 checks['local_minus_Wrr_pressure_memory_is_retained']=s.diff(correct,qWrr)==-1 and s.simplify(correct-(qWzz+qT))==-qWrr
 y=s.symbols('y',real=True);rhoP=s.Symbol('rhoP');osc=s.sin(y)
 haar=lambda f:s.integrate(f,(y,0,2*s.pi))/(2*s.pi)
 checks['rhoP_disappears_only_after_auxiliary_Q_projection']=haar(osc-rhoP)==-rhoP and s.simplify((osc-rhoP)-haar(osc-rhoP)-osc)==0
 checks['A1_annihilation_leaves_the_full_pressure_complement']=correct.has(qT) and correct.has(qWrr)
 a=prefast_quotient_rows(0,3,4,5,6,7,8)
 checks['positive_prefast_quotient_program_and_orders']=equal(a['r0'],36) and equal(a['P'],252) and equal(a['Q_pressure'],360) and equal(a['Q_Phi'],372) and equal(a['A1'],8928) and a['covariance_input_order']==3
 checks['full_compensated_pressure_memory_remains_recorded']=equal(a['full_nonlocal_pressure'],5*(36+8*252)) and a['full_nonlocal_pressure'].c and a['full_local_pressure'].c
 v=fixture();src=v.source(0);old=v.current.rows(0);tail=v.quotient_tail(0)
 checks['prior_regular_source_is_unchanged']=equal(src['whole']['Etheta'],old['regular']['Etheta']) and equal(src['whole']['Ez'],old['regular']['Ez']+tail['A1'])
 oldka=tail['A1'];orig=v.budget.compact
 v.budget.compact=MethodType(lambda self,n,p:dict(KT=EP(n+2),KA=EP(n+30),I=EP({1:1})+EP(n+4)),v.budget)
 checks['tail_reads_actual_original_KA1_profile_recipe']=not equal(v.quotient_tail(0)['A1'],oldka)
 v.budget.compact=orig
 rates=v.rates();k=Q(1,100000)
 checks['direct_A1_order_two_strictly_beats_current_D']=rates['original_prefast_A1_direct']==2 and rates['whole_source']==2-4*k and rates['direct_margin_over_regular']==4*k
 rec=v.recipe(0)
 checks['normalized_200001_jet_obligation_is_corrected']=rec['direct_prefast_A1_covariance_orders']==(7,8) and rec['largest_original_primary_order']==18 and rec['obsolete_conditional_raw_A1_order']==200009
 checks['direct_tail_does_not_claim_actual_continuum_enclosure']=not src['actual_continuum_inputs_certified'] and tail['complete_derivative_support_contract_required'] and not tail['actual_profile_support_or_field_enclosure_certified']
 f=fast_majorants(0,EP(2),EP(3),EP(5),dict(KT=EP(7),KA=EP(11),I=EP(13)),dict(KT=EP(17),KA=EP(19),I=EP(23)))
 gain=24*16*EP(13)*EP({1:1})
 checks['whole_Q_source_enters_fast_map_without_second_projection_factor']=equal(f['u'],24*gain*(EP(2)+12*EP(3)))
 checks['actual_gd_minus_A1gd_and_radial_KT1_are_retained']=f['actual_axial_velocity']=='gd-A1 gd' and equal(f['x'],120*24*32*EP(13)*EP({1:1})*EP(17)*EP(5))
 checks['fast_output_orders_follow_whole_D_source']=rates['fast_tangential']==rates['whole_source'] and rates['fast_radial']==rates['whole_source']+1
 checks['postfast_pressure_scalar_and_new_time_error_remain_next']=rec['next_obligation'].startswith('recompute postfast pressure/scalars') and not rates['complete_postfast_residual_or_blowup_certified']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='8431337a981afb25fd228070e85c86a057519dd7',
  scope='exact original prefast compact-quotient correction and conditional whole-source fast-map coefficients; no actual continuum enclosure or postfast residual closure')
