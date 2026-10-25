from pathlib import Path
import importlib.util,sys
from bounds import (EP,Q,a1_phi_recipe,final_row_envelope,complete_final_pressure,
 rate_ledger,strict_fast_closure,InitialFinalRowObstruction,previous,relaxation)

def equal(a,b):return {k:v.terms for k,v in EP(a).c.items()}=={k:v.terms for k,v in EP(b).c.items()}

def fixture():
 path=Path(__file__).resolve().parent.parent/'core_initial_postfast_refresh/controls.py'
 spec=importlib.util.spec_from_file_location('initial_final_fixture',path);mod=importlib.util.module_from_spec(spec)
 old=sys.modules.get('bounds');sys.modules['bounds']=previous
 try:spec.loader.exec_module(mod)
 finally:
  if old is None:sys.modules.pop('bounds',None)
  else:sys.modules['bounds']=old
 return previous.InitialPostfastRefresh(mod.fixture())

def review():
 checks={};n=0
 old={k:EP(i+2) for i,k in enumerate(('x','u','z'))};fast={k:EP(i+1) for i,k in enumerate(('x','u','z'))};slow={k:EP(2*i+1) for i,k in enumerate(('x','u','z'))}
 fixed={key:EP(i+1) for i,key in enumerate(('b','V','G','graph','c','r','rho','Kp','Kmom'))};haar=dict(Etheta=EP(11),Ez=EP(13));A1=EP(17)
 out=final_row_envelope(n,haar,A1,old,fast,slow,fixed)
 direct_fast=previous.previous.prior.fast.difference(n,old,fast,fixed);after={k:old[k]+fast[k] for k in old};direct_slow=previous.previous.prior.fast.difference(n,after,slow,fixed)
 checks['exact_final_row_adds_Haar_A1_fast_and_slow_once']=equal(out['Etheta'],haar['Etheta']+direct_fast['Etheta']+direct_slow['Etheta']) and equal(out['Ez'],haar['Ez']+A1+direct_fast['Ez']+direct_slow['Ez'])
 checks['sharp_slow_difference_uses_old_plus_actual_fast_state']=all(equal(out['slow_difference'][k],direct_slow[k]) for k in ('Etheta','Ez','pm'))
 checks['both_differences_retain_pressure_and_viscosity']=all(out[group][key].c for group in ('fast_difference','slow_difference') for key in ('pm','pressure_gradient','viscosity'))
 checks['full_differences_receive_no_extra_Qcirc_factor']=not out['extra_Qcirc_factor_on_full_differences']
 post=dict(parts={name:EP(i+1) for i,name in enumerate(('prefast','initial_fast','initial_slow','signed_wave','new_fast'))},pm=EP(15))
 fp=complete_final_pressure(post,direct_slow['pm'])
 checks['final_pressure_keeps_all_five_postfast_parts_and_sharp_slow']=tuple(fp['parts'])==('prefast','initial_fast','initial_slow','signed_wave','new_fast','sharp_slow') and equal(fp['pm'],EP(15)+direct_slow['pm'])
 geometry=dict(L='4',width='1/2',gap='1');a=a1_phi_recipe(0,geometry=geometry,pole=12)
 checks['A1_Phi_recipe_has_least_strict_count_and_one_extra_Phi_jet']=a['integration_order']==99995 and a['raw_Dz_Phi_order']==199993 and a['raw_Phi_order']==199994 and a['epsilon_order']=='59999/20000'
 checks['A1_recipe_retains_raw_profiles_floor_edge_and_same_source_pressure']=all(x in a['factors'] for x in ('N_0(k_1 chi_U)','N_199993(f_1)','Dz Qcirc Phi_postsigned')) and a['source_stage']=='post-signed/pre-new-fast' and not a['postfast_or_postslow_pressure_substituted']
 checks['A1_recipe_is_symbolic_not_actual_high_jet_certificate']=not a['materialized_actual_input_bound'] and not a['actual_Phi_profile_support_and_high_jet_bounds_materialized']
 r=rate_ledger();k=Q(1,100000)
 checks['complete_signed_covariance_is_unique_source_rate_bottleneck']=r['whole_source']==r['D'] and r['source_zero_margin_groups']==('complete_signed_covariance_response',) and r['source_margins']['corrected_prefast_A1']==4*k
 checks['radial_viscosity_is_unique_zero_margin_fast_group']=r['fast_difference']==r['nu'] and r['unique_zero_margin_fast_groups']==('radial_viscosity',) and r['fast_margins']['base_radial_flux']==k
 checks['every_sharp_slow_group_has_positive_margin']=r['slow_difference']==r['alpha']+1 and min(r['slow_margins'].values())==3*k
 inherited=relaxation.exponent_review();checks['rate_ledger_matches_independent_frozen_first_real_budget']=Q(inherited['fast_remainder_orders']['radial_viscosity'])==r['fast']['radial_viscosity'] and inherited['checks']['fast_viscosity_is_limiting']
 block=strict_fast_closure();checks['strict_closure_requires_positive_source_or_viscosity_gain']=not block['strictly_closes'] and block['strict_margin']==0 and block['required_inequality']=='delta_source+delta_viscosity>0'
 checks['either_positive_gain_would_resolve_this_specific_equality']=strict_fast_closure(k,0)['strictly_closes'] and strict_fast_closure(0,k)['strictly_closes']
 checks['conditional_nu_row_is_not_complete_residual_or_blowup_proof']=a['strict_margin']==k and not block['strictly_closes']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='3a966043f2b68ad1ab8f6b190c409df5734fc9da',
  scope='exact conditional first-real final-row and full pressure assembly, explicit same-state A1 raw recipe, and sharp radial-viscosity equality obstruction; no actual high-jet enclosure or complete residual closure')
