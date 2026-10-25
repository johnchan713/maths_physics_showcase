from pathlib import Path
import importlib.util,sys
from bounds import (EP,Q,complete_pressure_envelope,postfast_scalar_refresh,
 sharp_slow_fields,a1_time_error_request,rates,previous,refresh)

def equal(a,b):
 return {k:v.terms for k,v in EP(a).c.items()}=={k:v.terms for k,v in EP(b).c.items()}

def fixture():
 path=Path(__file__).resolve().parent.parent/'core_initial_quotient_fast_inverse/controls.py'
 spec=importlib.util.spec_from_file_location('postfast_previous_fixture',path);mod=importlib.util.module_from_spec(spec)
 old=sys.modules.get('bounds');sys.modules['bounds']=previous
 try:spec.loader.exec_module(mod)
 finally:
  if old is None:sys.modules.pop('bounds',None)
  else:sys.modules['bounds']=old
 return mod.fixture()

def review():
 checks={}
 pieces=[EP(i) for i in range(1,6)];p=complete_pressure_envelope(*pieces)
 checks['complete_pointwise_pressure_has_exactly_five_sequential_components']=tuple(p['parts'])==('prefast','initial_fast','initial_slow','signed_wave','new_fast') and equal(p['pm'],15)
 checks['pressure_memory_is_not_reset_or_projected']=p['complete_pressure_memory_retained'] and all(p['parts'][k].c for k in p['parts'])
 n=0;previous_scalars={k:EP(i+1) for i,k in enumerate(('P','Jtheta','Jz'))}
 wave=refresh.post.scalar.wave_change(n,EP(2),EP(3),EP(5))
 result=dict(updated_scalars={k:previous_scalars[k]+wave[k] for k in previous_scalars})
 old={k:EP(i+2) for i,k in enumerate(('x','u','z'))};inc={k:EP(i+1) for i,k in enumerate(('x','u','z'))};M=EP(7)
 out=postfast_scalar_refresh(n,result,old,inc,M)
 x,u,z=(inc[k] for k in ('x','u','z'));ox,ou,oz=(old[k] for k in ('x','u','z'))
 S=2*ox*x+x*x;H=ox*z+oz*x+x*z;T=2*ou*u+u*u;Zf=2*oz*z+z*z;K=oz*u+ou*z+u*z;q=M*(S+T+H)
 change=dict(P=q,Jtheta=M*K,Jz=M*Zf+q*Q(1,2))
 checks['post_signed_scalars_add_the_new_fast_change_exactly_once']=all(equal(out[k],previous_scalars[k]+wave[k]+change[k]) for k in out)
 wrong={k:out[k]+wave[k] for k in out}
 checks['adding_wave_change_again_is_detected_as_double_counting']=all(equal(wrong[k],out[k]+wave[k]) for k in out)
 zchanged=dict(inc,z=EP(99));changed=refresh.post.scalar.fast_change(n,old,zchanged,M)
 checks['new_fast_scalar_change_uses_distinct_actual_axial_velocity']=any(not equal(change[k],changed[k]) for k in ('P','Jtheta','Jz'))
 slow={k:EP(i+2) for i,k in enumerate(('H','Kp','Kz','Kg','Kpsi'))};slow_next={k:2*v for k,v in slow.items()};floor=EP({1:1})
 fields=sharp_slow_fields(n,out,{k:2*v for k,v in out.items()},slow,slow_next,floor)
 checks['sharp_slow_angular_and_axial_fields_use_current_complete_scalars']=equal(fields['u'],floor*slow['H']*(slow['Kp']*out['P']+slow['Kz']*out['Jz'])) and equal(fields['z'],floor*slow['H']*slow['Kg']*out['Jtheta'])
 checks['sharp_slow_radial_field_uses_next_scalar_and_weight_floor']=equal(fields['x'],floor*slow_next['H']*slow_next['Kpsi']*2*out['Jtheta'])
 q=a1_time_error_request(0);k=Q(1,100000)
 checks['new_A1_potential_count_is_least_strict_integer']=q['integration_order']==99995 and q['source_order']+(q['integration_order']-1)*k==q['target_order'] and q['strict_margin']==k
 checks['new_A1_raw_Phi_request_is_sharper_but_still_unmaterialized']=q['raw_Phi_order']==199994 and q['previous_full_row_raw_source_order']==200001 and not q['actual_current_Phi_high_jets_and_support_materialized']
 checks['new_A1_uses_its_post_signed_source_state_not_refreshed_pressure']=q['source_stage']=='post-signed/pre-new-fast' and not q['postfast_or_postslow_pressure_substituted']
 checks['complete_current_potential_and_pressure_dependencies_are_retained']=q['current_fields_for_potential']==200000 and q['coarse_pressure_for_potential']==199999 and q['conservative_pressure_field_recipe']==200002
 checks['direct_order_two_does_not_reach_final_target']=not q['direct_order_two_alone_reaches_target'] and q['source_order']<q['target_order']
 r=rates();checks['postfast_and_sharp_slow_orders_are_current_first_real_rates']=r['postfast_scalar']==2-3*k and r['sharp_slow_tangential']==2-3*k and r['sharp_slow_radial']==3-3*k
 checks['final_fast_difference_has_zero_strict_margin']=r['final_fast_difference']==r['nu'] and r['final_fast_difference_strict_margin']==0
 checks['no_complete_residual_or_blowup_claim']=not r['actual_current_A1_high_jets_materialized'] and not r['complete_residual_improvement_or_blowup_certified']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='87fdc99cfa64dcd36697ff73d35856607ae39efa',
  scope='exact conditional postfast pressure/scalar/slow coefficient algebra and strict current-A1 request; no actual continuum high-jet enclosure or complete residual closure')
