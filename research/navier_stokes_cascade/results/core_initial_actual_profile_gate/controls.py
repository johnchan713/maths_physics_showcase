from fractions import Fraction as Q
from profile_gate import (actual_profile_implication,exponents,
                          fourier_injectivity,mode_divisors,rate_gate)


def review():
    checks={}
    modes={(1,0):Q(3),(0,1):Q(-2),(-1,0):Q(3),(0,-1):Q(-2),(2,-1):Q(1),(-2,1):Q(1)}
    f=fourier_injectivity(modes)
    checks['every_nonzero_integer_mode_has_nonzero_radial_divisor']=f['all_radial_divisors_nonzero']
    checks['every_nonzero_integer_mode_has_nonzero_temporal_divisor']=f['all_temporal_divisors_nonzero']
    checks['L_squared_Nt_inverse_L_is_injective_on_zero_Haar_modes']=f['transfer_is_injective'] and f['input_modes']==f['output_modes']
    for k in ((1,0),(0,1),(13,-8),(239,-169)):
        d=mode_divisors(*k)
        checks['Pell_norm_nonzero_'+str(k)]=d['pell_norm']!=0 and d['source_to_viscosity_multiplier_nonzero']
    a=actual_profile_implication()
    checks['auxiliary_independent_signed_multiplier_does_not_make_physical_profile_constant']=a['multiplier_auxiliary_independent'] and not a['physical_profile_auxiliary_independent']
    checks['nonzero_selected_Haar_target_and_strict_transverse_support_force_nonconstant_profile']=a['profile_nonzero'] and a['profile_nonconstant']
    checks['nonconstant_actual_profile_excludes_internal_graph_and_viscosity_cancellation']=a['graph_source_nonzero'] and a['exact_profile_cancellation_excluded']
    zero=actual_profile_implication(selected_target_nonzero=False)
    checks['zero_target_branch_is_not_misreported_as_nonzero']=not zero['profile_nonzero'] and not zero['exact_profile_cancellation_excluded']
    r=rate_gate();o=r['orders'];k=o['kappa']
    checks['source_and_viscosity_orders_match_inherited_D_and_nu']=o['D']==2-4*k and o['nu']==3-6*k
    checks['exact_nonvanishing_is_not_promoted_to_quantitative_epsilon_lower_bound']=not a['quantitative_epsilon_gain_proved'] and not r['strict_gain_from_current_evidence']
    checks['remaining_strict_gate_is_explicit_positive_order_inequality']='delta>0' in r['strict_sufficient_inequality']
    return dict(checks={name:bool(value) for name,value in checks.items()},fixture=f,actual_profile=a,rate_gate=r,
      parent_commit='c7519993b5e45ec4e6f3fe08ebfdc32280d91b9b',
      scope=('conditional actual-profile Fourier injectivity and exact remaining rate gate; '
             'not an actual epsilon lower bound, high-jet enclosure, complete residual estimate or blow-up proof'))
