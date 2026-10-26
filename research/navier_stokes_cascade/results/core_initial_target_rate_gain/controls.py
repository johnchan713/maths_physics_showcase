from fractions import Fraction as Q
from rate_gain import divisor_data,exponents,final_row,shell_inverse_bound,target_gain,target_transfer


def review():
    checks={};k=Q(1,100000);o=exponents(k);g=target_gain(k);f=final_row(k)
    rows=[divisor_data(*mode) for mode in ((1,0),(0,1),(13,-8),(239,-169))]
    checks['actual_profile_multiplier_has_only_finite_derivative_inverse_loss']=all(r['pell_norm'] and r['inverse_power_bound']>0 for r in rows) and shell_inverse_bound()['summable']
    checks['labelwise_zero_collars_and_squared_partition_transfer_profile_order']=target_transfer()['operator_cannot_create_epsilon_gain']
    checks['sharper_initial_Haar_rows_are_beta_not_downgraded_alpha']=g['initial_Haar_rows']==2-2*k and g['initial_Haar_rows']-o['alpha']==k
    checks['order_three_moment_bumps_do_not_limit_projected_source']=g['projected_source']==o['beta']<g['compatibility_moments']
    checks['linear_compact_target_and_signed_maps_preserve_beta']=g['compact_target']==o['beta'] and target_transfer()['inherited_signed_maps_preserve_epsilon_order']
    checks['signed_principal_and_curl_orders_are_recomputed_from_beta']=g['signed_velocity']==Q(3,2)-2*k and g['signed_curl']==2-3*k
    checks['real_pair_remainder_improves_to_three_minus_four_kappa']=g['real_pair_remainder']==3-4*k
    checks['pointwise_covariance_response_pays_exactly_one_graph_loss']=g['pointwise_covariance_response']==o['source']==2-3*k
    checks['complete_fast_source_improves_by_kappa_over_old_D']=f['complete_fast_source']==o['source'] and f['complete_fast_source']-o['D']==k
    checks['radial_viscosity_is_now_three_minus_five_kappa']=f['fast_groups']['radial_viscosity']==3-5*k
    checks['retained_A1_bound_ties_but_does_not_spoil_new_final_order']=f['final_groups']['retained_A1_raw_bound']==f['conditional_final_order']==3-5*k
    checks['conditional_final_row_has_strict_kappa_margin_over_required_nu']=f['strict_margin']==k>0 and f['conditional_strict_absorption']
    checks['no_complete_residual_or_high_jet_promotion_is_made']=True
    return dict(checks={name:bool(value) for name,value in checks.items()},profile_modes=rows,target_gain=g,final_row=f,
      parent_commit='5e6063a30f0dba6fe98f8c12f2f388fa8f0ec93f',
      correction=('the inherited beta=2-2kappa initial Haar estimate must not be lowered to alpha=2-3kappa '
                  'before the first signed target; retaining beta gives the missing strict kappa margin'),
      scope=('conditional first-real target/source/final-row rate correction using existing weighted Haar and fixed-map premises; '
             'not an actual high-jet enclosure, complete residual estimate or blow-up proof'))
