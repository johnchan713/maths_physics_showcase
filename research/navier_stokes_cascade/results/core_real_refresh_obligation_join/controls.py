from bounds import rates,requests,current_postfast,slow_response,final_tangential,post,EP,Q

def review():
 checks={};r=rates();B,a,D,nu=(r[k] for k in ('B','alpha','D','nu'))
 checks['real_complete_change_is_not_assigned_fine_error_order']=a<nu and D==a-Q(1,100000)
 checks['current_postfast_scalar_keeps_signed_change']=r['postfast_scalar']==a and r['fast_scalar_change']>a
 checks['fresh_Haar_slow_reset_and_whole_pressure_orders']=r['postslow_scalar']==a+1 and r['whole_postslow_axial']==a+2
 checks['fast_difference_has_exactly_zero_strict_margin']=r['strict_fast_difference_margin']==0 and r['conditional_mean_order']==nu
 checks['slow_difference_has_positive_margin']=r['tangential_groups']['slow_difference']>nu
 for label,src,p in (('A1',D,r['A1_integrations']),('A0',B,r['A0_integrations'])):
  checks[label+'_least_strict_compact_count']=src+p*Q(1,100000)>nu and src+(p-1)*Q(1,100000)<=nu
 checks['compact_counts_differ_from_later_generic_stage']=r['A1_integrations']==99999 and r['A0_integrations']==199996
 q=requests(0,old_pole=1,target_allowance=12,retained_wave_pole=30)
 trace=q['trace'];checks['all_new_requests_traced_without_materialization']=trace==post.dependency_trace(0,99999,199996) and not trace['actual_jet_or_profile_values_evaluated']
 checks['exact_real_high_jet_obligations']=trace['A1_source_order']==200001 and trace['A0_source_order']==399995 and trace['maxima']['post_signed_row']==400005 and trace['maxima']['new_covariance']==400008 and trace['maxima']['next_signed_target_input']==400009
 checks['full_pressure_domain_exceeds_parent_low_output']=q['q']==2400066 and q['operator_domains']['Kp']==(4800132,9600264) and not q['low_order_parent_inputs_sufficient']
 for n in (0,1,7):
  z=requests(n,old_pole=2,target_allowance=24,retained_wave_pole=30)
  checks[f'order{n}_whole_new_graph_covered']=z['q']>=6*z['trace']['maxima']['new_covariance']+24+6 and z['operator_domains']['KT1']==2*z['q']
 # Marked positive norm diagnostics: no actual fields or fine rows certified.
 old={k:EP(i+2) for i,k in enumerate(('x','u','z'))};inc={k:EP(i+1) for i,k in enumerate(('x','u','z'))}
 previous={k:EP(i+1) for i,k in enumerate(('P','Jtheta','Jz'))};wave=post.scalar.wave_change(0,EP(2),EP(3),EP(3))
 result=dict(updated_scalars={k:previous[k]+wave[k] for k in previous})
 out=current_postfast(0,result,old,inc,EP(3));change=post.scalar.fast_change(0,old,inc,EP(3))
 checks['current_signed_scalar_change_added_exactly_once']=all(out[k].record()==(previous[k]+wave[k]+change[k]).record() for k in out)
 wrong=post.postfast_scalar_envelope(0,result['updated_scalars'],EP(2),EP(3),old,inc,EP(3))
 checks['generic_old_scalar_interface_would_duplicate_wave_change']=all(wrong[k].record()==(out[k]+wave[k]).record() for k in out)
 reset=current_postfast(0,result,{k:EP(0) for k in old},inc,EP(3))
 checks['full_old_fast_correlations_not_replaced_by_zero_Haar']=any(out[k].record()!=reset[k].record() for k in out)
 k={key:EP(i+1) for i,key in enumerate(('H','Kp','Kz','Kg','Kpsi'))}
 fields=slow_response(0,out,{key:2*value for key,value in out.items()},k,k)
 checks['slow_response_uses_current_pressure_and_axial_defect']=fields['u'].record()==(k['H']*(k['Kp']*out['P']+k['Kz']*out['Jz'])).record()
 checks['slow_radial_reads_next_scalar_order']=fields['x'].record()==(k['H']*k['Kpsi']*2*out['Jtheta']).record()
 fixed={key:EP(i+1) for i,key in enumerate(('b','V','G','graph','c','r','rho','Kp','Kmom'))}
 haar=dict(Etheta=EP(1),Ez=EP(2));slow={key:2*value for key,value in inc.items()}
 f=final_tangential(0,haar,EP(3),old,inc,slow,fixed)
 checks['final_rows_retain_both_viscosity_and_pressure_differences']=all(f[group][key].record() for group in ('fast_difference','slow_difference') for key in ('viscosity','pressure_gradient','pm'))
 after={key:old[key]+inc[key] for key in old}
 direct=post.fast.difference(0,after,slow,fixed)
 checks['second_difference_uses_old_plus_actual_fast']=all(f['slow_difference'][key].record()==direct[key].record() for key in ('Etheta','Ez','pm'))
 checks['whole_full_row_addition_has_no_extra_Qcirc_factor']=f['Ez'].record()==(haar['Ez']+3+f['fast_difference']['Ez']+f['slow_difference']['Ez']).record()
 checks['wave_feedback_does_not_inherit_mean_row_rate']=r['complete_wave_mean_feedback_order']==D and D<nu and not r['actual_mean_block_or_full_residual_certified']
 return dict(checks={k:bool(v) for k,v in checks.items()},parent_commit='81319e5c2341f1be017a595b347216f0a80611da',scope='exact first-real conditional rates, deferred new-layer requests and manufactured scalar/row norm algebra only')
