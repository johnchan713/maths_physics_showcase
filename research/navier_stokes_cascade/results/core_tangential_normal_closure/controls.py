"""Exact low-order sum and full tangential transport controls; no PDE samples."""
from fractions import Fraction as Q
from functools import lru_cache
from budgets import diagonal_tail,comparison,tangential_error,normal_gate,length_compatibility,head,HERE


def fixture():
    return comparison(3,1,2,2,3,1,1,2,6)


def sum_review():
    checks={};rows=[]
    for t in (Q(1,10000),Q(1,16),Q(1,2),Q(1)):
        total=diagonal_tail(t)
        for N in (1,2,7,19):
            partial=sum((t/2)**n for n in range(1,N+1))
            checks[f't{t}_N{N}_partial_plus_exact_remainder']=partial+diagonal_tail(t,start=N+1)==total
        checks[f't{t}_both_low_orders_bounded_by_t']=diagonal_tail(t,0)==total<=t
        rows.append(dict(t=str(t),majorant=str(total),ratio_to_t=str(total/t)))
    checks['unit_constant_is_sharp_for_this_majorant_at_t1']=diagonal_tail(1)==1
    checks['same_majorant_cannot_give_uniform_t_squared_rate']=all(diagonal_tail(Q(1,2**j))/Q(1,2**(2*j))>=2**(j-1) for j in (2,8,24))
    return dict(checks=checks,rows=rows,scope='exact infinite geometric majorant identity; not actual coefficient measurements or an infinite schedule certificate')


def transport_review():
    r=fixture();checks={};rows=[]
    for S in (Q(1),Q(2),Q(16),Q(1000)):
        # Saturate the slow budget in turn with every error contribution.
        for name,(e,e2,ki) in dict(axial=(1/S**2,0,0),rounding=(0,0,1/S**2),epsilon2=(0,1/S**2,0),mixed=(Q(1,3)/S**2,Q(1,3)/S**2,Q(1,3)/S**2)).items():
            g=normal_gate(r,1,S,e,e2,ki)
            checks[f'S{S}_{name}_original_slow_and_component_sum']=g['gates']['same_whole_path_slow_gate'] and g['error_upper']<=r['Ctan']/S
            rows.append(dict(S=str(S),case=name,error=str(g['error_upper']),bound=str(r['Ctan']/S)))
    S=Q(100)
    base=tangential_error(r,S,0,0,0)
    checks['full_axial_transport_is_not_dropped']=tangential_error(r,S,Q(1,100000),0,0)-base==r['groups']['axial_transport']/1000>0
    checks['nonzero_rounding_cost_is_one_not_half']=tangential_error(r,S,0,0,Q(1,100000))-base==r['groups']['angular_rounding']/100000
    checks['both_tilt_and_box_costs_survive_zero_transport']=base==r['groups']['tilt']/S+r['groups']['angular_box']/S**3>0
    checks['epsilon_squared_remains_in_the_original_slow_gate']=not normal_gate(r,1,S,0,Q(2)/S**2,0)['gates']['same_whole_path_slow_gate']
    checks['slow_failure_is_reported_not_assumed']=not normal_gate(r,1,S,Q(2)/S**2,0,0)['gates']['same_whole_path_slow_gate']
    # Consistent physical relations: epsilon2=epsilon^2 and
    # k=ceil(epsilon^(-1/2))=4*S^4, not independently selected errors.
    e=Q(1,16)/S**8;k=4*S**4
    physical=normal_gate(r,1,S,e,e*e,1/k)
    checks['consistent_epsilon_squared_and_ceiling_carrier_satisfy_slow_bound']=e*k*k==1 and physical['gates']['same_whole_path_slow_gate'] and physical['error_upper']<=r['Ctan']/S
    return dict(checks=checks,rows=rows)


def old_recipe_review():
    r=fixture();checks={};rows=[]
    for M in (Q(0),Q(1),Q(1000000)):
        for K in (Q(1),Q(4928)):
            old=head.BackgroundBudget(Q(1,200),lambda m:Q(m+2),lambda n,m:M,lambda m:K,1,2,'euler')
            o=old.comparison(U=2,lambda_max=3,gmin=1,Llow=1,Lhigh=2,Dbox=6)
            T=2*(r['carriers']['theta']+r['carriers']['axial'])*(r['chart']['F']+r['chart']['G'])
            Bcorr=old.bound(1)['head']['correction']
            checks[f'M{M}_K{K}_old_full_Cphi_strictly_larger']=o['Cphi']>r['Ctan']
            checks[f'M{M}_K{K}_old_axial_excess_exact']=o['groups']['axial_transport']-r['groups']['axial_transport']==T*(Bcorr-1)>0
            rows.append(dict(M=str(M),K=str(K),old_Cphi=str(o['Cphi']),new_Ctan=str(r['Ctan'])))
    # Parameter-space separation only: choose a supplied rational S exactly
    # on the new margin, then check both recipes on THAT same S.
    r= comparison(3,1,2,0,1,1,1,2,6)
    S=4*r['Ctan'] # beta=1/4 gives exact equality in 4*Ctan^2=beta*S^2.
    old=head.BackgroundBudget(Q(1,200),lambda m:Q(m+2),lambda n,m:Q(1),lambda m:Q(4928),1,2,'euler')
    o=old.comparison(U=0,lambda_max=1,gmin=1,Llow=1,Lhigh=2,Dbox=6)
    g=normal_gate(r,1,S,Q(1,3)/S**2,Q(1,3)/S**2,Q(1,3)/S**2)
    checks['same_supplied_scale_new_margin_passes_old_margin_fails']=all(g['gates'].values()) and 36*o['Cphi']**2>g['beta']*S*S
    checks['new_exact_margin_retains_inverse_normal_square']=g['inverse_normal_square']==16 and 4*r['Ctan']**2==g['beta']*S*S
    return dict(checks=checks,rows=rows,scope='manufactured recipe comparison, not a fluid counterexample or a newly selected actual band')


def join_review():
    checks={}
    native=head.load('tangential_native',HERE.parent/'core_native_chart_cover/budgets.py')
    c=native.cover_existing(16,1,64,1,8,Q(1,1000))
    rows=native.leading_global_rows(Q(1,200),1,16,c['native_Xhi'],7,128,13,101,303)
    lead=native.supported_reconstruct(Q(1,200),0,1,c['native_Xhi'],c['native_Rlo'],c['native_Rhi'],rows['F'][:2],rows['U'][:2],16)
    r=comparison(lead['through'][1],c['band_Rmin'],c['band_Rmax'],2,3,1,1,2,6)
    checks['existing_support_aware_leading_first_rows_feed_without_first_repair']=r['Ctan']>0 and not r['actual_normal_floor_verified']
    r0=Q(1,6057072);beta=Q(1,100);C=fixture()['Ctan']
    # U=2, lambda_min=1 => beta=1/100. L=40*r0*C is equality in (9).
    L=40*r0*C;S=L/(2*r0);out=length_compatibility(fixture(),1,r0,L)
    checks['original_radius_and_length_constants_match_direct_gate']=out['margin']==0 and beta*S*S==4*C*C
    checks['insufficient_existing_length_is_reported']=not length_compatibility(fixture(),1,r0,L/2)['sufficient_gate']
    checks['length_gate_does_not_select_scale_or_certify_actual_selector']=not out['scales_or_length_selected_or_changed'] and not out['actual_existing_selector_verified']
    checks['no_independent_first_K2_or_leading_second_input']=r['independent_inputs_removed']==['first coefficient-one tuple norm','Borel cutoff K2','leading second jet']
    checks['original_infinite_diagonal_remains_an_explicit_unverified_premise']=not r['actual_infinite_schedule_verified'] and 'every n>=1' in r['retained_premise']
    return dict(checks=checks,native_leading_first=str(lead['through'][1]),Ctan=str(r['Ctan']),scope='old native norm recipes on manufactured profile budgets; no actual chart/gap/schedule promoted')


@lru_cache(None)
def review():
    parts=dict(sum=sum_review(),transport=transport_review(),comparison=old_recipe_review(),join=join_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
