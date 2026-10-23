"""Exact null-moment stress tensor and mixed-rate radial/pressure controls."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
from budgets import HERE,load,fast,scalar,fresh,prefix,P,haar_rows,flux_remainder,rates,HaarStressJoin

old=sys.modules.get('budgets');sys.modules['budgets']=prefix.haar
try:initial=load('haar_stress_initial_controls',HERE.parent/'core_initial_haar_row_budget/controls.py')
finally:
    if old is None:sys.modules.pop('budgets',None)
    else:sys.modules['budgets']=old
means=initial.means;R,Z,T,Y=means.R,means.Z,means.T,means.U
dynamic=load('haar_stress_dynamic',HERE.parent/'core_dynamic_mean_recheck/proof_checks.py')


def fixture(trace=None):
    old=sys.modules.get('budgets');sys.modules['budgets']=fresh
    try:c=load('haar_stress_fresh_fixture',HERE.parent/'core_reweighted_fresh_mean_closure/controls.py')
    finally:
        if old is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=old
    return HaarStressJoin(c.fixture(trace))


def witness_review():
    eps=means.EPS;eta=Q(1,7)
    g=(R-1)**6*(2-R)**6;f=R**4*g*g
    plus=(means.Poly(1),R,means.Poly(-1));minus=(means.Poly(1),-R,means.Poly(1))
    tensor=lambda a,b:[[eps*f*(a*plus[i]*plus[j]+b*minus[i]*minus[j]) for j in range(3)] for i in range(3)]
    W0=tensor(Q(3,2),Q(1,2));W1=tensor(Q(3,2)+eta/2,Q(1,2)-eta/2)
    D=[[W1[i][j]-W0[i][j] for j in range(3)] for i in range(3)]
    normal=(means.Poly(),R.shift(-2),means.Poly(1))
    base=(eps*R*f,R**4*(2+Z),means.Poly());state=(means.Poly(),)*3
    stress=(W0[0][1],W0[0][2]);mass=f.integral().scalar();center=(R*f).integral().scalar()/mass
    rho=f/mass*(1+Z*(R-center))
    a=initial.full_rows(base,state,W0,stress,rho);b=initial.full_rows(base,state,W1,stress,rho)
    jt=lambda W:(R**2*W[1][2]).integral()
    jz=lambda W,gr:(R*W[2][2]).integral()-(R**2*gr).integral()/2
    et=b[3]-a[3];ez=b[4]-a[4]
    checks=dict(
        all_tensors_symmetric=all(W[i][j]==W[j][i] for W in (W0,W1) for i in range(3) for j in range(3)),
        both_covariances_transverse_to_same_normal=all(sum((W[i][j]*normal[j] for j in range(3)),means.Poly())==0 for W in (W0,W1) for i in range(3)),
        positive_rank_one_weights_give_PSD=Q(3,2)>0 and Q(1,2)>0 and Q(3,2)+eta/2>0 and Q(1,2)-eta/2>0,
        full_energy_trace_unchanged=sum((D[i][i] for i in range(3)),means.Poly())==0,
        complete_radial_source_and_pressure_unchanged=a[:3]==b[:3],
        all_three_scalar_defects_unchanged=a[1]==b[1] and jt(W0)==jt(W1) and jz(W0,a[0])==jz(W1,b[0]),
        reference_selected_stress_matched_exactly=W0[0][1]-stress[0]==W0[0][2]-stress[1]==0 and a[3]==0,
        reference_axial_full_moving_pressure_memory_nonzero=bool(a[4]),
        moving_density_has_exact_unit_mass=rho.integral()==1 and bool(rho.diff(1)),
        changed_selected_theta_flux=D[0][1]==eps*eta*R*f,
        changed_selected_axial_flux=D[0][2]==-eps*eta*f,
        nonzero_full_Haar_theta_row=bool(et) and et==means.radial(eps*eta*R*f,2),
        nonzero_full_Haar_axial_row=bool(ez) and ez==means.radial(-eps*eta*f),
        both_residual_global_moments_zero=(R**2*et).integral()==(R*ez).integral()==0,
        both_rows_are_auxiliary_independent=et.haar()==et and ez.haar()==ez,
        every_covariance_entry_smoothly_compact_to_fixture_order=all(w.at(0,1)==w.at(0,2)==0 for W in (W0,W1) for row in W for w in row),
        full_pressure_compact_at_both_ends=all(a[2].at(0,x)==b[2].at(0,x)==0 for x in (1,2)),
    )
    return dict(checks=checks,scope='exact compact PSD transverse covariance-map input, not the selected pulse/root family or its full curl; real compact curl adds a separate higher-order square')


def rate_review():
    checks={};rows=[]
    for k in (Q(1,100000),Q(1,10),Q(249,1000)):
        r=rates(k);tag=str(k)
        checks[tag+':fine_flux_minimum_is_C']=r['fine_flux']==r['C'] and all(a>=r['C'] for a in r['flux_remainder_terms'].values())
        checks[tag+':all_nonflux_rows_stronger_than_C']=r['nonflux_Haar_row']>r['C']>r['D']
        checks[tag+':radial_postslow_gain_is_B_plus_C']=r['postslow_scalar']==r['radial_Haar_compensator']==r['B']+r['C']>r['D']
        checks[tag+':normalized_principal_gain_is_positive']=r['principal_flux_required_gain']==Q(1,2)-2*k>0
        checks[tag+':actual_matching_not_promoted']=not r['actual_leading_stress_identity_certified'] and not r['full_Haar_from_coarse_or_scalar_bounds_alone']
        rows.append({a:str(r[a]) for a in ('B','C','D','principal_flux_required_gain','postslow_scalar')})
    # eps=t^10, B=9/10, C=14/10. Compare to the independent
    # unreduced dynamic scalar budget, with ALL base/old-mean terms.
    old={k:P(i+2) for i,k in enumerate(prefix.STATE)};inc={k:P(i+5) for i,k in enumerate(prefix.STATE)}
    for n in (0,1,2):
        upper=scalar.slow_reset(n,old,inc,P(2),P(3),P(8))
        for t in (Q(1,2),Q(1,8)):
            norms={k:Q(0) for k in dynamic.KEYS}
            norms.update(b=2*t**10,G=3,old_B=2*t**19,old_v=3*t**9,old_g=4*t**9,B=5*t**24,u=6*t**14,g=7*t**14)
            raw=dynamic.budget(n,t**10,1,2,norms)
            for key in ('P','Qrem','Jtheta','Jz'):
                checks[f'n{n}_t{t}:independent_mixed_postslow_'+key]=raw[key]<=t**23*upper[key].evaluate(1)
    for L in (2,8,32):
        eta=Q(1,L);a0,b0=Q(3,2),Q(1,2);a1,b1=a0+eta/2,b0-eta/2
        checks[f'L{L}:positive_close_weights_preserve_energy_not_flux']=a1>0 and b1>0 and a1+b1==a0+b0 and (a1-b1)-(a0-b0)==eta
    checks['algebraic_closeness_cannot_supply_positive_epsilon_power']=rates()['principal_flux_required_gain']>0
    return dict(checks=checks,rows=rows,scope='exact mixed-rate inequalities and rational sequence controls; exponential versus polynomial non-implication is analytical')


@lru_cache(None)
def composition_review():
    trace=set();h=fixture(trace);s=h.postslow_scalars(0);rad=h.radial_Haar_compensator(0)
    envelope=h.current_row_envelope(0,principal_defect=P(7),fine_flux=P(11))
    checks={}
    for key in ('P','Jtheta','Jz'):
        checks['fresh_postslow_complete_'+key]=s[key].evaluate(1)>0
    checks.update(
        radial_Haar_obstruction_uses_fresh_postslow_P=rad.evaluate(1)==(h.means.fixed(0)['rho']*s['P']).evaluate(1),
        whole_pressure_product_not_frozen=h.whole_axial_moment(0,P(5)).c==(s['Jz']+5*s['P']).c,
        current_Haar_pressure_full_and_not_order_promoted=envelope['pm'].evaluate(1)>0 and envelope['leading_order']==1,
        unmatched_leading_stress_retained_in_both_rows=envelope['leading_theta'].evaluate(1)>0 and envelope['leading_z'].evaluate(1)>0,
        fine_row_uses_all_current_means_and_wave=envelope['fine_theta'].evaluate(1)>0 and envelope['fine_z'].evaluate(1)>0,
        no_graph_or_fast_time_coefficient_in_Haar_grammar=True,
        common_higher_pole_recipes_unchanged=h.means.reweighting['covered_operator_poles']==h.scalars.means.reweighting['covered_operator_poles'],
    )
    fixed=h.means.fixed(3);state=h.means.final(3);w=h.means.wave(3)['U']
    altered=dict(fixed,graph=999*fixed['graph'],c=777*fixed['c'])
    a=haar_rows(0,state,fixed,w*w,P(3),P(5));b=haar_rows(0,state,altered,w*w,P(3),P(5))
    checks['no_graph_or_fast_time_coefficient_in_Haar_grammar']=all(a[k].c==b[k].c for k in ('gr','pm','fine_theta','fine_z'))
    checks['all_complete_covariance_differences_accumulate']=flux_remainder(P(2),P(3),P(5),P(7)).evaluate(1)==17
    try:h.certify_actual_principal_matching()
    except NotImplementedError:checks['actual_principal_matching_cannot_be_asserted_by_a_flag']=True
    else:checks['actual_principal_matching_cannot_be_asserted_by_a_flag']=False
    return dict(checks=checks,highest_original_initial_jet=max(n for name,n in trace if name=='initial'),
                radial_degree=max(rad.c),scope='existing actual-map bound recipe composition and an explicit two-rate Haar envelope, not actual principal tensor identification')


@lru_cache(None)
def review():
    parts=dict(null_moment_stress=witness_review(),mixed_rates=rate_review(),composition=composition_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
