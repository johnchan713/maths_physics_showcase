"""Exact original-length implications; all numerical data are manufactured."""
from fractions import Fraction as Q
from functools import lru_cache
from budgets import HERE,load,frame,selector,normal_length_certificate,tangential_error,verify_existing


def fixture():
    p=frame.inputs(Dbox=Q(6),Mz=Q(1));c=normal_length_certificate(p)
    sp=selector.fixture();L=c['length_rule_floor'];sp['length']=Q(-(-L.numerator//L.denominator))
    old=selector.select_tail(sp) # MANUFACTURED fixture only; verifier never selects.
    return dict(frame=p,selector=sp,selected=old,checked=verify_existing(p,sp,old['ell_tail']))


def length_review():
    checks={};rows=[]
    for f,g,bmin,u,Db in ((Q(1,10),Q(1,5),Q(1,100),Q(1,10000),Q(0)),
                         (Q(2),Q(3),Q(1,2),Q(3,4),Q(1,10**6)),
                         (Q(100),Q(200),Q(1,10),Q(1000),Q(10**12))):
        p=frame.inputs(Fmax=f,gmax=g,gmin=g/2,bmin=bmin,u=u,Db=Db,Dbox=6,
                       lam_min=min(f,g)/2,lam_max=max(Q(2),min(f,g)))
        c=normal_length_certificate(p);b=c['original_budgets']
        name=f'f{f}_g{g}_u{u}'
        checks[name+'_physical_consistency_and_uniform_length_implication']=all(c['gates'].values())
        checks[name+'_CE_is_primitive_positive_budget_before_any_ODE']=b['CE']>=b['DH']>=2*p['lam_min']*b['Da']
        checks[name+'_both_strict_domination_factors_retained']=c['proof_margins']['normal_length']>=0 and c['proof_margins']['radial_tilt_length']>=0
        checks[name+'_original_normal_threshold_dominated']=c['length_rule_floor']>=max(1,b['Da'],b['Ds'])
        rows.append(dict(case=name,normal_threshold=str(c['normal_threshold']),original_Riccati_length=str(b['Riccati_acceptance']),proof_margins={k:str(v) for k,v in c['proof_margins'].items()}))
    impossible=frame.inputs(lam_min=10**6,lam_max=10**6)
    checks['positive_unphysical_reference_inputs_are_not_promoted']=not normal_length_certificate(impossible)['gates']['physical_reference_consistency']
    return dict(checks=checks,rows=rows)


def growth_review():
    leading=load('normal_length_reference',HERE.parent/'leading_frame_foundation/frame_floors.py')
    checks={};count=0
    for f in (Q(1,10),Q(1),Q(3)):
        for a,b in ((Q(3),Q(0)),(Q(1),Q(4)),(Q(2),Q(3))):
            r=leading.leading_identities(f,a,b);x=r['lambda2']
            if x<=0:continue
            G=f*(a+abs(b));lower=x/(1+x)
            checks[f'f{f}_a{a}_b{b}_actual_reference_identity_implies_consistency']=x<=2*f*G and lower*lower<=x and lower*lower<=2*f*G
            count+=1
    return dict(checks=checks,cases=count,scope='exact reference algebra controls, not selected-profile measurements')


def error_review():
    p=frame.inputs(Dbox=6,Mz=1);checks={}
    for S,L in ((Q(1),Q(1)),(Q(4),Q(1)),(Q(4),Q(4)),(Q(100),Q(25))):
        for name,(e,e2,ki) in dict(axial=(1/S**2,0,0),rounding=(0,0,1/S**2),squared=(0,1/S**2,0),mixed=(Q(1,3)/S**2,Q(1,3)/S**2,Q(1,3)/S**2)).items():
            r=tangential_error(p,S,L,e,e2,ki)
            checks[f'S{S}_L{L}_{name}_full_error_bounded_by_original_Dt_over_L']=r['original_slow_gate'] and r['error_upper']<=r['Dt_over_L']
    S=Q(100);L=Q(25);e=Q(1,16)/S**8;k=4*S**4
    r=tangential_error(p,S,L,e,e*e,1/k)
    checks['consistent_epsilon_square_and_ceiling_carrier_relation']=e*k*k==1 and r['original_slow_gate'] and r['implication_valid']
    base=tangential_error(p,S,L,0,0,0)
    r=tangential_error(p,S,L,Q(1,1000000),0,0)
    checks['complete_axial_transport_keeps_full_Hz_and_actual_L']=r['error_upper']-base['error_upper']==L*frame.primitive(p)['Zc']/1000000>0
    checks['both_tangential_tilt_terms_and_box_displacement_survive']=all(base['terms'][k]>0 for k in ('angular_box','angular_tilt','axial_tilt'))
    checks['slow_failure_is_reported']=not tangential_error(p,S,L,1,0,0)['original_slow_gate']
    return dict(checks=checks)


def selector_review():
    f=fixture();p,sp,s=f['frame'],f['selector'],f['selected'];checks={};rows=[]
    for ell in (s['ell_tail'],s['ell_tail']+1,s['ell_tail']+10000):
        r=verify_existing(p,sp,ell)
        checks[f'ell{ell}_same_original_selector_and_length_close_normal_nonvanishing']=r['conditional_normal_floor_implication_passes']
        checks[f'ell{ell}_actual_radius_factor_and_supplied_band_retained']=r['actual_pulse_length_lower']==2*sp['radius']*ell**2 and r['existing_ell']==ell and not r['bands_lengths_supports_or_representatives_selected_or_changed']
        rows.append(dict(existing_ell=str(ell),pulse_length_lower=str(r['actual_pulse_length_lower']),length_budget=str(sp['length']),inverse_normal_square=str(r['inverse_normal_square'])))
    bad=dict(sp,length=Q(1));r=verify_existing(p,bad,s['ell_tail'])
    checks['old_length_without_its_declared_semantics_is_rejected_not_reselected']=not r['gates']['original_length_has_its_declared_Riccati_and_25_meaning'] and r['existing_length_budget']==1 and not r['conditional_normal_floor_implication_passes']
    r=verify_existing(p,sp,16)
    checks['unknown_early_band_is_not_identified_with_accepted_tail']=not r['conditional_normal_floor_implication_passes'] and r['existing_ell']==16
    larger=dict(p,Mz=p['Mz']*10**12);r=verify_existing(larger,sp,s['ell_tail'])
    checks['new_conservative_Mz_cannot_be_silently_substituted_into_old_length']=not r['gates']['original_length_has_its_declared_Riccati_and_25_meaning'] and r['existing_length_budget']==sp['length']
    distinct=dict(sp,dbox=2*sp['dbox']);r=verify_existing(p,distinct,s['ell_tail'])
    checks['independently_valid_box_enclosures_need_no_artificial_domination_gate']=r['conditional_normal_floor_implication_passes'] and p['Dbox']<distinct['dbox']
    return dict(checks=checks,rows=rows,scope='rule-consistent manufactured selector data; actual primitive enclosures and unknown historical-band membership not certified')


def denominator_review():
    p=frame.inputs(bmin=Q(1,100));c=normal_length_certificate(p)
    phase_beta=p['lam_min']/(4*(1+p['u']**2)**2)
    checks={
        'original_inverse_bound_uses_Bs_floor_squared':c['inverse_normal_square']==4/p['bmin']**2,
        'frame_linear_beta_is_not_phase_squared_beta':frame.budgets(p)['beta']==p['bmin']/2 and frame.budgets(p)['beta']!=phase_beta,
        'coarse_Bs_floor_does_not_prove_stronger_phase_beta_inverse_bound':p['bmin']**2<phase_beta and c['inverse_normal_square']>4/phase_beta,
        'no_stronger_phase_beta_or_actual_profile_promotion':not fixture()['checked']['stronger_phase_beta_bound_certified'] and not fixture()['checked']['actual_input_enclosures_or_arbitrary_historical_band_certified'],
    }
    return dict(checks=checks,original_inverse=str(c['inverse_normal_square']),stronger_optional_phase_inverse=str(4/phase_beta))


@lru_cache(None)
def review():
    parts=dict(length=length_review(),growth=growth_review(),error=error_review(),selector=selector_review(),denominators=denominator_review())
    return dict(checks={name+':'+k:v for name,p in parts.items() for k,v in p['checks'].items()},parts=parts)
