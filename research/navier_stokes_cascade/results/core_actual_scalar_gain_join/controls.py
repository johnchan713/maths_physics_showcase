"""Independent full-current identities, exact subunit rates and recipe composition."""
from fractions import Fraction as Q
from functools import lru_cache
import sys
from budgets import HERE,load,fresh,prefix,P,SCALARS,slow_reset,fast_change,rates,ScalarGainBudget

coupled=load('scalar_gain_coupled',HERE.parent/'core_coupled_defect_budget/proof_checks.py')
dynamic=coupled.dynamic


def fixture(trace=None):
    # Parent controls import their own budgets by this historical short name.
    previous=sys.modules.get('budgets');sys.modules['budgets']=fresh
    try:parent=load('scalar_gain_parent_controls',HERE.parent/'core_reweighted_fresh_mean_closure/controls.py')
    finally:
        if previous is None:sys.modules.pop('budgets',None)
        else:sys.modules['budgets']=previous
    return ScalarGainBudget(parent.fixture(trace))


def rate_review():
    checks={};rows=[]
    # eps=t^4, B=3/4: every power here is rational and exact.
    B=Q(3,4);n=2;M=P(8)
    old={k:P(i+2) for i,k in enumerate(prefix.STATE)}
    inc={k:P(i+5) for i,k in enumerate(prefix.STATE)}
    bound=slow_reset(n,old,inc,P(2),P(3),M)
    for t in (Q(1,2),Q(1,8),Q(1,64)):
        eps=t**4;norms={k:Q(0) for k in dynamic.KEYS}
        norms.update(b=2*eps,G=3,old_B=2*t**7,old_v=3*t**3,old_g=4*t**3,
                     B=5*t**7,u=6*t**3,g=7*t**3)
        raw=dynamic.budget(n,eps,Q(1),Q(2),norms)
        for k in ('P','Qrem','Jtheta','Jz'):
            checks[str(t)+':independent_general_dynamic_'+k]=raw[k]<=t**6*bound[k].evaluate(1)
        # Positive swirl-square moment cannot satisfy the invalid B+1 order.
        ratio=t**6/t**7
        checks[str(t)+':subunit_square_loses_B_plus_one']=ratio==1/t
        rows.append(dict(t=str(t),epsilon=str(eps),quadratic_over_false_B_plus_one=str(ratio)))
    for k in (Q(1,100000),Q(1,8),Q(249,1000)):
        r=rates(k);tag=str(k)
        checks[tag+':coarse_scalar_margin_strict_on_original_range']=r['scalar_margin']==Q(1,2)-k>0
        checks[tag+':current_scalars_reach_C_without_fast_D']=r['postfast']==r['C'] and r['fresh_slow']==r['C']
        checks[tag+':fast_feedback_still_blocks_D']=r['complete_feedback']==r['B']<r['D']
    # Every discarded epsilon exponent in the slow reset is nonnegative for
    # the whole interval 0<B<=2; this is exact affine endpoint arithmetic.
    for name,a,b in (('time_and_base',Q(2),Q(-1)),('square',Q(0),Q(0)),
                     ('radial_square',Q(2),Q(0))):
        checks['whole_interval_'+name]=min(a,a+2*b)>=0
    return dict(checks=checks,rows=rows,scope='exact rate arithmetic and independent unreduced budget; not continuum input enclosures')


def identity_review():
    means=dynamic.means;R,Z,T,U=means.R,means.Z,means.T,means.U
    f=(R-1)**6*(2-R)**6
    u=(R**4*f).diff(0).shift(-2)*(1+Z)
    psi=R*f*(1+Z+Z**2);g=means.radial(psi)
    G=R*f*(1+Z)
    leak_t=(R**2*G*u).integral();leak_z=(2*R*G*g).integral()
    cr=1+Z+Z**2;p=1+2*Z+Z**3;jz=2+Z**2
    checks=dict(
        support_contract_violation_leaks_linear_angular_scalar=bool(leak_t),
        support_contract_violation_leaks_linear_axial_scalar=bool(leak_z),
        witness_both_reserved_velocity_moments_zero=(R**2*u).integral()==(R*g).integral()==0,
        whole_moving_pressure_product_derivative=(jz+cr*p).diff(1)==jz.diff(1)+cr*p.diff(1)+cr.diff(1)*p,
        frozen_compensator_derivative_is_detectably_wrong=bool(cr.diff(1)*p),
        compact_swirl_square_pressure_strictly_positive=(u*u).shift(-1).integral().at(1,0).scalar()>0,
    )
    # Compare the new full three-component fast formula to independently
    # assembled physical products; Haar linear cancellation alone is not
    # enough to delete old/new correlations.
    old={k:P(i+2) for i,k in enumerate(prefix.STATE)};inc={k:P(i+5) for i,k in enumerate(prefix.STATE)}
    a=fast_change(1,old,inc,P(8))
    sf=2*2*5+5*5;hf=2*7+4*5+5*7;tf=2*3*6+6*6;zf=2*4*7+7*7;kf=4*6+3*7+6*7
    expected=dict(P=8*(sf+tf+2*hf),Qrem=8*(sf+tf+2*hf),Jtheta=8*kf,Jz=8*zf+4*(sf+tf+2*hf))
    checks.update({'independent_asymmetric_fast_'+k:a[k].evaluate(1)==v for k,v in expected.items()})
    return dict(checks=checks,scope='exact compact polynomial controls, including a linear support-leak witness; not selected profiles')


@lru_cache(None)
def composed_review():
    trace=set();s=fixture(trace);n=0
    old=s.old_postslow(n);wave=s.sourced_change(n);fast=s.fresh_fast_change(n);now=s.postfast(n)
    slow=s.slow_increment(n);an=s.postfast(1);k=s.means.slow(0);kn=s.means.slow(1)
    checks={}
    for key in SCALARS:
        checks['actual_old_reset_wave_fast_all_retained_'+key]=now[key].evaluate(1)==(old[key]+wave[key]+fast[key]).evaluate(1)
        checks['omitting_old_reset_changes_'+key]=old[key].evaluate(1)>0
        checks['omitting_complete_covariance_changes_'+key]=wave[key].evaluate(1)>0
        checks['omitting_actual_fast_correlations_changes_'+key]=fast[key].evaluate(1)>0
    expected=dict(u=k['H']*(k['Kp']*now['P']+k['Kz']*now['Jz']),z=k['H']*k['Kg']*now['Jtheta'],
                  x=kn['H']*kn['Kpsi']*an['Jtheta'])
    checks['fresh_slow_uses_derived_current_scalars_and_extra_angular_jet']=all(slow[key].c==expected[key].c for key in prefix.STATE)
    checks['whole_axial_moment_keeps_full_current_pressure']=s.whole_axial_moment(0,P(7)).c==(now['Jz']+7*now['P']).c
    checks['same_reweighted_original_operator_recipes']=s.means.reweighting['covered_operator_poles']==fresh.shift_recipe(1,4,Q(1,2),1)['covered_operator_poles']
    checks['long_shell_lower_pole_embedding_is_retained']=s.embedding==4**12>1
    checks['same_complete_fast_A1_field_budget']=s.means.KA1(0).evaluate(1)>0 and s.means.fast(0)['z'].evaluate(1)>0
    checks['original_actual_input_enclosures_not_promoted']=not s.generated.certificate['actual_primitive_enclosures_certified']
    try:s.unsupported_gain(0)
    except NotImplementedError:checks['improved_full_source_interface_still_blocked']=True
    else:checks['improved_full_source_interface_still_blocked']=False
    return dict(checks=checks,highest_original_initial_jet=max(n for name,n in trace if name=='initial'),
                postfast_degrees={key:max(now[key].c) for key in SCALARS},
                scope='positive-polynomial composition of the existing actual-map bound recipes on manufactured original inputs')


@lru_cache(None)
def review():
    parts=dict(rates=rate_review(),identities=identity_review(),composition=composed_review())
    return dict(checks={name+':'+key:value for name,part in parts.items() for key,value in part['checks'].items()},parts=parts)
