"""Exact algebra and manufactured-budget checks; no fluid norm sampling."""
from fractions import Fraction as Q
from functools import lru_cache
import sympy as s
from bounds import (normal,tangent,transfer,P,relative_error_floor,certificate,
    forced_from_reference,coefficientwise_leq,certificate_scope,load)

@lru_cache(None)
def fixture():
    p=normal.frame.inputs(bmin=Q(1,20),bmax=1,Fmax=1,gmax=2,gmin=1,u=Q(1,100),
        cmin=1,cmax=1,Dbox=6,Mz=1)
    sp=normal.selector.fixture();L=normal.normal_length_certificate(p)['length_rule_floor']
    sp['length']=Q(-(-L.numerator//L.denominator))
    ell=normal.selector.select_tail(sp)['ell_tail'] # MANUFACTURED control only
    kw=dict(F=lambda n:P(n+2),G=lambda n:P(n+3),btilde=lambda n:P(n+1),
        path=lambda n:P({1:1}) if n==0 else P(int(n==1)),source=lambda n:P(n+1))
    return p,sp,ell,kw

def symbolic_review():
    B,b,t=s.symbols('B b theta',positive=True)
    checks=dict(relative_triangle_slack_is_theta_times_reference_gap=(
        s.expand(B-t*b-(1-t)*B-t*(B-b))==0),
        original_Riccati_factor_leaves_one_over_32=(
            s.cancel(s.Symbol('Dt',positive=True)/(16*(1+s.Symbol('sr',positive=True))*2*s.Symbol('Dt',positive=True)/b*b))
            ==1/(32*(1+s.Symbol('sr',positive=True)))),
        universal_inverse_factor_below_four=Q(1024,961)<4,
        original_rounding_upper_gap_factorization=s.expand(4-(1+t)**2)==s.expand((1-t)*(3+t)))
    g,k,m=s.symbols('g k m',nonzero=True,real=True)
    n=s.Matrix([0,0,1]);f=g*n;Pi=s.eye(3)-n*n.T
    checks['pure_normal_source_velocity_projection_zero']=Pi*f==s.zeros(3,1)
    checks['pure_normal_source_pressure_not_zero']=s.simplify((n.dot(f))/(k*m)-g/(k*m))==0 and g/(k*m)!=0
    # Both moving-normal and full source terms in the pressure numerator.
    a,r,v,ff=s.symbols('nK nprime velocity nf')
    full=(a-r)*v+ff
    checks['moving_normal_pressure_cannot_be_dropped']=s.expand(full-(a*v+ff))==-r*v
    checks['projecting_pressure_source_erases_nonzero_term']=s.expand(full-(a-r)*v)==ff
    return checks,dict(scope='universal algebra and exact pure-normal pressure control')

def relative_review():
    checks={};rows=[]
    K=(Q(3,5),Q(4,5));Kperp=(-K[1],K[0])
    for i,(b,B,theta) in enumerate(((Q(1,100),Q(1,2),Q(1,2)),
        (Q(1,20),Q(1),Q(1,32)),(Q(1),Q(1),Q(1,32)))):
        e=theta*b;r=relative_error_floor(b,B,e)
        checks[f'case{i}_reference_relative_triangle_preserved']=r['comparison_margin']>=0 and r['relative_lower']==(1-theta)*B
        for j,direction in enumerate((K,Kperp,(-K[0],-K[1]))):
            nt=tuple(B*ki+e*di for ki,di in zip(K,direction))
            checks[f'case{i}_vector{j}_exact_squared_normal_floor']=sum(x*x for x in nt)>=(1-theta)**2*B*B
        rows.append(r)
    r=relative_error_floor(1,1,Q(1,32))
    checks['relative_31_over_32_bound_is_attained_by_antiparallel_error']=r['triangle_lower']==r['relative_lower']==Q(31,32)
    b,B,beta=Q(1,100),Q(1,2),Q(1,10);fake_n=b/2
    checks['coarse_uniform_floor_and_calibration_alone_do_not_prove_phase_floor']=(fake_n>=b/2 and B*B>=beta and fake_n*fake_n<beta/4)
    checks['coarse_countermodel_violates_retained_full_error_bound']=B-fake_n>b/2
    # Rational square-root charts: 1+u^2=w^2.
    for i,(u,w) in enumerate(((Q(3,4),Q(5,4)),(Q(4,3),Q(5,3)))):
        for a in (Q(1),Q(2),Q(4)):
            lo=Q(1,8);lam=3*lo;Bs2=lam/(a*w**3);beta=lo/(4*(1+u*u)**2)
            checks[f'calibration{i}_a{a}_original_three_halves_dominates_phase_beta']=w*w==1+u*u and Bs2>=beta
    u,w=Q(3,4),Q(5,4);Bs2=Q(1)/(4*w**3)
    checks['calibration_factor_four_cannot_be_dropped']=Bs2<1/(1+u*u)**2
    for i,t in enumerate((Q(1),Q(2,3),Q(1,7))):
        kinv=1/t;k=-(-kinv.numerator//kinv.denominator)
        checks[f'rounding{i}_same_ceiling_k_retains_both_calibration_bounds']=1<=t*t*k*k<=4
    return checks,dict(rows=rows,scope='manufactured rational vectors and calibration controls, not actual profiles')

def original_route_review():
    p,sp,ell,_=fixture();checks={};rows=[]
    for j,band in enumerate((ell,ell+1,ell+10000)):
        c=certificate(p,sp,band);theta=c['relative_error_fraction'];beta=c['phase_beta']
        checks[f'band{j}_original_Riccati_slack_closes_relative_reference_floor']=all(c['gates'].values()) and theta<=c['relative_error_cap']<=Q(1,32)
        checks[f'band{j}_requested_phase_inverse_needs_no_extra_native_margin']=c['reference_scaled_inverse']<4/beta and not c['extra_native_margin_gate_required_for_this_original_frame_route']
        checks[f'band{j}_coarse_bmin_square_can_be_below_calibrated_beta']=p['bmin']**2<beta
        checks[f'band{j}_existing_band_length_and_primitives_retained']=c['existing_ell']==band and c['existing_Lreq']==sp['length'] and not c['bands_lengths_norm_bounds_roots_or_fields_selected_or_changed']
        checks[f'band{j}_uniform_profile_constant_not_S_power']=theta==certificate(p,sp,ell)['relative_error_fraction']
        rows.append(dict(theta=theta,beta=beta,old_q=c['old_inverse_normal_square'],new_q=c['inverse_normal_square'],reference_q=c['reference_scaled_inverse']))
    for name,pp,ss,ee in (
        ('unknown_early_band',p,sp,16),
        ('length_semantics_lost',p,dict(sp,length=Q(1)),ell),
        ('new_Mz_cannot_replace_old_norm',dict(p,Mz=Q(10**30)),sp,ell)):
        try:certificate(pp,ss,ee)
        except ValueError:checks[name+'_rejected_without_reselection']=True
        else:checks[name+'_rejected_without_reselection']=False
    c=certificate(p,sp,ell)
    checks['calibrated_bound_not_inferred_from_coarse_uniform_floor_alone']=not c['phase_beta_bound_inferred_from_coarse_floor_alone']
    checks['actual_primitive_and_historical_band_claim_not_promoted']=not c['actual_primitive_enclosures_or_unknown_historical_band_certified']
    good_old=certificate(dict(p,bmin=Q(1)),sp,ell)
    checks['minimum_retains_old_bound_when_it_is_better']=good_old['inverse_normal_square']==good_old['old_inverse_normal_square']<good_old['reference_scaled_inverse']
    # A later conservative native envelope need not replace the ORIGINAL Mz
    # or invalidate a length already defined from valid original bounds.
    native=load('reference_floor_optional_native','core_leading_direction_enclosure/bounds.py')
    enc=native.leading_enclosures(Q(1,200),1,Q(1,4),2,10**12)
    comparison=native.comparison(enc,Q(1,200),Rmin=p['Rmin'],Rmax=p['Rmax'],
        tilt_upper=p['u'],lambda_max=p['lam_max'],gmin=p['gmin'],
        Llow=2*sp['radius'],Lhigh=12*sp['radius'],Dbox=p['Dbox'])
    optional=native.selector_margin(comparison,p['lam_min'],sp['length'],sp['radius'],sp['slow_cost'],sp['slow_tolerance'])
    S=Q(ell**2);cc=comparison['coefficients']
    same_band_error=cc['A']/S+cc['B']/S**2+cc['C']/S**3
    checks['optional_native_comparison_can_fail_while_original_reference_route_passes']=(
        not optional['normal_margin'] and 4*same_band_error**2>optional['beta'] and
        all(c['gates'].values()) and c['existing_Lreq']==optional['supplied_Lreq'])
    checks['optional_failure_does_not_replace_original_Mz_Db_or_length']=p['Mz']==1 and p['Db']==Q(1,10**6) and c['existing_Lreq']==sp['length']
    return checks,dict(rows=rows,optional_native_control=dict(error=optional['retained_power_error'],
        same_existing_band_error=same_band_error,beta=optional['beta'],margin_passes=optional['normal_margin']),
        scope='same-selector manufactured budget controls, not a fluid counterexample; analytical conclusion uses valid ORIGINAL real primitive bounds')

def inverse_review():
    p,sp,ell,kw=fixture();phase,new,old,c,tc=forced_from_reference(p,sp,ell,**kw)
    checks=dict(calibrated_reciprocal_is_strictly_smaller_in_this_coarse_floor_control=new.inverse_normal_square<old.inverse_normal_square,
        original_tangent_and_length_not_recomputed_from_new_reciprocal=new.tangent is old.tangent and new.length is old.length,
        full_source_callback_retained=new.source is old.source is kw['source'],
        all_full_normal_matrix_and_profile_callbacks_retained=all(getattr(new,k) is getattr(old,k) for k in ('normal','moving_normal','matrix','profile_rate')),
        original_tangent_certificate_retains_all_gates=all(tc['gates'].values()),
        old_majorant_not_mistaken_for_an_actual_field=not c['actual_primitive_enclosures_or_unknown_historical_band_certified'])
    cut=lambda n:P(n+1)
    kinds=('reciprocal','lift','undamped','normal_projector','generator','projected_source','velocity','pressure','potential')
    for j in range(3):
        for kind in kinds:
            a,b=getattr(new,kind)(j),getattr(old,kind)(j)
            checks[f'order{j}_{kind}_coefficientwise_monotone_at_same_other_inputs']=coefficientwise_leq(a,b)
        for kind in ('velocity','pressure','potential'):
            checks[f'order{j}_{kind}_same_cutoff_preserves_monotonicity']=coefficientwise_leq(new.cut(kind,cut,j),old.cut(kind,cut,j))
    checks['ambient_majorant_also_monotone']=coefficientwise_leq(new.ambient(),old.ambient())
    checks['new_reciprocal_strictly_improves_zero_order_pressure_majorant']=new.pressure(0).evaluate(1)<old.pressure(0).evaluate(1)
    return checks,dict(new_q=c['inverse_normal_square'],old_q=c['old_inverse_normal_square'],
        tangent_integer_bits=old.tangent.c[0].numerator.bit_length(),
        scope='exact complete recurrence controls through order2; all-order monotonicity follows from positivity induction, actual sources unverified')

@lru_cache(None)
def review():
    checks={};parts={}
    for name,func in (('algebra',symbolic_review),('relative',relative_review),('original_route',original_route_review),('inverse',inverse_review)):
        c,d=func();checks.update({name+':'+k:bool(v) for k,v in c.items()});parts[name]=d
    return dict(checks=checks,parts=parts,certificate=certificate_scope())
