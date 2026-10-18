"""Exact regional amplitude-floor algebra; no actual giant profile is sampled."""
from fractions import Fraction as F
from itertools import product
from math import exp

def rat(x):
    if isinstance(x,(bool,float)):raise TypeError('Exact rational inputs required')
    return F(x)

def serial(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serial(v) for v in x]
    return x

def compact_coordinates(fv,a,b,p1,p2):
    fv,a,b,p1,p2=map(rat,(fv,a,b,p1,p2))
    if min(fv,a)<=0:raise ValueError('Positive swirl and angular shear required')
    t=-b/a;v=a+b*b/a;D=p1+t*p2-v;J=p2-t*p1
    return dict(t=t,v=v,D=D,J=J,
        direct_stress2=fv*fv*((p1-a)**2+(p2+b)**2),
        rotated_stress2=fv*fv*(D*D+J*J)/(1+t*t))

def compact_lower(f,kappa,R):
    f,kappa,R=map(rat,(f,kappa,R))
    if min(f,kappa)<=0 or R<1:raise ValueError('Positive floors and R>=1 required')
    return f*kappa/(4*R*R)

def common_floor(f,G):
    """Only the amplitude formula. The actual regional hypotheses remain required."""
    f,G=map(rat,(f,G))
    if f<=0 or G<1:raise ValueError('Positive swirl floor and G>=1 required')
    return f/(64*G**10)

def first_relative_bound(f,G,c,Htheta,Hz):
    """Xc=c^2; propagate the previous sufficient Euclidean sum bound."""
    f,G,c,Htheta,Hz=map(rat,(f,G,c,Htheta,Hz))
    if c<=0 or min(Htheta,Hz)<0:raise ValueError('Positive compact radius and nonnegative source bounds required')
    return (c*c*Htheta/2+c*Hz)/common_floor(f,G)

def log_floor(log_f,log_G):
    """Return log(tau) with log(64) split off, without exponentiating."""
    log_f,log_G=map(rat,(log_f,log_G))
    if log_G<0:raise ValueError('log G>=0 required')
    return dict(logarithm_before_minus_log64=log_f-10*log_G,subtracted_constant='log(64)')

def inner_hierarchy_review():
    # y_inner=tc/16 and tc/t1=1/(10sqrt(log C)).
    inner_exp=100*16**2+2;G_power=4096;common_exp=10*G_power
    # log C<=C-1; C^13>=4096*C for C>=2; log G=C^131072.
    return dict(inner_start='Xa*exp(tc/16)',inner_stop='Xa*exp(tc/2)',
        s_inner='1/(160*sqrt(log C))',flat_exponent=inner_exp-2,
        inner_stress_lower='(9/16)*f_star*C^-25602',
        hierarchy='Aenv=C^4096; log G=Aenv^32=C^131072; C>=2',
        common_floor='f_star/(64*G^10)',
        inner_over_common_lower='36*C^15358',
        exponent_difference=common_exp-inner_exp,
        checks=dict(actual_offset_gives_exact_flat_exponent=100*16**2==25600,
            polynomial_log_comparison_coefficient=2**12==4096,
            log_G_dominates_4096_log_C=131072>13,
            inner_power_dominated_by_same_G=common_exp>inner_exp,
            comparison_prefactor_correct=F(9,16)*64==36,
            fixed_interval_begins_strictly_after_flat_zero=F(1,16)>0,
            inner_interval_ends_before_protected_collar_end=F(1,2)<1))

def compact_review():
    f=F(1,5);R=F(16);kappa=F(1,10);lower=compact_lower(f,kappa,R);rows=[]
    # New local coordinate controls, not values of a finalized global profile.
    for a,b,dt,dz in product((F(3,2),F(3),F(5)),(F(-1),F(0),F(1)),
                             (F(1),F(2),F(3)),(F(-1,4),F(0),F(1,4))):
        fv=F(2,5);p1=a+dt;p2=-b+dz;r=compact_coordinates(fv,a,b,p1,p2)
        inputs=(fv>=f and a>=1/R and max(abs(x) for x in (a,b,p1,p2))<=R and r['D']>=kappa/2)
        checks=dict(input_bounds=inputs,exact_rotated_norm=r['direct_stress2']==r['rotated_stress2'],
            slope_ratio=abs(r['t'])<=R*R,
            square_root_denominator=(1+r['t']**2)<=4*R**4,
            D_only_stress_lower=r['direct_stress2']>=fv*fv*r['D']**2/(4*R**4),
            compact_floor=r['direct_stress2']>=lower*lower)
        rows.append(dict(a=str(a),b=str(b),theta_excess=str(dt),axial_excess=str(dz),checks=checks))
    # This identity is valid for every positive G, independently of test grids.
    hierarchy=[]
    for G in (F(1),F(2),F(5),F(17)):
        hierarchy.append(dict(G=str(G),checks=dict(
            middle_floor_exactly_common=compact_lower(f,G**-4,4*G**3)==common_floor(f,G),
            common_floor_below_exterior_swirl=common_floor(f,G)<=f)))
    return dict(model='rational local stress-coordinate controls, not actual fields',rows=rows,
        illustrative_compact_lower=str(lower),hierarchy_identity_controls=hierarchy,
        checks=dict(all_81_exact_coordinate_controls=all(all(r['checks'].values()) for r in rows),
            exact_middle_hierarchy_identity=all(all(r['checks'].values()) for r in hierarchy)))

def exterior_review():
    rows=[]
    for h,epsilon,L,Xratio in product((F(1,200),F(1,1000)),(F(1,10**8),F(1,10**12)),
                                     (F(49,50),F(1)),tuple(map(F,(1,2,3,5,10,100)))):
        hXstar=F(4000)*10**12;X=Xratio*hXstar/h;Q=h/F(4000)
        p=X*Q/L;amax=4+h*epsilon
        rows.append(dict(h=str(h),relative_error=str(epsilon),L=str(L),X_over_Xstar=str(Xratio),
            checks=dict(edited_shear_below_five=amax<5,finite_radius_p_floor=p>=10**12,
                angular_stress_excess_exceeds_one=p-amax>1,
                endpoint_L_range=0<L<=1)))
    return dict(model='rational use of the retained hXstar, Q and edited-shear budgets; not a new heat derivation',rows=rows,
        assumptions=['Q>=h/4000','hXstar/4000>10^12','reference a<=4','|Delta a|<=h epsilon_Q<1',
                     'same forward stress by exact functional five-moment restoration'],
        conclusion='T0theta=F(p1-a)>F>=f_star on [Xrep,Xc]',
        checks=dict(all_48_exterior_algebra_controls=all(all(r['checks'].values()) for r in rows)))

def compose_review():
    rows=[]
    for f,G,c,ht,hz in product((F(1,10),F(1)),(F(1),F(2)),tuple(map(F,(1,2,3,5))),(F(0),F(5)),(F(0),F(7))):
        tau=common_floor(f,G);AW=first_relative_bound(f,G,c,ht,hz)
        rows.append(dict(checks=dict(reciprocal_stress_cost_exact=AW*tau==c*c*ht/2+c*hz,
            substituted_bound_exact=AW==64*G**10*(c*c*ht/2+c*hz)/f)))
    return dict(model='exact algebraic substitution in the previous first-coefficient estimate',rows=rows,
        formula='A_W<=64*G^10/f_star*(Xc*Htheta/2+sqrt(Xc)*Hz)',
        checks=dict(all_64_composed_controls=all(all(r['checks'].values()) for r in rows)))

def limits_review():
    neutral=compact_coordinates(1,3,1,3,-1)
    # A positive growing shear is compatible with zero target stress if its
    # D margin is omitted. This is a local algebra control, NOT the construction.
    # The actual flat edge has |T0| -> 0; the compact interval must start later.
    flat_log_sequence=[-n*n for n in (10,20,40,80)]
    logtau=log_floor(-10000,1000);badcone=compact_coordinates(1,3,0,4,100)
    badratio=(badcone['v']-2)*badcone['J']**2/(2*badcone['D']**2)
    return dict(neutral_target=serial(neutral),flat_edge_log_stress_sequence=flat_log_sequence,
        illustrative_log_floor=serial(logtau),positive_large_stress_but_bad_cone_ratio2=str(badratio),
        excluded_domains=['Xa itself','the full open annulus with Xa as an accumulation point','Xb itself'],
        checks=dict(swirl_floor_alone_not_stress_floor=neutral['v']>2 and neutral['direct_stress2']==0,
            flat_endpoint_has_no_positive_unweighted_floor=all(b<a for a,b in zip(flat_log_sequence,flat_log_sequence[1:])),
            no_floating_underflow_interpreted_as_zero=exp(float(logtau['logarithm_before_minus_log64']))==0 and logtau['logarithm_before_minus_log64']==-20000,
            same_tilt_cone_not_inferred_from_amplitude=badcone['direct_stress2']>1 and badratio>1))

def review():
    i=inner_hierarchy_review();c=compact_review();e=exterior_review();s=compose_review();l=limits_review()
    return dict(inner_offset_and_hierarchy=i,compact_coordinate_floor=c,matched_compensated_exterior=e,
        first_coefficient_substitution=s,scope_and_zero_edge_controls=l,
        checks=dict(inner_actual_offset_and_polynomial_hierarchy=all(i['checks'].values()),
            compact_exact_amplitude_estimate=all(c['checks'].values()),
            same_profile_exterior_amplitude_estimate=all(e['checks'].values()),
            first_weighted_constant_substitution=all(s['checks'].values()),
            amplitude_and_domain_scope_controls=all(l['checks'].values())))
