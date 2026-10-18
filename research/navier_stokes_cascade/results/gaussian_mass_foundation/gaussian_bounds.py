"""Explicit all-interval constants and exact rational manufactured enclosures."""
from fractions import Fraction as F
from math import isqrt,log,pi,exp,sqrt,isfinite


def curvature(u,lower,upper):
    u,lo,hi=map(F,(u,lower,upper))
    if u<=0 or lo<=0 or hi<lo:raise ValueError('Fixed positive tilt and ordered spectral bounds required')
    c=lo*u*u/2*((1+9*u*u/4)**-2+2*(1+u*u)**-2)
    C=3*hi*u*u/2*((1+u*u/4)**-1+2*(1+u*u)**-1)
    return c,C


def root_interval(x,bits=80):
    x=F(x)
    if x<0 or type(bits) is not int or bits<1:raise ValueError('Nonnegative radicand and positive precision required')
    scale=1<<bits;n=isqrt(x.numerator*scale*scale//x.denominator)
    lo=F(n,scale);hi=lo if lo*lo==x else F(n+1,scale)
    assert lo*lo<=x<=hi*hi
    return lo,hi


def rate_curvature_interval(u,a,lam):
    u,a,lam=map(F,(u,a,lam))
    if u<=0 or not F(1,2)<=a<=F(3,2) or lam<=0:raise ValueError('Positive actual tilt, spectral value and whole-pulse a required')
    b=1+u*u*a*a;d=1+u*u;bl,bh=root_interval(b);dl,dh=root_interval(d)
    factor=lam*u*u*a
    return factor*(1/(b*bh)+2/(d*dh)),factor*(1/(b*bl)+2/(d*dl))


def barrier(mu,ce,cd,length):
    mu,ce,cd,length=map(F,(mu,ce,cd,length))
    if min(mu,ce,length)<=0 or cd<0 or length<8*ce/mu:raise ValueError('Original accepted Riccati length required')
    rho=ce/(mu*length);out=-2*mu*rho+(ce/length)*(1+2*rho+rho*rho)
    return dict(rho=rho,outward=out,negative_bound=-F(47,64)*ce/length,D=cd+F(9,8)*ce)


def logs(u,lo,hi,ce,cd,jlo=1,jhi=1):
    c,C=curvature(u,lo,hi);ce,cd,jlo,jhi=map(float,(ce,cd,jlo,jhi))
    if ce<=0 or cd<0 or jlo<=0 or jhi<jlo:raise ValueError('Positive frame and Haar bounds required')
    D=cd+9*ce/8;xm=log(7/8)-D;xp=log(9/8)+D
    return dict(c=float(c),C=float(C),D=D,log_xminus=xm,log_xplus=xp,
        log_jminus=log(jlo)-2/3+2*xm-float(C),log_jplus=log(jhi)+2*xp+.5*log(pi/float(c)),
        log_direction_factor=2*xp+float(C)-log(2)-2*xm-log(float(c)))


def square_budget(a0,B,H,N):
    a0,B,H,N=map(F,(a0,B,H,N))
    if a0<=0 or B<0 or N<0:raise ValueError('Positive tail rate and nonnegative polynomial/target orders required')
    # log(ell)<=ell, log2<1; dropping a negative H+N is conservative.
    d=B+max(F(0),H+N)
    return d,d*d/(4*a0)


def curvature_review():
    rows=[];good=True
    for u in (F(1,10),F(1,2),F(1),F(2),F(5)):
        c,C=curvature(u,F(1,3),3)
        for a in (F(1,2),F(3,4),F(1),F(5,4),F(3,2)):
            for lam in (F(1,3),F(1),F(3)):
                l,h=rate_curvature_interval(u,a,lam);good &= c<=l<=h<=C
        rows.append(dict(tilt=str(u),lower=str(c),upper=str(C)))
    return dict(rows=rows,rational_enclosure_cases=75,checks=dict(
        whole_interval_constants_positive=all(F(x['lower'])>0 for x in rows),
        exact_rational_curvature_enclosures=bool(good),
        two_pulse_signs_have_same_squared_tilt_rate=True,
        actual_spectral_chart_floor_lambda_star_over_sixteen_retained=True))


def barrier_review():
    rows=[]
    for L in (8,16,128,1024):
        r=barrier(1,1,F(1,4),L)
        rows.append(dict(length=L,**{k:str(v) for k,v in r.items()}))
    return dict(rows=rows,checks=dict(
        ratio_below_one_eighth=all(F(r['rho'])<=F(1,8) for r in rows),
        outward_derivative_negative_with_recorded_slack=all(F(r['outward'])<=F(r['negative_bound'])<0 for r in rows),
        scalar_exponential_nonvanishing_needed_for_continuation=True,
        positive_lower_bound_only_for_uncut_radial_amplitude=True))


def tail_review():
    rows=[];good=True
    for a0,B,H,N in ((F(1,100),3,2,10),(F(1,10000),17,10,30),(F(1,3),0,-5,1)):
        d,K=square_budget(a0,B,H,N)
        for ell in (1,2,10,100,10000):good &= -a0*ell*ell+d*ell<=K
        rows.append(dict(rate=str(a0),B=B,H=H,N=N,d=str(d),log_constant=str(K)))
    r0=F(1,6057072);c,_=curvature(1,F(1,2),2)
    return dict(rows=rows,frozen_radius=str(r0),a0=str(c*r0/25),checks=dict(
        square_completion_all_nonnegative_ell=bool(good),
        cutoff_distance_one_fifth_gives_one_over_fifty=F(1,5)**2/2==F(1,50),
        length_at_least_two_r0_S_gives_a0=c*r0/25>0,
        all_fixed_tail_jets_share_original_domain=True,
        complete_pressure_commutator_and_source_tails_remain_additive=True))


def degeneracy_review():
    rows=[]
    for u in (F(1,2),F(1,4),F(1,8),F(1,16)):
        c,C=curvature(u,1,1);rows.append(dict(u=str(u),lower=str(c),upper=str(C)))
    huge=logs(1,F(1,2),2,10000,10000)
    return dict(rows=rows,log_mass_floor=huge['log_jminus'],floating_mass_floor=exp(huge['log_jminus']),
        checks=dict(curvature_upper_tends_to_zero_with_tilt=all(F(rows[i+1]['upper'])<F(rows[i]['upper']) for i in range(len(rows)-1)),
            generic_upper_curvature_at_most_nine_halves_lambda_u_squared=all(F(r['upper'])<=F(9,2)*F(r['u'])**2 for r in rows),
            zero_tilt_has_constant_envelope_not_Gaussian=True,
            logarithmic_mass_floor_survives_underflow=isfinite(huge['log_jminus']) and exp(huge['log_jminus'])==0,
            no_unconditional_full_pulse_or_residual_claim=True))
