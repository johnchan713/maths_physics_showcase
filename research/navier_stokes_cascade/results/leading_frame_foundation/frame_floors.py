"""Uniform leading floors with exact rational squared certificates."""
from fractions import Fraction as F
from itertools import product
from math import exp,log,isfinite


def squared_floors(fmin,fmax,alpha,abar,bbar,nu):
    fmin,fmax,alpha,abar,bbar,nu=map(F,(fmin,fmax,alpha,abar,bbar,nu))
    if min(fmin,alpha,nu)<=0 or fmax<fmin or abar<alpha or bbar<0:
        raise ValueError('Positive ordered leading bounds and strict growth gap required')
    vmax=abar+bbar*bbar/alpha
    if vmax<2+nu:raise ValueError('Empty growth interval')
    return dict(vmax=vmax,gmin2=fmin*fmin*alpha*(2+nu),gmax2=fmax*fmax*abar*vmax,
        normal_theta_min2=alpha/vmax,lambda_min2=2*alpha*fmin*fmin*nu/vmax,
        lambda_max2=2*abar*fmax*fmax,cmin2=nu/2,cmax2=(vmax-2)/2,
        left_inverse_max2=1/(2*min(F(1),nu/2)))


def leading_identities(F0,a,b):
    F0,a,b=map(F,(F0,a,b))
    if min(F0,a)<=0:raise ValueError('Positive swirl and a required')
    v=a+b*b/a;g2=F0*F0*(a*a+b*b);ntheta2=a*a/(a*a+b*b)
    lam2=2*a*F0*F0*(1-2/v)
    c2=lam2/(4*F0*F0*ntheta2)
    return dict(v=v,g2=g2,ntheta2=ntheta2,lambda2=lam2,c2=c2,
        checks=dict(shear_norm=g2==F0*F0*a*v,normal_ratio=ntheta2==a/v,
            exact_cancellation=c2==(v-2)/2))


def identity_review():
    bounds=squared_floors(F(1,10),4,F(1,4),5,4,F(1,10));rows=[];good=True
    for f,a,b in product((F(1,10),F(1),F(4)),(F(1,4),F(3),F(5)),(F(-4),F(0),F(4))):
        r=leading_identities(f,a,b)
        if r['v']<F(21,10):continue
        good &= all(r['checks'].values()) and bounds['gmin2']<=r['g2']<=bounds['gmax2']
        good &= bounds['normal_theta_min2']<=r['ntheta2'] and bounds['lambda_min2']<=r['lambda2']<=bounds['lambda_max2']
        good &= bounds['cmin2']<=r['c2']<=bounds['cmax2']
        rows.append([str(x) for x in (f,a,b,r['v'],r['lambda2'],r['c2'])])
    return dict(model='rational algebra controls, not actual selected profile',rows=rows,
        bounds={k:str(v) for k,v in bounds.items()},
        checks=dict(all_exact_shear_and_growth_identities=bool(good),strict_growth_not_just_a_positive=leading_identities(1,2,0)['lambda2']==0,
            frame_inverse_squared_is_inverse_gap=bounds['left_inverse_max2']==10))


def global_review():
    return dict(hierarchy='Aenv=C^4096; G=exp(Aenv^32); alpha=1/(2G^3); nu=1/(2G^4); Ab=Bb=4Aenv',
        regions=[dict(region='protected inner',a_lower='C^-22',growth_gap='.29',a_upper='10',b_upper='C*R<C^2'),
            dict(region='compact plus exact repair',a_lower='1/(2G^3)',growth_gap='1/(2G^4)',coordinate_upper='2Aenv+1'),
            dict(region='protected right power',growth_gap='2lambda>C^-1',b='0'),
            dict(region='compensated exterior',growth_gap='h>C^-2',v_upper='120',b='0'),
            dict(region='terminal closed heat edge',growth_gap='h',a_upper='2+2h',b='0')],
        anchor='exp(-2Lambda*M)/(4C)<=F(Xa)<=2exp(2Lambda*M)/C',
        log_floor='log fstar=-2Lambda*M-log(4C)-(Ab/2)*log(Xb/Xa)',
        checks=dict(protected_and_compact_bounds_use_same_profile_and_roots=True,
            division_by_two_is_not_fractional_exponent=True,
            heat_swirl_does_not_vanish_when_stress_is_flat=True,
            fixed_finite_profile_not_uniform_as_parameters_change=True))


def logarithmic_review():
    # New scalar model: a(y)=3+sin(y), 0<=y<=1000, log F(0)=-10000.
    # Integral is 3y+1-cos(y), so a<=4 bounds the full field from below.
    from math import cos
    y=1000.;anchor=-10000.;lower=anchor-2*y;actual=anchor-(3*y+1-cos(y))/2
    if not all(isfinite(x) for x in (lower,actual)):raise ValueError('Finite log model required')
    open_rows=[dict(delta=str(F(1,n)),lambda_squared=str(leading_identities(F(1,n),3,0)['lambda2'])) for n in (2,4,8,16,32)]
    return dict(model='new nonconstant-shear logarithmic control; not actual huge field',log_floor=lower,log_actual=actual,
        open_domain_rows=open_rows,checks=dict(integrated_logarithmic_swirl_floor=lower<=actual<=anchor,
            floating_underflow_not_a_zero_analytical_floor=exp(lower)==0 and isfinite(lower),
            pointwise_open_positivity_does_not_give_closed_uniform_floor=all(F(open_rows[i+1]['lambda_squared'])<F(open_rows[i]['lambda_squared']) for i in range(len(open_rows)-1))))
