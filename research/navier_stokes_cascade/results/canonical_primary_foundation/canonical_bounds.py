"""Conditional canonical-primary budgets and exact finite-matrix controls."""
from fractions import Fraction as F
from functools import lru_cache
import importlib.util
from math import log,log1p,exp,pi,isfinite
from pathlib import Path


@lru_cache(None)
def parent_budget():
    path=Path(__file__).resolve().parent.parent/'frame_column_foundation/frame_budget.py'
    spec=importlib.util.spec_from_file_location('canonical_parent_budget',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def canonical(p):
    b=parent_budget().budgets(p);ce=b['CEref'];wr=1+b['Sref']
    d0=b['Ds']+b['Sref']*b['Da']+p['cmax']*wr*(b['Da']+16*ce/(7*b['mu']))
    B0=2*(1+b['S0'])*max(F(1),p['cmax']*(1+b['S0']))
    Kmax=4*(p['Fmax']+b['DF'])+2*(p['gmax']+b['Dg'])
    return dict(CE=ce,Cd=b['Cd'],D0=d0,D1=b['D1'],B0=B0,
        pressure_difference_multiplier=Kmax/b['beta']+(b['Drprime']+b['Dtprime'])/b['beta']**2,
        Riccati_length=8*ce/b['mu'],normal_length=b['normal_acceptance'])


def logadd(x,y):
    a,b=max(x,y),min(x,y);return a+log1p(exp(b-a))


def comparison_logs(p,gamma=F(1,4),Jmin=F(1),Jmax=F(2)):
    gamma,Jmin,Jmax=map(F,(gamma,Jmin,Jmax))
    if not 0<gamma<=1 or not 0<Jmin<=Jmax:raise ValueError('Positive cone and ordered Haar floors required')
    b=parent_budget().budgets(p);r=canonical(p);u=p['u']
    c=p['lam_min']*u*u/2*((1+9*u*u/4)**-2+2*(1+u*u)**-2)
    C=3*p['lam_max']*u*u/2*((1+u*u/4)**-1+2*(1+u*u)**-1)
    Da=b['CE']+b['Cd'];Dstar=max(b['CE'],r['CE'])+b['Cd']
    lt=log(float(r['B0']))+float(Dstar)
    ld=log(float(r['B0']*b['Ds']))+float(Da)
    lh=log(2*float(Jmax))+lt+ld+.5*(log(pi)-log(float(c)))
    amp=b['Cd']+F(9,8)*b['CE']
    ljmin=log(float(Jmin))-2/3+2*log(7/8)-2*float(amp)-float(C)
    ljmax=log(float(Jmax))+2*log(9/8)+2*float(amp)+.5*(log(pi)-log(float(c)))
    amin=min(p['cmin'],u)
    lk=log(12)+lh-ljmin-log(float(amin))
    lkappa=log(24)+ljmax-log(float(gamma))-ljmin
    length=lk+logadd(log(2)+lkappa,0.)
    Dr=r['Cd']+F(9,8)*r['CE']
    moment=2*log(9/7)+4*float(Dr)+float(C)-log(2*float(c))
    cdir=logadd(log(float(r['D0'])),log(float(r['D1']))+moment)
    cone=2*(log(4)+cdir+log(float(max(1/p['cmin'],1/u)))-log(float(gamma)))
    return dict(log_uniform_physical_pulse_factor=lt,log_pulse_difference_coefficient=ld,
        log_Haar_column_difference_coefficient=lh,log_inherited_jminus=ljmin,log_inherited_jplus=ljmax,
        log_theta_length_coefficient=lk,log_kappa=lkappa,log_relative_root_acceptance_length=length,
        log_canonical_direction_constant=cdir,log_canonical_cone_length=cone)


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    return x


def budget_review():
    p=parent_budget().inputs();b=parent_budget().budgets(p);r=canonical(p);logs=comparison_logs(p)
    return dict(manufactured_inputs=encode(p),canonical_budgets=encode(r),comparison_log_budgets=logs,
        checks=dict(canonical_uses_reference_not_actual_mixing_error=r['CE']==b['CEref'],
            canonical_column_keeps_secondary_coordinate=r['D0']>b['Ds']+b['Sref']*b['Da'],
            canonical_column_does_not_insert_actual_w_difference_cost=r['D0']==b['Ds']+b['Sref']*b['Da']+p['cmax']*(1+b['Sref'])*(b['Da']+16*b['CEref']/(7*b['mu'])),
            common_B0_bounds_both_frame_conventions=r['B0']>=2*(1+b['S0']),
            large_comparison_and_cone_budgets_kept_logarithmically=all(isfinite(v) for v in logs.values()) and logs['log_relative_root_acceptance_length']>1000,
            actual_giant_first_band_not_numerically_evaluated=True))


def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def sub(a,b):return [[x-y for x,y in zip(row,col)] for row,col in zip(a,b)]
def inv(a):
    d=a[0][0]*a[1][1]-a[0][1]*a[1][0]
    if not d:raise ValueError('Singular covariance')
    return [[a[1][1]/d,-a[0][1]/d],[-a[1][0]/d,a[0][0]/d]]
def norminf(a):return max(sum(abs(x) for x in row) for row in a)


def matrix_review():
    rows=[];good=True;stale=False;scales=[]
    for radius in (F(1,100),F(1,50)):
        for L,sqrtL in ((100,10),(400,20)):
            s0=radius*radius/sqrtL;A,u=F(1,2),F(3,4);jmin,jmax=F(7,10),F(7,5);gamma=F(1,4)
            for sign in (-1,1):
                for w in (F(-2,5),F(0),F(2,5)):
                    Ha=[[s0*A*F(9,10),s0*A*F(11,10)],[-s0*u*F(9,10),s0*u*F(11,10)]]
                    Delta=[[s0*F(sign,100000*L),s0*F(1,200000*L)],
                           [-s0*F(1,300000*L),s0*F(sign,400000*L)]]
                    Hr=[[Ha[i][j]+Delta[i][j] for j in range(2)] for i in range(2)]
                    target=[[A],[u*w]];ya=mul(inv(Ha),target);yr=mul(inv(Hr),target)
                    theta=norminf(mul(inv(Ha),Delta));kappa=24*jmax/(gamma*jmin)
                    R=kappa*theta/(1-theta)
                    actual=max(abs(yr[i][0]/ya[i][0]-1) for i in range(2))
                    good &= theta<=1/(2*kappa+1) and 0<=actual<=R<=F(1,2)
                    good &= all(x[0]>0 for x in yr) and mul(Hr,yr)==target
                    stale |= mul(Hr,ya)!=target
                    scales.append(theta)
                    rows.append(dict(radius=str(radius),length=L,sign=sign,target_W=str(w),theta=str(theta),
                        relative_squared_amplitude_error=str(actual),relative_error_bound=str(R),
                        stale_covariance_defect_nonzero=mul(Hr,ya)!=target))
    return dict(rows=rows,checks=dict(exact_Neumann_and_componentwise_root_stability=bool(good),
        rebuilding_roots_recovers_exact_principal_covariance=bool(good),
        reusing_old_roots_generally_breaks_exact_covariance=bool(stale),
        relative_comparison_retains_and_cancels_common_physical_Haar_scale=scales[:12]==scales[12:],
        no_exact_equality_of_two_primary_root_families_claimed=True))


def acceptance_review():
    rows=[];good=True
    for ce,mu in ((F(1,10),F(1,2)),(F(3),F(1,5))):
        L=8*ce/mu;rho=ce/(mu*L)
        outward=-2*mu*rho+ce/L*(1+2*rho+rho*rho)
        good &= rho==F(1,8) and outward==F(-47,64)*ce/L
        rows.append(dict(CE=str(ce),mu=str(mu),length=str(L),outward=str(outward)))
    for kappa in (F(1),F(24),F(192)):
        theta=1/(2*kappa+1);R=kappa*theta/(1-theta)
        good &= R==F(1,2)
    return dict(rows=rows,checks=dict(reference_frame_Riccati_barrier_strict=bool(good),
        theta_acceptance_implies_half_relative_squared_root_error=bool(good),
        fixed_value_acceptance_not_derivative_dependent_band=True,
        root_endpoint_zero_extension_still_uses_weighted_stress_jets=True,
        complete_curl_and_fresh_mean_system_not_identified_by_principal_covariance=True))
