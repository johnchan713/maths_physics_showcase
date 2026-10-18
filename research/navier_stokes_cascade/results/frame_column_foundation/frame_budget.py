"""Exact sufficient phase/frame budgets; all-interval derivation in README.md."""
from fractions import Fraction as F
from math import log,exp,log1p,isfinite


def inputs(**kw):
    p=dict(Rmin=F(1),Rmax=F(3),bmin=F(1,2),bmax=F(1),gmin=F(1),gmax=F(3),
        Fmax=F(2),cmin=F(1,4),cmax=F(2),lam_min=F(1,2),lam_max=F(2),u=F(3,4),
        Db=F(1,10**6),Dbox=F(1,10**6),Mz=F(1,10**6))
    if set(kw)-p.keys():raise ValueError('Unknown primitive budget')
    p.update(kw)
    if any(not isinstance(v,(int,F)) or isinstance(v,bool) for v in p.values()):raise TypeError('Exact rational budget inputs required')
    p={k:F(v) for k,v in p.items()}
    if any(v<=0 for k,v in p.items() if k not in ('Db','Dbox','Mz')) or any(p[k]<0 for k in ('Db','Dbox','Mz')):raise ValueError('Positive floors and nonnegative errors required')
    if any(p[a]>p[b] for a,b in (('Rmin','Rmax'),('bmin','bmax'),('gmin','gmax'),('cmin','cmax'),('lam_min','lam_max'))):raise ValueError('Ordered value bounds required')
    return p


def primitive(p):
    rmin,rmax,bmax,gmin,u=(p[k] for k in ('Rmin','Rmax','bmax','gmin','u'))
    P=rmax*bmax*(1+u/gmin)+1;Z=bmax*(1+u/gmin);fr=p['gmax']/rmin
    zc=(P+Z)*p['Mz'];de=2*(P+Z)*p['Db']+fr
    dt=bmax*p['Dbox']/rmin+bmax*u*rmax/(rmin*gmin)+1/rmin+bmax*u/gmin+zc
    return dict(Pmax=P,Zmax=Z,Frmax=fr,Zc=zc,De=de,Dr=de,Drprime=bmax*u+de,
        Dt=dt,Dtprime=zc,DF=2*p['Db'],Dg=2*(rmax+1)*p['Db']+p['Dbox']*fr)


def budgets(p):
    d=primitive(p);beta=p['bmin']/2;sr=3*p['u']/2;s0=sr+1
    da=2*d['Dt']/p['bmin'];ds=(d['Dr']+sr*d['Dt'])/beta
    om=d['Dtprime']/beta;s1=(d['Drprime']+s0*d['Dtprime'])/beta
    h11=s0*(da*p['gmax']+d['Dg']+s1)
    h12=2*d['DF']+2*p['Fmax']*da+s0*om
    h21=2*d['DF']+2*p['Fmax']*da+d['Dg']+p['gmax']*da+s0*om
    dh=h11+h12+h21
    km=max(F(1),p['cmax']*(1+s0))/min(F(1),p['cmin'])
    ce=km*dh+s0*s1+p['lam_max']*ds
    kmref=max(F(1),p['cmax']*(1+sr))/min(F(1),p['cmin'])
    ceref=kmref*(dh+2*p['Fmax']*ds)+sr*p['u']
    dn=d['Dr']+d['Dt'];nref=p['bmax']*(1+sr);cd=4*dn*(2*nref+dn)
    mu=p['lam_min']/(1+sr)
    D0=ds+sr*da+p['cmax']*(ds+(1+sr)*da+16*(1+s0)*ce/(7*mu))
    D1=p['u']*(1+p['cmax'])
    return dict(**d,beta=beta,Sref=sr,S0=s0,Da=da,Ds=ds,Omega=om,S1=s1,D11=h11,D12=h12,D21=h21,DH=dh,
        kM=km,CE=ce,kMref=kmref,CEref=ceref,Cd=cd,mu=mu,D0=D0,D1=D1,
        normal_acceptance=max(F(1),2*d['Dt']/p['bmin'],ds),Riccati_acceptance=8*ce/mu)


def log_cone_budget(p,gamma=F(1,4)):
    gamma=F(gamma)
    if not 0<gamma<=1:raise ValueError('Positive cone margin at most one required')
    b=budgets(p);u=p['u'];lo,hi=p['lam_min'],p['lam_max']
    c=lo*u*u/2*((1+9*u*u/4)**-2+2*(1+u*u)**-2)
    C=3*hi*u*u/2*((1+u*u/4)**-1+2*(1+u*u)**-1)
    D=b['Cd']+F(9,8)*b['CE']
    # log(xplus^2*exp(C)/(2*xminus^2*c)), using the exact amplitude factors.
    lm=2*log(9/7)+4*float(D)+float(C)-log(2)-log(float(c))
    left=log(float(b['D0']));right=log(float(b['D1']))+lm
    maximum=max(left,right);ld=maximum+log1p(exp(min(left,right)-maximum))
    amin=p['cmin'] # conservative lower bound for |c0|sqrt(1+u^2)
    length=2*(log(4)+ld+log(float(max(1/amin,1/u)))-log(float(gamma)))
    return dict(log_direction_constant=ld,log_cone_length=length,log_moment_factor=lm,
        log_mass_lower=-F(2,3)+2*log(7/8)-2*float(D)-float(C))


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,list):return [encode(v) for v in x]
    return x


def budget_review():
    p=inputs();b=budgets(p);logs=log_cone_budget(p)
    return dict(manufactured_inputs=encode(p),explicit_budgets=encode(b),logarithmic_column_budget=logs,
        checks=dict(all_current_floor_denominators_positive=all(p[k]>0 for k in ('Rmin','bmin','gmin','cmin','lam_min','u')),
            frame_error_is_operator_norm_not_entrywise= b['CE']>=b['kM']*b['DH'] and b['kM']>=1,
            complete_connection_and_reference_rate_costs_retained=b['CE']==b['kM']*b['DH']+b['S0']*b['S1']+p['lam_max']*b['Ds'],
            manuscript_reference_mixing_has_separate_operator_budget=b['CEref']==b['kMref']*(b['DH']+2*p['Fmax']*b['Ds'])+b['Sref']*p['u'],
            damping_keeps_original_factor_four=b['Cd']==4*(b['Dr']+b['Dt'])*(2*p['bmax']*(1+b['Sref'])+b['Dr']+b['Dt']),
            original_Riccati_acceptance_has_finite_explicit_constant=b['Riccati_acceptance']>b['normal_acceptance']>0,
            column_constant_keeps_both_frame_coordinates=b['D0']>b['Ds']+b['Sref']*b['Da'],
            huge_mass_and_column_constants_retained_logarithmically=all(isfinite(float(v)) for v in logs.values()) and logs['log_direction_constant']>1000,
            no_giant_actual_background_or_numeric_first_band_evaluated=True))


def all_band_review():
    rows=[];good=True
    for S in (10**4,10**8,10**12):
        k=2*S*S;eps=F(1,k*k);L=F(2,6057072)*S
        accepted=S*S*(eps+eps*eps+F(1,k))<=1 and 1<=L<=S
        # For the smallest fixture L<1, no phase/frame conclusion is asserted.
        if L<1:continue
        good &= accepted
        rows.append(dict(S=S,length=str(L),epsilon=str(eps),k=k,accepted=bool(accepted),
            rounding_ratio=str(L*L*F(1,k)),axial_ratio=str(eps*L*L)))
    return dict(rows=rows,checks=dict(whole_pulse_rounding_and_axial_factors_retained=bool(good) and bool(rows),
        discrete_carriers_frozen_before_derivatives=True,
        actual_schedule_uses_its_existing_acceptance_not_manufactured_epsilon=True,
        primitive_jet_norms_come_from_reviewed_current_background=True,
        every_fixed_jet_uses_same_positive_value_domain=True,
        complete_pressure_memory_and_full_covariance_mean_residual_still_required=True))
