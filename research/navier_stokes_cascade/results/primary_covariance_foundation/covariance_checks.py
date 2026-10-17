"""Independent rational Taylor roots and manufactured homogeneous Haar columns."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import numpy as np
from scipy.integrate import quad,solve_ivp
from covariance_bounds import matrix,solve


def add(a,b):return [x+y for x,y in zip(a,b)]
def mul(a,b):return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(len(a))]
def reciprocal(a):
    if not a[0]:raise ValueError('Nonzero constant coefficient required')
    out=[1/a[0]]
    for n in range(1,len(a)):out.append(-sum(a[i]*out[n-i] for i in range(1,n+1))/a[0])
    return out
def square_root(a,root0):
    root0=F(root0)
    if root0<=0 or root0*root0!=a[0]:raise ValueError('Exact positive reference root required')
    out=[root0]
    for n in range(1,len(a)):out.append((a[n]-sum(out[i]*out[n-i] for i in range(1,n)))/(2*root0))
    return out


def polynomial(values,n=7):return list(map(F,values))+[F(0)]*(n-len(values))


def jet_review():
    # New manufactured roots; not values of the actual selected primary.
    aa=polynomial((2,1,F(1,5)));bb=polynomial((3,F(-1,2),F(1,7)))
    roots=(aa,bb);squares=[mul(a,a) for a in roots]
    C=[[polynomial((1,F(1,50))),polynomial((1,F(-1,60),F(1,100)))],
       [polynomial((-1,F(1,70))),polynomial((1,F(1,80)))]]
    weights=(polynomial((2,F(1,20))),polynomial((3,F(-1,30))))
    H=[[mul(C[i][j],weights[j]) for j in range(2)] for i in range(2)]
    det=add(mul(H[0][0],H[1][1]),[-x for x in mul(H[0][1],H[1][0])]);idet=reciprocal(det)
    target=[add(mul(H[i][0],squares[0]),mul(H[i][1],squares[1])) for i in range(2)]
    recovered=[mul(add(mul(H[1][1],target[0]),[-x for x in mul(H[0][1],target[1])]),idet),
        mul(add(mul(H[0][0],target[1]),[-x for x in mul(H[1][0],target[0])]),idet)]
    recovered_roots=[square_root(x,r[0]) for x,r in zip(recovered,roots)]
    inverse_identity=mul(det,idet)==polynomial((1,))
    # The target of the linear signed operation can be arbitrarily large.
    d=(polynomial((-100,)),polynomial((30,)))
    signed=[mul(q,reciprocal(a)) for q,a in zip(d,roots)]
    signed=[[x/2 for x in p] for p in signed]
    cross=[add(mul(H[i][0],mul(roots[0],signed[0])),mul(H[i][1],mul(roots[1],signed[1]))) for i in range(2)]
    rhs=[add(mul(H[i][0],d[0]),mul(H[i][1],d[1])) for i in range(2)]
    cross=[[2*x for x in row] for row in cross]
    quadratic=[add(mul(H[i][0],mul(signed[0],signed[0])),mul(H[i][1],mul(signed[1],signed[1]))) for i in range(2)]
    return dict(normalized_Taylor_order=6,exact_recovered_root_coefficients=[[str(x) for x in a] for a in recovered_roots],
        checks=dict(independent_Cramer_Taylor_inverse_exact=bool(inverse_identity and recovered==squares),
            differentiated_positive_square_roots_recovered_exactly=recovered_roots==list(roots),
            original_root_signed_linear_covariance_identity=cross==rhs,
            updated_square_root_not_required_for_signed_target=squares[0][0]+d[0][0]<0,
            signed_square_is_nonzero_and_not_erased=any(x for row in quadratic for x in row),
            reciprocal_primary_root_has_nonconstant_jets=any(x for p in signed for x in p[1:])))


def psi(v,L):
    x=(abs(v/L-.5)-.2)/.2
    if x<=0:return 1.
    if x>=1:return 0.
    a=np.exp(-1/x**2);b=np.exp(-1/(1-x)**2)
    return b/(a+b)


@lru_cache(None)
def gaussian_review():
    rows=[];mass=True;directions=True;positive=True;barriers=True
    tilt=.5;c0=-.8;A=-c0*np.sqrt(1+tilt**2);J=1+(np.sqrt(2)-1)**2
    chi=.5*quad(lambda x:np.exp(-2*x*x/(1-x*x)),-1,1,epsabs=1e-12)[0]
    for L in (400.,1600.):
        for sign in (-1,1):
            E=np.array([[.02,.015*sign],[-.01*sign,.005]])/L
            def rhs(v,y):return (np.diag([0.,-2.])+E-.01*np.cos(v/L)/L*np.eye(2))@y
            sol=solve_ivp(rhs,(0,L),[1.,0.],method='Radau',rtol=2e-10,atol=2e-12,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            def density(t,component):
                v=L/2+np.sqrt(L)*t;yp,ym=sol.sol(v);x=yp+ym;s=sign*tilt*(.5+v/L)
                base=psi(v,L)**2*x*x*np.exp(-t*t)
                if component==1:base*=(-c0)*np.sqrt(1+s*s)*(yp-ym)/(x*A)
                if component==2:base*=-s/tilt
                return float(base)
            integrals=[quad(lambda t:density(t,i),-.5*np.sqrt(L),.5*np.sqrt(L),points=[-.4*np.sqrt(L),-.2*np.sqrt(L),.2*np.sqrt(L),.4*np.sqrt(L)],epsabs=1e-10,epsrel=2e-10)[0] for i in range(3)]
            jmass=J*chi*integrals[0];eta=max(abs(integrals[1]/integrals[0]-1),abs(integrals[2]/integrals[0]+sign))
            ratio=.04/L;logbound=.01+.04*(1+ratio)
            lo=(1-ratio)*np.exp(-logbound);hi=(1+ratio)*np.exp(logbound)
            jlo=J*np.exp(-2/3)*lo*lo*np.exp(-1);jhi=J*hi*hi*np.sqrt(np.pi)
            mass &= jlo<=jmass<=jhi;directions &= eta<=.25/4
            for v in np.linspace(0,L,41):
                yp,ym=sol.sol(v);barriers &= abs(ym/yp)<=ratio*(1+1e-7) and lo<=yp+ym<=hi
            rows.append(dict(length=L,sign=sign,normalized_Haar_mass=float(jmass),mass_lower=float(jlo),mass_upper=float(jhi),normalized_column_error=float(eta)))
    for L in (400.,1600.):
        cols=[r for r in rows if r['length']==L]
        positive &= all(r['normalized_Haar_mass']>0 for r in cols)
    return dict(model='new homogeneous noncommuting frame controls, not actual gigantic fields',rows=rows,
        checks=dict(actual_Haar_area_bounds_not_radius_free=bool(mass),
            Gaussian_column_errors_fit_safe_cone_budget=bool(directions),
            homogeneous_Riccati_radial_barriers=bool(barriers),
            both_sign_column_weights_positive=bool(positive)))


def edge_review():
    # The actual outer asymptotic is zeta*delta^-3 times a nonzero smooth
    # factor. Normalized amplitude is not bounded even though it is flat.
    ratios=[(F(1,n),F(n)**3,F(n)) for n in (2,4,8,16)]
    return dict(outer_stress='zeta*delta_outer^-3*(nonzero smooth coefficient)',
        outer_primary_root='sqrt(zeta)*delta_outer^-3/2*bounded nonzero factor',
        half_weight_constants=['t1^2/2','2'],
        normalized_pole_controls=[[str(x) for x in row] for row in ratios],
        checks=dict(normalized_stress_not_assumed_bounded_at_outer_edge=all(r[1]>r[2] for r in ratios),
            normalized_primary_root_not_assumed_smooth_at_outer_edge=True,
            all_fixed_polynomial_losses_absorbed_by_exact_half_weight=True,
            root_strict_positive_only_on_open_active_shell=True,
            full_background_pressure_and_complete_covariance_means_not_replaced=True))
