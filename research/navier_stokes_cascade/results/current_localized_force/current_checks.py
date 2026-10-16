"""Manufactured actual-q localization and independent forced-pressure controls."""
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np
from forcing import TERMS

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'current_endpoint_jets'))
from endpoint_checks import field,curl
sys.path.pop(0)
from jets import Jet,arrays,cross,dot,variables
import waves


def chi(s):
    if s.value<=.5:return Jet(1)
    if s.value>=1:return Jet(0)
    left=(-1/(s-.5)).exp();right=(-1/(1-s)).exp()
    return right/(left+right)


def cutoff(coords):
    x,y,z,t=coords
    return chi((x*x+y*y)/.55**2)*chi(z*z/.4**2)*chi((1-t)/.0015)


def gradient(f):return [f.diff(i) for i in range(3)]
def laplacian(f):return sum(f.diff(i).diff(i) for i in range(3))
def relative(a,b):return float(np.linalg.norm(a-b)/(1+np.linalg.norm(a)+np.linalg.norm(b)))


def momentum(U,p,nu):
    return U[1][:,3]+waves.advection(U,U,None,True)-nu*waves.laplacian(U,None,True)+waves.pressure_gradient(p,None,True)


@lru_cache(None)
def localization_stats():
    z0=.08**.495
    points=((.4,.1,z0,1.),(.4,.1,z0,.9989),(.39,.13,-z0,.999),(.43,-.03,.25,.9991))
    error=product=div=orth=0.;naive=angular=0.;sizes={k:0. for k in TERMS}
    for point in points:
        coords=variables(point);data=field(coords)
        A,B,p=data['potential'],data['direct'],data['pressure'];u0=data['velocity']
        c=cutoff(coords);gc=gradient(c);d=cross(gc,A)
        u=[a+c*b for a,b in zip(curl([c*a for a in A]),B)]
        U,D,V=arrays(u0),arrays(d),arrays(u)
        product=max(product,*(relative(a,b) for a,b in zip(V,arrays([c*a+b for a,b in zip(u0,d)]))))
        div=max(div,abs(waves.divergence(V,None,True))/(1+np.linalg.norm(V[1])))
        orth=max(orth,abs(dot(d,gc).value)/(1+np.linalg.norm(arrays(d)[0])*np.linalg.norm(arrays(gc)[0])))
        naive=max(naive,abs(waves.divergence(arrays([c*a for a in u0]),None,True)))
        bad=c*(1+.1*coords[0])
        bad_u=[a+bad*b for a,b in zip(curl([bad*a for a in A]),B)]
        angular=max(angular,abs(waves.divergence(arrays(bad_u),None,True)))
        for nu in (.6,1.,1.7):
            terms=dict(base_residual=c.value*momentum(U,p,nu),
                cutoff_time=c.diff(3).value*U[0],curl_time=D[1][:,3],
                cutoff_laplacian=-nu*laplacian(c).value*U[0],
                gradient_viscosity=-2*nu*waves.advection(arrays(gc),U,None,True),
                curl_viscosity=-nu*waves.laplacian(D,None,True),
                pressure_cutoff=p.value*arrays(gc)[0],
                quadratic_rescaling=c.value*(c.value-1)*waves.advection(U,U,None,True),
                velocity_cutoff_transport=c.value*np.dot(U[0],arrays(gc)[0])*U[0],
                old_curl_transport=c.value*waves.advection(U,D,None,True),
                curl_old_transport=c.value*waves.advection(D,U,None,True),
                curl_square=waves.advection(D,D,None,True))
            assert tuple(terms)==TERMS
            direct=momentum(V,c*p,nu)
            error=max(error,relative(direct,sum(terms.values())))
            for name,value in terms.items():sizes[name]=max(sizes[name],float(np.linalg.norm(value)))
    return dict(identity_error=error,velocity_product_error=product,divergence_error=div,
        orthogonality_error=orth,naive_divergence=naive,nonaxisymmetric_direct_swirl_divergence=angular,
        retained_term_sizes=sizes,points=len(points),positive_viscosities=[.6,1.,1.7])


@lru_cache(None)
def forced_pressure_stats():
    # This polynomial potential has total degree <=3, so the Taylor algebra
    # contains every derivative needed for the value of div(f), exactly.
    coords=variables((.3,.2,.1,.6));x,y,z,t=coords
    A=[Jet(0),Jet(0),x*x*y-y*y*y/3+.2*t*y]
    u=curl(A);p=.3*x*z+.1*x*x+.2*t*y
    f=[v.diff(3)+sum(u[j]*v.diff(j) for j in range(3))-laplacian(v)+p.diff(i)
        for i,v in enumerate(u)]
    divf=sum(f[i].diff(i) for i in range(3)).value
    source=sum((u[i]*u[j]).diff(i).diff(j) for i in range(3) for j in range(3)).value
    left=-laplacian(p).value
    return dict(forced_Poisson_error=abs(left-source+divf),omitted_divergence_gap=abs(left-source),
        force_divergence=divf,velocity_divergence=sum(u[i].diff(i) for i in range(3)).value)


def diagnostic_review():
    a,b=localization_stats(),forced_pressure_stats()
    r=dict(manufactured_complete_fields_not_actual_global_flow=True,actual_implicit_q_and_terminal_points_used=True,
        localized_point_count=a['points'],positive_viscosities=a['positive_viscosities'],
        retained_term_names=list(a['retained_term_sizes']),diagnostic_algebra_viscosity_other_than_one_not_the_actual_construction=True,
        tolerances=dict(relative_identity=2e-10,relative_divergence=2e-12,forced_Poisson=2e-12),
        raw_roundoff_scores_tested_not_refrozen=True,
        checks=dict(full_current_q_localization_identity=a['identity_error']<2e-10,
            two_velocity_derivative_product_rule=a['velocity_product_error']<2e-10,
            exact_localized_divergence=a['divergence_error']<2e-12,
            all_twelve_terms_have_nonzero_omission_controls=min(a['retained_term_sizes'].values())>1e-7,
            excluded_d_dot_grad_c_is_identically_zero=a['orthogonality_error']<2e-12,
            naive_velocity_cutoff_fails=a['naive_divergence']>1e-4,
            theta_dependent_direct_swirl_cutoff_fails=a['nonaxisymmetric_direct_swirl_divergence']>1e-4,
            independent_forced_Poisson_identity=b['forced_Poisson_error']<2e-12,
            dropping_force_divergence_is_detected=b['omitted_divergence_gap']>1e-3,
            divergence_free_velocity_does_not_require_divergence_free_force=abs(b['velocity_divergence'])<2e-12 and abs(b['force_divergence'])>1e-3))
    r['checks']={k:bool(v) for k,v in r['checks'].items()}
    return r
